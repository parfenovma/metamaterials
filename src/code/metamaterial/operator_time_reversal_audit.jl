module OperatorTimeReversalAudit

using LinearAlgebra
using Statistics

export hann_burst_spectrum,
       normalized_pulse_weights,
       weighted_operator,
       singular_summary,
       observable_mode_count,
       time_reversal_kernel,
       coherence_kernel,
       normalized_focus_intensity,
       focus_metrics,
       subspace_novelty

"""Amplitude spectrum of the project-standard Hann-windowed tone burst."""
function hann_burst_spectrum(
    frequency_ratios::AbstractVector{<:Real};
    cycles::Real=5.0,
    sample_count::Integer=12_001,
)
    cycles > 0 || throw(ArgumentError("cycles must be positive"))
    sample_count >= 3 || throw(ArgumentError("sample_count must be at least three"))
    times = collect(range(0.0, Float64(cycles); length=sample_count))
    step = times[2] - times[1]
    waveform = @. sin(2pi * times) * 0.5 * (1 - cos(2pi * times / cycles))
    amplitudes = Float64[
        abs(sum(waveform .* exp.(-2pi * im * frequency * times)) * step)
        for frequency in frequency_ratios
    ]
    maximum(amplitudes) > 0 || error("pulse spectrum vanished")
    amplitudes
end

function normalized_pulse_weights(frequency_ratios; cycles=5.0)
    raw = hann_burst_spectrum(frequency_ratios; cycles)
    raw ./ norm(raw)
end

function validate_transfer_tensor(transfers, weights, state_indices)
    ndims(transfers) == 4 || throw(DimensionMismatch(
        "transfers must have output x port x state x frequency dimensions",
    ))
    size(transfers, 1) >= 1 || throw(DimensionMismatch("no output bins"))
    size(transfers, 2) >= 1 || throw(DimensionMismatch("no source ports"))
    length(weights) == size(transfers, 4) || throw(DimensionMismatch(
        "one pulse weight is required per frequency",
    ))
    isempty(state_indices) && throw(ArgumentError("at least one state is required"))
    all(index -> index in axes(transfers, 3), state_indices) ||
        throw(BoundsError(transfers, (:, :, state_indices, :)))
    all(isfinite, transfers) || throw(ArgumentError("transfers contain non-finite values"))
    all(isfinite, weights) || throw(ArgumentError("weights contain non-finite values"))
    nothing
end

"""
Concatenate the pulse-weighted output fields over states and frequencies.

Each selected state receives `1/Q` of the incident energy, hence its field is
scaled by `1/sqrt(Q)`. The resulting `A*A'` is the discrete broadband
time-reversal kernel at the output aperture.
"""
function weighted_operator(
    transfers::AbstractArray{<:Number,4},
    state_indices,
    weights::AbstractVector{<:Real},
)
    states = collect(state_indices)
    validate_transfer_tensor(transfers, weights, states)
    output_count, port_count, _, frequency_count = size(transfers)
    result = Matrix{ComplexF64}(undef,
        output_count, port_count * frequency_count * length(states),
    )
    energy_scale = inv(sqrt(length(states)))
    destination = 1
    for state in states, frequency in 1:frequency_count
        columns = destination:(destination + port_count - 1)
        result[:, columns] .= energy_scale * weights[frequency] .*
            transfers[:, :, state, frequency]
        destination += port_count
    end
    result
end

function observable_mode_count(singular_values, dynamic_range_db::Real)
    dynamic_range_db > 0 || throw(ArgumentError("dynamic range must be positive"))
    isempty(singular_values) && return 0
    threshold = maximum(singular_values) * 10.0^(-dynamic_range_db / 20)
    count(value -> value + eps(Float64) >= threshold, singular_values)
end

function singular_summary(operator, dynamic_ranges_db=(20.0, 30.0, 40.0, 60.0))
    singular_values = svdvals(operator)
    isempty(singular_values) && error("operator has no singular values")
    leading = first(singular_values)
    leading > 0 || error("operator is identically zero")
    power = abs2.(singular_values)
    probabilities = power ./ sum(power)
    entropy = -sum(value > 0 ? value * log(value) : 0.0 for value in probabilities)
    (
        singular_values,
        singular_value_db=20 .* log10.(max.(singular_values ./ leading, eps(Float64))),
        stable_rank=sum(power) / abs2(leading),
        entropy_rank=exp(entropy),
        observable_modes=Dict(
            Float64(range_db) => observable_mode_count(singular_values, range_db)
            for range_db in dynamic_ranges_db
        ),
    )
end

function time_reversal_kernel(operator)
    kernel = operator * operator'
    Matrix(Hermitian((kernel + kernel') / 2))
end

function coherence_kernel(kernel)
    size(kernel, 1) == size(kernel, 2) || throw(DimensionMismatch("kernel must be square"))
    diagonal = max.(real.(diag(kernel)), 0.0)
    result = zeros(Float64, size(kernel))
    for column in axes(kernel, 2), row in axes(kernel, 1)
        denominator = sqrt(diagonal[row] * diagonal[column])
        result[row, column] = denominator > 0 ?
            min(1.0, abs(kernel[row, column]) / denominator) : 0.0
    end
    result
end

function normalized_focus_intensity(kernel, focal_bin::Integer)
    focal_bin in axes(kernel, 2) || throw(BoundsError(kernel, (:, focal_bin)))
    target = abs2(kernel[focal_bin, focal_bin])
    target > 0 || error("zero time-reversal response at focal bin $focal_bin")
    abs2.(kernel[:, focal_bin]) ./ target
end

function contiguous_half_maximum(profile, focal_bin)
    left = focal_bin
    right = focal_bin
    while left > firstindex(profile) && profile[left - 1] >= 0.5
        left -= 1
    end
    while right < lastindex(profile) && profile[right + 1] >= 0.5
        right += 1
    end
    left:right
end

function focus_metrics(kernel, focal_bin::Integer; bin_spacing=1.0)
    bin_spacing > 0 || throw(ArgumentError("bin_spacing must be positive"))
    profile = normalized_focus_intensity(kernel, focal_bin)
    main_lobe = contiguous_half_maximum(profile, focal_bin)
    other_bins = [index for index in eachindex(profile) if index != focal_bin]
    outside_main_lobe = [index for index in eachindex(profile) if index ∉ main_lobe]
    mean_off_target = isempty(other_bins) ? 0.0 : mean(profile[other_bins])
    maximum_sidelobe = isempty(outside_main_lobe) ? 0.0 : maximum(profile[outside_main_lobe])
    (
        focal_bin,
        peak_bin=argmax(profile),
        target_is_global_maximum=profile[focal_bin] + 1.0e-12 >= maximum(profile),
        fwhm_bins=length(main_lobe),
        fwhm=length(main_lobe) * bin_spacing,
        focus_contrast_db=-10log10(max(mean_off_target, eps(Float64))),
        peak_outside_fwhm_db=10log10(max(maximum_sidelobe, eps(Float64))),
        profile,
    )
end

"""Energy fraction in `candidate` outside the observable reference subspace."""
function subspace_novelty(reference, candidate, dynamic_range_db::Real)
    decomposition = svd(reference; full=false)
    retained = observable_mode_count(decomposition.S, dynamic_range_db)
    retained == 0 && return 1.0
    basis = decomposition.U[:, 1:retained]
    residual = candidate - basis * (basis' * candidate)
    sum(abs2, residual) / max(sum(abs2, candidate), eps(Float64))
end

end
