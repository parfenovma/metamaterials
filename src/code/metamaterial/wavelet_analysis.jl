module WaveletAnalysis

using FFTW

export MorletConfig,
       CWTResult,
       WaveletComparison,
       continuous_wavelet_transform,
       compare_wavelets,
       morlet_scale_s,
       recommended_morlet_config,
       recommended_frequency_grid_hz,
       wavelet_value_at_frequency,
       wavelet_power

Base.@kwdef struct MorletConfig
    frequencies_hz::Vector{Float64} = collect(140.0e3:2.0e3:300.0e3)
    omega0::Float64 = 6.0
    truncate_sigma::Float64 = 5.0
    cone_of_influence_factor::Float64 = sqrt(2.0)
    reference_energy_floor_relative::Float64 = 1.0e-2
    amplitude_ratio_floor::Float64 = 5.0e-2
    correlation_floor::Float64 = 0.2
end

struct CWTResult
    time_s::Vector{Float64}
    frequency_hz::Vector{Float64}
    coefficients::Matrix{ComplexF64}
    cone_of_influence_valid::BitMatrix
end

struct WaveletComparison
    frequency_hz::Vector{Float64}
    reference_centroid_s::Vector{Float64}
    sample_centroid_s::Vector{Float64}
    centroid_delay_s::Vector{Float64}
    reference_peak_s::Vector{Float64}
    sample_peak_s::Vector{Float64}
    peak_delay_s::Vector{Float64}
    correlation_delay_s::Vector{Float64}
    amplitude_ratio::Vector{Float64}
    envelope_correlation::Vector{Float64}
    reference_energy_relative::Vector{Float64}
    valid::BitVector
end

function validate_config(config::MorletConfig)
    isempty(config.frequencies_hz) && throw(ArgumentError("frequency grid must not be empty"))
    all(isfinite, config.frequencies_hz) || throw(ArgumentError("frequencies must be finite"))
    all(>(0.0), config.frequencies_hz) || throw(ArgumentError("frequencies must be positive"))
    issorted(config.frequencies_hz) || throw(ArgumentError("frequencies must be sorted"))
    allunique(config.frequencies_hz) || throw(ArgumentError("frequencies must be unique"))
    config.omega0 > 0.0 || throw(ArgumentError("omega0 must be positive"))
    config.truncate_sigma >= 3.0 || throw(ArgumentError("truncate_sigma must be at least three"))
    config.cone_of_influence_factor > 0.0 ||
        throw(ArgumentError("cone_of_influence_factor must be positive"))
    0.0 <= config.reference_energy_floor_relative < 1.0 ||
        throw(ArgumentError("reference energy floor must be in [0, 1)"))
    0.0 <= config.amplitude_ratio_floor < 1.0 ||
        throw(ArgumentError("amplitude ratio floor must be in [0, 1)"))
    0.0 <= config.correlation_floor <= 1.0 ||
        throw(ArgumentError("correlation floor must be in [0, 1]"))
    nothing
end

function sample_interval(time_s::AbstractVector)
    length(time_s) >= 3 || throw(ArgumentError("at least three time samples are required"))
    all(isfinite, time_s) || throw(ArgumentError("time samples must be finite"))
    differences = diff(time_s)
    all(>(0.0), differences) || throw(ArgumentError("time samples must be strictly increasing"))
    dt = sum(differences) / length(differences)
    maximum(abs.(differences .- dt)) <= max(1.0e-10 * dt, eps(dt) * 32) ||
        throw(ArgumentError("CWT requires uniformly spaced samples"))
    dt
end

morlet_scale_s(frequency_hz::Real, omega0::Real=6.0) =
    Float64(omega0) / (2.0 * pi * Float64(frequency_hz))

"""
Build a signal-aware CWT frequency grid.

The lower bound is the first grid frequency for which the interval outside the
cone of influence occupies at least `minimum_coi_valid_fraction` of the time
record. The upper bound is limited both by `maximum_frequency_hz` and by a
conservative fraction of the Nyquist frequency.
"""
function recommended_frequency_grid_hz(
    time_s::AbstractVector;
    frequency_step_hz::Real=5.0e3,
    maximum_frequency_hz::Real=1.0e6,
    minimum_coi_valid_fraction::Real=0.5,
    omega0::Real=6.0,
    cone_of_influence_factor::Real=sqrt(2.0),
    nyquist_fraction::Real=0.9,
)
    dt = sample_interval(time_s)
    duration_s = last(time_s) - first(time_s)
    duration_s > 0.0 || throw(ArgumentError("time record must have positive duration"))
    frequency_step_hz > 0.0 || throw(ArgumentError("frequency step must be positive"))
    maximum_frequency_hz > 0.0 || throw(ArgumentError("maximum frequency must be positive"))
    0.0 <= minimum_coi_valid_fraction < 1.0 ||
        throw(ArgumentError("minimum COI-valid fraction must be in [0, 1)"))
    omega0 > 0.0 || throw(ArgumentError("omega0 must be positive"))
    cone_of_influence_factor > 0.0 ||
        throw(ArgumentError("cone-of-influence factor must be positive"))
    0.0 < nyquist_fraction < 1.0 ||
        throw(ArgumentError("Nyquist fraction must be in (0, 1)"))

    # The COI removes one margin c * omega0 / (2πf) from each edge. Requiring
    # the remaining interval to occupy fraction q of the full record gives
    # f >= c * omega0 / (π * T * (1 - q)).
    minimum_frequency_hz = cone_of_influence_factor * omega0 /
                           (pi * duration_s * (1.0 - minimum_coi_valid_fraction))
    nyquist_hz = 0.5 / dt
    upper_limit_hz = min(Float64(maximum_frequency_hz), nyquist_fraction * nyquist_hz)
    first_frequency_hz = ceil(minimum_frequency_hz / frequency_step_hz) * frequency_step_hz
    last_frequency_hz = floor(upper_limit_hz / frequency_step_hz) * frequency_step_hz
    first_frequency_hz <= last_frequency_hz || throw(ArgumentError(
        "time sampling and COI constraint leave no usable CWT frequencies",
    ))
    collect(Float64(first_frequency_hz):Float64(frequency_step_hz):Float64(last_frequency_hz))
end

function recommended_morlet_config(time_s::AbstractVector; kwargs...)
    frequencies_hz = recommended_frequency_grid_hz(time_s; kwargs...)
    MorletConfig(frequencies_hz=frequencies_hz)
end

function morlet_convolution_kernel(dt, frequency_hz, config::MorletConfig)
    scale_s = morlet_scale_s(frequency_hz, config.omega0)
    half_width = ceil(Int, config.truncate_sigma * scale_s / dt)
    offsets_s = collect((-half_width):half_width) .* dt
    normalization = pi^(-0.25) / sqrt(scale_s)
    wavelet = normalization .* exp.(2.0im * pi * frequency_hz .* offsets_s) .*
              exp.(-0.5 .* (offsets_s ./ scale_s) .^ 2)
    reverse(conj.(ComplexF64.(wavelet))) .* dt, scale_s
end

function fft_convolve_same(signal::AbstractVector, kernel::AbstractVector)
    signal_count = length(signal)
    kernel_count = length(kernel)
    isodd(kernel_count) || throw(ArgumentError("centered convolution kernel must have odd length"))
    full_count = signal_count + kernel_count - 1
    padded_count = nextpow(2, full_count)
    padded_signal = zeros(ComplexF64, padded_count)
    padded_kernel = zeros(ComplexF64, padded_count)
    padded_signal[1:signal_count] .= signal
    padded_kernel[1:kernel_count] .= kernel
    full = ifft(fft(padded_signal) .* fft(padded_kernel))[1:full_count]
    first_index = (kernel_count + 1) ÷ 2
    ComplexF64.(full[first_index:(first_index + signal_count - 1)])
end

function continuous_wavelet_transform(
    time_s::AbstractVector,
    signal::AbstractVector;
    config::MorletConfig=MorletConfig(),
)
    validate_config(config)
    length(time_s) == length(signal) ||
        throw(DimensionMismatch("time and signal must have equal lengths"))
    all(isfinite, signal) || throw(ArgumentError("signal samples must be finite"))
    dt = sample_interval(time_s)
    nyquist_hz = 0.5 / dt
    maximum(config.frequencies_hz) < nyquist_hz ||
        throw(ArgumentError("CWT frequency grid must remain below Nyquist frequency"))

    time = Float64.(time_s)
    values = Float64.(signal)
    coefficients = zeros(ComplexF64, length(config.frequencies_hz), length(time))
    cone_valid = falses(size(coefficients))
    edge_distance_s = min.(time .- first(time), last(time) .- time)

    for (frequency_index, frequency_hz) in enumerate(config.frequencies_hz)
        kernel, scale_s = morlet_convolution_kernel(dt, frequency_hz, config)
        coefficients[frequency_index, :] .= fft_convolve_same(values, kernel)
        cone_valid[frequency_index, :] .=
            edge_distance_s .>= config.cone_of_influence_factor * scale_s
    end

    CWTResult(time, copy(config.frequencies_hz), coefficients, BitMatrix(cone_valid))
end

wavelet_power(result::CWTResult) = abs2.(result.coefficients)

function full_cross_correlation(output::AbstractVector, input::AbstractVector)
    length(output) == length(input) ||
        throw(DimensionMismatch("cross-correlation signals must have equal lengths"))
    sample_count = length(input)
    padded_count = nextpow(2, 2 * sample_count - 1)
    convolution = ifft(
        fft(vcat(Float64.(output), zeros(padded_count - sample_count))) .*
        fft(vcat(reverse(Float64.(input)), zeros(padded_count - sample_count))),
    )
    real.(convolution[1:(2 * sample_count - 1)])
end

function energy_observables(time_s, coefficients, mask)
    power = abs2.(coefficients)
    power[.!mask] .= 0.0
    energy = sum(power)
    if !(energy > 0.0)
        return (energy=0.0, centroid_s=NaN, peak_s=NaN, envelope=sqrt.(power))
    end
    centroid_s = sum(time_s .* power) / energy
    peak_s = time_s[argmax(power)]
    (energy, centroid_s, peak_s, envelope=sqrt.(power))
end

function compare_wavelets(
    reference::CWTResult,
    sample::CWTResult;
    config::MorletConfig=MorletConfig(frequencies_hz=copy(reference.frequency_hz)),
)
    validate_config(config)
    reference.frequency_hz == sample.frequency_hz ||
        throw(ArgumentError("reference and sample must use the same CWT frequencies"))
    reference.frequency_hz == config.frequencies_hz ||
        throw(ArgumentError("comparison config must match the CWT frequencies"))
    reference.time_s == sample.time_s ||
        throw(ArgumentError("reference and sample must use the same time grid"))
    dt = sample_interval(reference.time_s)
    frequency_count = length(reference.frequency_hz)

    reference_centroid_s = fill(NaN, frequency_count)
    sample_centroid_s = fill(NaN, frequency_count)
    centroid_delay_s = fill(NaN, frequency_count)
    reference_peak_s = fill(NaN, frequency_count)
    sample_peak_s = fill(NaN, frequency_count)
    peak_delay_s = fill(NaN, frequency_count)
    correlation_delay_s = fill(NaN, frequency_count)
    amplitude_ratio = fill(NaN, frequency_count)
    envelope_correlation = fill(NaN, frequency_count)
    reference_energy = zeros(Float64, frequency_count)

    for frequency_index in 1:frequency_count
        mask = reference.cone_of_influence_valid[frequency_index, :] .&
               sample.cone_of_influence_valid[frequency_index, :]
        reference_values = energy_observables(
            reference.time_s,
            view(reference.coefficients, frequency_index, :),
            mask,
        )
        sample_values = energy_observables(
            sample.time_s,
            view(sample.coefficients, frequency_index, :),
            mask,
        )
        reference_energy[frequency_index] = reference_values.energy
        reference_centroid_s[frequency_index] = reference_values.centroid_s
        sample_centroid_s[frequency_index] = sample_values.centroid_s
        centroid_delay_s[frequency_index] = sample_values.centroid_s - reference_values.centroid_s
        reference_peak_s[frequency_index] = reference_values.peak_s
        sample_peak_s[frequency_index] = sample_values.peak_s
        peak_delay_s[frequency_index] = sample_values.peak_s - reference_values.peak_s

        if reference_values.energy > 0.0 && sample_values.energy > 0.0
            amplitude_ratio[frequency_index] =
                sqrt(sample_values.energy / reference_values.energy)
            correlation = full_cross_correlation(
                sample_values.envelope,
                reference_values.envelope,
            )
            normalization = sqrt(
                sum(abs2, sample_values.envelope) *
                sum(abs2, reference_values.envelope),
            )
            envelope_correlation[frequency_index] = maximum(correlation) / normalization
            lag_samples = argmax(correlation) - length(reference.time_s)
            correlation_delay_s[frequency_index] = lag_samples * dt
        end
    end

    maximum_reference_energy = maximum(reference_energy)
    reference_energy_relative = maximum_reference_energy > 0.0 ?
                                reference_energy ./ maximum_reference_energy :
                                zeros(Float64, frequency_count)
    valid = BitVector(
        (reference_energy_relative .>= config.reference_energy_floor_relative) .&
        (amplitude_ratio .>= config.amplitude_ratio_floor) .&
        (envelope_correlation .>= config.correlation_floor) .&
        isfinite.(centroid_delay_s) .&
        isfinite.(correlation_delay_s),
    )

    WaveletComparison(
        copy(reference.frequency_hz),
        reference_centroid_s,
        sample_centroid_s,
        centroid_delay_s,
        reference_peak_s,
        sample_peak_s,
        peak_delay_s,
        correlation_delay_s,
        amplitude_ratio,
        envelope_correlation,
        reference_energy_relative,
        valid,
    )
end

function wavelet_value_at_frequency(comparison::WaveletComparison, frequency_hz::Real)
    index = argmin(abs.(comparison.frequency_hz .- frequency_hz))
    (
        index,
        frequency_hz=comparison.frequency_hz[index],
        reference_centroid_s=comparison.reference_centroid_s[index],
        sample_centroid_s=comparison.sample_centroid_s[index],
        centroid_delay_s=comparison.centroid_delay_s[index],
        reference_peak_s=comparison.reference_peak_s[index],
        sample_peak_s=comparison.sample_peak_s[index],
        peak_delay_s=comparison.peak_delay_s[index],
        correlation_delay_s=comparison.correlation_delay_s[index],
        amplitude_ratio=comparison.amplitude_ratio[index],
        envelope_correlation=comparison.envelope_correlation[index],
        reference_energy_relative=comparison.reference_energy_relative[index],
        valid=comparison.valid[index],
    )
end

end
