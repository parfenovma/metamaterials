module SpectralAnalysis

using FFTW
using Statistics

export SpectrumConfig,
       SpectrumResult,
       TransferResult,
       one_sided_spectrum,
       analyze_transfer,
       analytic_envelope,
       envelope_delay,
       relative_transfer,
       value_at_frequency

Base.@kwdef struct SpectrumConfig
    window::Symbol = :rectangular
    zero_padding_factor::Int = 4
    input_floor_relative::Float64 = 1.0e-3
    regularization_relative::Float64 = 1.0e-12
end

struct SpectrumResult
    frequency_hz::Vector{Float64}
    values::Vector{ComplexF64}
end

struct TransferResult
    frequency_hz::Vector{Float64}
    transfer::Vector{ComplexF64}
    amplitude::Vector{Float64}
    magnitude_squared::Vector{Float64}
    phase_rad::Vector{Float64}
    group_delay_s::Vector{Float64}
    valid::BitVector
end

function validate_config(config::SpectrumConfig)
    config.window in (:hann, :rectangular) ||
        throw(ArgumentError("window must be :hann or :rectangular"))
    config.zero_padding_factor >= 1 ||
        throw(ArgumentError("zero_padding_factor must be at least one"))
    0.0 <= config.input_floor_relative < 1.0 ||
        throw(ArgumentError("input_floor_relative must be in [0, 1)"))
    config.regularization_relative >= 0.0 ||
        throw(ArgumentError("regularization_relative must be non-negative"))
end

function sample_interval(time_s::AbstractVector)
    length(time_s) >= 2 || throw(ArgumentError("at least two time samples are required"))
    all(isfinite, time_s) || throw(ArgumentError("time samples must be finite"))
    differences = diff(time_s)
    all(>(0.0), differences) || throw(ArgumentError("time samples must be strictly increasing"))
    dt = sum(differences) / length(differences)
    maximum(abs.(differences .- dt)) <= max(1.0e-10 * dt, eps(dt) * 32) ||
        throw(ArgumentError("FFT analysis requires uniformly spaced time samples"))
    dt
end

function window_values(kind::Symbol, sample_count::Integer)
    kind == :rectangular && return ones(Float64, sample_count)
    sample_count == 1 && return ones(Float64, 1)
    [0.5 * (1.0 - cos(2.0 * pi * n / (sample_count - 1))) for n in 0:(sample_count - 1)]
end

function one_sided_spectrum(
    time_s::AbstractVector,
    signal::AbstractVector;
    config::SpectrumConfig=SpectrumConfig(),
)
    validate_config(config)
    length(time_s) == length(signal) ||
        throw(DimensionMismatch("time and signal must have equal lengths"))
    all(isfinite, signal) || throw(ArgumentError("signal samples must be finite"))

    dt = sample_interval(time_s)
    sample_count = length(signal)
    padded_count = nextpow(2, sample_count * config.zero_padding_factor)
    window = window_values(config.window, sample_count)
    normalization = sum(window)

    padded = zeros(Float64, padded_count)
    padded[1:sample_count] .= Float64.(signal) .* window
    transformed = fft(padded) ./ normalization
    last_index = padded_count ÷ 2 + 1
    frequency_hz = collect(0:(last_index - 1)) ./ (padded_count * dt)
    SpectrumResult(frequency_hz, ComplexF64.(transformed[1:last_index]))
end

function analytic_envelope(signal::AbstractVector)
    all(isfinite, signal) || throw(ArgumentError("signal samples must be finite"))
    sample_count = length(signal)
    sample_count > 0 || return Float64[]
    multiplier = zeros(Float64, sample_count)
    multiplier[1] = 1.0
    if iseven(sample_count)
        multiplier[sample_count ÷ 2 + 1] = 1.0
        multiplier[2:(sample_count ÷ 2)] .= 2.0
    else
        multiplier[2:((sample_count + 1) ÷ 2)] .= 2.0
    end
    abs.(ifft(fft(Float64.(signal)) .* multiplier))
end

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

"""Legacy envelope cross-correlation delay retained for comparison with phase delay."""
function envelope_delay(time_s, input_signal, output_signal)
    length(time_s) == length(input_signal) == length(output_signal) ||
        throw(DimensionMismatch("time and signals must have equal lengths"))
    dt = sample_interval(time_s)
    input_envelope = analytic_envelope(input_signal)
    output_envelope = analytic_envelope(output_signal)
    correlation = full_cross_correlation(output_envelope, input_envelope)
    lag_samples = argmax(correlation) - length(input_signal)
    lag_samples * dt
end

function unwrap_phase(values::AbstractVector{<:Complex}, valid::AbstractVector{Bool})
    length(values) == length(valid) || throw(DimensionMismatch("values and mask must match"))
    result = fill(NaN, length(values))
    previous_raw = 0.0
    previous_unwrapped = 0.0
    active_segment = false

    for i in eachindex(values)
        if !valid[i]
            active_segment = false
            continue
        end

        raw = angle(values[i])
        if !active_segment
            result[i] = raw
            previous_raw = raw
            previous_unwrapped = raw
            active_segment = true
            continue
        end

        delta = mod(raw - previous_raw + pi, 2.0 * pi) - pi
        previous_unwrapped += delta
        result[i] = previous_unwrapped
        previous_raw = raw
    end
    result
end

function group_delay(frequency_hz, phase_rad, valid)
    result = fill(NaN, length(frequency_hz))
    for i in 2:(length(frequency_hz) - 1)
        if valid[i - 1] && valid[i] && valid[i + 1]
            result[i] = -(phase_rad[i + 1] - phase_rad[i - 1]) /
                        (2.0 * pi * (frequency_hz[i + 1] - frequency_hz[i - 1]))
        end
    end
    result
end

function transfer_from_spectra(
    input::SpectrumResult,
    output::SpectrumResult;
    config::SpectrumConfig=SpectrumConfig(),
)
    input.frequency_hz == output.frequency_hz ||
        throw(ArgumentError("input and output spectra must use the same frequency grid"))

    input_power = abs2.(input.values)
    maximum_power = maximum(input_power)
    maximum_power > 0.0 || throw(ArgumentError("input spectrum is identically zero"))
    valid = BitVector(input_power .>= (config.input_floor_relative^2 * maximum_power))
    regularization = config.regularization_relative * maximum_power
    transfer = output.values .* conj.(input.values) ./ (input_power .+ regularization)
    transfer[.!valid] .= ComplexF64(NaN, NaN)

    phase = unwrap_phase(transfer, valid)
    delay = group_delay(input.frequency_hz, phase, valid)
    TransferResult(
        input.frequency_hz,
        transfer,
        abs.(transfer),
        abs2.(transfer),
        phase,
        delay,
        valid,
    )
end

function analyze_transfer(
    time_s::AbstractVector,
    input_signal::AbstractVector,
    output_signal::AbstractVector;
    config::SpectrumConfig=SpectrumConfig(),
)
    input = one_sided_spectrum(time_s, input_signal; config)
    output = one_sided_spectrum(time_s, output_signal; config)
    transfer_from_spectra(input, output; config)
end

"""
Return the response of `sample` relative to an equal-length reference block.

This removes the source spectrum, common propagation phase and common port
response. It is a calibrated transfer response, not a complete S-matrix.
"""
function relative_transfer(sample::TransferResult, reference::TransferResult)
    sample.frequency_hz == reference.frequency_hz ||
        throw(ArgumentError("sample and reference must use the same frequency grid"))
    valid = sample.valid .& reference.valid .& (reference.amplitude .> 0.0)
    transfer = fill(ComplexF64(NaN, NaN), length(sample.transfer))
    transfer[valid] .= sample.transfer[valid] ./ reference.transfer[valid]
    phase = unwrap_phase(transfer, valid)
    delay = group_delay(sample.frequency_hz, phase, valid)
    TransferResult(
        sample.frequency_hz,
        transfer,
        abs.(transfer),
        abs2.(transfer),
        phase,
        delay,
        BitVector(valid),
    )
end

function value_at_frequency(result::TransferResult, frequency_hz::Real)
    index = argmin(abs.(result.frequency_hz .- frequency_hz))
    (
        index=index,
        frequency_hz=result.frequency_hz[index],
        transfer=result.transfer[index],
        amplitude=result.amplitude[index],
        magnitude_squared=result.magnitude_squared[index],
        phase_rad=result.phase_rad[index],
        group_delay_s=result.group_delay_s[index],
        valid=result.valid[index],
    )
end

end
