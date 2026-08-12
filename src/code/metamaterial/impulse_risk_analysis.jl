module ImpulseRiskAnalysis

using FFTW
using LinearAlgebra
using Statistics

export PulseConfig,
       ApertureConfig,
       PulseSpectrum,
       pulse_signal,
       pulse_spectrum,
       spectral_band,
       spectral_energy_fraction,
       aperture_centers_m,
       ideal_delays_s,
       quantize_delays,
       focus_waveform,
       peak_profile,
       contiguous_width,
       pulse_metrics

"""Provisional finite-band input used by the impulse-lens risk experiments."""
Base.@kwdef struct PulseConfig
    center_frequency_hz::Float64 = 242.0e3
    cycles::Float64 = 5.0
    samples_per_period::Int = 80
    fft_length::Int = 65536
end

"""Scalar aperture used only as an optimistic delay-bandwidth feasibility bound."""
Base.@kwdef struct ApertureConfig
    element_count::Int = 15
    element_width_m::Float64 = 4.2e-3
    pitch_m::Float64 = 4.8e-3
    focal_distance_m::Float64 = 35.0e-3
    propagation_speed_m_per_s::Float64 = 2340.0
    rayleigh_alpha_per_s::Float64 = 79560.0
    rayleigh_beta_s::Float64 = 2.5e-9
end

struct PulseSpectrum
    config::PulseConfig
    time_s::Vector{Float64}
    signal::Vector{Float64}
    frequency_hz::Vector{Float64}
    spectrum::Vector{ComplexF64}
    relative_amplitude::Vector{Float64}
end

function validate(config::PulseConfig)
    config.center_frequency_hz > 0 || throw(ArgumentError("center frequency must be positive"))
    config.cycles > 0 || throw(ArgumentError("pulse cycles must be positive"))
    config.samples_per_period >= 8 || throw(ArgumentError("at least eight samples per period are required"))
    config.fft_length > config.cycles * config.samples_per_period ||
        throw(ArgumentError("FFT length must exceed the driven pulse length"))
    ispow2(config.fft_length) || throw(ArgumentError("FFT length must be a power of two"))
    nothing
end

function validate(config::ApertureConfig)
    isodd(config.element_count) || throw(ArgumentError("element count must be odd"))
    config.element_count >= 3 || throw(ArgumentError("at least three elements are required"))
    config.element_width_m > 0 || throw(ArgumentError("element width must be positive"))
    config.pitch_m >= config.element_width_m || throw(ArgumentError("pitch must not be smaller than element width"))
    config.focal_distance_m > 0 || throw(ArgumentError("focal distance must be positive"))
    config.propagation_speed_m_per_s > 0 || throw(ArgumentError("propagation speed must be positive"))
    nothing
end

function pulse_signal(t_s::Real, config::PulseConfig=PulseConfig())
    duration_s = config.cycles / config.center_frequency_hz
    0 <= t_s < duration_s || return 0.0
    0.5 * (1 - cospi(2t_s / duration_s)) *
    sinpi(2config.center_frequency_hz * t_s)
end

function pulse_spectrum(config::PulseConfig=PulseConfig())
    validate(config)
    dt_s = 1 / (config.samples_per_period * config.center_frequency_hz)
    time_s = collect(0:(config.fft_length - 1)) .* dt_s
    signal = pulse_signal.(time_s, Ref(config))
    spectrum = rfft(signal)
    frequency_hz = collect(0:(length(spectrum) - 1)) ./ (config.fft_length * dt_s)
    scale = maximum(abs, spectrum)
    relative_amplitude = abs.(spectrum) ./ scale
    PulseSpectrum(config, time_s, signal, frequency_hz, spectrum, relative_amplitude)
end

"""Contiguous spectral interval around the carrier peak above an amplitude-dB threshold."""
function spectral_band(spectrum::PulseSpectrum, threshold_db::Real)
    threshold_db < 0 || throw(ArgumentError("spectral threshold must be negative"))
    threshold = 10.0^(threshold_db / 20)
    peak_index = argmax(spectrum.relative_amplitude)
    left = peak_index
    right = peak_index
    while left > firstindex(spectrum.relative_amplitude) &&
          spectrum.relative_amplitude[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(spectrum.relative_amplitude) &&
          spectrum.relative_amplitude[right + 1] >= threshold
        right += 1
    end
    (
        lower_hz=spectrum.frequency_hz[left],
        upper_hz=spectrum.frequency_hz[right],
        peak_hz=spectrum.frequency_hz[peak_index],
        lower_index=left,
        upper_index=right,
    )
end

function spectral_energy_fraction(
    spectrum::PulseSpectrum;
    lower_hz::Real=-Inf,
    upper_hz::Real=Inf,
)
    weights = abs2.(spectrum.spectrum)
    selected = (spectrum.frequency_hz .>= lower_hz) .&
               (spectrum.frequency_hz .<= upper_hz)
    sum(weights[selected]) / sum(weights)
end

function aperture_centers_m(config::ApertureConfig=ApertureConfig())
    validate(config)
    half = (config.element_count - 1) ÷ 2
    collect((-half):half) .* config.pitch_m
end

"""Non-negative delays that equalize geometric travel times at the target focus."""
function ideal_delays_s(config::ApertureConfig=ApertureConfig())
    centers = aperture_centers_m(config)
    distances = hypot.(config.focal_distance_m, centers)
    (maximum(distances) .- distances) ./ config.propagation_speed_m_per_s
end

function quantize_delays(delays_s::AbstractVector{<:Real}, state_count::Integer)
    state_count >= 2 || throw(ArgumentError("at least two delay states are required"))
    maximum_delay = maximum(delays_s)
    iszero(maximum_delay) && return zeros(Float64, length(delays_s))
    levels = collect(range(0.0, maximum_delay; length=state_count))
    [levels[argmin(abs.(levels .- delay))] for delay in delays_s]
end

function focus_transfer(
    frequencies_hz::AbstractVector{<:Real},
    config::ApertureConfig,
    delays_s::AbstractVector{<:Real};
    observation_x_m::Real=config.focal_distance_m,
    observation_y_m::Real=0.0,
    loss_scale::Real=0.0,
)
    validate(config)
    centers = aperture_centers_m(config)
    length(delays_s) == length(centers) || throw(DimensionMismatch("one delay per aperture element is required"))
    distances = hypot.(observation_x_m, observation_y_m .- centers)
    result = zeros(ComplexF64, length(frequencies_hz))
    for (frequency_index, frequency_hz) in enumerate(frequencies_hz)
        frequency_hz <= 0 && continue
        omega = 2pi * frequency_hz
        wavelength_m = config.propagation_speed_m_per_s / frequency_hz
        decay_rate_per_s = loss_scale * (
            config.rayleigh_alpha_per_s + config.rayleigh_beta_s * omega^2
        ) / 2
        result[frequency_index] = sum(eachindex(centers)) do index
            propagation = config.element_width_m /
                          sqrt(wavelength_m * distances[index]) *
                          cis(-omega * distances[index] / config.propagation_speed_m_per_s)
            cell = exp(-decay_rate_per_s * delays_s[index]) *
                   cis(-omega * delays_s[index])
            propagation * cell
        end
    end
    result
end

function waveform_from_transfer(spectrum::PulseSpectrum, transfer)
    length(transfer) == length(spectrum.spectrum) ||
        throw(DimensionMismatch("transfer and pulse spectra must have equal lengths"))
    irfft(spectrum.spectrum .* transfer, length(spectrum.signal))
end

"""
Return the focused waveform for a specified delay mask and the matching uniform
aperture waveform. Loss is applied only to the additional cell delay, making
this an optimistic upper bound for a real structure.
"""
function focus_waveform(
    spectrum::PulseSpectrum,
    config::ApertureConfig,
    delays_s::AbstractVector{<:Real}=ideal_delays_s(config);
    loss_scale::Real=0.0,
    observation_x_m::Real=config.focal_distance_m,
    observation_y_m::Real=0.0,
)
    focus = focus_transfer(
        spectrum.frequency_hz,
        config,
        delays_s;
        observation_x_m,
        observation_y_m,
        loss_scale,
    )
    uniform = focus_transfer(
        spectrum.frequency_hz,
        config,
        zeros(length(delays_s));
        observation_x_m,
        observation_y_m,
        loss_scale=0.0,
    )
    (
        focused=waveform_from_transfer(spectrum, focus),
        uniform=waveform_from_transfer(spectrum, uniform),
        focused_transfer=focus,
        uniform_transfer=uniform,
    )
end

function peak_profile(
    spectrum::PulseSpectrum,
    config::ApertureConfig,
    delays_s::AbstractVector{<:Real},
    coordinate_m::AbstractVector{<:Real};
    direction::Symbol,
    loss_scale::Real=0.0,
)
    direction in (:axial, :transverse) ||
        throw(ArgumentError("profile direction must be :axial or :transverse"))
    [begin
        observation_x_m = direction == :axial ? coordinate : config.focal_distance_m
        observation_y_m = direction == :transverse ? coordinate : 0.0
        response = focus_waveform(
            spectrum,
            config,
            delays_s;
            loss_scale,
            observation_x_m,
            observation_y_m,
        )
        maximum(abs, response.focused)
    end for coordinate in coordinate_m]
end

function contiguous_width(coordinate, amplitude; relative_level=inv(sqrt(2.0)))
    length(coordinate) == length(amplitude) ||
        throw(DimensionMismatch("coordinate and amplitude must have equal lengths"))
    peak_index = argmax(amplitude)
    threshold = relative_level * amplitude[peak_index]
    left = peak_index
    right = peak_index
    while left > firstindex(amplitude) && amplitude[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(amplitude) && amplitude[right + 1] >= threshold
        right += 1
    end
    (
        width=coordinate[right] - coordinate[left],
        start=coordinate[left],
        stop=coordinate[right],
        peak_coordinate=coordinate[peak_index],
        peak_amplitude=amplitude[peak_index],
    )
end

function analytic_signal(signal::AbstractVector{<:Real})
    count = length(signal)
    spectrum = fft(signal)
    multiplier = zeros(Float64, count)
    multiplier[1] = 1.0
    if iseven(count)
        multiplier[2:(count ÷ 2)] .= 2.0
        multiplier[count ÷ 2 + 1] = 1.0
    else
        multiplier[2:((count + 1) ÷ 2)] .= 2.0
    end
    ifft(spectrum .* multiplier)
end

function minimum_energy_interval(signal, fraction::Real=0.9)
    0 < fraction <= 1 || throw(ArgumentError("energy fraction must lie in (0, 1]"))
    energy = abs2.(analytic_signal(signal))
    target = fraction * sum(energy)
    best_left, best_right = firstindex(energy), lastindex(energy)
    right = firstindex(energy) - 1
    accumulated = 0.0
    for left in eachindex(energy)
        while right < lastindex(energy) && accumulated < target
            right += 1
            accumulated += energy[right]
        end
        if accumulated >= target && right - left < best_right - best_left
            best_left, best_right = left, right
        end
        accumulated -= energy[left]
    end
    (left=best_left, right=best_right, energy=energy)
end

function maximum_normalized_correlation(first_signal, second_signal)
    length(first_signal) == length(second_signal) ||
        throw(DimensionMismatch("correlated signals must have equal lengths"))
    first = analytic_signal(first_signal)
    second = analytic_signal(second_signal)
    padded_count = nextpow(2, 2length(first) - 1)
    padded_first = zeros(ComplexF64, padded_count)
    padded_second = zeros(ComplexF64, padded_count)
    padded_first[eachindex(first)] .= first
    padded_second[eachindex(second)] .= second
    correlation = ifft(conj.(fft(padded_first)) .* fft(padded_second))
    denominator = norm(first) * norm(second)
    denominator > 0 ? maximum(abs, correlation) / denominator : NaN
end

function pulse_metrics(
    focused::AbstractVector{<:Real},
    uniform::AbstractVector{<:Real},
    ideal::AbstractVector{<:Real},
    dt_s::Real,
)
    focus_interval = minimum_energy_interval(focused)
    uniform_interval = minimum_energy_interval(uniform)
    focus_width_s = (focus_interval.right - focus_interval.left + 1) * dt_s
    uniform_width_s = (uniform_interval.right - uniform_interval.left + 1) * dt_s
    main_energy = sum(focus_interval.energy[focus_interval.left:focus_interval.right])
    postcursor_energy = focus_interval.right < lastindex(focus_interval.energy) ?
        sum(focus_interval.energy[(focus_interval.right + 1):end]) : 0.0
    (
        gain_peak=maximum(abs, focused) / maximum(abs, uniform),
        focused_peak=maximum(abs, focused),
        uniform_peak=maximum(abs, uniform),
        duration90_s=focus_width_s,
        broadening_ratio=focus_width_s / uniform_width_s,
        pulse_correlation=maximum_normalized_correlation(focused, ideal),
        postcursor_ratio=postcursor_energy / main_energy,
    )
end

end
