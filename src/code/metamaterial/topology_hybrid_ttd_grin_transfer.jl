module TopologyHybridTTDGRINTransfer

using FFTW
using LinearAlgebra: I

if !isdefined(parentmodule(@__MODULE__), :TopologyHybridTTDGRINDesign)
    include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
    using .TopologyHybridTTDGRINDesign
else
    using ..TopologyHybridTTDGRINDesign
end
include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .ImpulseRiskAnalysis

export channel_segments,
       segment_matrix,
       channel_scattering,
       transfer_sweep,
       group_delay_s,
       reconstruct_pulse,
       weighted_band_metrics,
       run_transfer_gate

"""Force/particle-velocity transfer matrix of a nondispersive rod segment."""
function segment_matrix(frequency_hz, length_mm, density, speed, width_mm)
    frequency_hz >= 0 || throw(ArgumentError("frequency must be non-negative"))
    length_mm >= 0 || throw(ArgumentError("segment length must be non-negative"))
    width_mm > 0 || throw(ArgumentError("segment width must be positive"))
    phase = 2pi * Float64(frequency_hz) * Float64(length_mm) * 1e-3 / speed
    impedance = density * speed * width_mm * 1e-3 # unit out-of-plane depth
    ComplexF64[
        cos(phase) im * impedance * sin(phase)
        im * sin(phase) / impedance cos(phase)
    ]
end

function channel_segments(
    delay_s::Real,
    config::HybridTTDGRINConfig=HybridTTDGRINConfig();
    transition_steps::Integer=5,
)
    validate(config)
    delay = Float64(delay_s)
    0 <= delay <= maximum_delay_s(config) * (1 + 1e-10) ||
        throw(ArgumentError("channel delay lies outside the aperture budget"))
    extra_path_mm = config.path_delay_fraction * delay *
                    config.calibrated_path_delay_speed_m_s * 1e3
    total_path_mm = config.device_depth_mm + extra_path_mm
    insert = slow_insert_profile(delay, config; transition_steps)
    remaining_aluminium_mm = total_path_mm - insert.total_length_mm
    remaining_aluminium_mm >= 0 || error("slow insert is longer than the guide")
    al = guided_material(0.0, config)
    segments = NamedTuple[(
        copper_fraction=0.0,
        length_mm=remaining_aluminium_mm / 2,
        width_mm=config.guide_width_mm,
        density_kg_m3=al.density_kg_m3,
        speed_m_s=al.guided_speed_m_s,
    )]
    for segment in insert.segments
        material = guided_material(segment.copper_fraction, config)
        push!(segments, (
            copper_fraction=segment.copper_fraction,
            length_mm=segment.length_mm,
            width_mm=segment.width_mm,
            density_kg_m3=material.density_kg_m3,
            speed_m_s=material.guided_speed_m_s,
        ))
    end
    push!(segments, (
        copper_fraction=0.0,
        length_mm=remaining_aluminium_mm / 2,
        width_mm=config.guide_width_mm,
        density_kg_m3=al.density_kg_m3,
        speed_m_s=al.guided_speed_m_s,
    ))
    segments
end

function channel_scattering(
    frequency_hz::Real,
    delay_s::Real,
    config::HybridTTDGRINConfig=HybridTTDGRINConfig();
    transition_steps::Integer=5,
    width_error_mm::Real=0.0,
)
    segments = channel_segments(delay_s, config; transition_steps)
    matrix = Matrix{ComplexF64}(I, 2, 2)
    for segment in segments
        width = segment.width_mm + Float64(width_error_mm)
        width > 0 || throw(ArgumentError("width error closes a guide segment"))
        matrix *= segment_matrix(
            frequency_hz,
            segment.length_mm,
            segment.density_kg_m3,
            segment.speed_m_s,
            width,
        )
    end
    al = guided_material(0.0, config)
    reference_impedance = al.density_kg_m3 * al.guided_speed_m_s *
                          config.guide_width_mm * 1e-3
    a, b, c, d = matrix[1, 1], matrix[1, 2], matrix[2, 1], matrix[2, 2]
    denominator = a + b / reference_impedance + c * reference_impedance + d
    s21 = 2 / denominator
    s11 = (a + b / reference_impedance - c * reference_impedance - d) / denominator
    (
        frequency_hz=Float64(frequency_hz),
        delay_s=Float64(delay_s),
        S21=ComplexF64(s21),
        S11=ComplexF64(s11),
        transmission=abs2(s21),
        reflection=abs2(s11),
        power_sum=abs2(s21) + abs2(s11),
    )
end

function transfer_sweep(
    delay_s::Real=maximum_delay_s(),
    config::HybridTTDGRINConfig=HybridTTDGRINConfig();
    frequencies_hz=collect(162.2e3:2.0e3:322.2e3),
    transition_steps::Integer=5,
    width_error_mm::Real=0.0,
)
    [
        channel_scattering(
            frequency,
            delay_s,
            config;
            transition_steps,
            width_error_mm,
        )
        for frequency in frequencies_hz
    ]
end

function unwrap_phase(values)
    result = Float64[first(values)]
    for value in Iterators.drop(values, 1)
        push!(result, last(result) + mod(value - last(result) + pi, 2pi) - pi)
    end
    result
end

function group_delay_s(rows)
    frequencies = getproperty.(rows, :frequency_hz)
    phase = unwrap_phase(angle.(getproperty.(rows, :S21)))
    delay = similar(phase)
    delay[1] = -(phase[2] - phase[1]) / (2pi * (frequencies[2] - frequencies[1]))
    for index in 2:(length(rows) - 1)
        delay[index] = -(phase[index + 1] - phase[index - 1]) /
                       (2pi * (frequencies[index + 1] - frequencies[index - 1]))
    end
    delay[end] = -(phase[end] - phase[end - 1]) /
                 (2pi * (frequencies[end] - frequencies[end - 1]))
    delay
end

function reconstruct_pulse(
    delay_s::Real=maximum_delay_s(),
    config::HybridTTDGRINConfig=HybridTTDGRINConfig();
    transition_steps::Integer=5,
    width_error_mm::Real=0.0,
    fft_length::Integer=16384,
)
    pulse = pulse_spectrum(PulseConfig(
        center_frequency_hz=config.frequency_hz,
        cycles=5.0,
        samples_per_period=80,
        fft_length=Int(fft_length),
    ))
    hybrid_transfer = ComplexF64[
        channel_scattering(
            frequency,
            delay_s,
            config;
            transition_steps,
            width_error_mm,
        ).S21
        for frequency in pulse.frequency_hz
    ]
    common_delay_s = config.device_depth_mm * 1e-3 /
                     config.calibrated_aluminium_guide_speed_m_s + Float64(delay_s)
    ideal_transfer = cis.(-2pi .* pulse.frequency_hz .* common_delay_s)
    hybrid = irfft(pulse.spectrum .* hybrid_transfer, length(pulse.signal))
    ideal = irfft(pulse.spectrum .* ideal_transfer, length(pulse.signal))
    dt = pulse.time_s[2] - pulse.time_s[1]
    metrics = pulse_metrics(hybrid, ideal, ideal, dt)
    (; pulse, hybrid_transfer, ideal_transfer, hybrid, ideal, metrics)
end

function linear_interpolate(x, values, query)
    query <= first(x) && return first(values)
    query >= last(x) && return last(values)
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function weighted_band_metrics(rows, config::HybridTTDGRINConfig=HybridTTDGRINConfig())
    pulse = pulse_spectrum(PulseConfig(
        center_frequency_hz=config.frequency_hz,
        cycles=5.0,
        samples_per_period=80,
        fft_length=16384,
    ))
    frequencies = getproperty.(rows, :frequency_hz)
    spectral_power = abs2.(ComplexF64[
        linear_interpolate(pulse.frequency_hz, pulse.spectrum, frequency)
        for frequency in frequencies
    ])
    function weighted(values)
        sum(values .* spectral_power) / sum(spectral_power)
    end
    delays = group_delay_s(rows)
    (
        weighted_transmission=weighted(getproperty.(rows, :transmission)),
        weighted_reflection=weighted(getproperty.(rows, :reflection)),
        minimum_transmission=minimum(getproperty.(rows, :transmission)),
        maximum_reflection=maximum(getproperty.(rows, :reflection)),
        mean_group_delay_us=weighted(delays) * 1e6,
        group_delay_ripple_us=(maximum(delays) - minimum(delays)) * 1e6,
        maximum_power_error=maximum(abs.(getproperty.(rows, :power_sum) .- 1)),
    )
end

function run_transfer_gate(config::HybridTTDGRINConfig=HybridTTDGRINConfig())
    rows = transfer_sweep(maximum_delay_s(config), config)
    band = weighted_band_metrics(rows, config)
    pulse = reconstruct_pulse(maximum_delay_s(config), config)
    metrics = pulse.metrics
    passed = band.weighted_transmission >= 0.80 &&
             metrics.broadening_ratio <= 1.25 &&
             metrics.pulse_correlation >= 0.90 &&
             metrics.postcursor_ratio <= 0.10
    (; rows, band, pulse, passed)
end

end
