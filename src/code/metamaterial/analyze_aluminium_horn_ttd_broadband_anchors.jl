module AnalyseAluminiumHornTTDBroadbandAnchors

ENV["GKSwstype"] = "100"

using FFTW
using JLD2
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const LENS_ROOT = joinpath(
    PROJECT_ROOT, "tmp", "aluminium_horn_ttd_full_diffuser_aperture_242khz",
)
const REFERENCE_ROOT = joinpath(
    PROJECT_ROOT, "tmp", "aluminium_horn_ttd_full_aperture_242khz",
)
const FREQUENCIES_HZ = Float64[193.8e3, 242.0e3, 290.1e3]

include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
include(joinpath(@__DIR__, "spectral_analysis.jl"))
using .ImpulseRiskAnalysis
using .SpectralAnalysis

frequency_suffix(frequency_hz) = isapprox(frequency_hz, 242.0e3; atol=1e-6) ? "" :
    "_f$(round(Int, frequency_hz))hz"
lens_path(frequency_hz) = joinpath(
    LENS_ROOT, "lens_half_weighted_harmonic$(frequency_suffix(frequency_hz)).jld2",
)
reference_path(frequency_hz) = joinpath(
    REFERENCE_ROOT, "uniform_half_harmonic$(frequency_suffix(frequency_hz)).jld2",
)

function unwrap_phase(phases)
    result = Float64[first(phases)]
    for phase in Iterators.drop(phases, 1)
        push!(result, last(result) + mod(phase - last(result) + pi, 2pi) - pi)
    end
    result
end

function linear_interpolate(x, y, query)
    right = searchsortedfirst(x, query)
    right == firstindex(x) && return y[right]
    right > lastindex(x) && return y[end]
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * y[left] + fraction * y[right]
end

function analytic_envelope(signal)
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
    abs.(ifft(spectrum .* multiplier))
end

function main()
    lens_data = JLD2.load.(lens_path.(FREQUENCIES_HZ))
    reference_data = JLD2.load.(reference_path.(FREQUENCIES_HZ))
    lens_focus = ComplexF64[item["focus_displacement_m"][1] for item in lens_data]
    reference_focus = ComplexF64[item["focus_displacement_m"][1] for item in reference_data]
    lens_power = abs.(Float64[item["active_input_power_w"] for item in lens_data])
    reference_power = abs.(Float64[item["active_input_power_w"] for item in reference_data])
    relative = lens_focus ./ reference_focus
    relative_phase = unwrap_phase(angle.(relative))
    phase_steps_deg = rad2deg.(diff(relative_phase))

    pulse = pulse_spectrum(PulseConfig(
        center_frequency_hz=242.0e3,
        cycles=5.0,
        samples_per_period=80,
        fft_length=65536,
    ))
    active = (pulse.frequency_hz .>= first(FREQUENCIES_HZ)) .&
             (pulse.frequency_hz .<= last(FREQUENCIES_HZ))
    lens_transfer = zeros(ComplexF64, length(pulse.frequency_hz))
    reference_transfer = zeros(ComplexF64, length(pulse.frequency_hz))
    interpolated_lens_power = zeros(Float64, length(pulse.frequency_hz))
    interpolated_reference_power = zeros(Float64, length(pulse.frequency_hz))
    log_amplitude = log.(abs.(relative))
    for index in findall(active)
        frequency_hz = pulse.frequency_hz[index]
        lens_transfer[index] = exp(linear_interpolate(
            FREQUENCIES_HZ, log_amplitude, frequency_hz,
        )) * cis(linear_interpolate(FREQUENCIES_HZ, relative_phase, frequency_hz))
        reference_transfer[index] = 1.0 + 0im
        interpolated_lens_power[index] = exp(linear_interpolate(
            FREQUENCIES_HZ, log.(lens_power), frequency_hz,
        ))
        interpolated_reference_power[index] = exp(linear_interpolate(
            FREQUENCIES_HZ, log.(reference_power), frequency_hz,
        ))
    end
    spectral_weight = abs2.(pulse.spectrum)
    lens_energy_proxy = sum(spectral_weight[active] .* interpolated_lens_power[active])
    reference_energy_proxy = sum(
        spectral_weight[active] .* interpolated_reference_power[active],
    )
    equal_energy_scale = sqrt(reference_energy_proxy / lens_energy_proxy)
    lens_transfer .*= equal_energy_scale
    lens_waveform = irfft(pulse.spectrum .* lens_transfer, length(pulse.signal))
    reference_waveform = irfft(
        pulse.spectrum .* reference_transfer, length(pulse.signal),
    )
    metrics = pulse_metrics(
        lens_waveform, reference_waveform, reference_waveform,
        pulse.time_s[2] - pulse.time_s[1],
    )
    energy_fraction = spectral_energy_fraction(
        pulse; lower_hz=first(FREQUENCIES_HZ), upper_hz=last(FREQUENCIES_HZ),
    )
    carrier_index = 2
    carrier_equal_power_gain = abs(relative[carrier_index]) *
        sqrt(reference_power[carrier_index] / lens_power[carrier_index])
    interpolation_reliable = maximum(abs, phase_steps_deg) <= 175.0
    passed_to_dense = interpolation_reliable && metrics.gain_peak >= 1.90 &&
                      metrics.broadening_ratio <= 1.25 &&
                      metrics.pulse_correlation >= 0.90 &&
                      metrics.postcursor_ratio <= 0.12

    mkpath(LENS_ROOT)
    frequency_path = joinpath(LENS_ROOT, "15_weighted_broadband_anchors.csv")
    open(frequency_path, "w") do io
        println(io, "frequency_hz,lens_focus_ux_abs_nm,reference_focus_ux_abs_nm,pressure_gain,lens_active_power_w,reference_active_power_w,equal_power_gain,relative_phase_deg")
        for index in eachindex(FREQUENCIES_HZ)
            equal_power_gain = abs(relative[index]) *
                sqrt(reference_power[index] / lens_power[index])
            println(io, join((
                FREQUENCIES_HZ[index], abs(lens_focus[index]) * 1e9,
                abs(reference_focus[index]) * 1e9, abs(relative[index]),
                lens_power[index], reference_power[index], equal_power_gain,
                rad2deg(relative_phase[index]),
            ), ','))
        end
    end
    summary_path = joinpath(LENS_ROOT, "15_weighted_broadband_anchor_summary.csv")
    open(summary_path, "w") do io
        println(io, "lower_frequency_hz,carrier_frequency_hz,upper_frequency_hz,band_spectral_energy_fraction,equal_energy_scale,provisional_impulse_peak_gain,carrier_equal_power_gain,broadening_ratio,pulse_correlation,postcursor_ratio,lower_phase_step_deg,upper_phase_step_deg,interpolation_reliable,passed_to_dense")
        println(io, join((
            FREQUENCIES_HZ[1], FREQUENCIES_HZ[2], FREQUENCIES_HZ[3],
            energy_fraction, equal_energy_scale, metrics.gain_peak,
            carrier_equal_power_gain, metrics.broadening_ratio,
            metrics.pulse_correlation, metrics.postcursor_ratio,
            phase_steps_deg[1], phase_steps_deg[2], interpolation_reliable,
            passed_to_dense,
        ), ','))
    end
    response_panel = plot(
        FREQUENCIES_HZ ./ 1e3,
        abs.(relative) .* sqrt.(reference_power ./ lens_power);
        marker=:circle, linewidth=2.5, xlabel="frequency, kHz",
        ylabel="equal-power focus gain", title="Weighted Al lens: broadband anchors",
        label="lens / abrupt straight", gridalpha=0.25,
    )
    hline!(response_panel, [2.0]; linestyle=:dash, color=:black, label="target 2")
    time_us = pulse.time_s .* 1e6
    peak_index = argmax(analytic_envelope(lens_waveform))
    left = max(firstindex(time_us), peak_index - 700)
    right = min(lastindex(time_us), peak_index + 900)
    waveform_panel = plot(
        time_us[left:right], analytic_envelope(lens_waveform)[left:right];
        linewidth=2.5, xlabel="relative time, μs", ylabel="normalized envelope",
        title="Three-anchor provisional 5-cycle pulse", label="weighted lens",
        gridalpha=0.25,
    )
    plot!(waveform_panel, time_us[left:right],
          analytic_envelope(reference_waveform)[left:right];
          linewidth=2, linestyle=:dash, label="abrupt straight")
    figure_path = joinpath(LENS_ROOT, "15_weighted_broadband_anchors.png")
    savefig(plot(response_panel, waveform_panel; layout=(2, 1), size=(1050, 1000),
                 margin=5Plots.mm), figure_path)
    data_path = joinpath(LENS_ROOT, "15_weighted_broadband_anchors.jld2")
    JLD2.jldsave(
        data_path;
        format_version=1,
        frequencies_hz=FREQUENCIES_HZ,
        lens_focus_m=lens_focus,
        reference_focus_m=reference_focus,
        lens_power_w=lens_power,
        reference_power_w=reference_power,
        relative_phase_rad=relative_phase,
        pulse_time_s=pulse.time_s,
        lens_waveform,
        reference_waveform,
        equal_energy_scale,
        metrics,
        energy_fraction,
        interpolation_reliable,
        passed_to_dense,
    )
    println("[+] anchor equal-power gains=" * string(
        abs.(relative) .* sqrt.(reference_power ./ lens_power),
    ))
    println("[+] relative phase steps=$phase_steps_deg deg")
    println("[+] provisional G_peak=$(metrics.gain_peak), Bt=$(metrics.broadening_ratio), " *
            "rho=$(metrics.pulse_correlation), post=$(metrics.postcursor_ratio)")
    println("[+] band energy fraction=$energy_fraction, passed_to_dense=$passed_to_dense")
    println("[+] $summary_path")
    println("[+] $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
