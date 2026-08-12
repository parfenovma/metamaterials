module ImpulseRiskPilot

ENV["GKSwstype"] = "100"

using Plots
using Printf
using Statistics

include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .ImpulseRiskAnalysis

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_IMPULSE_RISK_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "impulse_risk_pilot"),
)

const APERTURES = (
    (name="current_h4p2", config=ApertureConfig()),
    (name="single_mode_h3p2", config=ApertureConfig(
        element_count=19,
        element_width_m=3.2e-3,
        pitch_m=3.8e-3,
    )),
)

function write_rows(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function pulse_summary(spectrum)
    band6 = spectral_band(spectrum, -6.0)
    band20 = spectral_band(spectrum, -20.0)
    config = spectrum.config
    row = (
        center_frequency_hz=config.center_frequency_hz,
        cycles=config.cycles,
        duration_s=config.cycles / config.center_frequency_hz,
        dt_s=spectrum.time_s[2] - spectrum.time_s[1],
        fft_length=config.fft_length,
        frequency_resolution_hz=spectrum.frequency_hz[2] - spectrum.frequency_hz[1],
        peak_frequency_hz=band6.peak_hz,
        band_minus6_lower_hz=band6.lower_hz,
        band_minus6_upper_hz=band6.upper_hz,
        band_minus20_lower_hz=band20.lower_hz,
        band_minus20_upper_hz=band20.upper_hz,
        energy_above_253khz=spectral_energy_fraction(spectrum; lower_hz=253.0e3),
    )
    write_rows(joinpath(OUTPUT_ROOT, "pulse_config.csv"), [row])
    row
end

function sensitivity_cases(config)
    exact = ideal_delays_s(config)
    cases = [(label="exact_lossless", delays=exact, loss_scale=0.0, states=0)]
    for scale in (0.25, 0.5, 1.0)
        push!(cases, (
            label="exact_loss_$(replace(string(scale), '.' => 'p'))",
            delays=exact,
            loss_scale=scale,
            states=0,
        ))
    end
    for states in (2, 3, 4, 6, 8)
        quantized = quantize_delays(exact, states)
        push!(cases, (
            label="quantized_$(states)_lossless",
            delays=quantized,
            loss_scale=0.0,
            states,
        ))
        push!(cases, (
            label="quantized_$(states)_loss_1p0",
            delays=quantized,
            loss_scale=1.0,
            states,
        ))
    end
    cases
end

function aperture_rows(spectrum)
    dt_s = spectrum.time_s[2] - spectrum.time_s[1]
    rows = NamedTuple[]
    waveforms = Dict{Tuple{String, String}, Any}()
    for aperture in APERTURES
        config = aperture.config
        ideal_delays = ideal_delays_s(config)
        ideal_response = focus_waveform(spectrum, config, ideal_delays; loss_scale=0.0)
        for case in sensitivity_cases(config)
            response = focus_waveform(
                spectrum,
                config,
                case.delays;
                loss_scale=case.loss_scale,
            )
            metrics = pulse_metrics(
                response.focused,
                response.uniform,
                ideal_response.focused,
                dt_s,
            )
            delay_error = case.delays .- ideal_delays
            push!(rows, (
                aperture=aperture.name,
                case=case.label,
                element_count=config.element_count,
                element_width_mm=config.element_width_m * 1e3,
                pitch_mm=config.pitch_m * 1e3,
                aperture_mm=((config.element_count - 1) * config.pitch_m + config.element_width_m) * 1e3,
                focal_distance_mm=config.focal_distance_m * 1e3,
                delay_states=case.states,
                loss_scale=case.loss_scale,
                maximum_required_delay_us=maximum(ideal_delays) * 1e6,
                maximum_delay_error_us=maximum(abs, delay_error) * 1e6,
                rms_delay_error_us=sqrt(mean(abs2, delay_error)) * 1e6,
                gain_peak=metrics.gain_peak,
                broadening_ratio=metrics.broadening_ratio,
                pulse_correlation=metrics.pulse_correlation,
                postcursor_ratio=metrics.postcursor_ratio,
                pass_gain=metrics.gain_peak >= 2.0,
                pass_broadening=metrics.broadening_ratio <= 1.25,
                pass_correlation=metrics.pulse_correlation >= 0.90,
                pass_postcursor=metrics.postcursor_ratio <= 0.10,
            ))
            waveforms[(aperture.name, case.label)] = response
        end
    end
    rows, waveforms
end

function save_pulse_plot(spectrum, summary)
    frequency_khz = spectrum.frequency_hz ./ 1e3
    db = 20 .* log10.(max.(spectrum.relative_amplitude, 1.0e-8))
    panel = plot(
        frequency_khz,
        db;
        xlim=(120, 360),
        ylim=(-60, 1),
        xlabel="frequency, kHz",
        ylabel="relative amplitude, dB",
        title="5-cycle Hann pulse spectrum",
        label=false,
        linewidth=2,
    )
    hline!(panel, [-6, -20]; color=[:darkorange :firebrick], linestyle=:dash, label=["-6 dB" "-20 dB"])
    vline!(panel, [253.0]; color=:black, linestyle=:dot, label="extra symmetric cut-on")
    annotate!(
        panel,
        128,
        -52,
        text(@sprintf("B-20: %.1f--%.1f kHz", summary.band_minus20_lower_hz / 1e3,
                      summary.band_minus20_upper_hz / 1e3), 9, :left),
    )
    savefig(panel, joinpath(OUTPUT_ROOT, "pulse_spectrum.png"))
end

function save_sensitivity_plot(rows)
    panels = Any[]
    for aperture in getproperty.(APERTURES, :name)
        selected = filter(row -> row.aperture == aperture && occursin("quantized", row.case), rows)
        lossless = sort(filter(row -> row.loss_scale == 0.0, selected); by=row -> row.delay_states)
        lossy = sort(filter(row -> row.loss_scale == 1.0, selected); by=row -> row.delay_states)
        panel = plot(
            getproperty.(lossless, :delay_states),
            getproperty.(lossless, :gain_peak);
            marker=:circle,
            linewidth=2,
            label="lossless",
            xlabel="delay states",
            ylabel="G_peak over uniform",
            title=aperture,
            ylim=(0, max(3.5, maximum(getproperty.(lossless, :gain_peak)) * 1.1)),
        )
        plot!(panel, getproperty.(lossy, :delay_states), getproperty.(lossy, :gain_peak);
              marker=:diamond, linewidth=2, label="Rayleigh loss x1")
        hline!(panel, [2.0]; color=:black, linestyle=:dash, label="G=2")
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(1, length(panels)), size=(1200, 500)),
            joinpath(OUTPUT_ROOT, "delay_state_sensitivity.png"))
end

function save_waveform_plot(spectrum, waveforms)
    time_us = spectrum.time_s .* 1e6
    panels = Any[]
    for aperture in getproperty.(APERTURES, :name)
        exact = waveforms[(aperture, "exact_lossless")]
        lossy = waveforms[(aperture, "exact_loss_1p0")]
        peak_index = argmax(abs.(exact.focused))
        time_min = max(0.0, time_us[peak_index] - 20.0)
        time_max = time_us[peak_index] + 35.0
        panel = plot(
            time_us,
            exact.uniform;
            xlim=(time_min, time_max),
            xlabel="time, us",
            ylabel="focus observable, a.u.",
            title=aperture,
            label="uniform",
            linewidth=1.5,
        )
        plot!(panel, time_us, exact.focused; label="ideal TTD", linewidth=2)
        plot!(panel, time_us, lossy.focused; label="TTD + Rayleigh loss", linewidth=2)
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(length(panels), 1), size=(1100, 750)),
            joinpath(OUTPUT_ROOT, "ideal_focus_waveforms.png"))
end

function spatial_rows_and_plot(spectrum)
    rows = NamedTuple[]
    panels = Any[]
    for aperture in APERTURES
        config = aperture.config
        delays = ideal_delays_s(config)
        axial_m = collect(5.0e-3:0.5e-3:75.0e-3)
        transverse_m = collect(-35.0e-3:0.5e-3:35.0e-3)
        axial = peak_profile(spectrum, config, delays, axial_m; direction=:axial)
        transverse = peak_profile(spectrum, config, delays, transverse_m; direction=:transverse)
        axial_width = contiguous_width(axial_m, axial)
        transverse_width = contiguous_width(transverse_m, transverse)
        push!(rows, (
            aperture=aperture.name,
            axial_peak_x_mm=axial_width.peak_coordinate * 1e3,
            dof_minus3db_mm=axial_width.width * 1e3,
            transverse_peak_y_mm=transverse_width.peak_coordinate * 1e3,
            transverse_fwhm_mm=transverse_width.width * 1e3,
        ))
        axial_panel = plot(
            axial_m .* 1e3,
            axial ./ maximum(axial);
            xlabel="x, mm",
            ylabel="normalized pulse peak",
            title="$(aperture.name): axial",
            label=false,
            linewidth=2,
        )
        hline!(axial_panel, [inv(sqrt(2.0))]; color=:black, linestyle=:dash, label="-3 dB")
        transverse_panel = plot(
            transverse_m .* 1e3,
            transverse ./ maximum(transverse);
            xlabel="y, mm",
            ylabel="normalized pulse peak",
            title="$(aperture.name): transverse",
            label=false,
            linewidth=2,
        )
        hline!(transverse_panel, [inv(sqrt(2.0))]; color=:black, linestyle=:dash, label="-3 dB")
        push!(panels, axial_panel, transverse_panel)
    end
    write_rows(joinpath(OUTPUT_ROOT, "ideal_aperture_spatial_metrics.csv"), rows)
    savefig(plot(panels...; layout=(length(APERTURES), 2), size=(1200, 800)),
            joinpath(OUTPUT_ROOT, "ideal_aperture_spatial_profiles.png"))
    rows
end

function run()
    mkpath(OUTPUT_ROOT)
    spectrum = pulse_spectrum()
    summary = pulse_summary(spectrum)
    rows, waveforms = aperture_rows(spectrum)
    write_rows(joinpath(OUTPUT_ROOT, "ideal_aperture_sensitivity.csv"), rows)
    save_pulse_plot(spectrum, summary)
    save_sensitivity_plot(rows)
    save_waveform_plot(spectrum, waveforms)
    spatial_rows = spatial_rows_and_plot(spectrum)

    @printf(
        "[+] B-6 = %.2f--%.2f kHz; B-20 = %.2f--%.2f kHz\n",
        summary.band_minus6_lower_hz / 1e3,
        summary.band_minus6_upper_hz / 1e3,
        summary.band_minus20_lower_hz / 1e3,
        summary.band_minus20_upper_hz / 1e3,
    )
    @printf("[+] spectral energy above 253 kHz = %.3f%%\n", 100summary.energy_above_253khz)
    for aperture in getproperty.(APERTURES, :name)
        exact = only(filter(row -> row.aperture == aperture && row.case == "exact_lossless", rows))
        lossy = only(filter(row -> row.aperture == aperture && row.case == "exact_loss_1p0", rows))
        @printf(
            "[+] %s: max delay %.3f us; ideal G=%.3f; optimistic lossy G=%.3f\n",
            aperture,
            exact.maximum_required_delay_us,
            exact.gain_peak,
            lossy.gain_peak,
        )
    end
    for row in spatial_rows
        @printf(
            "[+] %s: broadband DOF=%.2f mm; transverse FWHM=%.2f mm\n",
            row.aperture,
            row.dof_minus3db_mm,
            row.transverse_fwhm_mm,
        )
    end
    println("[+] $OUTPUT_ROOT")
    (; spectrum, summary, rows, waveforms, spatial_rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
