module DispersiveDelayProxy

ENV["GKSwstype"] = "100"

using FFTW
using Plots
using Printf

include(joinpath(@__DIR__, "port_mode_solver.jl"))
using .ElasticPortModes
include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .ImpulseRiskAnalysis

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_DISPERSIVE_DELAY_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "dispersive_delay_proxy"),
)
const PORT_HEIGHT_M = 3.2e-3
const PORT_ELEMENTS = 60
const APERTURE = ApertureConfig(
    element_count=19,
    element_width_m=3.2e-3,
    pitch_m=3.8e-3,
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

function branch_data(upper_hz)
    frequencies_hz = collect(10.0e3:2.0e3:ceil(upper_hz / 2e3) * 2e3)
    branch = track_fundamental_branch(
        PortModeConfig(
            height_m=PORT_HEIGHT_M,
            frequency_hz=first(frequencies_hz),
            element_count=PORT_ELEMENTS,
        ),
        frequencies_hz,
    )
    wavenumber = real.(getproperty.(branch.modes, :wavenumber_per_m))
    group_velocity = similar(wavenumber)
    group_velocity[1] = 2pi * (frequencies_hz[2] - frequencies_hz[1]) /
                        (wavenumber[2] - wavenumber[1])
    for index in 2:(length(frequencies_hz) - 1)
        group_velocity[index] = 2pi *
            (frequencies_hz[index + 1] - frequencies_hz[index - 1]) /
            (wavenumber[index + 1] - wavenumber[index - 1])
    end
    group_velocity[end] = 2pi * (frequencies_hz[end] - frequencies_hz[end - 1]) /
                          (wavenumber[end] - wavenumber[end - 1])
    rows = [(
        frequency_hz=frequencies_hz[index],
        wavenumber_per_m=wavenumber[index],
        phase_velocity_m_per_s=2pi * frequencies_hz[index] / wavenumber[index],
        group_velocity_m_per_s=group_velocity[index],
        p_fraction=branch.modes[index].p_fraction,
        axial_fraction=branch.modes[index].axial_displacement_fraction,
        overlap=branch.overlaps[index],
    ) for index in eachindex(frequencies_hz)]
    (; frequencies_hz, wavenumber, group_velocity, rows)
end

function interpolate_linear(x, y, value)
    first(x) <= value <= last(x) || return NaN
    upper = searchsortedfirst(x, value)
    upper == firstindex(x) && return y[upper]
    upper > lastindex(x) && return y[end]
    lower = upper - 1
    weight = (value - x[lower]) / (x[upper] - x[lower])
    (1 - weight) * y[lower] + weight * y[upper]
end

function truncated_spectrum(spectrum, band)
    values = copy(spectrum.spectrum)
    active = (spectrum.frequency_hz .>= band.lower_hz) .&
             (spectrum.frequency_hz .<= band.upper_hz)
    values[.!active] .= 0.0
    relative = abs.(values) ./ maximum(abs, values)
    PulseSpectrum(
        spectrum.config,
        spectrum.time_s,
        irfft(values, length(spectrum.signal)),
        spectrum.frequency_hz,
        values,
        relative,
    )
end

function proxy_transfer(
    spectrum,
    branch,
    target_delays_s;
    model::Symbol,
    loss_scale::Real,
)
    model in (:ideal_ttd, :phase_matched_path, :group_matched_trim) ||
        throw(ArgumentError("unknown delay proxy model"))
    config = APERTURE
    centers = aperture_centers_m(config)
    distances = hypot.(config.focal_distance_m, centers)
    omega0 = 2pi * spectrum.config.center_frequency_hz
    k0 = interpolate_linear(branch.frequencies_hz, branch.wavenumber,
                            spectrum.config.center_frequency_hz)
    vg0 = interpolate_linear(branch.frequencies_hz, branch.group_velocity,
                             spectrum.config.center_frequency_hz)
    path_lengths_m = if model == :phase_matched_path
        omega0 .* target_delays_s ./ k0
    elseif model == :group_matched_trim
        vg0 .* target_delays_s
    else
        zeros(length(target_delays_s))
    end
    phase_trim_rad = model == :group_matched_trim ?
        k0 .* path_lengths_m .- omega0 .* target_delays_s :
        zeros(length(target_delays_s))

    transfer = zeros(ComplexF64, length(spectrum.frequency_hz))
    for (frequency_index, frequency_hz) in enumerate(spectrum.frequency_hz)
        iszero(spectrum.spectrum[frequency_index]) && continue
        omega = 2pi * frequency_hz
        k = interpolate_linear(branch.frequencies_hz, branch.wavenumber, frequency_hz)
        vg = interpolate_linear(branch.frequencies_hz, branch.group_velocity, frequency_hz)
        isfinite(k) && isfinite(vg) || continue
        wavelength_m = config.propagation_speed_m_per_s / frequency_hz
        decay_rate = loss_scale * (
            config.rayleigh_alpha_per_s + config.rayleigh_beta_s * omega^2
        ) / 2
        transfer[frequency_index] = sum(eachindex(centers)) do index
            propagation = config.element_width_m /
                          sqrt(wavelength_m * distances[index]) *
                          cis(-omega * distances[index] / config.propagation_speed_m_per_s)
            if model == :ideal_ttd
                delay = target_delays_s[index]
                cell = exp(-decay_rate * delay) * cis(-omega * delay)
            else
                modal_delay = path_lengths_m[index] / vg
                cell = exp(-decay_rate * modal_delay) *
                       cis(-k * path_lengths_m[index] + phase_trim_rad[index])
            end
            propagation * cell
        end
    end
    (; transfer, path_lengths_m, phase_trim_rad, k0, vg0)
end

function uniform_transfer(spectrum)
    config = APERTURE
    centers = aperture_centers_m(config)
    distances = hypot.(config.focal_distance_m, centers)
    result = zeros(ComplexF64, length(spectrum.frequency_hz))
    for (index, frequency_hz) in enumerate(spectrum.frequency_hz)
        iszero(spectrum.spectrum[index]) && continue
        omega = 2pi * frequency_hz
        wavelength_m = config.propagation_speed_m_per_s / frequency_hz
        result[index] = sum(eachindex(centers)) do cell_index
            config.element_width_m /
            sqrt(wavelength_m * distances[cell_index]) *
            cis(-omega * distances[cell_index] / config.propagation_speed_m_per_s)
        end
    end
    result
end

function run_cases(spectrum, branch)
    exact_delays = ideal_delays_s(APERTURE)
    delay_masks = (
        (name="exact", states=0, delays=exact_delays),
        (name="four_state", states=4, delays=quantize_delays(exact_delays, 4)),
        (name="six_state", states=6, delays=quantize_delays(exact_delays, 6)),
        (name="eight_state", states=8, delays=quantize_delays(exact_delays, 8)),
    )
    uniform = irfft(spectrum.spectrum .* uniform_transfer(spectrum), length(spectrum.signal))
    ideal_proxy = proxy_transfer(
        spectrum,
        branch,
        exact_delays;
        model=:ideal_ttd,
        loss_scale=0.0,
    )
    ideal_waveform = irfft(spectrum.spectrum .* ideal_proxy.transfer, length(spectrum.signal))
    dt_s = spectrum.time_s[2] - spectrum.time_s[1]
    rows = NamedTuple[]
    waveforms = Dict{String, Vector{Float64}}("uniform" => uniform, "ideal_reference" => ideal_waveform)
    for mask in delay_masks
        for model in (:ideal_ttd, :phase_matched_path, :group_matched_trim)
            for loss_scale in (0.0, 0.5, 1.0, 1.5, 2.0)
                proxy = proxy_transfer(
                    spectrum,
                    branch,
                    mask.delays;
                    model,
                    loss_scale,
                )
                waveform = irfft(spectrum.spectrum .* proxy.transfer, length(spectrum.signal))
                metrics = pulse_metrics(waveform, uniform, ideal_waveform, dt_s)
                loss_slug = replace(string(loss_scale), '.' => 'p')
                label = "$(mask.name)_$(model)_loss_$(loss_slug)"
                waveforms[label] = waveform
                push!(rows, (
                    delay_mask=mask.name,
                    delay_states=mask.states,
                    model=string(model),
                    loss_scale,
                    phase_velocity_at_f0_m_per_s=omega0(spectrum) / proxy.k0,
                    group_velocity_at_f0_m_per_s=proxy.vg0,
                    maximum_target_delay_us=maximum(mask.delays) * 1e6,
                    maximum_extra_path_mm=maximum(proxy.path_lengths_m) * 1e3,
                    maximum_phase_trim_rad=maximum(abs, proxy.phase_trim_rad),
                    gain_peak=metrics.gain_peak,
                    gain_if_common_T08=metrics.gain_peak * sqrt(0.8),
                    minimum_common_power_transmission_for_G2=(2.0 / metrics.gain_peak)^2,
                    broadening_ratio=metrics.broadening_ratio,
                    pulse_correlation=metrics.pulse_correlation,
                    postcursor_ratio=metrics.postcursor_ratio,
                    pass_all=metrics.gain_peak >= 2.0 &&
                             metrics.broadening_ratio <= 1.25 &&
                             metrics.pulse_correlation >= 0.90 &&
                             metrics.postcursor_ratio <= 0.10,
                ))
            end
        end
    end
    rows, waveforms
end

omega0(spectrum) = 2pi * spectrum.config.center_frequency_hz

function save_branch_plot(branch, band)
    selected = filter(row -> band.lower_hz <= row.frequency_hz <= band.upper_hz, branch.rows)
    panel = plot(
        getproperty.(selected, :frequency_hz) ./ 1e3,
        getproperty.(selected, :phase_velocity_m_per_s);
        xlabel="frequency, kHz",
        ylabel="velocity, m/s",
        linewidth=2,
        label="phase velocity",
        title="3.2 mm single-symmetric-mode branch",
    )
    plot!(panel, getproperty.(selected, :frequency_hz) ./ 1e3,
          getproperty.(selected, :group_velocity_m_per_s);
          linewidth=2, label="group velocity")
    savefig(panel, joinpath(OUTPUT_ROOT, "single_mode_branch_dispersion.png"))
end

function save_waveform_plot(spectrum, waveforms)
    time_us = spectrum.time_s .* 1e6
    labels = (
        ("ideal_reference", "ideal TTD"),
        ("exact_phase_matched_path_loss_0p0", "phase-matched path"),
        ("exact_group_matched_trim_loss_0p0", "group-matched + phase trim"),
        ("six_state_phase_matched_path_loss_1p0", "6-state path + loss"),
    )
    ideal = waveforms["ideal_reference"]
    peak_index = argmax(abs.(ideal))
    panel = plot(
        time_us,
        waveforms["uniform"];
        xlim=(max(0.0, time_us[peak_index] - 20), time_us[peak_index] + 35),
        xlabel="time, us",
        ylabel="focus observable, a.u.",
        linewidth=1.5,
        label="uniform",
        title="Dispersive geometric-delay proxy",
    )
    for (key, label) in labels
        plot!(panel, time_us, waveforms[key]; linewidth=2, label)
    end
    savefig(panel, joinpath(OUTPUT_ROOT, "dispersive_delay_waveforms.png"))
end

function run()
    mkpath(OUTPUT_ROOT)
    raw_spectrum = pulse_spectrum()
    band20 = spectral_band(raw_spectrum, -20.0)
    spectrum = truncated_spectrum(raw_spectrum, band20)
    branch = branch_data(band20.upper_hz)
    rows, waveforms = run_cases(spectrum, branch)
    write_rows(joinpath(OUTPUT_ROOT, "single_mode_branch.csv"), branch.rows)
    write_rows(joinpath(OUTPUT_ROOT, "dispersive_delay_proxy.csv"), rows)
    save_branch_plot(branch, band20)
    save_waveform_plot(spectrum, waveforms)

    for row in rows
        row.loss_scale == 0.0 || continue
        @printf(
            "[+] %s / %s: G=%.3f, Bt=%.3f, rho=%.4f, path=%.2f mm, trim=%.2f rad, pass=%s\n",
            row.delay_mask,
            row.model,
            row.gain_peak,
            row.broadening_ratio,
            row.pulse_correlation,
            row.maximum_extra_path_mm,
            row.maximum_phase_trim_rad,
            row.pass_all,
        )
    end
    println("[+] $OUTPUT_ROOT")
    (; spectrum, band20, branch, rows, waveforms)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
