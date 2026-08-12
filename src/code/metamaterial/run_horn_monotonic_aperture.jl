module HornMonotonicAperture

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BINARY_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_long_smooth_aperture")
const PILOT_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_pilot")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_MONOTONIC_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_aperture"),
)

ENV["GKSwstype"] = "100"
using FFTW
using JLD2
using Plots

include(joinpath(@__DIR__, "monotonic_horn_lens.jl"))
include(joinpath(@__DIR__, "run_horn_long_smooth_aperture.jl"))
using .MonotonicHornLens
using .HornLongSmoothAperture
using .HornLongSmoothAperture.ImpulseRiskAnalysis
using .HornLongSmoothAperture.SpectralAnalysis

const F0_HZ = HornLongSmoothAperture.F0_HZ
const REFERENCE_EXTRA_PATH_MM = 4.42
const SCALE_MIN = 0.5
const SCALE_MAX = 3.2
const SCALE_STEP = 0.01
const CONFIG = MonotonicHornLensConfig()

function unwrap_near(phase_rad::Real, reference_rad::Real)
    Float64(phase_rad) + 2pi * round((Float64(reference_rad) - phase_rad) / (2pi))
end

function transfer_log_on_pulse_grid(result, frequencies_hz)
    segment = HornLongSmoothAperture.carrier_segment(result)
    source_frequency = result.frequency_hz[segment]
    source_log_amplitude = log.(result.amplitude[segment])
    source_phase = result.phase_rad[segment]
    active = findall((frequencies_hz .>= first(source_frequency)) .&
                     (frequencies_hz .<= last(source_frequency)))
    log_transfer = zeros(ComplexF64, length(frequencies_hz))
    for index in active
        frequency_hz = frequencies_hz[index]
        log_transfer[index] = HornLongSmoothAperture.interpolate_linear(
            source_frequency, source_log_amplitude, frequency_hz,
        ) + im * HornLongSmoothAperture.interpolate_linear(
            source_frequency, source_phase, frequency_hz,
        )
    end
    carrier_log = HornLongSmoothAperture.interpolate_linear(
        source_frequency, source_log_amplitude, F0_HZ,
    ) + im * HornLongSmoothAperture.interpolate_linear(
        source_frequency, source_phase, F0_HZ,
    )
    (; log_transfer, active, carrier_log,
       lower_hz=first(source_frequency), upper_hz=last(source_frequency))
end

function context_log_correction(reference_carrier_log)
    library = HornLongSmoothAperture.load_context_library()
    context = library.interior["SSS"]
    phase = unwrap_near(angle(context), imag(reference_carrier_log))
    log(abs(context)) + im * phase - reference_carrier_log
end

function pilot_relative_transfer()
    maximum_path = joinpath(PILOT_ROOT, "signals", "monotonic_max.jld2")
    device_path = joinpath(PILOT_ROOT, "signals", "monotonic_device.jld2")
    isfile(maximum_path) && isfile(device_path) || return nothing
    maximum_data = JLD2.load(maximum_path)
    device_data = JLD2.load(device_path)
    probe_name = "near_radiator"
    maximum_index = findfirst(==(probe_name), maximum_data["probe_names"])
    device_index = findfirst(==(probe_name), device_data["probe_names"])
    spectral_config = SpectrumConfig(
        window=:rectangular,
        zero_padding_factor=16,
        input_floor_relative=1e-3,
    )
    maximum_transfer = analyze_transfer(
        maximum_data["time_s"],
        maximum_data["source_drive_mpa"],
        vec(maximum_data["probe_velocity_x_m_per_s"][maximum_index, :]);
        config=spectral_config,
    )
    device_transfer = analyze_transfer(
        device_data["time_s"],
        device_data["source_drive_mpa"],
        vec(device_data["probe_velocity_x_m_per_s"][device_index, :]);
        config=spectral_config,
    )
    relative_transfer(maximum_transfer, device_transfer)
end

function transfer_from_endpoint_log(endpoint_log, fractions, kernels, active)
    result = zeros(ComplexF64, size(kernels, 1))
    for frequency_index in active
        result[frequency_index] = sum(eachindex(fractions)) do element_index
            exp(fractions[element_index] * endpoint_log[frequency_index]) *
            kernels[frequency_index, element_index]
        end
    end
    result
end

function optimize_continuation(spectrum, reference, correction_log, fractions, kernels)
    best_scale, best_peak = NaN, -Inf
    best_transfer = ComplexF64[]
    scale_grid = collect(SCALE_MIN:SCALE_STEP:SCALE_MAX)
    peaks = similar(scale_grid)
    for (scale_index, scale) in enumerate(scale_grid)
        endpoint_log = scale .* (reference.log_transfer .+ correction_log)
        transfer = transfer_from_endpoint_log(endpoint_log, fractions, kernels, reference.active)
        peak = maximum(analytic_envelope(HornLongSmoothAperture.waveform(spectrum, transfer)))
        peaks[scale_index] = peak
        if peak > best_peak
            best_scale, best_peak, best_transfer = scale, peak, transfer
        end
    end
    coarse_index = argmax(peaks)
    lower = scale_grid[max(firstindex(scale_grid), coarse_index - 2)]
    upper = scale_grid[min(lastindex(scale_grid), coarse_index + 2)]
    fine_grid = collect(lower:0.00025:upper)
    fine_peaks = similar(fine_grid)
    for (scale_index, scale) in enumerate(fine_grid)
        endpoint_log = scale .* (reference.log_transfer .+ correction_log)
        transfer = transfer_from_endpoint_log(endpoint_log, fractions, kernels, reference.active)
        peak = maximum(analytic_envelope(HornLongSmoothAperture.waveform(spectrum, transfer)))
        fine_peaks[scale_index] = peak
        if peak > best_peak
            best_scale, best_peak, best_transfer = scale, peak, transfer
        end
    end
    (; scale=best_scale, peak=best_peak, transfer=best_transfer,
       scale_grid=vcat(scale_grid, fine_grid), peak_grid=vcat(peaks, fine_peaks))
end

function uniform_transfer(kernels, active)
    result = zeros(ComplexF64, size(kernels, 1))
    result[active] .= vec(sum(kernels[active, :]; dims=2))
    result
end

function profile(endpoint_log, fractions, spectrum, active, coordinate; direction)
    [begin
        x_m, y_m = direction == :transverse ?
            (HornLongSmoothAperture.FOCUS_M, value) : (value, 0.0)
        kernels = HornLongSmoothAperture.kernel_matrix(
            spectrum.frequency_hz, active, x_m, y_m,
        )
        transfer = transfer_from_endpoint_log(endpoint_log, fractions, kernels, active)
        maximum(analytic_envelope(HornLongSmoothAperture.waveform(spectrum, transfer)))
    end for value in coordinate]
end

function axial_metrics(x_mm, axial, target_x_mm)
    target_index = argmin(abs.(x_mm .- target_x_mm))
    threshold = axial[target_index] / sqrt(2)
    left, right = target_index, target_index
    while left > firstindex(axial) && axial[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(axial) && axial[right + 1] >= threshold
        right += 1
    end
    local_peak = left - 1 + argmax(axial[left:right])
    (; width_mm=x_mm[right] - x_mm[left], peak_x_mm=x_mm[local_peak])
end

function save_geometry_figure(path, geometry)
    profile_panel = plot(
        xlabel="axial coordinate in delay section, mm",
        ylabel="out-of-plane bend z, mm",
        title="Eight monotonic sin⁴ delay states",
        gridalpha=0.2,
        legend=:topleft,
    )
    x_mm = range(0.0, geometry.axial_length_mm; length=401)
    unique_indices = 8:15
    palette = cgrad(:viridis, length(unique_indices), categorical=true)
    for (color_index, element_index) in enumerate(unique_indices)
        z_mm = smooth_centerline_mm.(
            x_mm, geometry.axial_length_mm, geometry.amplitude_mm[element_index],
        )
        plot!(
            profile_panel,
            x_mm,
            z_mm;
            linewidth=2.2,
            color=palette[color_index],
            label="y=$(round(geometry.centers_mm[element_index]; digits=1)) mm",
        )
    end

    delay_us = geometry.delay_s .* 1e6
    delay_panel = plot(
        geometry.centers_mm,
        delay_us;
        marker=:circle,
        linewidth=2.5,
        label="required delay",
        xlabel="aperture y, mm",
        ylabel="delay, μs",
        title="Unwrapped true-time-delay law",
        gridalpha=0.2,
    )
    plot!(
        twinx(delay_panel),
        geometry.centers_mm,
        geometry.extra_path_mm;
        marker=:diamond,
        linewidth=2,
        color=:darkorange,
        label="extra path",
        ylabel="extra path, mm",
    )

    layout = @layout [a{0.63w} b]
    savefig(plot(profile_panel, delay_panel; layout, size=(1450, 650), margin=6Plots.mm), path)
end

function run()
    mkpath(OUTPUT_ROOT)
    spectrum = pulse_spectrum(PulseConfig(center_frequency_hz=F0_HZ, cycles=5.0))
    reference_result = HornLongSmoothAperture.transient_relative_transfer()
    reference = transfer_log_on_pulse_grid(reference_result, spectrum.frequency_hz)
    correction_log = context_log_correction(reference.carrier_log)
    fractions = normalized_delays(CONFIG)
    focus_kernels = HornLongSmoothAperture.kernel_matrix(
        spectrum.frequency_hz, reference.active, HornLongSmoothAperture.FOCUS_M, 0.0,
    )
    provisional = optimize_continuation(
        spectrum, reference, correction_log, fractions, focus_kernels,
    )

    physical = pilot_relative_transfer()
    if isnothing(physical)
        source_kind = "calibrated_continuation"
        maximum_scale = provisional.scale
        endpoint_log = maximum_scale .* (reference.log_transfer .+ correction_log)
        focused_transfer = provisional.transfer
    else
        source_kind = "maximum_state_transient"
        maximum_scale = provisional.scale
        physical_grid = transfer_log_on_pulse_grid(physical, spectrum.frequency_hz)
        expected_carrier_phase = maximum_scale * imag(reference.carrier_log)
        phase_shift = 2pi * round(
            (expected_carrier_phase - imag(physical_grid.carrier_log)) / (2pi),
        )
        aligned_log_transfer = copy(physical_grid.log_transfer)
        aligned_log_transfer[physical_grid.active] .+= im * phase_shift
        physical_correction = maximum_scale .* correction_log
        endpoint_log = aligned_log_transfer .+ physical_correction
        focus_kernels = HornLongSmoothAperture.kernel_matrix(
            spectrum.frequency_hz,
            physical_grid.active,
            HornLongSmoothAperture.FOCUS_M,
            0.0,
        )
        focused_transfer = transfer_from_endpoint_log(
            endpoint_log, fractions, focus_kernels, physical_grid.active,
        )
        reference = physical_grid
    end

    maximum_extra_path_mm = maximum_scale * REFERENCE_EXTRA_PATH_MM
    geometry = synthesize_monotonic_geometry(maximum_extra_path_mm; config=CONFIG)
    focused = HornLongSmoothAperture.waveform(spectrum, focused_transfer)
    uniform_focus_transfer = uniform_transfer(focus_kernels, reference.active)
    uniform = HornLongSmoothAperture.waveform(spectrum, uniform_focus_transfer)
    ideal_transfer = HornLongSmoothAperture.ideal_transfer(
        spectrum.frequency_hz, reference.active, focus_kernels,
    )
    ideal = HornLongSmoothAperture.waveform(spectrum, ideal_transfer)
    dt_s = spectrum.time_s[2] - spectrum.time_s[1]
    metrics = pulse_metrics(focused, uniform, ideal, dt_s)

    transverse_y_m = collect(-20.0e-3:0.25e-3:20.0e-3)
    axial_x_m = collect(15.0e-3:0.5e-3:55.0e-3)
    transverse = profile(
        endpoint_log, fractions, spectrum, reference.active, transverse_y_m;
        direction=:transverse,
    )
    axial = profile(
        endpoint_log, fractions, spectrum, reference.active, axial_x_m;
        direction=:axial,
    )
    transverse_width = contiguous_width(transverse_y_m .* 1e3, transverse)
    axial_result = axial_metrics(axial_x_m .* 1e3, axial, CONFIG.focal_distance_mm)
    side_ratio = HornLongSmoothAperture.sidelobe_ratio(transverse)
    carrier_index = argmin(abs.(spectrum.frequency_hz .- F0_HZ))
    weights = exp.(fractions .* endpoint_log[carrier_index])
    carrier_gain = abs(focused_transfer[carrier_index]) /
                   abs(uniform_focus_transfer[carrier_index])
    effective_delay_us = -imag(endpoint_log[carrier_index]) / (2pi * F0_HZ) * 1e6

    selection_path = joinpath(OUTPUT_ROOT, "monotonic_lens_geometry.csv")
    open(selection_path, "w") do io
        println(io, "element_index,center_y_mm,delay_fraction,required_delay_us,extra_path_mm,path_length_mm,bend_amplitude_mm,minimum_inner_radius_mm,weight_real,weight_imag,weight_amplitude,weight_phase_deg")
        for index in eachindex(fractions)
            println(io, join((
                index,
                geometry.centers_mm[index],
                fractions[index],
                geometry.delay_s[index] * 1e6,
                geometry.extra_path_mm[index],
                geometry.path_length_mm[index],
                geometry.amplitude_mm[index],
                geometry.inner_radius_mm[index],
                real(weights[index]),
                imag(weights[index]),
                abs(weights[index]),
                rad2deg(angle(weights[index])),
            ), ','))
        end
    end

    summary_path = joinpath(OUTPUT_ROOT, "monotonic_aperture_impulse_summary.csv")
    open(summary_path, "w") do io
        println(io, "transfer_source,maximum_scale,maximum_extra_path_mm,common_axial_length_mm,maximum_bend_amplitude_mm,minimum_inner_radius_mm,required_maximum_delay_us,effective_carrier_delay_us,calibrated_lower_hz,calibrated_upper_hz,impulse_peak_gain,carrier_gain,broadening_ratio,pulse_correlation,postcursor_ratio,transverse_fwhm_mm,axial_dof_mm,peak_x_mm,peak_y_mm,sidelobe_amplitude_ratio")
        println(io, join((
            source_kind,
            maximum_scale,
            maximum_extra_path_mm,
            geometry.axial_length_mm,
            maximum(geometry.amplitude_mm),
            minimum(geometry.inner_radius_mm),
            maximum(geometry.delay_s) * 1e6,
            effective_delay_us,
            reference.lower_hz,
            reference.upper_hz,
            metrics.gain_peak,
            carrier_gain,
            metrics.broadening_ratio,
            metrics.pulse_correlation,
            metrics.postcursor_ratio,
            transverse_width.width,
            axial_result.width_mm,
            axial_result.peak_x_mm,
            transverse_width.peak_coordinate,
            side_ratio,
        ), ','))
    end

    result_path = joinpath(OUTPUT_ROOT, "monotonic_aperture_impulse.jld2")
    jldsave(
        result_path;
        format_version=1,
        transfer_source=source_kind,
        maximum_scale,
        maximum_extra_path_mm,
        common_axial_length_mm=geometry.axial_length_mm,
        maximum_bend_amplitude_mm=maximum(geometry.amplitude_mm),
        centers_mm=geometry.centers_mm,
        delay_fraction=fractions,
        required_delay_s=geometry.delay_s,
        extra_path_mm=geometry.extra_path_mm,
        path_length_mm=geometry.path_length_mm,
        bend_amplitude_mm=geometry.amplitude_mm,
        inner_radius_mm=geometry.inner_radius_mm,
        weights,
        endpoint_log_transfer=endpoint_log,
        time_s=spectrum.time_s,
        focused,
        uniform,
        ideal,
        frequencies_hz=spectrum.frequency_hz,
        focused_transfer,
        uniform_transfer=uniform_focus_transfer,
        transverse_y_mm=transverse_y_m .* 1e3,
        transverse_peak=transverse,
        axial_x_mm=axial_x_m .* 1e3,
        axial_peak=axial,
        provisional_scale_grid=provisional.scale_grid,
        provisional_peak_grid=provisional.peak_grid,
    )

    waveform_panel = plot(
        spectrum.time_s .* 1e6,
        analytic_envelope(focused);
        label="monotonic TTD",
        linewidth=2.5,
        xlabel="time, μs",
        ylabel="envelope, a.u.",
        title="Five-cycle focus waveform",
        gridalpha=0.25,
        xlims=(5, 70),
    )
    plot!(waveform_panel, spectrum.time_s .* 1e6, analytic_envelope(uniform);
          label="uniform aperture", linewidth=2)
    plot!(waveform_panel, spectrum.time_s .* 1e6, analytic_envelope(ideal);
          label="lossless ideal TTD", linewidth=1.8, linestyle=:dash)
    transverse_panel = plot(
        transverse_y_m .* 1e3,
        transverse ./ maximum(transverse);
        label="monotonic TTD",
        linewidth=2.5,
        xlabel="y at x=35 mm, mm",
        ylabel="normalized peak envelope",
        title="Impulse transverse profile",
        gridalpha=0.25,
    )
    axial_panel = plot(
        axial_x_m .* 1e3,
        axial ./ maximum(axial);
        label="monotonic TTD",
        linewidth=2.5,
        xlabel="x, mm",
        ylabel="normalized peak envelope",
        title="Impulse axial profile",
        gridalpha=0.25,
    )
    delay_panel = plot(
        geometry.centers_mm,
        geometry.extra_path_mm;
        marker=:circle,
        linewidth=2.5,
        label="extra path",
        xlabel="radiator y, mm",
        ylabel="extra path, mm",
        title="Monotonic physical delay profile",
        gridalpha=0.25,
    )
    figure_path = joinpath(OUTPUT_ROOT, "monotonic_aperture_impulse.png")
    savefig(plot(
        waveform_panel, transverse_panel, axial_panel, delay_panel;
        layout=(2, 2), size=(1450, 950), margin=5Plots.mm,
    ), figure_path)
    geometry_path = joinpath(OUTPUT_ROOT, "monotonic_lens_geometry.png")
    save_geometry_figure(geometry_path, geometry)

    println("[+] transfer source=$source_kind")
    println("[+] monotonic scale=$maximum_scale, max extra path=$maximum_extra_path_mm mm")
    println("[+] common axial=$(geometry.axial_length_mm) mm, max bend=$(maximum(geometry.amplitude_mm)) mm")
    println("[+] required delay=$(maximum(geometry.delay_s) * 1e6) μs, carrier delay=$effective_delay_us μs")
    println("[+] impulse gain=$(metrics.gain_peak), carrier gain=$carrier_gain")
    println("[+] Bt=$(metrics.broadening_ratio), rho=$(metrics.pulse_correlation), postcursor=$(metrics.postcursor_ratio)")
    println("[+] FWHM=$(transverse_width.width) mm, DOF=$(axial_result.width_mm) mm, sidelobe=$side_ratio")
    println("[+] $summary_path")
    println("[+] $selection_path")
    println("[+] $figure_path")
    println("[+] $geometry_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
