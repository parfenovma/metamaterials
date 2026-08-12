module HornLongSmoothAperture

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TRANSIENT_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_long_smooth_pilot", "signals")
const NEIGHBOR_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_neighbor_3d_242khz")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_LONG_APERTURE_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_long_smooth_aperture"),
)

ENV["GKSwstype"] = "100"
using FFTW
using JLD2
using Plots

include(joinpath(@__DIR__, "transient_model_solver.jl"))
include(joinpath(@__DIR__, "spectral_analysis.jl"))
include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .SpectralAnalysis
using .ImpulseRiskAnalysis

const F0_HZ = 242.0e3
const C_P_M_S = 2340.0
const ELEMENT_COUNT = 15
const ELEMENT_WIDTH_M = 1.6e-3
const PITCH_M = 4.8e-3
const FOCUS_M = 35.0e-3
const QUADRATURE_POINTS = 7
const STATES = ('D', 'S')

signal_path(name) = joinpath(TRANSIENT_ROOT, "$(name).jld2")
neighbor_path(name) = joinpath(NEIGHBOR_ROOT, "$(name)_harmonic.jld2")

function probe(data, name)
    index = findfirst(==(name), data["probe_names"])
    isnothing(index) && error("probe $name was not found")
    vec(data["probe_displacement_m"][index, :])
end

function preout_ratio(data, control, index)
    probe(data, "preout_$index")[1] / probe(control, "preout_$index")[1]
end

function context_case(triple::AbstractString)
    direct = Dict(
        "DDD" => "long_DDD",
        "DSD" => "long_DSD",
        "DSS" => "long_DSS",
        "SDD" => "long_SDD",
        "SDS" => "long_SDS",
        "SSS" => "long_SSS",
    )
    haskey(direct, triple) && return direct[triple]
    reversed = reverse(triple)
    haskey(direct, reversed) || error("missing binary context $triple")
    direct[reversed]
end

function load_context_library()
    names = ("long_DDD", "long_DSD", "long_DSS", "long_SDD", "long_SDS", "long_SSS")
    data = Dict(name => JLD2.load(neighbor_path(name)) for name in names)
    control = data["long_DDD"]
    interior = Dict{String, ComplexF64}()
    for triple in ("DDD", "DSD", "DSS", "SDD", "SDS", "SSS", "DDS", "SSD")
        interior[triple] = preout_ratio(data[context_case(triple)], control, 2)
    end
    edge_cases = Dict(
        "DD" => ("long_DDD", 1),
        "DS" => ("long_DSD", 1),
        "SD" => ("long_SDD", 1),
        "SS" => ("long_SSS", 1),
    )
    edge = Dict(
        pair => preout_ratio(data[case_name], control, index)
        for (pair, (case_name, index)) in edge_cases
    )
    (; data, interior, edge)
end

function context_coefficients(mask, library)
    length(mask) == ELEMENT_COUNT || throw(DimensionMismatch("wrong aperture mask length"))
    coefficients = Vector{ComplexF64}(undef, ELEMENT_COUNT)
    for index in eachindex(mask)
        if index == firstindex(mask)
            coefficients[index] = library.edge[string(mask[index], mask[index + 1])]
        elseif index == lastindex(mask)
            coefficients[index] = library.edge[string(mask[index], mask[index - 1])]
        else
            triple = string(mask[index - 1], mask[index], mask[index + 1])
            coefficients[index] = library.interior[triple]
        end
    end
    coefficients
end

function transient_relative_transfer()
    device = JLD2.load(signal_path("long_device"))
    smooth = JLD2.load(signal_path("long_smooth"))
    probe_name = "near_radiator"
    device_index = findfirst(==(probe_name), device["probe_names"])
    smooth_index = findfirst(==(probe_name), smooth["probe_names"])
    config = SpectrumConfig(window=:rectangular, zero_padding_factor=16, input_floor_relative=1e-3)
    device_transfer = analyze_transfer(
        device["time_s"],
        device["source_drive_mpa"],
        vec(device["probe_velocity_x_m_per_s"][device_index, :]);
        config,
    )
    smooth_transfer = analyze_transfer(
        smooth["time_s"],
        smooth["source_drive_mpa"],
        vec(smooth["probe_velocity_x_m_per_s"][smooth_index, :]);
        config,
    )
    relative_transfer(smooth_transfer, device_transfer)
end

function carrier_segment(result)
    carrier = argmin(abs.(result.frequency_hz .- F0_HZ))
    result.valid[carrier] || error("measured transfer is invalid at the carrier")
    left, right = carrier, carrier
    while left > firstindex(result.valid) && result.valid[left - 1]
        left -= 1
    end
    while right < lastindex(result.valid) && result.valid[right + 1]
        right += 1
    end
    left:right
end

function interpolate_linear(x, y, target)
    target <= first(x) && return first(y)
    target >= last(x) && return last(y)
    right = searchsortedfirst(x, target)
    left = right - 1
    fraction = (target - x[left]) / (x[right] - x[left])
    (1 - fraction) * y[left] + fraction * y[right]
end

function measured_ratio_on_pulse_grid(result, frequencies_hz)
    segment = carrier_segment(result)
    source_frequency = result.frequency_hz[segment]
    source_amplitude = result.amplitude[segment]
    source_phase = result.phase_rad[segment]
    ratio = zeros(ComplexF64, length(frequencies_hz))
    active = findall((frequencies_hz .>= first(source_frequency)) .&
                     (frequencies_hz .<= last(source_frequency)))
    for index in active
        frequency = frequencies_hz[index]
        amplitude = interpolate_linear(source_frequency, source_amplitude, frequency)
        phase = interpolate_linear(source_frequency, source_phase, frequency)
        ratio[index] = amplitude * cis(phase)
    end
    carrier_ratio = interpolate_linear(source_frequency, source_amplitude, F0_HZ) *
                    cis(interpolate_linear(source_frequency, source_phase, F0_HZ))
    (; ratio, active, carrier_ratio,
       lower_hz=first(source_frequency), upper_hz=last(source_frequency))
end

aperture_centers_m() = collect(-((ELEMENT_COUNT - 1) ÷ 2):((ELEMENT_COUNT - 1) ÷ 2)) .* PITCH_M

function element_kernel(frequency_hz, center_y_m, x_m, y_m)
    frequency_hz > 0 || return 0.0 + 0.0im
    dy = ELEMENT_WIDTH_M / QUADRATURE_POINTS
    first_y = center_y_m - ELEMENT_WIDTH_M / 2 + dy / 2
    wavelength = C_P_M_S / frequency_hz
    k = 2pi / wavelength
    dy * sum(1:QUADRATURE_POINTS) do quadrature_index
        source_y = first_y + (quadrature_index - 1) * dy
        distance = hypot(x_m, y_m - source_y)
        cis(-k * distance + pi / 4) / sqrt(wavelength * distance)
    end
end

function kernel_matrix(frequencies_hz, active, x_m, y_m)
    centers = aperture_centers_m()
    result = zeros(ComplexF64, length(frequencies_hz), ELEMENT_COUNT)
    for frequency_index in active, element_index in eachindex(centers)
        result[frequency_index, element_index] = element_kernel(
            frequencies_hz[frequency_index], centers[element_index], x_m, y_m,
        )
    end
    result
end

function symmetric_mask(bits)
    centers = aperture_centers_m()
    groups = [((bits >> group) & 1) == 1 ? 'S' : 'D' for group in 0:7]
    [groups[round(Int, abs(center) / PITCH_M) + 1] for center in centers]
end

function aperture_transfer(mask, context, measured, kernels, active)
    scale = ComplexF64[
        mask[index] == 'S' ? context[index] / measured.carrier_ratio : context[index]
        for index in eachindex(mask)
    ]
    result = zeros(ComplexF64, size(kernels, 1))
    for frequency_index in active
        result[frequency_index] = sum(eachindex(mask)) do element_index
            state = mask[element_index] == 'S' ? measured.ratio[frequency_index] : 1.0 + 0.0im
            scale[element_index] * state * kernels[frequency_index, element_index]
        end
    end
    result
end

function ideal_transfer(frequencies_hz, active, kernels)
    centers = aperture_centers_m()
    distance = hypot.(FOCUS_M, centers)
    delay = (maximum(distance) .- distance) ./ C_P_M_S
    result = zeros(ComplexF64, length(frequencies_hz))
    for frequency_index in active
        omega = 2pi * frequencies_hz[frequency_index]
        result[frequency_index] = sum(eachindex(centers)) do element_index
            cis(-omega * delay[element_index]) * kernels[frequency_index, element_index]
        end
    end
    result
end

waveform(spectrum, transfer) = irfft(spectrum.spectrum .* transfer, length(spectrum.signal))

function optimize_mask(spectrum, measured, library, kernels)
    best_peak = -Inf
    best_mask = Char[]
    best_context = ComplexF64[]
    best_transfer = ComplexF64[]
    for bits in 0:255
        mask = symmetric_mask(bits)
        context = context_coefficients(mask, library)
        transfer = aperture_transfer(mask, context, measured, kernels, measured.active)
        peak = maximum(analytic_envelope(waveform(spectrum, transfer)))
        if peak > best_peak
            best_peak = peak
            best_mask = mask
            best_context = context
            best_transfer = transfer
        end
    end
    (; mask=best_mask, context=best_context, transfer=best_transfer, peak=best_peak)
end

function sidelobe_ratio(profile)
    center_peak = argmax(profile)
    left_minimum = center_peak
    while left_minimum > firstindex(profile) + 1
        left_minimum -= 1
        profile[left_minimum] <= profile[left_minimum - 1] &&
            profile[left_minimum] <= profile[left_minimum + 1] && break
    end
    right_minimum = center_peak
    while right_minimum < lastindex(profile) - 1
        right_minimum += 1
        profile[right_minimum] <= profile[right_minimum - 1] &&
            profile[right_minimum] <= profile[right_minimum + 1] && break
    end
    outside = vcat(collect(firstindex(profile):(left_minimum - 1)),
                   collect((right_minimum + 1):lastindex(profile)))
    isempty(outside) ? 0.0 : maximum(profile[outside]) / profile[center_peak]
end

function peak_profile(mask, context, spectrum, measured, coordinate; direction)
    [begin
        x_m, y_m = direction == :transverse ? (FOCUS_M, value) : (value, 0.0)
        kernels = kernel_matrix(spectrum.frequency_hz, measured.active, x_m, y_m)
        transfer = aperture_transfer(mask, context, measured, kernels, measured.active)
        maximum(analytic_envelope(waveform(spectrum, transfer)))
    end for value in coordinate]
end

function run()
    mkpath(OUTPUT_ROOT)
    library = load_context_library()
    relative = transient_relative_transfer()
    spectrum = pulse_spectrum(PulseConfig(center_frequency_hz=F0_HZ, cycles=5.0))
    measured = measured_ratio_on_pulse_grid(relative, spectrum.frequency_hz)
    focus_kernels = kernel_matrix(spectrum.frequency_hz, measured.active, FOCUS_M, 0.0)
    optimized = optimize_mask(spectrum, measured, library, focus_kernels)
    focused = waveform(spectrum, optimized.transfer)

    uniform_mask = fill('D', ELEMENT_COUNT)
    uniform_context = context_coefficients(uniform_mask, library)
    uniform_transfer = aperture_transfer(
        uniform_mask, uniform_context, measured, focus_kernels, measured.active,
    )
    uniform = waveform(spectrum, uniform_transfer)
    ideal = waveform(spectrum, ideal_transfer(spectrum.frequency_hz, measured.active, focus_kernels))
    dt_s = spectrum.time_s[2] - spectrum.time_s[1]
    metrics = pulse_metrics(focused, uniform, ideal, dt_s)

    transverse_y_m = collect(-20.0e-3:0.25e-3:20.0e-3)
    axial_x_m = collect(15.0e-3:0.5e-3:55.0e-3)
    transverse = peak_profile(
        optimized.mask, optimized.context, spectrum, measured, transverse_y_m;
        direction=:transverse,
    )
    axial = peak_profile(
        optimized.mask, optimized.context, spectrum, measured, axial_x_m;
        direction=:axial,
    )
    transverse_width = contiguous_width(transverse_y_m .* 1e3, transverse)
    target_index = argmin(abs.(axial_x_m .- FOCUS_M))
    threshold = axial[target_index] / sqrt(2)
    left, right = target_index, target_index
    while left > firstindex(axial) && axial[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(axial) && axial[right + 1] >= threshold
        right += 1
    end
    axial_dof_mm = (axial_x_m[right] - axial_x_m[left]) * 1e3
    local_peak = left - 1 + argmax(axial[left:right])
    carrier_index = argmin(abs.(spectrum.frequency_hz .- F0_HZ))
    carrier_gain = abs(optimized.transfer[carrier_index]) / abs(uniform_transfer[carrier_index])
    side_ratio = sidelobe_ratio(transverse)

    selection_path = joinpath(OUTPUT_ROOT, "long_smooth_binary_selection.csv")
    open(selection_path, "w") do io
        println(io, "element_index,center_y_mm,state,context_amplitude,context_phase_deg")
        for index in eachindex(optimized.mask)
            println(io, join((
                index,
                aperture_centers_m()[index] * 1e3,
                optimized.mask[index],
                abs(optimized.context[index]),
                rad2deg(angle(optimized.context[index])),
            ), ','))
        end
    end
    context_path = joinpath(OUTPUT_ROOT, "long_smooth_context_library.csv")
    open(context_path, "w") do io
        println(io, "kind,context,amplitude,phase_deg")
        for (name, value) in sort(collect(library.interior); by=first)
            println(io, "interior,$name,$(abs(value)),$(rad2deg(angle(value)))")
        end
        for (name, value) in sort(collect(library.edge); by=first)
            println(io, "edge,$name,$(abs(value)),$(rad2deg(angle(value)))")
        end
    end
    summary_path = joinpath(OUTPUT_ROOT, "long_smooth_aperture_impulse_summary.csv")
    open(summary_path, "w") do io
        println(io, "mask,smooth_count,calibrated_lower_hz,calibrated_upper_hz,impulse_peak_gain,carrier_gain,broadening_ratio,pulse_correlation,postcursor_ratio,transverse_fwhm_mm,axial_dof_mm,peak_x_mm,peak_y_mm,sidelobe_amplitude_ratio")
        println(io, join((
            String(optimized.mask),
            count(==('S'), optimized.mask),
            measured.lower_hz,
            measured.upper_hz,
            metrics.gain_peak,
            carrier_gain,
            metrics.broadening_ratio,
            metrics.pulse_correlation,
            metrics.postcursor_ratio,
            transverse_width.width,
            axial_dof_mm,
            axial_x_m[local_peak] * 1e3,
            transverse_width.peak_coordinate,
            side_ratio,
        ), ','))
    end
    result_path = joinpath(OUTPUT_ROOT, "long_smooth_aperture_impulse.jld2")
    jldsave(
        result_path;
        format_version=1,
        mask=String(optimized.mask),
        context_coefficients=optimized.context,
        time_s=spectrum.time_s,
        focused,
        uniform,
        ideal,
        frequencies_hz=spectrum.frequency_hz,
        focused_transfer=optimized.transfer,
        uniform_transfer,
        transverse_y_mm=transverse_y_m .* 1e3,
        transverse_peak=transverse,
        axial_x_mm=axial_x_m .* 1e3,
        axial_peak=axial,
    )

    waveform_panel = plot(
        spectrum.time_s .* 1e6, analytic_envelope(focused);
        label="context-calibrated binary", linewidth=2.5,
        xlabel="time, μs", ylabel="envelope, a.u.",
        title="Five-cycle focus waveform", gridalpha=0.25,
    )
    plot!(waveform_panel, spectrum.time_s .* 1e6, analytic_envelope(uniform);
          label="uniform device aperture", linewidth=2)
    plot!(waveform_panel, spectrum.time_s .* 1e6, analytic_envelope(ideal);
          label="ideal true-time-delay", linewidth=1.8, linestyle=:dash, xlims=(5, 65))
    transverse_panel = plot(
        transverse_y_m .* 1e3, transverse ./ maximum(transverse);
        label="binary", linewidth=2.5, xlabel="y at x=35 mm, mm",
        ylabel="normalized peak envelope", title="Impulse transverse profile", gridalpha=0.25,
    )
    axial_panel = plot(
        axial_x_m .* 1e3, axial ./ maximum(axial);
        label="binary", linewidth=2.5, xlabel="x, mm", ylabel="normalized peak envelope",
        title="Impulse axial profile", gridalpha=0.25,
    )
    mask_panel = bar(
        aperture_centers_m() .* 1e3,
        [state == 'S' ? 1.0 : 0.0 for state in optimized.mask];
        label=false, xlabel="radiator y, mm", ylabel="state (S=1)",
        title="Optimized context-aware binary mask", ylims=(0, 1.15),
    )
    figure_path = joinpath(OUTPUT_ROOT, "long_smooth_aperture_impulse.png")
    savefig(plot(
        waveform_panel, transverse_panel, axial_panel, mask_panel;
        layout=(2, 2), size=(1450, 950), margin=5Plots.mm,
    ), figure_path)

    println("[+] mask=$(String(optimized.mask))")
    println("[+] impulse gain=$(metrics.gain_peak), carrier gain=$carrier_gain")
    println("[+] Bt=$(metrics.broadening_ratio), rho=$(metrics.pulse_correlation), postcursor=$(metrics.postcursor_ratio)")
    println("[+] FWHM=$(transverse_width.width) mm, DOF=$axial_dof_mm mm, sidelobe=$side_ratio")
    println("[+] $summary_path")
    println("[+] $selection_path")
    println("[+] $context_path")
    println("[+] $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
