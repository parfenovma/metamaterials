module HornMonotonicV1CarrierAnalysis

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const V1_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_V1_FIELD_ROOT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_array_3d_v1"),
)
const V0_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_array_3d_half_symmetry")
const TRANSIENT_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_pilot_v1")
const V1_DESIGN_PATH = joinpath(
    PROJECT_ROOT,
    "results",
    "horn_monotonic_group_delay_193p8_290p1khz",
    "corrected_monotonic_lens_v1.jld2",
)
const V0_DESIGN_PATH = joinpath(
    PROJECT_ROOT,
    "tmp",
    "horn_monotonic_aperture",
    "monotonic_aperture_impulse.jld2",
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_V1_CARRIER_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_v1_carrier"),
)
const LEVEL_ID = "phase3"

using JLD2
using Statistics: mean
ENV["GKSwstype"] = "100"
using Plots

if !isdefined(parentmodule(@__MODULE__), :HornMonotonicResponseMatrix)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "run_horn_monotonic_response_matrix.jl"),
    )
end
using ..HornMonotonicResponseMatrix: profile_metrics
using ..MonotonicHornLens: geometric_delays_s

export wrap_phase_rad, output_phase_error_deg, analyze_v1_carrier

wrap_phase_rad(phase_rad::Real) = mod(Float64(phase_rad) + pi, 2pi) - pi
wrap_phase_rad(phase_rad) = wrap_phase_rad.(phase_rad)

function output_phase_error_deg(channel_output, frequency_hz, target_delay_s)
    length(channel_output) == length(target_delay_s) ||
        throw(DimensionMismatch("one target delay per channel"))
    relative_phase_rad = angle.(channel_output ./ first(channel_output))
    ideal_relative_phase_rad = -2pi * Float64(frequency_hz) .* target_delay_s
    rad2deg.(wrap_phase_rad(relative_phase_rad .- ideal_relative_phase_rad))
end

function full_profile(data)
    nonnegative_y_mm = Float64.(data["scan_y_mm"])
    nonnegative_ux_m = vec(data["scan_ux_m"])
    y_mm = vcat(-reverse(nonnegative_y_mm[2:end]), nonnegative_y_mm)
    ux_m = vcat(reverse(nonnegative_ux_m[2:end]), nonnegative_ux_m)
    profile_metrics(y_mm, ux_m)
end

function analyze_v1_carrier()
    mkpath(OUTPUT_ROOT)
    v1_path = joinpath(V1_ROOT, "field_half_lens_$(LEVEL_ID).jld2")
    v0_lens_path = joinpath(V0_ROOT, "field_half_lens_$(LEVEL_ID).jld2")
    uniform_path = joinpath(V0_ROOT, "field_half_uniform_$(LEVEL_ID).jld2")
    mesh_info_path = joinpath(V1_ROOT, "mesh_half_lens_$(LEVEL_ID).jld2")
    all(isfile, (v1_path, v0_lens_path, uniform_path, mesh_info_path)) ||
        error("v1/v0 carrier inputs are incomplete")

    v1 = JLD2.load(v1_path)
    v0_lens = JLD2.load(v0_lens_path)
    uniform = JLD2.load(uniform_path)
    mesh_info = JLD2.load(mesh_info_path)
    v1_design = JLD2.load(V1_DESIGN_PATH)
    v0_design = JLD2.load(V0_DESIGN_PATH)
    frequency_hz = Float64(v1["frequency_hz"])
    frequency_hz == uniform["frequency_hz"] == v0_lens["frequency_hz"] ||
        error("carrier fields use different frequencies")
    matched_reference_reused = isapprox(
        v1_design["common_axial_length_mm"],
        v0_design["common_axial_length_mm"];
        atol=1e-10,
        rtol=0,
    )
    matched_reference_reused || error("v1 straight reference is not identical to v0")

    v1_profile = full_profile(v1)
    v0_profile = full_profile(v0_lens)
    uniform_profile = full_profile(uniform)
    v1_focus_nm = abs(v1["focus_displacement_m"][1]) * 1e9
    v0_focus_nm = abs(v0_lens["focus_displacement_m"][1]) * 1e9
    uniform_focus_nm = abs(uniform["focus_displacement_m"][1]) * 1e9
    matched_gain = v1_focus_nm / uniform_focus_nm
    v0_matched_gain = v0_focus_nm / uniform_focus_nm
    focus_change = v1_focus_nm / v0_focus_nm - 1

    v1_channels = ComplexF64.(v1["channel_output_ux_m"])
    v0_channels = ComplexF64.(v0_lens["channel_output_ux_m"])
    target_delay_s = geometric_delays_s()
    v1_phase_error_deg = output_phase_error_deg(v1_channels, frequency_hz, target_delay_s)
    v0_phase_error_deg = output_phase_error_deg(v0_channels, frequency_hz, target_delay_s)
    central_channel_ratio = abs(v1_channels[8]) / mean(abs.(v1_channels[[1, 15]]))
    maximum_phase_error_deg = maximum(abs, v1_phase_error_deg)
    mirror_amplitude_mismatch = maximum(
        abs(abs(v1_channels[index]) - abs(v1_channels[16 - index])) /
        max(abs(v1_channels[index]), abs(v1_channels[16 - index]))
        for index in 1:7
    )
    mirror_phase_mismatch_deg = maximum(
        abs(rad2deg(wrap_phase_rad(
            angle(v1_channels[index]) - angle(v1_channels[16 - index]),
        )))
        for index in 1:7
    )

    gain_passed = matched_gain >= 2.0
    central_amplitude_passed = central_channel_ratio >= 0.85
    phase_passed = maximum_phase_error_deg <= 10.0
    fwhm_passed = v1_profile.fwhm_mm <= 6.0
    sidelobe_passed = v1_profile.sidelobe_amplitude_ratio <= 0.40
    mirror_amplitude_passed = mirror_amplitude_mismatch <= 0.03
    mirror_phase_passed = mirror_phase_mismatch_deg <= 3.0
    carrier_gate_passed = gain_passed && central_amplitude_passed && phase_passed &&
                          fwhm_passed && sidelobe_passed &&
                          mirror_amplitude_passed && mirror_phase_passed
    full_transient_authorized = carrier_gate_passed
    phase4_partial_count = count(
        name -> startswith(name, "field_half_lens_phase4.crash-partial"),
        readdir(V1_ROOT),
    )

    summary_path = joinpath(OUTPUT_ROOT, "v1_carrier_gate_summary.csv")
    open(summary_path, "w") do io
        println(io, "level_id,frequency_khz,path_mesh_size_mm,focus_mesh_size_mm,receiver_mesh_size_mm,lens_nodes,lens_elements,v1_focus_nm,v0_focus_nm,uniform_focus_nm,v1_matched_gain,v0_matched_gain,v1_focus_change_from_v0,v1_fwhm_mm,v1_sidelobe_ratio,central_channel_over_outer_edge,maximum_output_phase_error_deg,mirror_amplitude_mismatch,mirror_phase_mismatch_deg,gain_passed,central_amplitude_passed,phase_passed,fwhm_passed,sidelobe_passed,mirror_amplitude_passed,mirror_phase_passed,carrier_gate_passed,full_transient_authorized,matched_reference_reused,phase4_allocator_crash_count,path_h_over_lambda_s,focus_h_over_lambda_s,receiver_h_over_lambda_s")
        println(io, join((
            LEVEL_ID,
            frequency_hz / 1e3,
            mesh_info["size_path_mm"],
            mesh_info["size_focus_mm"],
            mesh_info["size_receiver_mm"],
            mesh_info["node_count"],
            mesh_info["element_count"],
            v1_focus_nm,
            v0_focus_nm,
            uniform_focus_nm,
            matched_gain,
            v0_matched_gain,
            focus_change,
            v1_profile.fwhm_mm,
            v1_profile.sidelobe_amplitude_ratio,
            central_channel_ratio,
            maximum_phase_error_deg,
            mirror_amplitude_mismatch,
            mirror_phase_mismatch_deg,
            gain_passed,
            central_amplitude_passed,
            phase_passed,
            fwhm_passed,
            sidelobe_passed,
            mirror_amplitude_passed,
            mirror_phase_passed,
            carrier_gate_passed,
            full_transient_authorized,
            matched_reference_reused,
            phase4_partial_count,
            mesh_info["path_h_over_lambda_s"],
            mesh_info["focus_h_over_lambda_s"],
            mesh_info["receiver_h_over_lambda_s"],
        ), ','))
    end

    channels_path = joinpath(OUTPUT_ROOT, "v1_carrier_channels.csv")
    centers_mm = collect(-7:7) .* 4.8
    open(channels_path, "w") do io
        println(io, "channel_index,center_y_mm,v1_amplitude_nm,v1_phase_deg,v1_phase_error_deg,v0_amplitude_nm,v0_phase_deg,v0_phase_error_deg,target_delay_us,center_y_over_pitch,target_delay_cycles")
        for index in eachindex(centers_mm)
            println(io, join((
                index,
                centers_mm[index],
                abs(v1_channels[index]) * 1e9,
                rad2deg(angle(v1_channels[index])),
                v1_phase_error_deg[index],
                abs(v0_channels[index]) * 1e9,
                rad2deg(angle(v0_channels[index])),
                v0_phase_error_deg[index],
                target_delay_s[index] * 1e6,
                centers_mm[index] / 4.8,
                frequency_hz * target_delay_s[index],
            ), ','))
        end
    end

    result_path = joinpath(OUTPUT_ROOT, "v1_carrier_gate.jld2")
    jldsave(
        result_path;
        format_version=1,
        level_id=LEVEL_ID,
        frequency_hz,
        v1_design_path=V1_DESIGN_PATH,
        v1_field_path=v1_path,
        v0_lens_path,
        matched_uniform_path=uniform_path,
        matched_reference_reused,
        v1_focus_nm,
        v0_focus_nm,
        uniform_focus_nm,
        matched_gain,
        v0_matched_gain,
        focus_change,
        v1_fwhm_mm=v1_profile.fwhm_mm,
        v1_sidelobe_ratio=v1_profile.sidelobe_amplitude_ratio,
        central_channel_ratio,
        maximum_phase_error_deg,
        mirror_amplitude_mismatch,
        mirror_phase_mismatch_deg,
        carrier_gate_passed,
        full_transient_authorized,
        phase4_partial_count,
    )

    profile_panel = plot(
        v1_profile.coordinate_mm,
        v1_profile.amplitude_nm;
        linewidth=2.7,
        xlabel="y at x=35 mm, mm",
        ylabel="|ux|, nm",
        title="Carrier focal profile",
        label="path-corrected v1",
        gridalpha=0.25,
    )
    plot!(profile_panel, v0_profile.coordinate_mm, v0_profile.amplitude_nm;
          linewidth=2.2, label="v0 lens")
    plot!(profile_panel, uniform_profile.coordinate_mm, uniform_profile.amplitude_nm;
          linewidth=2.0, label="matched straight")

    channel_panel = plot(
        centers_mm,
        abs.(v1_channels) .* 1e9;
        marker=:square,
        linewidth=2.5,
        xlabel="channel centre y, mm",
        ylabel="preout |ux|, nm",
        title="Channel amplitudes",
        label="v1",
        gridalpha=0.25,
    )
    plot!(channel_panel, centers_mm, abs.(v0_channels) .* 1e9;
          marker=:circle, linewidth=2.0, label="v0")

    phase_panel = plot(
        centers_mm,
        v1_phase_error_deg;
        marker=:square,
        linewidth=2.5,
        xlabel="channel centre y, mm",
        ylabel="output phase error, deg",
        title="Phase error relative to unwrapped target",
        label="v1",
        gridalpha=0.25,
    )
    plot!(phase_panel, centers_mm, v0_phase_error_deg;
          marker=:circle, linewidth=2.0, label="v0")
    hline!(phase_panel, [-10.0, 10.0]; linestyle=:dash, color=:black, label="±10 deg gate")

    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_v1_carrier_gate.png")
    savefig(plot(
        profile_panel,
        channel_panel,
        phase_panel;
        layout=(1, 3),
        size=(1800, 600),
        margin=5Plots.mm,
    ), figure_path)

    println("[+] v1 phase3 matched gain=$matched_gain (gate >=2.0)")
    println("[+] v1/v0 focus change=$(100focus_change)%")
    println("[+] FWHM=$(v1_profile.fwhm_mm) mm, sidelobe=$(v1_profile.sidelobe_amplitude_ratio)")
    println("[+] center/edge=$central_channel_ratio, max phase error=$maximum_phase_error_deg deg")
    println("[+] carrier gate passed=$carrier_gate_passed; full transient authorized=$full_transient_authorized")
    println("[+] $summary_path")
    println("[+] $figure_path")
    (; summary_path, channels_path, result_path, figure_path, carrier_gate_passed)
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_v1_carrier()
end

end
