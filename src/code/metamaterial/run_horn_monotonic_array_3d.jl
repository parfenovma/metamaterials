module HornMonotonicArray3D

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_aperture")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_MONOTONIC_ARRAY_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_array_3d"),
)

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "report")
const VARIANT = argument("variant", "lens")
VARIANT in ("lens", "uniform") || error("--variant must be lens or uniform")

using JLD2
include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
include(joinpath(@__DIR__, "monotonic_horn_lens.jl"))
include(joinpath(@__DIR__, "horn_monotonic_array_3d_mesher.jl"))
using .HornMonotonicArray3DMesher

if STAGE == "mesh"
    using Gmsh: gmsh
elseif STAGE == "solve"
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
    include(joinpath(@__DIR__, "horn_monotonic_array_3d_harmonic_solver.jl"))
    using .SinusoidalMaterialLens
    using .HornMonotonicArray3DHarmonicSolver
elseif STAGE in ("report", "compare")
    ENV["GKSwstype"] = "100"
    using Plots
    using Statistics: mean
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    using .ImpulseRiskAnalysis: contiguous_width
end

const DESIGN = JLD2.load(joinpath(DESIGN_ROOT, "monotonic_aperture_impulse.jld2"))
const DESIGN_AMPLITUDES_MM = Float64.(DESIGN["bend_amplitude_mm"])
const CONFIG = HornMonotonicArray3DConfig(
    guide_axial_length_mm=Float64(DESIGN["common_axial_length_mm"]),
    bend_amplitude_mm=VARIANT == "lens" ? DESIGN_AMPLITUDES_MM : zeros(length(DESIGN_AMPLITUDES_MM)),
)

mesh_path(variant=VARIANT) = joinpath(
    OUTPUT_ROOT,
    variant == "lens" ? "mesh_horn_monotonic_array_3d.msh" :
    "mesh_horn_uniform_array_3d.msh",
)
result_path(variant=VARIANT) = joinpath(
    OUTPUT_ROOT,
    variant == "lens" ? "horn_monotonic_array_3d_harmonic.jld2" :
    "horn_uniform_array_3d_harmonic.jld2",
)

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        result = build_horn_monotonic_array_3d_mesh(
            mesh_path();
            config=CONFIG,
            size_path_mm=0.90,
            size_focus_mm=1.80,
            size_receiver_mm=2.60,
        )
        mkpath(OUTPUT_ROOT)
        open(joinpath(OUTPUT_ROOT, "mesh_summary_$(VARIANT).csv"), "w") do io
            println(io, "channels,nodes,elements,outlet_x_mm,focus_x_mm,receiver_half_width_mm,receiver_half_height_mm")
            println(io, join((
                length(CONFIG.bend_amplitude_mm),
                result.node_count,
                result.element_count,
                outlet_x_mm(CONFIG),
                focus_x_mm(CONFIG),
                CONFIG.receiver_half_width_mm,
                CONFIG.receiver_half_height_mm,
            ), ','))
        end
    finally
        gmsh.finalize()
    end
end

function run_solve_stage()
    probes = if VARIANT == "uniform"
        (names=["focus"], x_mm=[focus_x_mm(CONFIG)], y_mm=[0.0], z_mm=[0.0])
    else
        probe_points_mm(CONFIG)
    end
    result = solve_horn_monotonic_array_3d_harmonic(
        mesh_path(),
        length(CONFIG.bend_amplitude_mm),
        outlet_x_mm(CONFIG),
        probes,
        photopolymer();
        config=HornMonotonicArray3DHarmonicConfig(),
    )
    jldsave(
        result_path();
        format_version=1,
        variant=VARIANT,
        frequency_hz=242.0e3,
        channel_count=length(CONFIG.bend_amplitude_mm),
        outlet_x_mm=outlet_x_mm(CONFIG),
        focal_distance_mm=CONFIG.focal_distance_mm,
        centers_mm=channel_centers_mm(CONFIG),
        bend_amplitude_mm=CONFIG.bend_amplitude_mm,
        probe_names=result.probe_names,
        probe_x_mm=result.probe_x_mm,
        probe_y_mm=result.probe_y_mm,
        probe_z_mm=result.probe_z_mm,
        probe_displacement_m=result.probe_displacement_m,
        focus_displacement_m=result.focus_displacement_m,
        scan_x_mm=result.scan_x_mm,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=result.scan_ux_m,
        scan_uy_m=result.scan_uy_m,
        scan_uz_m=result.scan_uz_m,
        scan_total_amplitude_m=result.scan_total_amplitude_m,
    )
    println("[+] full-array focus |ux|=$(abs(result.focus_displacement_m[1])) m")
    println("[+] $(result_path())")
end

function run_compare_stage()
    lens = JLD2.load(result_path("lens"))
    uniform = JLD2.load(result_path("uniform"))
    lens_focus_nm = abs(lens["focus_displacement_m"][1]) * 1e9
    uniform_focus_nm = abs(uniform["focus_displacement_m"][1]) * 1e9
    focus_gain = lens_focus_nm / uniform_focus_nm
    x_mm = lens["scan_x_mm"]
    y_mm = lens["scan_y_mm"]
    x_mm == uniform["scan_x_mm"] || error("lens and uniform x scans differ")
    y_mm == uniform["scan_y_mm"] || error("lens and uniform y scans differ")
    target_x = argmin(abs.(x_mm .- lens["focal_distance_mm"]))
    lens_map_nm = abs.(lens["scan_ux_m"]) .* 1e9
    uniform_map_nm = abs.(uniform["scan_ux_m"]) .* 1e9
    lens_transverse = lens_map_nm[:, target_x]
    uniform_transverse = uniform_map_nm[:, target_x]
    lens_width = contiguous_width(y_mm, lens_transverse)
    uniform_width = contiguous_width(y_mm, uniform_transverse)
    common_limit = max(maximum(lens_map_nm), maximum(uniform_map_nm))

    summary_path = joinpath(OUTPUT_ROOT, "horn_monotonic_array_3d_gain_summary.csv")
    open(summary_path, "w") do io
        println(io, "lens_focus_ux_nm,uniform_focus_ux_nm,focus_gain,lens_fwhm_mm,uniform_fwhm_mm")
        println(io, join((
            lens_focus_nm,
            uniform_focus_nm,
            focus_gain,
            lens_width.width,
            uniform_width.width,
        ), ','))
    end

    lens_panel = heatmap(
        x_mm,
        y_mm,
        lens_map_nm;
        xlabel="distance from radiator plane, mm",
        ylabel="y, mm",
        title="Monotonic 15-channel array",
        color=:viridis,
        colorbar_title="|ux|, nm",
        clims=(0.0, common_limit),
        aspect_ratio=:equal,
        grid=false,
    )
    uniform_panel = heatmap(
        x_mm,
        y_mm,
        uniform_map_nm;
        xlabel="distance from radiator plane, mm",
        ylabel="y, mm",
        title="Uniform 15-channel reference",
        color=:viridis,
        colorbar_title="|ux|, nm",
        clims=(0.0, common_limit),
        aspect_ratio=:equal,
        grid=false,
    )
    for panel in (lens_panel, uniform_panel)
        scatter!(panel, [lens["focal_distance_mm"]], [0.0];
                 marker=:xcross, color=:white, markerstrokewidth=2, label=false)
    end
    profile_panel = plot(
        y_mm,
        lens_transverse;
        linewidth=2.7,
        label="monotonic lens",
        xlabel="y at x=35 mm, mm",
        ylabel="|ux|, nm",
        title="Matched full-3D focus comparison",
        gridalpha=0.25,
    )
    plot!(profile_panel, y_mm, uniform_transverse;
          linewidth=2.2, label="uniform reference")
    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_array_3d_gain.png")
    savefig(plot(
        lens_panel,
        uniform_panel,
        profile_panel;
        layout=(2, 2),
        size=(1450, 950),
        margin=5Plots.mm,
    ), figure_path)
    println("[+] matched full-3D gain=$focus_gain")
    println("[+] lens=$lens_focus_nm nm, uniform=$uniform_focus_nm nm")
    println("[+] lens FWHM=$(lens_width.width) mm, uniform FWHM=$(uniform_width.width) mm")
    println("[+] $summary_path")
    println("[+] $figure_path")
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
    outside = vcat(
        collect(firstindex(profile):(left_minimum - 1)),
        collect((right_minimum + 1):lastindex(profile)),
    )
    isempty(outside) ? 0.0 : maximum(profile[outside]) / profile[center_peak]
end

function run_report_stage()
    data = JLD2.load(result_path())
    x_mm = data["scan_x_mm"]
    y_mm = data["scan_y_mm"]
    ux_nm = abs.(data["scan_ux_m"]) .* 1e9
    total_nm = data["scan_total_amplitude_m"] .* 1e9
    target_x = argmin(abs.(x_mm .- data["focal_distance_mm"]))
    axis_y = argmin(abs.(y_mm))
    transverse = ux_nm[:, target_x]
    axial = ux_nm[axis_y, :]
    width = contiguous_width(y_mm, transverse)
    search_x = findall((x_mm .>= 25.0) .& (x_mm .<= 45.0))
    search_y = findall(abs.(y_mm) .<= 8.0)
    local_submap = ux_nm[search_y, search_x]
    local_linear = argmax(local_submap)
    local_y_sub, local_x_sub = Tuple(local_linear)
    peak_x_mm = x_mm[search_x[local_x_sub]]
    peak_y_mm = y_mm[search_y[local_y_sub]]
    focus_ux_nm = abs(data["focus_displacement_m"][1]) * 1e9
    outer = abs.(y_mm) .>= 12.0
    center_to_outer = transverse[axis_y] / mean(transverse[outer])
    concentration = [ux_nm[axis_y, index] / mean(ux_nm[outer, index]) for index in eachindex(x_mm)]
    widths_mm = [contiguous_width(y_mm, ux_nm[:, index]).width for index in eachindex(x_mm)]
    minimum_width = minimum(widths_mm[search_x])
    waist_candidates = [index for index in search_x if widths_mm[index] == minimum_width]
    waist_index = waist_candidates[argmax(concentration[waist_candidates])]
    waist_x_mm = x_mm[waist_index]
    focus_exists = width.width <= 10.0 && center_to_outer > 1.2
    side_ratio = sidelobe_ratio(transverse)

    preout_indices = [
        findfirst(==("preout_$index"), data["probe_names"])
        for index in 1:data["channel_count"]
    ]
    preout_ux = ComplexF64[data["probe_displacement_m"][index, 1] for index in preout_indices]
    edge_reference = (first(preout_ux) + last(preout_ux)) / 2
    actual_relative = preout_ux ./ edge_reference
    expected_weights = ComplexF64.(DESIGN["weights"])
    expected_relative = expected_weights ./ ((first(expected_weights) + last(expected_weights)) / 2)
    phase_error_deg = rad2deg.(angle.(actual_relative ./ expected_relative))
    frequency_hz = Float64(data["frequency_hz"])
    speed_m_s = 2340.0
    centers_m = data["centers_mm"] .* 1e-3
    function phase_fit_error(focal_distance_mm)
        focal_distance_m = focal_distance_mm * 1e-3
        distance_m = hypot.(focal_distance_m, centers_m)
        extra_distance_m = maximum(distance_m) .- distance_m
        ideal = cis.(-2pi * frequency_hz .* extra_distance_m ./ speed_m_s)
        unit_actual = actual_relative ./ abs.(actual_relative)
        weights = abs.(actual_relative) ./ maximum(abs.(actual_relative))
        sum(weights .* abs2.(unit_actual .- ideal)) / sum(weights)
    end
    focal_grid_mm = collect(20.0:0.05:60.0)
    phase_fit_errors = phase_fit_error.(focal_grid_mm)
    phase_fit_focal_mm = focal_grid_mm[argmin(phase_fit_errors)]

    channel_path = joinpath(OUTPUT_ROOT, "horn_monotonic_array_3d_channel_outputs.csv")
    open(channel_path, "w") do io
        println(io, "element_index,center_y_mm,actual_amplitude_over_edge,actual_phase_deg_over_edge,expected_amplitude_over_edge,expected_phase_deg_over_edge,phase_error_deg")
        for index in eachindex(actual_relative)
            println(io, join((
                index,
                data["centers_mm"][index],
                abs(actual_relative[index]),
                rad2deg(angle(actual_relative[index])),
                abs(expected_relative[index]),
                rad2deg(angle(expected_relative[index])),
                phase_error_deg[index],
            ), ','))
        end
    end

    summary_path = joinpath(OUTPUT_ROOT, "horn_monotonic_array_3d_summary.csv")
    open(summary_path, "w") do io
        println(io, "focus_exists,focus_ux_nm,local_peak_x_mm,local_peak_y_mm,target_transverse_fwhm_mm,beam_waist_x_mm,beam_waist_fwhm_mm,phase_fit_focal_mm,center_to_outer_amplitude,sidelobe_amplitude_ratio,maximum_output_phase_error_deg")
        println(io, join((
            focus_exists,
            focus_ux_nm,
            peak_x_mm,
            peak_y_mm,
            width.width,
            waist_x_mm,
            minimum_width,
            phase_fit_focal_mm,
            center_to_outer,
            side_ratio,
            maximum(abs, phase_error_deg),
        ), ','))
    end

    map_panel = heatmap(
        x_mm,
        y_mm,
        ux_nm;
        xlabel="distance from radiator plane, mm",
        ylabel="y, mm",
        title="Full 15-channel 3D lens: |ux| at z=0",
        color=:viridis,
        colorbar_title="|ux|, nm",
        aspect_ratio=:equal,
        grid=false,
    )
    scatter!(map_panel, [data["focal_distance_mm"]], [0.0];
             marker=:xcross, color=:white, markerstrokewidth=2, label="target")
    scatter!(map_panel, [peak_x_mm], [peak_y_mm];
             marker=:circle, color=:red, markerstrokewidth=1.5, label="local peak")
    profile_panel = plot(
        y_mm,
        transverse;
        linewidth=2.5,
        label="x=$(data["focal_distance_mm"]) mm",
        xlabel="y, mm",
        ylabel="|ux|, nm",
        title="Full-array transverse focus",
        gridalpha=0.25,
    )
    axial_panel = plot(
        x_mm,
        axial;
        linewidth=2.5,
        label="axis",
        xlabel="distance from radiator plane, mm",
        ylabel="|ux|, nm",
        title="Full-array axial profile",
        gridalpha=0.25,
    )
    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_array_3d_focus.png")
    savefig(plot(
        map_panel,
        profile_panel,
        axial_panel;
        layout=(2, 2),
        size=(1500, 900),
        margin=5Plots.mm,
    ), figure_path)
    println("[+] focus exists=$focus_exists, |ux(target)|=$focus_ux_nm nm")
    println("[+] local peak=($peak_x_mm, $peak_y_mm) mm, target FWHM=$(width.width) mm")
    println("[+] beam waist x=$waist_x_mm mm, width=$minimum_width mm, phase-fit F=$phase_fit_focal_mm mm")
    println("[+] center/outer=$center_to_outer, sidelobe=$side_ratio")
    println("[+] $summary_path")
    println("[+] $channel_path")
    println("[+] $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "solve" ? run_solve_stage() :
    STAGE == "report" ? run_report_stage() :
    STAGE == "compare" ? run_compare_stage() :
    error("unknown stage: $STAGE")
end

end
