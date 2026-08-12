module HornMonotonicArray3DConvergence

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_aperture")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_MONOTONIC_CONVERGENCE_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_array_3d_convergence"),
)

using JLD2
using Statistics: mean
using Gmsh: gmsh
ENV["GKSwstype"] = "100"
using Plots

if !isdefined(parentmodule(@__MODULE__), :HornPointRadiatorMesher)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
end
if !isdefined(parentmodule(@__MODULE__), :MonotonicHornLens)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "monotonic_horn_lens.jl"))
end
if !isdefined(parentmodule(@__MODULE__), :HornMonotonicArray3DMesher)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "horn_monotonic_array_3d_mesher.jl"),
    )
end
if !isdefined(parentmodule(@__MODULE__), :ImpulseRiskAnalysis)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "impulse_risk_analysis.jl"))
end
if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
end
if !isdefined(parentmodule(@__MODULE__), :DimensionlessWaveScaling)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "dimensionless_wave_scaling.jl"),
    )
end
if !isdefined(parentmodule(@__MODULE__), :HornMonotonicArray3DHarmonicSolver)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "horn_monotonic_array_3d_harmonic_solver.jl"),
    )
end

using ..MonotonicHornLens: smooth_slope
using ..HornMonotonicArray3DMesher
using ..ImpulseRiskAnalysis: contiguous_width
using ..SinusoidalMaterialLens: photopolymer
using ..DimensionlessWaveScaling:
    WaveScale,
    pressure_wavelength_mm,
    shear_wavelength_mm,
    length_over_pressure_wavelength,
    length_over_shear_wavelength,
    physical_shear_length_mm,
    mesh_resolution
using ..HornMonotonicArray3DHarmonicSolver:
    HornMonotonicArray3DHarmonicConfig,
    solve_horn_monotonic_array_3d_harmonic

export MeshLevel,
       LEVELS,
       CARRIER_WAVE_SCALE,
       CARRIER_LAMBDA_P_MM,
       CARRIER_LAMBDA_S_MM,
       parse_options,
       convergence_probe_points_mm,
       average_channel_outputs,
       mirror_pair_metrics,
       convergence_gate

const CARRIER_FREQUENCY_HZ = 242.0e3
const CARRIER_WAVE_SCALE = WaveScale(
    photopolymer();
    reference_frequency_hz=CARRIER_FREQUENCY_HZ,
)
const CARRIER_LAMBDA_P_MM = pressure_wavelength_mm(CARRIER_WAVE_SCALE)
const CARRIER_LAMBDA_S_MM = shear_wavelength_mm(CARRIER_WAVE_SCALE)

struct MeshLevel
    id::String
    path_h_over_lambda_s::Float64
    focus_h_over_lambda_s::Float64
    receiver_h_over_lambda_s::Float64
    size_path_mm::Float64
    size_focus_mm::Float64
    size_receiver_mm::Float64
end

function MeshLevel(;
    id::AbstractString,
    path_h_over_lambda_s::Real,
    focus_h_over_lambda_s::Real,
    receiver_h_over_lambda_s::Real,
    scale::WaveScale=CARRIER_WAVE_SCALE,
)
    ratios = Float64.(
        (path_h_over_lambda_s, focus_h_over_lambda_s, receiver_h_over_lambda_s),
    )
    all(>(0.0), ratios) || throw(ArgumentError("mesh wavelength ratios must be positive"))
    MeshLevel(
        String(id),
        ratios...,
        physical_shear_length_mm(ratios[1], scale),
        physical_shear_length_mm(ratios[2], scale),
        physical_shear_length_mm(ratios[3], scale),
    )
end

# The reference level is the prescription used for the first integrated v0.
# Coarse is a cheap sensitivity screen; the acceptance gate compares reference
# against refined, i.e. it never certifies v0 from a coarse-to-baseline pair.
const LEVELS = [
    MeshLevel(
        id="coarse",
        path_h_over_lambda_s=1.10 / CARRIER_LAMBDA_S_MM,
        focus_h_over_lambda_s=2.20 / CARRIER_LAMBDA_S_MM,
        receiver_h_over_lambda_s=3.20 / CARRIER_LAMBDA_S_MM,
    ),
    MeshLevel(
        id="reference",
        path_h_over_lambda_s=0.90 / CARRIER_LAMBDA_S_MM,
        focus_h_over_lambda_s=1.80 / CARRIER_LAMBDA_S_MM,
        receiver_h_over_lambda_s=2.60 / CARRIER_LAMBDA_S_MM,
    ),
    MeshLevel(
        id="refined",
        path_h_over_lambda_s=0.80 / CARRIER_LAMBDA_S_MM,
        focus_h_over_lambda_s=1.60 / CARRIER_LAMBDA_S_MM,
        receiver_h_over_lambda_s=2.30 / CARRIER_LAMBDA_S_MM,
    ),
]

function parse_options(args=ARGS)
    options = Dict{String, String}()
    for argument in args
        startswith(argument, "--") || throw(ArgumentError("unknown argument: $argument"))
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || throw(ArgumentError("expected --name=value: $argument"))
        options[fields[1]] = fields[2]
    end
    stage = get(options, "stage", "analyze")
    level = get(options, "level", "all")
    variant = get(options, "variant", "all")
    stage in ("mesh", "solve", "analyze") ||
        throw(ArgumentError("--stage must be mesh, solve, or analyze"))
    level in vcat(getproperty.(LEVELS, :id), ["all"]) ||
        throw(ArgumentError("--level must be coarse, reference, refined, or all"))
    variant in ("lens", "uniform", "all") ||
        throw(ArgumentError("--variant must be lens, uniform, or all"))
    stage == "analyze" && (level != "all" || variant != "all") &&
        throw(ArgumentError("analyze consumes every level and variant"))
    (; stage, level, variant)
end

level_by_id(id::AbstractString) = only(filter(level -> level.id == id, LEVELS))
selected_levels(id::AbstractString) = id == "all" ? LEVELS : [level_by_id(id)]
selected_variants(variant::AbstractString) =
    variant == "all" ? ["lens", "uniform"] : [variant]

function load_design()
    path = joinpath(DESIGN_ROOT, "monotonic_aperture_impulse.jld2")
    isfile(path) || error("monotonic aperture design not found: $path")
    JLD2.load(path)
end

function array_config(variant::AbstractString)
    variant in ("lens", "uniform") || throw(ArgumentError("unknown variant: $variant"))
    design = load_design()
    amplitudes_mm = Float64.(design["bend_amplitude_mm"])
    HornMonotonicArray3DConfig(
        guide_axial_length_mm=Float64(design["common_axial_length_mm"]),
        bend_amplitude_mm=variant == "lens" ? amplitudes_mm : zeros(length(amplitudes_mm)),
    )
end

mesh_path(level, variant) =
    joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant)_$(level.id).msh")
mesh_info_path(level, variant) =
    joinpath(OUTPUT_ROOT, "mesh_$(variant)_$(level.id).jld2")
result_path(level, variant) =
    joinpath(OUTPUT_ROOT, "field_$(variant)_$(level.id).jld2")

"""
Return nine interior samples per guide and one focus point.

The channel observable is a 3×3 arithmetic average across the local guide
cross-section. This is substantially less sensitive to tetrahedron placement
than the old single center-point probe while remaining inside the guide.
"""
function convergence_probe_points_mm(
    config::HornMonotonicArray3DConfig;
    channel_indices=eachindex(config.bend_amplitude_mm),
    nonnegative_half::Bool=false,
)
    names = String[]
    x_mm = Float64[]
    y_mm = Float64[]
    z_mm = Float64[]
    channel_groups = Vector{Vector{Int}}()
    local_x_mm = config.guide_axial_length_mm - 1.0
    center_x_mm = config.horn_length_mm + local_x_mm
    full_y_offsets_mm = (-0.45, 0.0, 0.45) .* (config.channel_depth_mm / 2)
    normal_offsets_mm = (-0.45, 0.0, 0.45) .* (config.throat_height_mm / 2)
    centers_mm = channel_centers_mm(config)
    for channel_index in channel_indices
        center_y_mm = centers_mm[channel_index]
        y_offsets_mm = if nonnegative_half && iszero(center_y_mm)
            (0.10, 0.45, 0.80) .* (config.channel_depth_mm / 2)
        else
            full_y_offsets_mm
        end
        center_z_mm = channel_center_z_mm(config, channel_index, local_x_mm)
        slope = smooth_slope(
            local_x_mm,
            config.guide_axial_length_mm,
            config.bend_amplitude_mm[channel_index],
        )
        normalization = hypot(1.0, slope)
        normal_x, normal_z = -slope / normalization, 1 / normalization
        group = Int[]
        for y_offset_mm in y_offsets_mm, normal_offset_mm in normal_offsets_mm
            push!(names, "preout_$(channel_index)_$(length(group) + 1)")
            push!(x_mm, center_x_mm + normal_offset_mm * normal_x)
            push!(y_mm, center_y_mm + y_offset_mm)
            push!(z_mm, center_z_mm + normal_offset_mm * normal_z)
            push!(group, length(names))
        end
        push!(channel_groups, group)
    end
    push!(names, "focus")
    push!(x_mm, focus_x_mm(config))
    push!(y_mm, 0.0)
    push!(z_mm, 0.0)
    (; names, x_mm, y_mm, z_mm, channel_groups, focus_index=length(names))
end

function average_channel_outputs(probe_displacement_m, channel_groups)
    size(probe_displacement_m, 2) == 3 ||
        throw(ArgumentError("probe displacement must have three components"))
    ComplexF64[
        mean(probe_displacement_m[group, 1])
        for group in channel_groups
    ]
end

function mirror_pair_metrics(channel_outputs)
    count = length(channel_outputs)
    isodd(count) || throw(ArgumentError("mirror metrics need an odd channel count"))
    half = (count - 1) ÷ 2
    left_index = collect(1:half)
    right_index = count .+ 1 .- left_index
    amplitude_mismatch = Float64[]
    phase_mismatch_deg = Float64[]
    for (left, right) in zip(left_index, right_index)
        left_value = channel_outputs[left]
        right_value = channel_outputs[right]
        denominator = (abs(left_value) + abs(right_value)) / 2
        push!(
            amplitude_mismatch,
            iszero(denominator) ? Inf : abs(abs(left_value) - abs(right_value)) / denominator,
        )
        push!(phase_mismatch_deg, abs(rad2deg(angle(left_value / right_value))))
    end
    (;
        left_index,
        right_index,
        amplitude_mismatch,
        phase_mismatch_deg,
        maximum_amplitude_mismatch=maximum(amplitude_mismatch),
        maximum_phase_mismatch_deg=maximum(phase_mismatch_deg),
    )
end

function convergence_gate(
    baseline_gain,
    refined_gain,
    baseline_fwhm_mm,
    refined_fwhm_mm,
    refined_amplitude_mismatch,
    refined_phase_mismatch_deg;
    scan_step_mm=1.0,
)
    relative_gain_change = abs(refined_gain - baseline_gain) / abs(refined_gain)
    fwhm_change_mm = abs(refined_fwhm_mm - baseline_fwhm_mm)
    gain_stable = relative_gain_change <= 0.05
    fwhm_stable = fwhm_change_mm <= scan_step_mm
    amplitude_symmetric = refined_amplitude_mismatch <= 0.03
    phase_symmetric = refined_phase_mismatch_deg <= 3.0
    (;
        relative_gain_change,
        fwhm_change_mm,
        gain_stable,
        fwhm_stable,
        amplitude_symmetric,
        phase_symmetric,
        passed=gain_stable && fwhm_stable && amplitude_symmetric && phase_symmetric,
    )
end

function run_mesh(level, variant)
    config = array_config(variant)
    result = build_horn_monotonic_array_3d_mesh(
        mesh_path(level, variant);
        config,
        size_path_mm=level.size_path_mm,
        size_focus_mm=level.size_focus_mm,
        size_receiver_mm=level.size_receiver_mm,
    )
    jldsave(
        mesh_info_path(level, variant);
        level_id=level.id,
        variant,
        size_path_mm=level.size_path_mm,
        size_focus_mm=level.size_focus_mm,
        size_receiver_mm=level.size_receiver_mm,
        path_h_over_lambda_s=level.path_h_over_lambda_s,
        focus_h_over_lambda_s=level.focus_h_over_lambda_s,
        receiver_h_over_lambda_s=level.receiver_h_over_lambda_s,
        path_k_s_h=2pi * level.path_h_over_lambda_s,
        focus_k_s_h=2pi * level.focus_h_over_lambda_s,
        receiver_k_s_h=2pi * level.receiver_h_over_lambda_s,
        carrier_lambda_p_mm=CARRIER_LAMBDA_P_MM,
        carrier_lambda_s_mm=CARRIER_LAMBDA_S_MM,
        node_count=result.node_count,
        element_count=result.element_count,
    )
end

function run_mesh_stage(levels, variants)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for level in levels, variant in variants
            println("[+] meshing $(level.id) $variant")
            run_mesh(level, variant)
        end
    finally
        gmsh.finalize()
    end
end

function run_solve(level, variant)
    config = array_config(variant)
    probes = convergence_probe_points_mm(config)
    solver_config = HornMonotonicArray3DHarmonicConfig(
        scan_x_min_mm=config.focal_distance_mm,
        scan_x_max_mm=config.focal_distance_mm,
        # FWHM is below 6 mm; ±12 mm keeps a generous background margin while
        # avoiding a known Gridap point-search partition at y=-18 mm on the
        # deliberately coarse receiver mesh.
        scan_y_half_width_mm=12.0,
        scan_step_mm=1.0,
    )
    result = solve_horn_monotonic_array_3d_harmonic(
        mesh_path(level, variant),
        length(config.bend_amplitude_mm),
        outlet_x_mm(config),
        probes,
        photopolymer();
        config=solver_config,
    )
    channel_output_ux_m = average_channel_outputs(
        result.probe_displacement_m,
        probes.channel_groups,
    )
    jldsave(
        result_path(level, variant);
        format_version=1,
        level_id=level.id,
        variant,
        frequency_hz=solver_config.frequency_hz,
        centers_mm=channel_centers_mm(config),
        channel_output_ux_m,
        focus_displacement_m=result.focus_displacement_m,
        scan_x_mm=result.scan_x_mm,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=result.scan_ux_m,
    )
    println(
        "[+] $(level.id) $variant focus |ux|=",
        abs(result.focus_displacement_m[1]) * 1e9,
        " nm",
    )
end

function run_solve_stage(levels, variants)
    for level in levels, variant in variants
        println("[+] solving $(level.id) $variant")
        run_solve(level, variant)
    end
end

function level_metrics(level)
    lens = JLD2.load(result_path(level, "lens"))
    uniform = JLD2.load(result_path(level, "uniform"))
    mesh_lens = JLD2.load(mesh_info_path(level, "lens"))
    mesh_uniform = JLD2.load(mesh_info_path(level, "uniform"))
    lens_focus_nm = abs(lens["focus_displacement_m"][1]) * 1e9
    uniform_focus_nm = abs(uniform["focus_displacement_m"][1]) * 1e9
    lens_profile_nm = abs.(vec(lens["scan_ux_m"])) .* 1e9
    uniform_profile_nm = abs.(vec(uniform["scan_ux_m"])) .* 1e9
    lens_width = contiguous_width(lens["scan_y_mm"], lens_profile_nm).width
    uniform_width = contiguous_width(uniform["scan_y_mm"], uniform_profile_nm).width
    lens_mirror = mirror_pair_metrics(lens["channel_output_ux_m"])
    uniform_mirror = mirror_pair_metrics(uniform["channel_output_ux_m"])
    (;
        level,
        lens_focus_nm,
        uniform_focus_nm,
        gain=lens_focus_nm / uniform_focus_nm,
        lens_fwhm_mm=lens_width,
        uniform_fwhm_mm=uniform_width,
        lens_mirror,
        uniform_mirror,
        lens_nodes=Int(mesh_lens["node_count"]),
        uniform_nodes=Int(mesh_uniform["node_count"]),
    )
end

function write_pair_rows(io, row, variant, metrics)
    for index in eachindex(metrics.left_index)
        println(io, join((
            row.level.id,
            variant,
            metrics.left_index[index],
            metrics.right_index[index],
            metrics.amplitude_mismatch[index],
            metrics.phase_mismatch_deg[index],
        ), ','))
    end
end

function write_scale_manifest()
    config = array_config("lens")
    design = load_design()
    rows = [
        ("carrier_frequency", CARRIER_FREQUENCY_HZ / 1e3, "kHz", 1.0, "f/f0"),
        ("polymer_P_wavelength", CARRIER_LAMBDA_P_MM, "mm", 1.0, "lambdaP(f0)"),
        ("polymer_S_wavelength", CARRIER_LAMBDA_S_MM, "mm", 1.0, "lambdaS(f0)"),
        (
            "throat_height",
            config.throat_height_mm,
            "mm",
            length_over_shear_wavelength(config.throat_height_mm, CARRIER_WAVE_SCALE),
            "L/lambdaS(f0)",
        ),
        (
            "channel_depth",
            config.channel_depth_mm,
            "mm",
            length_over_shear_wavelength(config.channel_depth_mm, CARRIER_WAVE_SCALE),
            "L/lambdaS(f0)",
        ),
        (
            "channel_pitch",
            config.channel_pitch_mm,
            "mm",
            length_over_shear_wavelength(config.channel_pitch_mm, CARRIER_WAVE_SCALE),
            "L/lambdaS(f0)",
        ),
        (
            "focal_distance",
            config.focal_distance_mm,
            "mm",
            length_over_pressure_wavelength(config.focal_distance_mm, CARRIER_WAVE_SCALE),
            "L/lambdaP(f0)",
        ),
        (
            "guide_axial_length",
            config.guide_axial_length_mm,
            "mm",
            length_over_pressure_wavelength(config.guide_axial_length_mm, CARRIER_WAVE_SCALE),
            "L/lambdaP(f0)",
        ),
    ]
    if haskey(design, "extra_path_mm")
        maximum_extra_path_mm = maximum(Float64.(design["extra_path_mm"]))
        push!(rows, (
            "maximum_extra_path",
            maximum_extra_path_mm,
            "mm",
            length_over_pressure_wavelength(maximum_extra_path_mm, CARRIER_WAVE_SCALE),
            "L/lambdaP(f0)",
        ))
    end
    path = joinpath(OUTPUT_ROOT, "physical_and_dimensionless_scales.csv")
    open(path, "w") do io
        println(io, "quantity,physical_value,physical_unit,dimensionless_value,normalization")
        for row in rows
            println(io, join(row, ','))
        end
    end
    path
end

function run_analysis_stage()
    rows = level_metrics.(LEVELS)
    coarse, reference, refined = rows
    gate = convergence_gate(
        reference.gain,
        refined.gain,
        reference.lens_fwhm_mm,
        refined.lens_fwhm_mm,
        refined.lens_mirror.maximum_amplitude_mismatch,
        refined.lens_mirror.maximum_phase_mismatch_deg,
    )
    mkpath(OUTPUT_ROOT)
    metrics_path = joinpath(OUTPUT_ROOT, "convergence_metrics.csv")
    open(metrics_path, "w") do io
        println(io, "level,size_path_mm,size_focus_mm,size_receiver_mm,path_h_over_lambda_s,focus_h_over_lambda_s,receiver_h_over_lambda_s,path_k_s_h,focus_k_s_h,receiver_k_s_h,lens_nodes,uniform_nodes,lens_focus_ux_nm,uniform_focus_ux_nm,focus_gain,lens_fwhm_mm,uniform_fwhm_mm,lens_max_pair_amplitude_mismatch,lens_max_pair_phase_mismatch_deg,uniform_max_pair_amplitude_mismatch,uniform_max_pair_phase_mismatch_deg")
        for row in rows
            println(io, join((
                row.level.id,
                row.level.size_path_mm,
                row.level.size_focus_mm,
                row.level.size_receiver_mm,
                row.level.path_h_over_lambda_s,
                row.level.focus_h_over_lambda_s,
                row.level.receiver_h_over_lambda_s,
                2pi * row.level.path_h_over_lambda_s,
                2pi * row.level.focus_h_over_lambda_s,
                2pi * row.level.receiver_h_over_lambda_s,
                row.lens_nodes,
                row.uniform_nodes,
                row.lens_focus_nm,
                row.uniform_focus_nm,
                row.gain,
                row.lens_fwhm_mm,
                row.uniform_fwhm_mm,
                row.lens_mirror.maximum_amplitude_mismatch,
                row.lens_mirror.maximum_phase_mismatch_deg,
                row.uniform_mirror.maximum_amplitude_mismatch,
                row.uniform_mirror.maximum_phase_mismatch_deg,
            ), ','))
        end
    end
    pairs_path = joinpath(OUTPUT_ROOT, "mirror_pair_metrics.csv")
    open(pairs_path, "w") do io
        println(io, "level,variant,left_index,right_index,amplitude_mismatch,phase_mismatch_deg")
        for row in rows
            write_pair_rows(io, row, "lens", row.lens_mirror)
            write_pair_rows(io, row, "uniform", row.uniform_mirror)
        end
    end
    summary_path = joinpath(OUTPUT_ROOT, "convergence_gate.csv")
    open(summary_path, "w") do io
        println(io, "baseline_gain,refined_gain,relative_gain_change,baseline_lens_fwhm_mm,refined_lens_fwhm_mm,fwhm_change_mm,refined_max_pair_amplitude_mismatch,refined_max_pair_phase_mismatch_deg,gain_stable,fwhm_stable,amplitude_symmetric,phase_symmetric,passed")
        println(io, join((
            reference.gain,
            refined.gain,
            gate.relative_gain_change,
            reference.lens_fwhm_mm,
            refined.lens_fwhm_mm,
            gate.fwhm_change_mm,
            refined.lens_mirror.maximum_amplitude_mismatch,
            refined.lens_mirror.maximum_phase_mismatch_deg,
            gate.gain_stable,
            gate.fwhm_stable,
            gate.amplitude_symmetric,
            gate.phase_symmetric,
            gate.passed,
        ), ','))
    end
    scales_path = write_scale_manifest()

    path_sizes = getproperty.(LEVELS, :size_path_mm)
    gain_panel = plot(
        path_sizes,
        getproperty.(rows, :gain);
        marker=:circle,
        linewidth=2.5,
        xlabel="path mesh size, mm",
        ylabel="matched focus gain",
        title="Matched gain convergence",
        label=false,
        xflip=true,
        gridalpha=0.25,
    )
    hline!(gain_panel, [2.0]; linestyle=:dash, color=:black, label="system target")
    focus_panel = plot(
        path_sizes,
        getproperty.(rows, :lens_focus_nm);
        marker=:circle,
        linewidth=2.5,
        label="lens",
        xlabel="path mesh size, mm",
        ylabel="|ux(target)|, nm",
        title="Absolute matched focus",
        xflip=true,
        gridalpha=0.25,
    )
    plot!(focus_panel, path_sizes, getproperty.(rows, :uniform_focus_nm);
          marker=:circle, linewidth=2.5, label="straight reference")
    amplitude_panel = plot(
        path_sizes,
        100 .* [row.lens_mirror.maximum_amplitude_mismatch for row in rows];
        marker=:circle,
        linewidth=2.5,
        label=false,
        xlabel="path mesh size, mm",
        ylabel="maximum pair mismatch, %",
        title="Lens mirror amplitude",
        xflip=true,
        gridalpha=0.25,
    )
    hline!(amplitude_panel, [3.0]; linestyle=:dash, color=:black, label="gate")
    phase_panel = plot(
        path_sizes,
        [row.lens_mirror.maximum_phase_mismatch_deg for row in rows];
        marker=:circle,
        linewidth=2.5,
        label=false,
        xlabel="path mesh size, mm",
        ylabel="maximum pair mismatch, deg",
        title="Lens mirror phase",
        xflip=true,
        gridalpha=0.25,
    )
    hline!(phase_panel, [3.0]; linestyle=:dash, color=:black, label="gate")
    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_array_3d_convergence.png")
    savefig(plot(
        gain_panel,
        focus_panel,
        amplitude_panel,
        phase_panel;
        layout=(2, 2),
        size=(1400, 950),
        margin=5Plots.mm,
    ), figure_path)
    println("[+] reference-to-refined matched-gain change=$(100gate.relative_gain_change)%")
    println("[+] lens FWHM change=$(gate.fwhm_change_mm) mm")
    println(
        "[+] refined mirror mismatch=",
        100refined.lens_mirror.maximum_amplitude_mismatch,
        "% / ",
        refined.lens_mirror.maximum_phase_mismatch_deg,
        " deg",
    )
    println("[+] P5.1A gate passed=$(gate.passed)")
    println("[+] $metrics_path")
    println("[+] $pairs_path")
    println("[+] $summary_path")
    println("[+] $scales_path")
    println("[+] $figure_path")
end

function main(args=ARGS)
    options = parse_options(args)
    levels = selected_levels(options.level)
    variants = selected_variants(options.variant)
    if options.stage == "mesh"
        run_mesh_stage(levels, variants)
    elseif options.stage == "solve"
        run_solve_stage(levels, variants)
    else
        run_analysis_stage()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
