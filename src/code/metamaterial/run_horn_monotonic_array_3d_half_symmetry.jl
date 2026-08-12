module HornMonotonicArray3DHalfSymmetry

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_aperture")
const DESIGN_PATH = get(
    ENV,
    "METAMATERIALS_HORN_MONOTONIC_DESIGN",
    joinpath(DESIGN_ROOT, "monotonic_aperture_impulse.jld2"),
)
const FULL_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_array_3d_convergence")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_MONOTONIC_HALF_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_array_3d_half_symmetry"),
)

using JLD2
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
if !isdefined(parentmodule(@__MODULE__), :HornMonotonicArray3DHarmonicSolver)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "horn_monotonic_array_3d_harmonic_solver.jl"),
    )
end
if !isdefined(parentmodule(@__MODULE__), :HornMonotonicArray3DConvergence)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "run_horn_monotonic_array_3d_convergence.jl"),
    )
end

using ..HornMonotonicArray3DMesher
using ..SinusoidalMaterialLens: photopolymer
using ..HornMonotonicArray3DHarmonicSolver:
    HornMonotonicArray3DHarmonicConfig,
    solve_horn_monotonic_array_3d_harmonic
using ..HornMonotonicArray3DConvergence:
    MeshLevel,
    LEVELS,
    CARRIER_LAMBDA_S_MM,
    convergence_probe_points_mm,
    average_channel_outputs,
    convergence_gate
using ..ImpulseRiskAnalysis: contiguous_width

export reconstruct_full_outputs, reconstruct_full_profile

const HALF_LEVELS = vcat(
    filter(level -> level.id in ("reference", "refined"), LEVELS),
    [MeshLevel(
        id="fine",
        path_h_over_lambda_s=0.70 / CARRIER_LAMBDA_S_MM,
        focus_h_over_lambda_s=1.40 / CARRIER_LAMBDA_S_MM,
        receiver_h_over_lambda_s=2.00 / CARRIER_LAMBDA_S_MM,
    ), MeshLevel(
        id="phase1",
        path_h_over_lambda_s=0.12,
        focus_h_over_lambda_s=0.20,
        receiver_h_over_lambda_s=0.30,
    ), MeshLevel(
        id="phase2",
        path_h_over_lambda_s=0.10,
        focus_h_over_lambda_s=0.16,
        receiver_h_over_lambda_s=0.25,
    ), MeshLevel(
        id="phase3",
        path_h_over_lambda_s=0.08,
        focus_h_over_lambda_s=0.16,
        receiver_h_over_lambda_s=0.25,
    ), MeshLevel(
        id="phase4",
        path_h_over_lambda_s=0.07,
        focus_h_over_lambda_s=0.16,
        receiver_h_over_lambda_s=0.25,
    )],
)

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
    level in ("reference", "refined", "fine", "phase1", "phase2", "phase3", "phase4", "all") ||
        throw(ArgumentError("--level must be reference, refined, fine, phase1, phase2, phase3, phase4, or all"))
    variant in ("lens", "uniform", "all") ||
        throw(ArgumentError("--variant must be lens, uniform, or all"))
    stage == "analyze" && (level != "all" || variant != "all") &&
        throw(ArgumentError("analyze consumes every half-domain level and variant"))
    (; stage, level, variant)
end

selected_levels(id) = id == "all" ? HALF_LEVELS : filter(level -> level.id == id, HALF_LEVELS)
selected_variants(id) = id == "all" ? ["lens", "uniform"] : [id]

function load_design()
    isfile(DESIGN_PATH) || error("monotonic aperture design not found: $DESIGN_PATH")
    JLD2.load(DESIGN_PATH)
end

function array_config(variant)
    design = load_design()
    amplitudes_mm = Float64.(design["bend_amplitude_mm"])
    HornMonotonicArray3DConfig(
        guide_axial_length_mm=Float64(design["common_axial_length_mm"]),
        bend_amplitude_mm=variant == "lens" ? amplitudes_mm : zeros(length(amplitudes_mm)),
    )
end

mesh_path(level, variant) =
    joinpath(OUTPUT_ROOT, "meshes", "mesh_half_$(variant)_$(level.id).msh")
mesh_info_path(level, variant) =
    joinpath(OUTPUT_ROOT, "mesh_half_$(variant)_$(level.id).jld2")
result_path(level, variant) =
    joinpath(OUTPUT_ROOT, "field_half_$(variant)_$(level.id).jld2")

function reconstruct_full_outputs(nonnegative_outputs)
    length(nonnegative_outputs) >= 2 ||
        throw(ArgumentError("need the center and at least one positive channel"))
    vcat(reverse(nonnegative_outputs[2:end]), nonnegative_outputs)
end

function reconstruct_full_profile(nonnegative_coordinate, nonnegative_profile)
    length(nonnegative_coordinate) == length(nonnegative_profile) ||
        throw(ArgumentError("coordinate and profile lengths differ"))
    first(nonnegative_coordinate) == 0 ||
        throw(ArgumentError("half-domain profile must start at y=0"))
    coordinate = vcat(-reverse(nonnegative_coordinate[2:end]), nonnegative_coordinate)
    profile = vcat(reverse(nonnegative_profile[2:end]), nonnegative_profile)
    (; coordinate, profile)
end

function run_mesh_stage(levels, variants)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for level in levels, variant in variants
            config = array_config(variant)
            result = build_horn_monotonic_array_3d_half_mesh(
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
                node_count=result.node_count,
                element_count=result.element_count,
                source_count=result.source_count,
                size_path_mm=level.size_path_mm,
                size_focus_mm=level.size_focus_mm,
                size_receiver_mm=level.size_receiver_mm,
                path_h_over_lambda_s=level.path_h_over_lambda_s,
                focus_h_over_lambda_s=level.focus_h_over_lambda_s,
                receiver_h_over_lambda_s=level.receiver_h_over_lambda_s,
                path_k_s_h=2pi * level.path_h_over_lambda_s,
                focus_k_s_h=2pi * level.focus_h_over_lambda_s,
                receiver_k_s_h=2pi * level.receiver_h_over_lambda_s,
            )
        end
    finally
        gmsh.finalize()
    end
end

function run_solve(level, variant)
    config = array_config(variant)
    positive_indices = 8:15
    probes = convergence_probe_points_mm(
        config;
        channel_indices=positive_indices,
        nonnegative_half=true,
    )
    solver_config = HornMonotonicArray3DHarmonicConfig(
        scan_x_min_mm=config.focal_distance_mm,
        scan_x_max_mm=config.focal_distance_mm,
        scan_y_half_width_mm=12.0,
        scan_y_min_mm=0.0,
        scan_step_mm=1.0,
        symmetry_y=true,
    )
    result = solve_horn_monotonic_array_3d_harmonic(
        mesh_path(level, variant),
        length(positive_indices),
        outlet_x_mm(config),
        probes,
        photopolymer();
        config=solver_config,
    )
    nonnegative_channel_output_ux_m = average_channel_outputs(
        result.probe_displacement_m,
        probes.channel_groups,
    )
    channel_output_ux_m = reconstruct_full_outputs(nonnegative_channel_output_ux_m)
    jldsave(
        result_path(level, variant);
        format_version=1,
        level_id=level.id,
        variant,
        frequency_hz=solver_config.frequency_hz,
        nonnegative_channel_indices=collect(positive_indices),
        nonnegative_channel_output_ux_m,
        channel_output_ux_m,
        focus_displacement_m=result.focus_displacement_m,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=result.scan_ux_m,
    )
    println(
        "[+] half $(level.id) $variant focus |ux|=",
        abs(result.focus_displacement_m[1]) * 1e9,
        " nm",
    )
end

function run_solve_stage(levels, variants)
    for level in levels, variant in variants
        run_solve(level, variant)
    end
end

function half_level_metrics(level)
    lens = JLD2.load(result_path(level, "lens"))
    uniform = JLD2.load(result_path(level, "uniform"))
    lens_focus_nm = abs(lens["focus_displacement_m"][1]) * 1e9
    uniform_focus_nm = abs(uniform["focus_displacement_m"][1]) * 1e9
    lens_half_profile_nm = abs.(vec(lens["scan_ux_m"])) .* 1e9
    uniform_half_profile_nm = abs.(vec(uniform["scan_ux_m"])) .* 1e9
    lens_profile = reconstruct_full_profile(lens["scan_y_mm"], lens_half_profile_nm)
    uniform_profile = reconstruct_full_profile(uniform["scan_y_mm"], uniform_half_profile_nm)
    focus_window = abs.(lens_profile.coordinate) .<= 3.0
    lens_peak_nm = maximum(lens_profile.profile)
    uniform_peak_nm = maximum(uniform_profile.profile)
    focus_window_energy_gain = sum(abs2, lens_profile.profile[focus_window]) /
                               sum(abs2, uniform_profile.profile[focus_window])
    (;
        level,
        lens_focus_nm,
        uniform_focus_nm,
        gain=lens_focus_nm / uniform_focus_nm,
        plane_peak_gain=lens_peak_nm / uniform_peak_nm,
        focus_window_amplitude_gain=sqrt(focus_window_energy_gain),
        focus_window_energy_gain,
        lens_fwhm_mm=contiguous_width(lens_profile.coordinate, lens_profile.profile).width,
        uniform_fwhm_mm=contiguous_width(uniform_profile.coordinate, uniform_profile.profile).width,
        lens_profile,
        uniform_profile,
    )
end

function run_analysis_stage()
    rows = half_level_metrics.(HALF_LEVELS)
    reference, refined, fine, phase1, phase2, phase3, phase4 = rows
    gate = convergence_gate(
        phase3.gain,
        phase4.gain,
        phase3.lens_fwhm_mm,
        phase4.lens_fwhm_mm,
        0.0,
        0.0,
    )
    full_refined = JLD2.load(joinpath(FULL_ROOT, "field_lens_refined.jld2"))
    full_uniform_refined = JLD2.load(joinpath(FULL_ROOT, "field_uniform_refined.jld2"))
    full_refined_gain = abs(full_refined["focus_displacement_m"][1]) /
                        abs(full_uniform_refined["focus_displacement_m"][1])
    half_to_full_gain_difference = abs(refined.gain - full_refined_gain) / abs(refined.gain)

    mkpath(OUTPUT_ROOT)
    metrics_path = joinpath(OUTPUT_ROOT, "half_symmetry_convergence.csv")
    open(metrics_path, "w") do io
        println(io, "level,size_path_mm,size_focus_mm,size_receiver_mm,path_h_over_lambda_s,focus_h_over_lambda_s,receiver_h_over_lambda_s,path_k_s_h,focus_k_s_h,receiver_k_s_h,lens_focus_ux_nm,uniform_focus_ux_nm,target_focus_gain,plane_peak_gain,focus_window_amplitude_gain,focus_window_energy_gain,lens_fwhm_mm,uniform_fwhm_mm")
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
                row.lens_focus_nm,
                row.uniform_focus_nm,
                row.gain,
                row.plane_peak_gain,
                row.focus_window_amplitude_gain,
                row.focus_window_energy_gain,
                row.lens_fwhm_mm,
                row.uniform_fwhm_mm,
            ), ','))
        end
    end
    gate_path = joinpath(OUTPUT_ROOT, "half_symmetry_gate.csv")
    open(gate_path, "w") do io
        println(io, "baseline_gain,phase_refined_gain,relative_gain_change,baseline_lens_fwhm_mm,phase_refined_lens_fwhm_mm,fwhm_change_mm,full_refined_gain,half_refined_to_full_refined_gain_difference,gain_stable,fwhm_stable,symmetry_exact,passed")
        println(io, join((
            phase3.gain,
            phase4.gain,
            gate.relative_gain_change,
            phase3.lens_fwhm_mm,
            phase4.lens_fwhm_mm,
            gate.fwhm_change_mm,
            full_refined_gain,
            half_to_full_gain_difference,
            gate.gain_stable,
            gate.fwhm_stable,
            true,
            gate.passed,
        ), ','))
    end

    gain_panel = plot(
        getproperty.(HALF_LEVELS, :size_path_mm),
        getproperty.(rows, :gain);
        marker=:circle,
        linewidth=2.5,
        label="half-domain",
        xlabel="path mesh size, mm",
        ylabel="matched focus gain",
        title="Symmetry-constrained convergence",
        xflip=true,
        gridalpha=0.25,
    )
    hline!(gain_panel, [full_refined_gain]; linestyle=:dash, label="full refined")
    profile_panel = plot(
        phase4.lens_profile.coordinate,
        phase4.lens_profile.profile;
        linewidth=2.7,
        label="half lens",
        xlabel="y at x=35 mm, mm",
        ylabel="|ux|, nm",
        title="Reconstructed phase-refined focus",
        gridalpha=0.25,
    )
    plot!(
        profile_panel,
        phase4.uniform_profile.coordinate,
        phase4.uniform_profile.profile;
        linewidth=2.2,
        label="half straight reference",
    )
    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_array_3d_half_symmetry.png")
    savefig(plot(
        gain_panel,
        profile_panel;
        layout=(1, 2),
        size=(1400, 600),
        margin=5Plots.mm,
    ), figure_path)
    println("[+] half phase3-to-phase4 gain change=$(100gate.relative_gain_change)%")
    println("[+] half FWHM change=$(gate.fwhm_change_mm) mm")
    println("[+] half phase4 gain=$(phase4.gain), full refined gain=$full_refined_gain")
    println("[+] half-domain P5.1A gate passed=$(gate.passed)")
    println("[+] $metrics_path")
    println("[+] $gate_path")
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
