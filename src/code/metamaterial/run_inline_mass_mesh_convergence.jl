module InlineMassMeshConvergence

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_INLINE_MASS_CONVERGENCE_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "inline_mass_mesh_convergence"),
)
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "inline_mass_component_mesher.jl"))
using .InlineMassComponentMesher

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif REQUESTED_STAGE == "solve" || REQUESTED_STAGE == "analyze"
    include(joinpath(@__DIR__, "modal_solver.jl"))
    using .ConservativeElasticModes
    using JLD2
end

if REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using Plots
end

const CONFIG = InlineMassConfig(
    mass_length_mm=1.4,
    mass_height_mm=2.4,
    neck_length_mm=0.35,
    neck_height_mm=0.95,
)

const LEVELS = [
    (id="coarse", size_min_mm=0.070, size_max_mm=0.220),
    (id="reference", size_min_mm=0.045, size_max_mm=0.150),
    (id="fine", size_min_mm=0.030, size_max_mm=0.100),
]

mesh_path(level) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(level.id).msh")
model_path(level) = joinpath(OUTPUT_ROOT, "models", "model_$(level.id).json")
result_path(level) = joinpath(OUTPUT_ROOT, "mode_$(level.id).jld2")

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for level in LEVELS
            build_inline_mass_component_mesh(
                mesh_path(level);
                config=CONFIG,
                size_min_mm=level.size_min_mm,
                size_max_mm=level.size_max_mm,
            )
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    for level in LEVELS
        convert_mesh(mesh_path(level); output_dir=dirname(model_path(level)))
    end
end

function solve_level(level)
    modes = solve_conservative_modes(
        model_path(level);
        config=ModeConfig(
            element_order=2,
            quadrature_degree=4,
            target_frequency_hz=250.0e3,
            mode_count=8,
            tolerance=1.0e-9,
            clamp_tags=["FixedInterface"],
        ),
    )
    metrics = [inline_mass_mode_metrics(mode, CONFIG) for mode in modes]
    selected_index = select_inline_mass_mode(modes, metrics)
    JLD2.jldsave(result_path(level); modes, metrics, selected_index, level)
    println("[+] $(level.id): $(modes[selected_index].frequency_hz / 1e3) kHz")
end

function run_solve_stage()
    mkpath(OUTPUT_ROOT)
    for level in LEVELS
        solve_level(level)
    end
end

function run_analysis_stage()
    rows = map(LEVELS) do level
        data = JLD2.load(result_path(level))
        index = data["selected_index"]
        mode = data["modes"][index]
        metric = data["metrics"][index]
        (; level, mode, metric)
    end
    fine_frequency = last(rows).mode.frequency_hz
    csv_path = joinpath(OUTPUT_ROOT, "inline_mass_mesh_convergence.csv")
    open(csv_path, "w") do io
        println(io, "level,size_min_mm,size_max_mm,frequency_hz,relative_to_fine,mass_fraction,translation_coherence,longitudinal_fraction,parity_y,residual")
        for row in rows
            println(io, join((
                row.level.id,
                row.level.size_min_mm,
                row.level.size_max_mm,
                row.mode.frequency_hz,
                (row.mode.frequency_hz - fine_frequency) / fine_frequency,
                row.metric.mass_fraction,
                row.metric.translation_coherence,
                row.mode.longitudinal_fraction,
                row.mode.parity_y,
                row.mode.relative_residual,
            ), ','))
        end
    end
    panel = plot(
        getproperty.(getproperty.(rows, :level), :size_max_mm),
        getproperty.(getproperty.(rows, :mode), :frequency_hz) ./ 1e3;
        marker=:circle,
        linewidth=2,
        xlabel="maximum mesh size, mm",
        ylabel="selected f_M, kHz",
        title="Mesh convergence of tuned fixed-interface M",
        label=false,
        gridalpha=0.25,
        xflip=true,
    )
    hline!(panel, [242.0]; color=:black, linestyle=:dash, label="242 kHz")
    png_path = joinpath(OUTPUT_ROOT, "inline_mass_mesh_convergence.png")
    savefig(panel, png_path)
    println("[+] $csv_path")
    println("[+] $png_path")
    for row in rows
        println(
            "  $(row.level.id): ", round(row.mode.frequency_hz / 1e3; digits=4),
            " kHz, relative to fine=",
            round(100 * (row.mode.frequency_hz - fine_frequency) / fine_frequency; digits=4),
            "%",
        )
    end
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_solve_stage()
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown convergence stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Inline mass convergence stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_inline_mass_mesh_convergence.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
