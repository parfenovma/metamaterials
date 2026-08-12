module StiffnessMeshConvergence

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_STIFFNESS_CONVERGENCE_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "stiffness_mesh_convergence"),
)
const TARGET_FREQUENCY_HZ = 242.0e3
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "stiffness_component_mesher.jl"))
using .StiffnessComponentMesher

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif REQUESTED_STAGE == "solve" || REQUESTED_STAGE == "analyze"
    include(joinpath(@__DIR__, "modal_solver.jl"))
    include(joinpath(@__DIR__, "dynamic_stiffness_solver.jl"))
    using .ConservativeElasticModes
    using .DynamicStiffnessSolver
    using JLD2
end

if REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using Plots
end

const CONFIG = StiffnessConfig(ligament_height_mm=0.60)
const FREQUENCIES_HZ = [0.0, 222.0e3, TARGET_FREQUENCY_HZ, 262.0e3]
const LEVELS = [
    (id="coarse", size_min_mm=0.070, size_max_mm=0.220),
    (id="reference", size_min_mm=0.045, size_max_mm=0.150),
    (id="fine", size_min_mm=0.030, size_max_mm=0.100),
]

mesh_path(level) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(level.id).msh")
model_path(level) = joinpath(OUTPUT_ROOT, "models", "model_$(level.id).json")
result_path(level) = joinpath(OUTPUT_ROOT, "result_$(level.id).jld2")

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for level in LEVELS
            build_stiffness_component_mesh(
                mesh_path(level);
                config=CONFIG,
                variant=:slotted,
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

function run_solve_stage()
    mkpath(OUTPUT_ROOT)
    for level in LEVELS
        modes = solve_conservative_modes(
            model_path(level);
            config=ModeConfig(
                element_order=2,
                quadrature_degree=4,
                target_frequency_hz=TARGET_FREQUENCY_HZ,
                mode_count=12,
                tolerance=1.0e-9,
                clamp_tags=["FixedInterface"],
            ),
        )
        points = solve_dynamic_stiffness_sweep(
            model_path(level),
            FREQUENCIES_HZ;
            component_length_m=CONFIG.length_mm * 1e-3,
            config=DynamicStiffnessConfig(),
        )
        symmetric = filter(mode -> mode.parity_y >= 0.8, modes)
        nearest_symmetric = symmetric[
            argmin(abs(mode.frequency_hz - TARGET_FREQUENCY_HZ) for mode in symmetric)
        ]
        JLD2.jldsave(result_path(level); modes, points, nearest_symmetric, level)
        println(
            "[+] $(level.id): f_sym=",
            round(nearest_symmetric.frequency_hz / 1e3; digits=3), " kHz",
        )
    end
end

function point_at(points, frequency_hz)
    only(filter(point -> point.frequency_hz == frequency_hz, points))
end

function run_analysis_stage()
    rows = map(LEVELS) do level
        data = JLD2.load(result_path(level))
        points = data["points"]
        (
            level,
            static_stiffness=real(point_at(points, 0.0).stiffness_n_per_m2),
            target_stiffness=real(point_at(points, TARGET_FREQUENCY_HZ).stiffness_n_per_m2),
            band_low_stiffness=real(point_at(points, 222.0e3).stiffness_n_per_m2),
            band_high_stiffness=real(point_at(points, 262.0e3).stiffness_n_per_m2),
            symmetric_frequency=data["nearest_symmetric"].frequency_hz,
        )
    end
    fine = last(rows)
    csv_path = joinpath(OUTPUT_ROOT, "stiffness_mesh_convergence.csv")
    open(csv_path, "w") do io
        println(io, "level,size_min_mm,size_max_mm,K0_n_per_m2,K242_n_per_m2,K222_n_per_m2,K262_n_per_m2,symmetric_mode_hz,K0_relative_to_fine,K242_relative_to_fine,mode_relative_to_fine")
        for row in rows
            println(io, join((
                row.level.id,
                row.level.size_min_mm,
                row.level.size_max_mm,
                row.static_stiffness,
                row.target_stiffness,
                row.band_low_stiffness,
                row.band_high_stiffness,
                row.symmetric_frequency,
                (row.static_stiffness - fine.static_stiffness) / fine.static_stiffness,
                (row.target_stiffness - fine.target_stiffness) / fine.target_stiffness,
                (row.symmetric_frequency - fine.symmetric_frequency) / fine.symmetric_frequency,
            ), ','))
        end
    end
    x = getproperty.(getproperty.(rows, :level), :size_max_mm)
    stiffness_panel = plot(
        x,
        getproperty.(rows, :static_stiffness) ./ 1e9;
        marker=:circle,
        linewidth=2,
        label="K(0)",
        xlabel="maximum mesh size, mm",
        ylabel="stiffness, GN/m²",
        title="Dynamic-stiffness convergence",
        xflip=true,
        gridalpha=0.25,
    )
    plot!(
        stiffness_panel,
        x,
        getproperty.(rows, :target_stiffness) ./ 1e9;
        marker=:diamond,
        linewidth=2,
        label="K(242 kHz)",
    )
    mode_panel = plot(
        x,
        getproperty.(rows, :symmetric_frequency) ./ 1e3;
        marker=:circle,
        linewidth=2,
        label=false,
        xlabel="maximum mesh size, mm",
        ylabel="nearest symmetric mode, kHz",
        title="Fixed-interface spectrum convergence",
        xflip=true,
        gridalpha=0.25,
    )
    hline!(mode_panel, [TARGET_FREQUENCY_HZ / 1e3]; color=:black, linestyle=:dash,
           label="242 kHz")
    figure = plot(
        stiffness_panel,
        mode_panel;
        layout=(1, 2),
        size=(1300, 500),
        margin=5Plots.mm,
        plot_title="Mesh convergence of selected K/2 (t=0.60 mm)",
    )
    png_path = joinpath(OUTPUT_ROOT, "stiffness_mesh_convergence.png")
    savefig(figure, png_path)
    println("[+] $csv_path")
    println("[+] $png_path")
    for row in rows
        println(
            "  $(row.level.id): K242=", round(row.target_stiffness / 1e9; digits=5),
            " GN/m², f_sym=", round(row.symmetric_frequency / 1e3; digits=4), " kHz",
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
        error("unknown stiffness-convergence stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Stiffness convergence stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_stiffness_mesh_convergence.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
