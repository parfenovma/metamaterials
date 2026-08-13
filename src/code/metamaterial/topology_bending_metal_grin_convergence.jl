module TopologyBendingMetalGRINConvergence

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV, "METAMATERIALS_BENDING_METAL_GRIN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "bending_metal_grin"),
)
const STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
              split(ARGS[1], '='; limit=2)[2] : error("use --stage=mesh|convert|solve")

include(joinpath(@__DIR__, "topology_bending_metal_grin_mesher.jl"))
using .TopologyBendingMetalGRINMesher
if STAGE == "mesh"
    using Gmsh: gmsh
elseif STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif STAGE == "solve"
    include(joinpath(@__DIR__, "multimode_port_harmonic_solver.jl"))
    using .MultimodePortHarmonicSolver
    using .SinusoidalMaterialLens: aluminium_6061
end

const LEVELS = (
    (name="coarse", ligament_mm=0.20, bulk_mm=0.40),
    (name="baseline", ligament_mm=0.16, bulk_mm=0.34),
    (name="fine", ligament_mm=0.11, bulk_mm=0.24),
)

mesh_path(level) = joinpath(OUTPUT_ROOT, "convergence", "mesh_$(level.name).msh")
model_path(level) = joinpath(OUTPUT_ROOT, "convergence", "model_$(level.name).json")

function mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for level in LEVELS
            build_bending_metal_grin_mesh(
                mesh_path(level);
                variant=:dense_c,
                size_ligament_mm=level.ligament_mm,
                size_bulk_mm=level.bulk_mm,
            )
            println("[+] convergence mesh $(level.name)")
        end
    finally
        gmsh.finalize()
    end
end

function convert_stage()
    for level in LEVELS
        convert_mesh(mesh_path(level); output_dir=dirname(model_path(level)))
    end
end

function solve_stage()
    rows = NamedTuple[]
    for level in LEVELS
        state = solve_multimode_port_harmonic(
            model_path(level), ["Microphone"], aluminium_6061();
            config=MultimodePortConfig(
                frequency_hz=242e3,
                element_order=2,
                quadrature_degree=4,
                port_element_count=120,
                evanescent_modes_per_port=8,
            ),
        )
        incident = state.incident_power_w_per_m
        transmission = only(state.transmitted_power_w_per_m) / incident
        reflection = state.reflected_power_w_per_m / incident
        push!(rows, (
            mesh_level=level.name,
            size_ligament_mm=level.ligament_mm,
            size_bulk_mm=level.bulk_mm,
            file_size_bytes=filesize(mesh_path(level)),
            T_all=transmission,
            R_all=reflection,
            power_sum=transmission + reflection,
        ))
        println("[+] $(level.name): T=$transmission R=$reflection")
    end
    baseline = only(filter(r -> r.mesh_level == "baseline", rows))
    enriched = [merge(r, (
        delta_T_vs_baseline=r.T_all - baseline.T_all,
        delta_R_vs_baseline=r.R_all - baseline.R_all,
    )) for r in rows]
    path = joinpath(OUTPUT_ROOT, "g1_mesh_convergence.csv")
    open(path, "w") do io
        columns = propertynames(first(enriched))
        println(io, join(string.(columns), ','))
        for row in enriched
            println(io, join((string(getproperty(row, c)) for c in columns), ','))
        end
    end
    println("[+] $path")
end

STAGE == "mesh" ? mesh_stage() :
STAGE == "convert" ? convert_stage() :
STAGE == "solve" ? solve_stage() : error("unknown stage $STAGE")

end
