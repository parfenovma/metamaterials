module TaperedShuntTrimPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_TAPERED_SHUNT_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "tapered_shunt_trim_pilot"),
)
const VARIANTS = (:solid, :mid, :max)
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "tapered_shunt_trim_mesher.jl"))
using .TaperedShuntTrimMesher

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
end

const CONFIG = TaperedShuntTrimConfig(
    lead_length_mm=parse(
        Float64,
        get(ENV, "METAMATERIALS_TAPERED_SHUNT_LEAD_MM", "12.0"),
    ),
)

if REQUESTED_STAGE in (
    "lossless", "fine-lossless", "fine-lossy", "analyze", "analyze-fine", "analyze-lossy",
    "validate",
)
    ENV["METAMATERIALS_SLOW_WAVE_OUTPUT"] = OUTPUT_ROOT
    ENV["METAMATERIALS_SLOW_WAVE_LEAD_MM"] = string(CONFIG.lead_length_mm)
    ENV["METAMATERIALS_THREE_STATE_TITLE"] =
        "Tapered direct path with symmetric shunt phase trim"
    ENV["METAMATERIALS_THREE_STATE_TOTAL_LENGTH_MM"] = string(total_length_mm(CONFIG))
    include(joinpath(@__DIR__, "run_distributed_slow_wave_pilot.jl"))
end

mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            build_tapered_shunt_trim_mesh(mesh_path(variant); config=CONFIG, variant)
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    for variant in VARIANTS
        convert_mesh(mesh_path(variant); output_dir=joinpath(OUTPUT_ROOT, "models"))
    end
end

function delegate_three_state_stage()
    DistributedSlowWavePilot.main(ARGS)
end

function main(args=ARGS)
    REQUESTED_STAGE == "mesh" && return run_mesh_stage()
    REQUESTED_STAGE == "convert" && return run_convert_stage()
    REQUESTED_STAGE in (
        "lossless", "fine-lossless", "fine-lossy", "analyze", "analyze-fine", "analyze-lossy",
        "validate",
    ) &&
        return delegate_three_state_stage()
    error("usage: julia run_tapered_shunt_trim_pilot.jl --stage=mesh|convert|lossless|fine-lossless|fine-lossy|analyze|analyze-fine|analyze-lossy|validate")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
