module RunVirtualNDTStand

include(joinpath(@__DIR__, "ndt_benchmark.jl"))
using .NDTBenchmark
include(joinpath(@__DIR__, "virtual_ndt_stand_mesher.jl"))
using .VirtualNDTStandMesher

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_NDT_STAND_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "virtual_ndt_stand"),
)

function selected_cases(args)
    isempty(args) && return required_defects()
    args == ["--required"] && return required_defects()
    args == ["--extended"] && return all_defects()
    if length(args) == 1 && startswith(only(args), "--case=")
        id = Symbol(split(only(args), '='; limit=2)[2])
        matches = filter(defect -> defect.id == id, all_defects())
        isempty(matches) && error("unknown defect case: $id")
        return [only(matches)]
    end
    error("usage: run_virtual_ndt_stand.jl [--required|--extended|--case=ID]")
end

function main(args=ARGS)
    config = NDTStandConfig()
    validate(config)
    println("[i] lambda_P = $(pressure_wavelength_mm(config)) mm")
    println("[i] lambda_S = $(shear_wavelength_mm(config)) mm")
    println("[i] burst duration = $(burst_duration_s(config) * 1e6) us")
    println("[i] nominal round trip = $(round_trip_time_s(config.nominal_focus_depth_mm, config) * 1e6) us")
    println("[i] burst-only dead-zone depth = $(nominal_dead_zone_depth_mm(config)) mm")
    for defect in selected_cases(args)
        path = joinpath(OUTPUT_ROOT, "meshes", "stand_$(defect.id).msh")
        info = build_virtual_ndt_stand_mesh(path; config, defect)
        println("[+] $(info.defect_id): $(info.mesh_path)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
