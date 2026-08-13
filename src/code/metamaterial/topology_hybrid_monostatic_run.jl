module TopologyHybridMonostaticRun

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HYBRID_MONOSTATIC_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_monostatic_g3"),
)
const TOPOLOGIES = (:bare, :straight_footprint, :pure_ttd, :hybrid)
const CASES = (:none, :void_r2_focus)

const REQUESTED_STAGE = begin
    match = findfirst(arg -> startswith(arg, "--stage="), ARGS)
    isnothing(match) ? :mesh : Symbol(split(ARGS[match], '='; limit=2)[2])
end

if REQUESTED_STAGE == :mesh
    include(joinpath(@__DIR__, "topology_hybrid_monostatic_mesher.jl"))
elseif REQUESTED_STAGE == :convert
    using Gridap
    using GridapGmsh
elseif REQUESTED_STAGE == :solve
    include(joinpath(@__DIR__, "topology_hybrid_monostatic_transient.jl"))
end

function option(args, name, default)
    prefix = "--$(name)="
    matches = filter(arg -> startswith(arg, prefix), args)
    length(matches) <= 1 || error("duplicate --$name")
    isempty(matches) ? default : split(only(matches), '='; limit=2)[2]
end

function selected_symbols(raw, allowed)
    raw == "all" && return collect(allowed)
    values = Symbol.(split(raw, ','))
    all(in(allowed), values) || error("unsupported selection: $raw")
    values
end

mesh_path(level, topology, case_id) =
    joinpath(OUTPUT_ROOT, "meshes", "mesh_$(level)_$(topology)_$(case_id).msh")
model_path(level, topology, case_id) =
    joinpath(OUTPUT_ROOT, "models", "model_$(level)_$(topology)_$(case_id).json")
signal_path(level, topology, case_id) =
    joinpath(OUTPUT_ROOT, "signals", "signal_$(level)_$(topology)_$(case_id).jld2")
material_path(level, topology, case_id) =
    joinpath(OUTPUT_ROOT, "manifests", "materials_$(level)_$(topology)_$(case_id).txt")
mesh_info_path(level, topology, case_id) =
    joinpath(OUTPUT_ROOT, "manifests", "mesh_$(level)_$(topology)_$(case_id).csv")

function level_mesh_parameters(level)
    level == :coarse && return (
        specimen_mesh_mm=1.35,
        guide_mesh_mm=0.45,
        defect_mesh_mm=0.30,
    )
    level == :refined && return (
        specimen_mesh_mm=1.00,
        guide_mesh_mm=0.32,
        defect_mesh_mm=0.22,
    )
    error("unknown level: $level")
end

function level_transient_parameters(level)
    level == :coarse && return (samples_per_period=18, final_time_s=80e-6)
    level == :refined && return (samples_per_period=24, final_time_s=80e-6)
    error("unknown level: $level")
end

function write_mesh_manifest(level, topology, case_id, info, mesh_parameters)
    mkpath(dirname(material_path(level, topology, case_id)))
    open(material_path(level, topology, case_id), "w") do io
        println(io, join(info.material_keys, ','))
    end
    path = mesh_info_path(level, topology, case_id)
    open(path, "w") do io
        println(io, "level,topology,case_id,node_count,element_count,source_curve_count,source_length_mm,pressure_scale,minimum_clearance_mm,specimen_mesh_mm,guide_mesh_mm,defect_mesh_mm")
        println(io, join((
            level, topology, case_id, info.node_count, info.element_count,
            info.source_curve_count, info.source_length_mm, info.pressure_scale,
            info.minimum_clearance_mm, mesh_parameters.specimen_mesh_mm,
            mesh_parameters.guide_mesh_mm, mesh_parameters.defect_mesh_mm,
        ), ','))
    end
end

function run_mesh(levels, topologies, cases)
    using_module = TopologyHybridMonostaticMesher
    for level in levels
        parameters = level_mesh_parameters(level)
        for topology in topologies, case_id in cases
            defect = only(filter(
                item -> item.id == case_id,
                using_module.NDTBenchmark.required_defects(),
            ))
            config = using_module.HybridMonostaticMeshConfig(
                ; topology, defect, parameters...,
            )
            started = time()
            info = using_module.build_hybrid_monostatic_mesh(
                mesh_path(level, topology, case_id); config,
            )
            write_mesh_manifest(level, topology, case_id, info, parameters)
            println("[+] mesh $(level)/$(topology)/$(case_id): $(info.node_count) nodes, $(info.element_count) cells, $(round(time()-started;digits=2)) s")
        end
    end
end

function run_convert(levels, topologies, cases)
    for level in levels, topology in topologies, case_id in cases
        input = mesh_path(level, topology, case_id)
        isfile(input) || error("mesh missing: $input")
        output = model_path(level, topology, case_id)
        mkpath(dirname(output))
        model = GridapGmsh.GmshDiscreteModel(input)
        Gridap.Io.to_json_file(model, output)
        println("[+] model $(level)/$(topology)/$(case_id) -> $output")
    end
end

function read_material_fractions(level, topology, case_id)
    path = material_path(level, topology, case_id)
    isfile(path) || error("material manifest missing: $path")
    raw = strip(read(path, String))
    isempty(raw) && error("empty material manifest: $path")
    parse.(Int, split(raw, ',')) ./ 1000
end

function run_solve(levels, topologies, cases)
    using_module = TopologyHybridMonostaticTransient
    for level in levels
        parameters = level_transient_parameters(level)
        config = using_module.HybridMonostaticTransientConfig(; parameters...)
        for topology in topologies, case_id in cases
            input = model_path(level, topology, case_id)
            fractions = read_material_fractions(level, topology, case_id)
            using_module.solve_hybrid_monostatic_transient(
                input,
                signal_path(level, topology, case_id);
                topology,
                defect_id=case_id,
                material_fractions=fractions,
                config,
            )
        end
    end
end

function main(args=ARGS)
    stage = Symbol(option(args, "stage", "mesh"))
    levels = selected_symbols(option(args, "level", "coarse"), (:coarse, :refined))
    default_topologies = levels == [:refined] ? "pure_ttd,hybrid" : "all"
    topologies = selected_symbols(option(args, "topology", default_topologies), TOPOLOGIES)
    cases = selected_symbols(option(args, "case", "all"), CASES)
    stage == :mesh && return run_mesh(levels, topologies, cases)
    stage == :convert && return run_convert(levels, topologies, cases)
    stage == :solve && return run_solve(levels, topologies, cases)
    error("usage: --stage=mesh|convert|solve [--level=coarse|refined] [--topology=all|...] [--case=all|...] ")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
