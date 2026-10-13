module RunOperatorO5AFreshTestScenes

using JLD2

include(joinpath(@__DIR__, "run_operator_o3_f1c_scenes.jl"))

using .RunOperatorO3F1CScenes

export build_scene, solve_one, run_campaign

const CONTEXT_ID = :o5a_test_1
const CONTEXT_SEED = UInt64(0x05a52026)
const TARGET_Y_MM = RunOperatorO3F1CScenes.TARGET_Y_MM

mesh_root() = RunOperatorO3F1CScenes.mesh_root(CONTEXT_ID)
response_path(role) = RunOperatorO3F1CScenes.response_path(CONTEXT_ID, role)
mesh_path(role) = RunOperatorO3F1CScenes.mesh_path(CONTEXT_ID, role)
roles() = vcat([:h0], [
    Symbol("target_y$(replace(string(y), '.' => 'p', '-' => 'm'))")
    for y in TARGET_Y_MM
])

function build_scene(; force=false)
    expected = [mesh_path(role) for role in roles()]
    if force || !all(isfile, expected)
        println("[+] building O5a fresh structural test scene")
        RunOperatorO3F1CScenes.NDTMultiscatterSceneMesher.build_multiscatter_scene_meshes(
            mesh_root(),
            RunOperatorO3F1CScenes.generate_backgrounds(CONTEXT_SEED);
            target_y_mm=TARGET_Y_MM,
        )
    end
    csv_path = joinpath(RunOperatorO3F1CScenes.context_root(CONTEXT_ID), "structural_context.csv")
    open(csv_path, "w") do io
        println(io, "context_id,split,seed,scatterer_id,x_mm,y_mm,radius_mm")
        for defect in RunOperatorO3F1CScenes.generate_backgrounds(CONTEXT_SEED)
            println(io, join((
                CONTEXT_ID, :fresh_test, CONTEXT_SEED, defect.id,
                defect.center_x_mm, defect.center_y_mm, defect.radius_mm,
            ), ','))
        end
    end
    expected
end

function solve_one(role; force=false)
    role in roles() || throw(ArgumentError("unknown O5a role: $role"))
    output = response_path(role)
    if isfile(output) && !force
        println("[=] Existing O5a response: $output")
        return output
    end
    config = RunOperatorO3F1CScenes.NDTFullWaveTransientConfig()
    println("[+] O5a fresh transient: $role")
    solved = RunOperatorO3F1CScenes.solve_ndt_element_transient_matrix(
        mesh_path(role); config,
    )
    RunOperatorO3F1CScenes.save_response(
        output, CONTEXT_ID, role, mesh_path(role), solved, config,
    )
    println("[+] Saved: $output")
    output
end

function child_command(role; force=false)
    command = `$(Base.julia_cmd()) --startup-file=no $(@__FILE__) --solve-one=$(role)`
    force ? `$command --force` : command
end

function run_campaign(; force_mesh=false, force_solve=false)
    build_scene(; force=force_mesh)
    for role in roles()
        if force_solve || !isfile(response_path(role))
            run(child_command(role; force=force_solve))
        end
    end
    println("[+] O5a fresh structural test scene complete")
end

function main(args=ARGS)
    if any(startswith(argument, "--solve-one=") for argument in args)
        argument = only(filter(argument -> startswith(argument, "--solve-one="), args))
        solve_one(Symbol(split(argument, '='; limit=2)[2]); force="--force" in args)
    elseif "--build-only" in args
        build_scene(; force="--force" in args)
    else
        run_campaign(
            ; force_mesh="--force-mesh" in args,
            force_solve="--force-solve" in args,
        )
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

end
