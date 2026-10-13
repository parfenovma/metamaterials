module RunOperatorF1DContextScenes

include(joinpath(@__DIR__, "run_operator_o3_f1c_scenes.jl"))

using .RunOperatorO3F1CScenes

export context_rows, build_context, solve_one, run_context, run_campaign

const CONTEXT_ROWS = [
    (id=:f1d_validation_1, split=:validation, seed=UInt64(0x051d2101)),
    (id=:f1d_validation_2, split=:validation, seed=UInt64(0x051d2102)),
    (id=:f1d_test_1, split=:test, seed=UInt64(0x051d2301)),
    (id=:f1d_test_2, split=:test, seed=UInt64(0x051d2302)),
]
const TARGET_Y_MM = RunOperatorO3F1CScenes.TARGET_Y_MM

context_rows() = copy(CONTEXT_ROWS)
roles() = vcat([:h0], [
    Symbol("target_y$(replace(string(y), '.' => 'p', '-' => 'm'))")
    for y in TARGET_Y_MM
])
row_for(context_id) = only(filter(row -> row.id == context_id, CONTEXT_ROWS))
mesh_root(context_id) = RunOperatorO3F1CScenes.mesh_root(context_id)
mesh_path(context_id, role) = RunOperatorO3F1CScenes.mesh_path(context_id, role)
response_path(context_id, role) = RunOperatorO3F1CScenes.response_path(context_id, role)

function build_context(context_id; force=false)
    row = row_for(context_id)
    expected = [mesh_path(context_id, role) for role in roles()]
    if force || !all(isfile, expected)
        println("[+] building F1d $(row.split) context $(row.id)")
        RunOperatorO3F1CScenes.NDTMultiscatterSceneMesher.build_multiscatter_scene_meshes(
            mesh_root(context_id),
            RunOperatorO3F1CScenes.generate_backgrounds(row.seed);
            target_y_mm=TARGET_Y_MM,
        )
    end
    root = RunOperatorO3F1CScenes.context_root(context_id)
    open(joinpath(root, "structural_context.csv"), "w") do io
        println(io, "context_id,split,seed,scatterer_id,x_mm,y_mm,radius_mm")
        for defect in RunOperatorO3F1CScenes.generate_backgrounds(row.seed)
            println(io, join((
                row.id, row.split, row.seed, defect.id,
                defect.center_x_mm, defect.center_y_mm, defect.radius_mm,
            ), ','))
        end
    end
    expected
end

function solve_one(context_id, role; force=false)
    row_for(context_id)
    role in roles() || throw(ArgumentError("unknown F1d role: $role"))
    output = response_path(context_id, role)
    if isfile(output) && !force
        println("[=] Existing F1d response: $output")
        return output
    end
    config = RunOperatorO3F1CScenes.NDTFullWaveTransientConfig()
    println("[+] F1d transient: $context_id / $role")
    solved = RunOperatorO3F1CScenes.solve_ndt_element_transient_matrix(
        mesh_path(context_id, role); config,
    )
    RunOperatorO3F1CScenes.save_response(
        output, context_id, role, mesh_path(context_id, role), solved, config,
    )
    println("[+] Saved: $output")
    output
end

function child_command(context_id, role; force=false)
    command = `$(Base.julia_cmd()) --startup-file=no $(@__FILE__) --solve-one=$(context_id):$(role)`
    force ? `$command --force` : command
end

function run_context(context_id; force_mesh=false, force_solve=false)
    build_context(context_id; force=force_mesh)
    for role in roles()
        if force_solve || !isfile(response_path(context_id, role))
            run(child_command(context_id, role; force=force_solve))
        end
    end
    println("[+] F1d context complete: $context_id")
end

function run_campaign(; force_mesh=false, force_solve=false)
    for row in CONTEXT_ROWS
        run_context(row.id; force_mesh, force_solve)
    end
    println("[+] F1d expanded physical-context scenes complete")
end

function main(args=ARGS)
    if any(startswith(argument, "--solve-one=") for argument in args)
        argument = only(filter(argument -> startswith(argument, "--solve-one="), args))
        context_text, role_text = split(split(argument, '='; limit=2)[2], ':'; limit=2)
        solve_one(Symbol(context_text), Symbol(role_text); force="--force" in args)
    elseif any(startswith(argument, "--context=") for argument in args)
        argument = only(filter(argument -> startswith(argument, "--context="), args))
        context_id = Symbol(split(argument, '='; limit=2)[2])
        run_context(
            context_id;
            force_mesh="--force-mesh" in args,
            force_solve="--force-solve" in args,
        )
    elseif "--build-only" in args
        for row in CONTEXT_ROWS
            build_context(row.id; force="--force" in args)
        end
    else
        run_campaign(
            ; force_mesh="--force-mesh" in args,
            force_solve="--force-solve" in args,
        )
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

end
