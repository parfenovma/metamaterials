module RunOperatorO3F1CScenes

using JLD2
using Random

include(joinpath(@__DIR__, "ndt_benchmark.jl"))
include(joinpath(@__DIR__, "ndt_multiscatter_scene_mesher.jl"))
include(joinpath(@__DIR__, "ndt_fullwave_transient_solver.jl"))

using .NDTBenchmark
using .NDTMultiscatterSceneMesher
using .NDTFullWaveTransientSolver

export context_rows,
       generate_backgrounds,
       build_scenes,
       solve_one,
       response_path,
       run_campaign

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_O3_F1C_SCENE_OUTPUT",
    joinpath(PROJECT_ROOT, "results", "operator_space", "fullwave_f1c_scenes_v0"),
)
const BASE_SEED = UInt64(0x051c2026)
const TARGET_Y_MM = [-28.0, -14.0, 0.0, 14.0, 28.0]

context_rows() = [
    (id=:train_1, split=:train, seed=BASE_SEED + 0x101),
    (id=:train_2, split=:train, seed=BASE_SEED + 0x102),
    (id=:validation_1, split=:validation, seed=BASE_SEED + 0x201),
    (id=:test_1, split=:test, seed=BASE_SEED + 0x301),
]

function valid_background(candidate, accepted)
    radius = candidate.radius_mm
    22.0 + radius < candidate.center_x_mm < 108.0 - radius || return false
    abs(candidate.center_y_mm) + radius < 32.0 || return false
    for target_y in TARGET_Y_MM
        hypot(candidate.center_x_mm - 60.0, candidate.center_y_mm - target_y) >
            radius + 1.0 + 3.0 || return false
    end
    for other in accepted
        hypot(
            candidate.center_x_mm - other.center_x_mm,
            candidate.center_y_mm - other.center_y_mm,
        ) > radius + other.radius_mm + 4.0 || return false
    end
    true
end

function generate_backgrounds(seed; count=5)
    rng = Xoshiro(seed)
    backgrounds = NDTDefectCase[]
    attempts = 0
    while length(backgrounds) < count
        attempts += 1
        attempts <= 10_000 || error("background rejection sampler stalled")
        candidate = NDTDefectCase(
            id=Symbol("background_$(length(backgrounds) + 1)"),
            center_x_mm=25.0 + 80.0rand(rng),
            center_y_mm=-27.0 + 54.0rand(rng),
            radius_mm=0.6 + 0.5rand(rng),
        )
        valid_background(candidate, backgrounds) && push!(backgrounds, candidate)
    end
    backgrounds
end

context_root(context_id) = joinpath(OUTPUT_ROOT, string(context_id))
mesh_root(context_id) = joinpath(context_root(context_id), "meshes")
response_root(context_id) = joinpath(context_root(context_id), "responses")
mesh_path(context_id, role) = role == :h0 ?
    joinpath(mesh_root(context_id), "h0_background.msh") :
    joinpath(mesh_root(context_id), "h1_$(role).msh")
response_path(context_id, role) = joinpath(response_root(context_id), "$(role)_transient.jld2")

function write_context_csv(path, rows)
    open(path, "w") do io
        println(io, "context_id,split,scatterer_id,x_mm,y_mm,radius_mm")
        for row in rows, defect in generate_backgrounds(row.seed)
            println(io, join((
                row.id, row.split, defect.id, defect.center_x_mm,
                defect.center_y_mm, defect.radius_mm,
            ), ','))
        end
    end
end

function build_scenes(; force=false)
    mkpath(OUTPUT_ROOT)
    for row in context_rows()
        root = mesh_root(row.id)
        expected = vcat(
            [joinpath(root, "h0_background.msh")],
            [joinpath(root, "h1_target_y$(replace(string(y), '.' => 'p', '-' => 'm')).msh") for y in TARGET_Y_MM],
        )
        if force || !all(isfile, expected)
            println("[+] building F1c structural scene $(row.id)")
            build_multiscatter_scene_meshes(
                root,
                generate_backgrounds(row.seed);
                target_y_mm=TARGET_Y_MM,
            )
        end
    end
    write_context_csv(joinpath(OUTPUT_ROOT, "structural_contexts.csv"), context_rows())
end

function save_response(path, context_id, role, mesh, solved, config)
    mkpath(dirname(path))
    jldsave(
        path;
        format_version=1,
        context_id,
        role,
        mesh_path=abspath(mesh),
        time_s=solved.time_s,
        source_drive=solved.source_drive,
        normal_response_m2=solved.normal_response_m2,
        tangent_response_m2=solved.tangent_response_m2,
        dt_s=solved.dt_s,
        matrix_size=solved.matrix_size,
        mass_nonzeros=solved.mass_nonzeros,
        damping_nonzeros=solved.damping_nonzeros,
        stiffness_nonzeros=solved.stiffness_nonzeros,
        effective_nonzeros=solved.effective_nonzeros,
        element_lengths_m=solved.element_lengths_m,
        solver_config=(
            frequency_hz=config.frequency_hz,
            pulse_cycles=config.pulse_cycles,
            pressure_amplitude_pa=config.pressure_amplitude_pa,
            final_time_s=config.final_time_s,
            samples_per_cycle=config.samples_per_cycle,
            density_kg_m3=config.density_kg_m3,
            pressure_wave_speed_m_s=config.pressure_wave_speed_m_s,
            shear_wave_speed_m_s=config.shear_wave_speed_m_s,
            element_count=config.element_count,
            element_order=config.element_order,
            quadrature_degree=config.quadrature_degree,
            newmark_beta=config.newmark_beta,
            newmark_gamma=config.newmark_gamma,
        ),
    )
end

function solve_one(context_id::Symbol, role::Symbol; force=false)
    context_id in getproperty.(context_rows(), :id) || throw(ArgumentError("unknown context"))
    valid_roles = vcat([:h0], [
        Symbol("target_y$(replace(string(y), '.' => 'p', '-' => 'm'))") for y in TARGET_Y_MM
    ])
    role in valid_roles || throw(ArgumentError("unknown scene role: $role"))
    output = response_path(context_id, role)
    if isfile(output) && !force
        println("[=] Existing F1c response: $output")
        return output
    end
    mesh = mesh_path(context_id, role)
    isfile(mesh) || error("F1c scene mesh not found: $mesh")
    config = NDTFullWaveTransientConfig()
    println("[+] F1c transient: $context_id / $role")
    solved = solve_ndt_element_transient_matrix(mesh; config)
    save_response(output, context_id, role, mesh, solved, config)
    println("[+] Saved: $output")
    output
end

function child_command(context_id, role; force=false)
    julia = Base.julia_cmd()
    command = `$julia --startup-file=no $(@__FILE__) --solve-one=$(context_id):$(role)`
    force ? `$command --force` : command
end

function run_campaign(; force_mesh=false, force_solve=false)
    build_scenes(; force=force_mesh)
    roles = vcat([:h0], [
        Symbol("target_y$(replace(string(y), '.' => 'p', '-' => 'm'))") for y in TARGET_Y_MM
    ])
    for row in context_rows(), role in roles
        if force_solve || !isfile(response_path(row.id, role))
            run(child_command(row.id, role; force=force_solve))
        end
    end
    println("[+] F1c structural transient scenes complete")
end

function main(args=ARGS)
    if any(startswith(argument, "--solve-one=") for argument in args)
        argument = only(filter(argument -> startswith(argument, "--solve-one="), args))
        context_text, role_text = split(split(argument, '='; limit=2)[2], ':'; limit=2)
        solve_one(Symbol(context_text), Symbol(role_text); force="--force" in args)
    elseif "--build-only" in args
        build_scenes(; force="--force" in args)
    else
        run_campaign(
            ; force_mesh="--force-mesh" in args,
            force_solve="--force-solve" in args,
        )
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

end
