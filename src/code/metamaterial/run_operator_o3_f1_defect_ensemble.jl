module RunOperatorO3F1DefectEnsemble

using JLD2

include(joinpath(@__DIR__, "ndt_benchmark.jl"))
include(joinpath(@__DIR__, "ndt_element_array_mesher.jl"))
include(joinpath(@__DIR__, "ndt_fullwave_transient_solver.jl"))

using .NDTBenchmark
using .NDTElementArrayMesher
using .NDTFullWaveTransientSolver

export f1_defect_cases, build_ensemble_meshes, solve_one, run_campaign

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_O3_F1_ENSEMBLE_OUTPUT",
    joinpath(PROJECT_ROOT, "results", "operator_space", "fullwave_transient_ensemble_v0"),
)
const CENTER_RESPONSE_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "fullwave_transient_v0", "responses",
)

function stand_config()
    NDTStandConfig(mesh_size_mm=1.0, defect_mesh_size_mm=0.25)
end

function f1_defect_cases()
    config = stand_config()
    Dict(defect.id => defect for defect in vcat(
        [required_defects(config)[2]],
        [
            NDTDefectCase(id=:void_r1_focus, center_x_mm=60.0, radius_mm=1.0),
            NDTDefectCase(id=:void_r2_near, center_x_mm=35.0, radius_mm=2.0),
            NDTDefectCase(id=:void_r2_deep, center_x_mm=85.0, radius_mm=2.0),
            NDTDefectCase(id=:void_r2_offset, center_x_mm=60.0, center_y_mm=10.0, radius_mm=2.0),
            NDTDefectCase(
                id=:slit_30_deg,
                kind=:slit,
                center_x_mm=60.0,
                major_radius_mm=4.0,
                minor_radius_mm=0.20,
                angle_deg=30.0,
            ),
        ],
    ))
end

case_root(case_id) = joinpath(OUTPUT_ROOT, string(case_id))
mesh_path(case_id, role) = joinpath(case_root(case_id), "meshes", "$(role).msh")
response_path(case_id, role) = if case_id == :void_r2_focus
    joinpath(CENTER_RESPONSE_ROOT, role == :free ? "none_transient.jld2" : "void_r2_focus_transient.jld2")
else
    joinpath(case_root(case_id), "responses", "$(role)_transient.jld2")
end

function build_ensemble_meshes(; force=false)
    config = stand_config()
    for (case_id, defect) in sort!(collect(f1_defect_cases()); by=first)
        case_id == :void_r2_focus && continue
        free_path = mesh_path(case_id, :free)
        defect_path = mesh_path(case_id, :defect)
        if force || !isfile(free_path) || !isfile(defect_path)
            println("[+] building common F1b mesh pair: $case_id")
            build_ndt_element_array_mesh_pair(
                free_path,
                defect_path;
                config,
                defect,
                element_count=17,
            )
        end
    end
end

function save_response(path, case_id, role, mesh, solved, config)
    mkpath(dirname(path))
    jldsave(
        path;
        format_version=1,
        case_id,
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
            rayleigh_alpha_s_inv=config.rayleigh_alpha_s_inv,
            rayleigh_beta_s=config.rayleigh_beta_s,
            element_count=config.element_count,
            element_order=config.element_order,
            quadrature_degree=config.quadrature_degree,
            newmark_beta=config.newmark_beta,
            newmark_gamma=config.newmark_gamma,
        ),
    )
    path
end

function solve_one(case_id::Symbol, role::Symbol; force=false)
    case_id in keys(f1_defect_cases()) || throw(ArgumentError("unknown defect case: $case_id"))
    role in (:free, :defect) || throw(ArgumentError("role must be free or defect"))
    case_id == :void_r2_focus && return response_path(case_id, role)
    path = response_path(case_id, role)
    if isfile(path) && !force
        println("[=] Existing F1b response: $path")
        return path
    end
    mesh = mesh_path(case_id, role)
    isfile(mesh) || error("ensemble mesh not found: $mesh")
    config = NDTFullWaveTransientConfig()
    println("[+] F1b transient: $case_id / $role")
    solved = solve_ndt_element_transient_matrix(mesh; config)
    save_response(path, case_id, role, mesh, solved, config)
    println("[+] Saved: $path")
    path
end

function child_command(case_id, role; force=false)
    julia = Base.julia_cmd()
    command = `$julia --startup-file=no --project=$PROJECT_ROOT $(@__FILE__) --solve-one=$(case_id):$(role)`
    force ? `$command --force` : command
end

function run_campaign(; force_mesh=false, force_solve=false)
    build_ensemble_meshes(; force=force_mesh)
    cases = sort!(collect(keys(f1_defect_cases())); by=string)
    for case_id in cases
        case_id == :void_r2_focus && continue
        for role in (:free, :defect)
            path = response_path(case_id, role)
            if force_solve || !isfile(path)
                run(child_command(case_id, role; force=force_solve))
            end
        end
    end
    println("[+] F1b finite-defect transient ensemble complete")
end

function main(args=ARGS)
    if any(startswith(argument, "--solve-one=") for argument in args)
        argument = only(filter(argument -> startswith(argument, "--solve-one="), args))
        case_text, role_text = split(split(argument, '='; limit=2)[2], ':'; limit=2)
        solve_one(Symbol(case_text), Symbol(role_text); force="--force" in args)
    elseif "--build-only" in args
        build_ensemble_meshes(; force="--force" in args)
    else
        run_campaign(
            ; force_mesh="--force-mesh" in args,
            force_solve="--force-solve" in args,
        )
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

end

