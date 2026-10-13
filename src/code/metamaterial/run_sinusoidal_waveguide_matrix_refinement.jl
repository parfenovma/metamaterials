module RunSinusoidalWaveguideMatrixRefinement

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PACKAGE_ROOT = joinpath(
    PROJECT_ROOT, "results", "sinusoidal_material_transfer_p6_242khz",
)
const OUTPUT_ROOT = joinpath(
    PACKAGE_ROOT, "time_reversal_broadband_matrix", "refinement_321700hz",
)
const DESIGN_ROOT = joinpath(
    PACKAGE_ROOT, "design_5cycle",
)
const FREQUENCY_HZ = 321_700.0
const OUTPUT_Y_MM = collect(-57.4:8.2:57.4)
const JOBS = parse(Int, get(ENV, "METAMATERIALS_SINUSOIDAL_REFINEMENT_JOBS", "2"))
const ROLE_MAP = Dict(:sinusoidal => :selected, :straight => :uniform)

ENV["METAMATERIALS_MATERIAL_LENS_DESIGN"] = DESIGN_ROOT
ENV["METAMATERIALS_MATERIAL_LENS_FEM_OUTPUT"] = OUTPUT_ROOT
ENV["METAMATERIALS_MATERIAL_LENS_FEM_ORDER"] = "1"
ENV["METAMATERIALS_MATERIAL_LENS_H_LENS_MM"] = "0.32"
ENV["METAMATERIALS_MATERIAL_LENS_H_OUTPUT_MM"] = "1.20"
ENV["METAMATERIALS_MATERIAL_LENS_SCAN"] = "false"

include(joinpath(@__DIR__, "run_material_lens_fem.jl"))
using .RunMaterialLensFEM

const Solver = RunMaterialLensFEM.MaterialLensHarmonicSolver

result_path(role) = joinpath(OUTPUT_ROOT, "$(role)_matrix.jld2")

function solve_one(role)
    role in keys(ROLE_MAP) || throw(ArgumentError("unknown role: $role"))
    path = result_path(role)
    isfile(path) && return println("[=] Existing refinement: $path")
    fem_role = ROLE_MAP[role]
    mesh_info = RunMaterialLensFEM.build_case_mesh(:aluminium, fem_role)
    mesh_path = mesh_info isa AbstractString ? mesh_info : mesh_info.mesh_path
    material = RunMaterialLensFEM.aluminium_6061()
    result = Solver.solve_material_lens_response_matrix(
        mesh_path, material, material;
        config=Solver.MaterialLensResponseConfig(
            frequency_hz=FREQUENCY_HZ,
            pressure_amplitude_pa=1.0e6,
            element_order=1,
            quadrature_degree=2,
            source_count=15,
            output_x_mm=33.0,
            output_y_mm=copy(OUTPUT_Y_MM),
        ),
    )
    jldsave(
        path;
        format_version=1,
        role=String(role),
        frequency_hz=FREQUENCY_HZ,
        lens_mesh_size_mm=0.32,
        output_mesh_size_mm=1.20,
        output_ux_m=result.output_ux_m,
        output_uy_m=result.output_uy_m,
        port_work_compliance=result.port_work_compliance,
        accepted_work=result.accepted_work,
        reciprocity_residual=result.reciprocity_residual,
        reflection_residual=result.reflection_residual,
        matrix_size=result.matrix_size,
        matrix_nonzeros=result.matrix_nonzeros,
        mesh_path=abspath(mesh_path),
    )
    println(
        "[+] refined $role: mirror=$(result.reflection_residual), " *
        "reciprocity=$(result.reciprocity_residual)",
    )
end

function child_command(role)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --solve-one=$(role)`
end

function solve_all()
    pending = [role for role in keys(ROLE_MAP) if !isfile(result_path(role))]
    isempty(pending) && return println("[=] Refinement matrices exist")
    jobs = min(max(JOBS, 1), length(pending), Sys.CPU_THREADS)
    semaphore = Base.Semaphore(jobs)
    @sync for role in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command(role))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function main(args=ARGS)
    if isempty(args) || args == ["--solve"]
        solve_all()
    elseif length(args) == 1 && startswith(only(args), "--solve-one=")
        solve_one(Symbol(split(only(args), '='; limit=2)[2]))
    else
        error("usage: run_sinusoidal_waveguide_matrix_refinement.jl --solve|--solve-one=ROLE")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

end
