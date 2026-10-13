module RunSinusoidalWaveguideBroadbandMatrix

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra

include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
using .SinusoidalMaterialLens
include(joinpath(@__DIR__, "material_lens_harmonic_solver.jl"))
using .MaterialLensHarmonicSolver

export run_pilot, solve_all

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PACKAGE_ROOT = joinpath(
    PROJECT_ROOT, "results", "sinusoidal_material_transfer_p6_242khz",
)
const TAGGED_MESH_ROOT = joinpath(PACKAGE_ROOT, "response_matrix", "base", "meshes")
const OUTPUT_ROOT = joinpath(PACKAGE_ROOT, "time_reversal_broadband_matrix")
const LEGACY_MATRIX_PATH = joinpath(
    PACKAGE_ROOT, "response_matrix", "base",
    "aluminium_strip_group_response_matrix.jld2",
)
const FREQUENCIES_HZ = Float64[
    162_200, 178_160, 194_120, 202_100, 210_080, 218_060, 226_040,
    234_020, 242_000, 249_970, 257_940, 265_910, 273_880, 281_850,
    289_820, 305_760, 321_700,
]
const ROLES = (:sinusoidal, :straight)
const ROLE_STEMS = Dict(:sinusoidal => :selected, :straight => :uniform)
const CARRIER_HZ = 242.0e3
const OUTPUT_Y_MM = collect(-57.4:8.2:57.4)
const JOBS = parse(Int, get(ENV, "METAMATERIALS_SINUSOIDAL_MATRIX_JOBS", "2"))

mesh_path(role) = joinpath(
    TAGGED_MESH_ROOT, "aluminium_$(ROLE_STEMS[role]).msh",
)

result_path(role, frequency_hz) = joinpath(
    OUTPUT_ROOT, "frequency_matrices",
    "$(role)_f$(round(Int, frequency_hz))hz.jld2",
)

function solve_one(role, frequency_hz)
    role in ROLES || throw(ArgumentError("unknown role: $role"))
    frequency_hz in FREQUENCIES_HZ || throw(ArgumentError("frequency is outside frozen grid"))
    path = result_path(role, frequency_hz)
    isfile(path) && return println("[=] Existing matrix: $path")
    mesh = mesh_path(role)
    isfile(mesh) || error("missing frozen broadband mesh: $mesh")
    material = aluminium_6061()
    result = solve_material_lens_response_matrix(
        mesh, material, material;
        config=MaterialLensResponseConfig(
            frequency_hz=frequency_hz,
            pressure_amplitude_pa=1.0e6,
            element_order=1,
            quadrature_degree=2,
            source_count=15,
            output_x_mm=33.0,
            output_y_mm=copy(OUTPUT_Y_MM),
        ),
    )
    mkpath(dirname(path))
    jldsave(
        path;
        format_version=1,
        role=String(role),
        frequency_hz=result.frequency_hz,
        pressure_amplitude_pa=result.pressure_amplitude_pa,
        output_x_mm=result.output_x_mm,
        output_y_mm=result.output_y_mm,
        output_ux_m=result.output_ux_m,
        output_uy_m=result.output_uy_m,
        port_work_compliance=result.port_work_compliance,
        accepted_work=result.accepted_work,
        source_lengths_m=result.source_lengths_m,
        reciprocity_residual=result.reciprocity_residual,
        port_reflection_residual=result.port_reflection_residual,
        reflection_residual=result.reflection_residual,
        matrix_size=result.matrix_size,
        matrix_nonzeros=result.matrix_nonzeros,
        mesh_path=abspath(mesh),
    )
    println(
        "[+] $role $(frequency_hz / 1e3) kHz: mirror=" *
        "$(result.reflection_residual), reciprocity=$(result.reciprocity_residual)",
    )
    println("[+] $path")
    path
end

function child_command(role, frequency_hz)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --solve-one=$(role):$(frequency_hz)`
end

function solve_specs(specifications)
    pending = filter(
        spec -> !isfile(result_path(spec.role, spec.frequency_hz)),
        specifications,
    )
    isempty(pending) && return println("[=] Requested broadband matrices exist")
    jobs = min(max(JOBS, 1), length(pending), Sys.CPU_THREADS)
    semaphore = Base.Semaphore(jobs)
    println("=== Full sinusoidal response matrices: $(length(pending)) factorizations, $jobs processes ===")
    @sync for spec in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command(spec.role, spec.frequency_hz))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function grouped_even_matrix(full_matrix)
    columns = Matrix{ComplexF64}[]
    push!(columns, full_matrix[:, 8:8])
    for offset in 1:7
        push!(columns, full_matrix[:, (8 - offset):(8 - offset)] .+
            full_matrix[:, (8 + offset):(8 + offset)])
    end
    hcat(columns...)
end

function pilot_rows()
    isfile(LEGACY_MATRIX_PATH) || error("missing legacy carrier matrix")
    legacy = JLD2.load(LEGACY_MATRIX_PATH)
    rows = NamedTuple[]
    for role in ROLES
        current = JLD2.load(result_path(role, CARRIER_HZ))
        legacy_role = role == :sinusoidal ? "selected" : "uniform"
        expected = ComplexF64.(legacy[legacy_role][:preout])
        grouped = grouped_even_matrix(ComplexF64.(current["output_ux_m"]))
        reconstructed = grouped[8:15, :]
        reconstruction_error = norm(reconstructed - expected) /
            max(norm(expected), eps(Float64))
        projected = 0.5 .* (
            current["output_ux_m"] +
            reverse(reverse(current["output_ux_m"]; dims=1); dims=2)
        )
        projected_reflection_residual = norm(
            reverse(reverse(projected; dims=1); dims=2) - projected,
        ) / max(norm(projected), eps(Float64))
        push!(rows, (
            role,
            reconstruction_error,
            reflection_residual=current["reflection_residual"],
            projected_reflection_residual,
            reciprocity_residual=current["reciprocity_residual"],
            minimum_accepted_work=minimum(current["accepted_work"]),
            passed=reconstruction_error <= 1.0e-10 &&
                current["reflection_residual"] <= 0.02 &&
                projected_reflection_residual <= 1.0e-12 &&
                current["reciprocity_residual"] <= 1.0e-8 &&
                minimum(current["accepted_work"]) > 0,
        ))
    end
    rows
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function run_pilot()
    solve_specs([(; role, frequency_hz=CARRIER_HZ) for role in ROLES])
    rows = pilot_rows()
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "carrier_pilot.csv"), rows)
    for row in rows
        println(
            "[+] carrier pilot $(row.role): reconstruction=" *
            "$(row.reconstruction_error), pass=$(row.passed)",
        )
    end
    all(getproperty.(rows, :passed)) || error("carrier matrix pilot failed")
    rows
end

function solve_all()
    rows = run_pilot()
    all(getproperty.(rows, :passed)) || error("carrier pilot did not pass")
    solve_specs([
        (; role, frequency_hz)
        for role in ROLES for frequency_hz in FREQUENCIES_HZ
    ])
end

function parse_solve_one(argument)
    value = split(argument, '='; limit=2)[2]
    role_text, frequency_text = split(value, ':'; limit=2)
    Symbol(role_text), parse(Float64, frequency_text)
end

function main(args=ARGS)
    if args == ["--pilot"]
        run_pilot()
    elseif args == ["--solve"] || isempty(args)
        solve_all()
    elseif length(args) == 1 && startswith(only(args), "--solve-one=")
        role, frequency_hz = parse_solve_one(only(args))
        solve_one(role, frequency_hz)
    else
        error("usage: run_sinusoidal_waveguide_broadband_matrix.jl --pilot|--solve|--solve-one=ROLE:FREQUENCY")
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end

end
