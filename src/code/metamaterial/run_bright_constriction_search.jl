module BrightConstrictionSearch

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "bright_constriction_search")
const LENGTH_MM = 24.0
const HEIGHT_MM = 7.0
const FREQUENCIES_HZ = collect(210.0e3:2.0e3:270.0e3)
const GAPS_MM = [1.8, 2.4, 3.0]
const WIDTHS_MM = [1.5, 2.5, 3.5]

const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "profiles.jl"))
using .MetamaterialProfiles

if REQUESTED_STAGE == "mesh"
    include(joinpath(@__DIR__, "step1_mesher.jl"))
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif REQUESTED_STAGE == "solve-reference" ||
       (!isnothing(REQUESTED_STAGE) && startswith(REQUESTED_STAGE, "solve-one-"))
    include(joinpath(@__DIR__, "harmonic_solver.jl"))
    using .HarmonicElasticity
    using JLD2
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
end

struct BrightCase
    id::String
    gap_mm::Float64
    width_mm::Float64
    profile::CoupledConstrictionProfile
    geometry::GeometryConfig
end

label(value) = replace(string(Float64(value)), "." => "p")

function bright_case(gap_mm, width_mm)
    gap = Float64(gap_mm)
    width = Float64(width_mm)
    BrightCase(
        "single_G_$(label(gap))_W_$(label(width))",
        gap,
        width,
        CoupledConstrictionProfile([gap], [width], [LENGTH_MM / 2]; height_mm=HEIGHT_MM),
        GeometryConfig(length_mm=LENGTH_MM, height_mm=HEIGHT_MM, samples=240),
    )
end

bright_cases() = [bright_case(gap, width) for gap in GAPS_MM for width in WIDTHS_MM]
reference_profile() = SinusoidalProfile(0.0; periods=1)
reference_geometry() = GeometryConfig(length_mm=LENGTH_MM, height_mm=HEIGHT_MM, samples=240)
case_dir(kind, case) = joinpath(OUTPUT_ROOT, kind, case.id)
reference_dir(kind) = joinpath(OUTPUT_ROOT, kind, "reference")
mesh_path(case) = joinpath(case_dir("meshes", case), "mesh_$(profile_slug(case.profile)).msh")
model_path(case) = joinpath(case_dir("models", case), "model_$(profile_slug(case.profile)).json")
reference_mesh_path() = joinpath(reference_dir("meshes"), "mesh_$(profile_slug(reference_profile())).msh")
reference_model_path() = joinpath(reference_dir("models"), "model_$(profile_slug(reference_profile())).json")
result_path(case) = joinpath(case_dir("harmonic", case), "response.jld2")
reference_result_path() = joinpath(reference_dir("harmonic"), "response.jld2")

function run_mesh_stage()
    if !isfile(reference_mesh_path())
        generate_meshes(
            [reference_profile()];
            geometry=reference_geometry(),
            mesh=MeshConfig(output_dir=reference_dir("meshes")),
        )
    end
    for case in bright_cases()
        isfile(mesh_path(case)) && continue
        generate_meshes(
            [case.profile];
            geometry=case.geometry,
            mesh=MeshConfig(output_dir=case_dir("meshes", case)),
        )
    end
end

function run_convert_stage()
    if !isfile(reference_model_path())
        convert_mesh(reference_mesh_path(); output_dir=reference_dir("models"))
    end
    for case in bright_cases()
        isfile(model_path(case)) && continue
        convert_mesh(mesh_path(case); output_dir=case_dir("models", case))
    end
end

function save_response(path, model_path)
    points = solve_harmonic_sweep(model_path, FREQUENCIES_HZ)
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        frequency_hz=getproperty.(points, :frequency_hz),
        left_displacement=getproperty.(points, :left_displacement),
        right_displacement=getproperty.(points, :right_displacement),
        left_traction_pa=getproperty.(points, :left_traction_pa),
        right_traction_pa=getproperty.(points, :right_traction_pa),
    )
    println("[+] $path")
end


function run_reference_stage()
    isfile(reference_result_path()) && return println("[=] Existing reference response")
    save_response(reference_result_path(), reference_model_path())
end

function run_single_stage(index::Integer)
    cases = bright_cases()
    1 <= index <= length(cases) || error("invalid bright-case index: $index")
    case = cases[index]
    isfile(result_path(case)) && return println("[=] Existing response: $(case.id)")
    save_response(result_path(case), model_path(case))
end

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_solve_stage()
    run(child_command("solve-reference"))
    cases = bright_cases()
    pending = findall(case -> !isfile(result_path(case)), cases)
    isempty(pending) && return println("[=] All bright-resonator responses exist")
    jobs = min(parse(Int, get(ENV, "METAMATERIALS_GEOMETRY_JOBS", "4")), length(pending), Sys.CPU_THREADS)
    semaphore = Base.Semaphore(jobs)
    println("=== Harmonic bright search: $(length(pending)) cases, $jobs processes ===")
    @sync for index in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command("solve-one-$index"))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function response_rows(case, reference)
    sample = JLD2.load(result_path(case))
    [
        let transfer = sample["right_displacement"][index] / reference["right_displacement"][index]
            (
                case_id=case.id,
                gap_mm=case.gap_mm,
                width_mm=case.width_mm,
                frequency_hz=frequency_hz,
                H_real=real(transfer),
                H_imag=imag(transfer),
                amplitude=abs(transfer),
                phase_rad=angle(transfer),
            )
        end
        for (index, frequency_hz) in enumerate(FREQUENCIES_HZ)
    ]
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function save_spectra(rows)
    panels = Plots.Plot[]
    for gap in GAPS_MM
        panel = plot(
            xlabel="frequency, kHz",
            ylabel="calibrated |H|",
            title="single constriction, gap=$(gap) mm",
            gridalpha=0.25,
            legend=:best,
        )
        for width in WIDTHS_MM
            current = filter(row -> row.gap_mm == gap && row.width_mm == width, rows)
            plot!(
                panel,
                getproperty.(current, :frequency_hz) ./ 1e3,
                getproperty.(current, :amplitude);
                marker=:circle,
                markersize=3,
                linewidth=2,
                label="w=$(width) mm",
            )
        end
        vline!(panel, [242.0]; color=:gray, linestyle=:dash, label=false)
        push!(panels, panel)
    end
    figure = plot(panels...; layout=(3, 1), size=(1000, 1250), margin=5Plots.mm)
    path = joinpath(OUTPUT_ROOT, "single_constriction_spectra.png")
    savefig(figure, path)
    path
end

function run_analysis_stage()
    reference = JLD2.load(reference_result_path())
    rows = reduce(vcat, response_rows(case, reference) for case in bright_cases())
    summary = [
        let current = filter(row -> row.case_id == case.id, rows),
            minimum_row = current[argmin(getproperty.(current, :amplitude))],
            target_row = current[argmin(abs.(getproperty.(current, :frequency_hz) .- 242e3))]
            (
                case_id=case.id,
                gap_mm=case.gap_mm,
                width_mm=case.width_mm,
                minimum_frequency_hz=minimum_row.frequency_hz,
                minimum_amplitude=minimum_row.amplitude,
                amplitude_at_242khz=target_row.amplitude,
                phase_at_242khz_rad=target_row.phase_rad,
            )
        end
        for case in bright_cases()
    ]
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "single_constriction_response.csv"), rows)
    write_csv(joinpath(OUTPUT_ROOT, "single_constriction_summary.csv"), summary)
    println("[+] $(save_spectra(rows))")
    ranked = sort(summary; by=row -> (abs(row.minimum_frequency_hz - 242e3), row.minimum_amplitude))
    println("Closest broad-dip candidate: $(first(ranked))")
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_solve_stage()
    elseif stage == "solve-reference"
        run_reference_stage()
    elseif startswith(stage, "solve-one-")
        run_single_stage(parse(Int, split(stage, '-')[3]))
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown bright-search stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Stage: $stage ===")
            run(child_command(stage))
        end
    elseif args == ["--list"]
        foreach(case -> println(case.id), bright_cases())
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_bright_constriction_search.jl [--list | --stage=...]" )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
