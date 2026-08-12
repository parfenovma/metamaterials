module CoupledPairSearch

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(ENV, "METAMATERIALS_PAIR_OUTPUT", joinpath(PROJECT_ROOT, "tmp", "coupled_pair_search"))
const BRIGHT_ROOT = joinpath(PROJECT_ROOT, "tmp", "bright_constriction_search")
const LENGTH_MM = 24.0
const HEIGHT_MM = 7.0
const GAP_MM = 1.8
const WIDTH_MM = 3.5
parse_values(name, default) = parse.(Float64, split(get(ENV, name, default), ','))
const SECOND_GAPS_MM = parse_values("METAMATERIALS_PAIR_SECOND_GAPS_MM", "1.8")
const DISTANCES_MM = parse_values("METAMATERIALS_PAIR_DISTANCES_MM", "4,5,6,7,8")
const FREQUENCIES_HZ = collect(
    parse(Float64, get(ENV, "METAMATERIALS_PAIR_FREQUENCY_START_HZ", "210000")):
    parse(Float64, get(ENV, "METAMATERIALS_PAIR_FREQUENCY_STEP_HZ", "1000")):
    parse(Float64, get(ENV, "METAMATERIALS_PAIR_FREQUENCY_STOP_HZ", "270000")),
)
const DAMPING_SCALE = parse(Float64, get(ENV, "METAMATERIALS_PAIR_DAMPING_SCALE", "1"))

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

struct PairCase
    id::String
    distance_mm::Float64
    second_gap_mm::Float64
    profile::CoupledConstrictionProfile
    geometry::GeometryConfig
end

label(value) = replace(string(Float64(value)), "." => "p")

function pair_case(distance_mm, second_gap_mm=GAP_MM)
    distance = Float64(distance_mm)
    second_gap = Float64(second_gap_mm)
    centre = LENGTH_MM / 2
    centers = [centre - distance / 2, centre + distance / 2]
    PairCase(
        second_gap == GAP_MM ?
        "pair_G_$(label(GAP_MM))_W_$(label(WIDTH_MM))_D_$(label(distance))" :
        "pair_G_$(label(GAP_MM))-$(label(second_gap))_W_$(label(WIDTH_MM))_D_$(label(distance))",
        distance,
        second_gap,
        CoupledConstrictionProfile(
            [GAP_MM, second_gap],
            [WIDTH_MM, WIDTH_MM],
            centers;
            height_mm=HEIGHT_MM,
        ),
        GeometryConfig(length_mm=LENGTH_MM, height_mm=HEIGHT_MM, samples=240),
    )
end

pair_cases() = [pair_case(distance, second_gap) for distance in DISTANCES_MM for second_gap in SECOND_GAPS_MM]
reference_profile() = SinusoidalProfile(0.0; periods=1)
reference_geometry() = GeometryConfig(length_mm=LENGTH_MM, height_mm=HEIGHT_MM, samples=240)
case_dir(kind, case) = joinpath(OUTPUT_ROOT, kind, case.id)
reference_dir(kind) = joinpath(OUTPUT_ROOT, kind, "reference")
mesh_path(case) = joinpath(case_dir("meshes", case), "mesh_$(profile_slug(case.profile)).msh")
model_path(case) = joinpath(case_dir("models", case), "model_$(profile_slug(case.profile)).json")
result_path(case) = joinpath(case_dir("harmonic", case), "response.jld2")
reference_mesh_path() = joinpath(reference_dir("meshes"), "mesh_$(profile_slug(reference_profile())).msh")
reference_model_path() = joinpath(reference_dir("models"), "model_$(profile_slug(reference_profile())).json")
reference_result_path() = joinpath(reference_dir("harmonic"), "response.jld2")

function run_mesh_stage()
    if !isfile(reference_mesh_path())
        generate_meshes([reference_profile()]; geometry=reference_geometry(), mesh=MeshConfig(output_dir=reference_dir("meshes")))
    end
    for case in pair_cases()
        isfile(mesh_path(case)) && continue
        generate_meshes([case.profile]; geometry=case.geometry, mesh=MeshConfig(output_dir=case_dir("meshes", case)))
    end
end

function run_convert_stage()
    isfile(reference_model_path()) || convert_mesh(reference_mesh_path(); output_dir=reference_dir("models"))
    for case in pair_cases()
        isfile(model_path(case)) || convert_mesh(mesh_path(case); output_dir=case_dir("models", case))
    end
end

function save_response(path, model)
    harmonic_config = HarmonicConfig(
        rayleigh_alpha=DAMPING_SCALE * 79560.0,
        rayleigh_beta=DAMPING_SCALE * 2.5e-9,
    )
    points = solve_harmonic_sweep(model, FREQUENCIES_HZ; config=harmonic_config)
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

run_reference_stage() = isfile(reference_result_path()) ? println("[=] Existing reference response") : save_response(reference_result_path(), reference_model_path())

function run_single_stage(index)
    cases = pair_cases()
    1 <= index <= length(cases) || error("invalid pair-case index: $index")
    case = cases[index]
    isfile(result_path(case)) ? println("[=] Existing response: $(case.id)") : save_response(result_path(case), model_path(case))
end

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_solve_stage()
    run(child_command("solve-reference"))
    cases = pair_cases()
    pending = findall(case -> !isfile(result_path(case)), cases)
    isempty(pending) && return println("[=] All coupled-pair responses exist")
    jobs = min(parse(Int, get(ENV, "METAMATERIALS_GEOMETRY_JOBS", "4")), length(pending), Sys.CPU_THREADS)
    semaphore = Base.Semaphore(jobs)
    println("=== Harmonic pair search: $(length(pending)) cases, $jobs processes ===")
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

function unwrap_phase(phases)
    result = Float64[first(phases)]
    for raw in Iterators.drop(phases, 1)
        push!(result, last(result) + mod(raw - last(result) + pi, 2pi) - pi)
    end
    result
end

function response_rows(case, reference)
    sample = JLD2.load(result_path(case))
    transfer = sample["right_displacement"] ./ reference["right_displacement"]
    phases = unwrap_phase(angle.(transfer))
    delays = similar(phases)
    omega = 2pi .* FREQUENCIES_HZ
    delays[1] = -(phases[2] - phases[1]) / (omega[2] - omega[1])
    for index in 2:(length(phases) - 1)
        delays[index] = -(phases[index + 1] - phases[index - 1]) / (omega[index + 1] - omega[index - 1])
    end
    delays[end] = -(phases[end] - phases[end - 1]) / (omega[end] - omega[end - 1])
    [
        (
            case_id=case.id,
            distance_mm=case.distance_mm,
            second_gap_mm=case.second_gap_mm,
            frequency_hz=frequency_hz,
            H_real=real(transfer[index]),
            H_imag=imag(transfer[index]),
            amplitude=abs(transfer[index]),
            phase_rad=phases[index],
            group_delay_s=delays[index],
        )
        for (index, frequency_hz) in enumerate(FREQUENCIES_HZ)
    ]
end

function local_peaks(rows)
    peaks = NamedTuple[]
    amplitudes = getproperty.(rows, :amplitude)
    for index in 2:(length(rows) - 1)
        amplitudes[index] > amplitudes[index - 1] || continue
        amplitudes[index] >= amplitudes[index + 1] || continue
        left_minimum = minimum(amplitudes[max(1, index - 10):index])
        right_minimum = minimum(amplitudes[index:min(length(rows), index + 10)])
        push!(peaks, (
            index,
            frequency_hz=rows[index].frequency_hz,
            amplitude=amplitudes[index],
            prominence=amplitudes[index] - max(left_minimum, right_minimum),
            group_delay_s=rows[index].group_delay_s,
        ))
    end
    peaks
end

function read_csv_rows(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1) if !isempty(strip(line))]
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

function save_plot(rows)
    amplitude_plot = plot(
        xlabel="frequency, kHz",
        ylabel="calibrated |H|",
        title="Coupled constrictions: transparency search (damping scale=$(DAMPING_SCALE))",
        gridalpha=0.25,
        legend=:outertopright,
    )
    delay_plot = plot(
        xlabel="frequency, kHz",
        ylabel="group delay, us",
        title="Phase-derived delay",
        gridalpha=0.25,
        legend=:outertopright,
    )
    bright_rows = read_csv_rows(joinpath(BRIGHT_ROOT, "single_constriction_response.csv"))
    bright = filter(row -> parse(Float64, row["gap_mm"]) == GAP_MM && parse(Float64, row["width_mm"]) == WIDTH_MM, bright_rows)
    plot!(
        amplitude_plot,
        parse.(Float64, getindex.(bright, "frequency_hz")) ./ 1e3,
        parse.(Float64, getindex.(bright, "amplitude"));
        color=:black,
        linestyle=:dash,
        linewidth=2,
        label="single",
    )
    for case in pair_cases()
        current = filter(row -> row.case_id == case.id, rows)
        curve_label = length(SECOND_GAPS_MM) == 1 ?
                      "d=$(case.distance_mm) mm" :
                      "d=$(case.distance_mm), g2=$(case.second_gap_mm) mm"
        plot!(amplitude_plot, getproperty.(current, :frequency_hz) ./ 1e3, getproperty.(current, :amplitude); linewidth=2, label=curve_label)
        reliable_delay = [row.amplitude >= 0.01 ? row.group_delay_s * 1e6 : NaN for row in current]
        plot!(delay_plot, getproperty.(current, :frequency_hz) ./ 1e3, reliable_delay; linewidth=2, label=curve_label)
    end
    vline!(amplitude_plot, [242.0]; color=:gray, linestyle=:dot, label=false)
    vline!(delay_plot, [242.0]; color=:gray, linestyle=:dot, label=false)
    figure = plot(amplitude_plot, delay_plot; layout=(2, 1), size=(1100, 950), margin=5Plots.mm)
    path = joinpath(OUTPUT_ROOT, "coupled_pair_spectra.png")
    savefig(figure, path)
    path
end

function run_analysis_stage()
    reference = JLD2.load(reference_result_path())
    rows = reduce(vcat, response_rows(case, reference) for case in pair_cases())
    peak_rows = NamedTuple[]
    for case in pair_cases()
        current = filter(row -> row.case_id == case.id, rows)
        for peak in local_peaks(current)
            push!(peak_rows, (
                case_id=case.id,
                distance_mm=case.distance_mm,
                second_gap_mm=case.second_gap_mm,
                frequency_hz=peak.frequency_hz,
                amplitude=peak.amplitude,
                prominence=peak.prominence,
                group_delay_s=peak.group_delay_s,
            ))
        end
    end
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "coupled_pair_response.csv"), rows)
    isempty(peak_rows) || write_csv(joinpath(OUTPUT_ROOT, "coupled_pair_peaks.csv"), peak_rows)
    println("[+] $(save_plot(rows))")
    if isempty(peak_rows)
        println("[-] No internal transparency peak found")
    else
        best = peak_rows[argmax(getproperty.(peak_rows, :prominence))]
        println("[+] Most prominent internal peak: $best")
    end
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
        error("unknown pair-search stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Stage: $stage ===")
            run(child_command(stage))
        end
    elseif args == ["--list"]
        foreach(case -> println(case.id), pair_cases())
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_coupled_pair_search.jl [--list | --stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
