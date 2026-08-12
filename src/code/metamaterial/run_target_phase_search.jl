module TargetPhaseSearch

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "target_phase_search_242khz")
const REFERENCE_ROOT = joinpath(PROJECT_ROOT, "tmp", "period_length_pilot_220khz")
const CARRIER_FREQUENCY_HZ = 220.0e3
const TARGET_FREQUENCY_HZ = 242.0e3
const FINAL_TIME_S = 120.0e-6
const HEIGHT_MM = 7.0
const LENGTHS_MM = [33.0, 34.0, 35.0]
const GAPS_MM = [3.6, 3.7, 3.8, 3.9]

const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "profiles.jl"))
using .MetamaterialProfiles

if REQUESTED_STAGE == "mesh"
    include(joinpath(@__DIR__, "step1_mesher.jl"))
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) && startswith(REQUESTED_STAGE, "solve-one-")
    include(joinpath(@__DIR__, "step2_solver.jl"))
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    include(joinpath(@__DIR__, "step3_analyzer.jl"))
    using JLD2
    using Plots
end

struct SearchCase
    id::String
    length_mm::Float64
    gap_mm::Float64
    profile::SinusoidalProfile
    geometry::GeometryConfig
end

number_label(value::Real) = replace(string(Float64(value)), "." => "p")

function search_case(length_mm::Real, gap_mm::Real)
    length_value = Float64(length_mm)
    gap_value = Float64(gap_mm)
    profile = SinusoidalProfile((HEIGHT_MM - gap_value) / 2; periods=1)
    SearchCase(
        "sin_N_1_L_$(number_label(length_value))_G_$(number_label(gap_value))",
        length_value,
        gap_value,
        profile,
        GeometryConfig(length_mm=length_value, height_mm=HEIGHT_MM),
    )
end

search_cases() = [search_case(length_mm, gap_mm) for length_mm in LENGTHS_MM for gap_mm in GAPS_MM]
case_dir(kind, case) = joinpath(OUTPUT_ROOT, kind, case.id)
mesh_path(case) = joinpath(case_dir("meshes", case), "mesh_$(profile_slug(case.profile)).msh")
model_path(case) = joinpath(case_dir("models", case), "model_$(profile_slug(case.profile)).json")
signal_path(case) = joinpath(case_dir("signals", case), "data_$(profile_slug(case.profile))_F_220.0.jld2")

function reference_signal_path(case)
    label = number_label(case.length_mm)
    joinpath(
        REFERENCE_ROOT,
        "signals",
        "reference_L_$label",
        "data_sin_A_0.0_N_2_F_220.0.jld2",
    )
end

function run_mesh_stage()
    for case in search_cases()
        isfile(mesh_path(case)) && continue
        generate_meshes(
            [case.profile];
            geometry=case.geometry,
            mesh=MeshConfig(output_dir=case_dir("meshes", case)),
        )
    end
end

function run_conversion_stage()
    for case in search_cases()
        isfile(model_path(case)) && continue
        convert_mesh(mesh_path(case); output_dir=case_dir("models", case))
    end
end

function simulation(case)
    SimulationConfig(
        frequencies_hz=[CARRIER_FREQUENCY_HZ],
        final_time_s=FINAL_TIME_S,
        samples_per_period=30,
        pulse_cycles=4.0,
        save_vtk=false,
        skip_existing=true,
        model_dir=case_dir("models", case),
        vtk_dir=case_dir("vtk", case),
        signal_dir=case_dir("signals", case),
    )
end

function run_single_solver_stage(case_index::Integer)
    cases = search_cases()
    1 <= case_index <= length(cases) || error("invalid search case index: $case_index")
    case = cases[case_index]
    run_acoustic_simulation(case.profile, CARRIER_FREQUENCY_HZ; simulation=simulation(case))
end

function child_command(stage; threads::Integer=1)
    `$((Base.julia_cmd())) --startup-file=no --project=$PROJECT_ROOT --threads=$threads $(@__FILE__) --stage=$stage`
end

function process_count()
    requested = parse(Int, get(ENV, "METAMATERIALS_GEOMETRY_JOBS", "4"))
    requested > 0 || error("METAMATERIALS_GEOMETRY_JOBS must be positive")
    min(requested, length(search_cases()), Sys.CPU_THREADS)
end

function run_solver_stage()
    cases = search_cases()
    pending = findall(case -> !isfile(signal_path(case)), cases)
    isempty(pending) && return println("[=] All targeted FEM results already exist")
    jobs = min(process_count(), length(pending))
    println("=== Target-phase FEM: $(length(pending)) cases, $jobs processes ===")
    semaphore = Base.Semaphore(jobs)
    @sync for case_index in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command("solve-one-$case_index"; threads=1))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function circular_phase_distance_deg(phase_deg, target_deg)
    abs(mod(phase_deg - target_deg + 180, 360) - 180)
end

function analyze_case(case)
    output_path = save_analysis(
        signal_path(case);
        reference_path=reference_signal_path(case),
        output_dir=case_dir("characteristics", case),
    )
    transfer = JLD2.load(output_path)["calibrated_transmission"]
    point = SpectralAnalysis.value_at_frequency(transfer, TARGET_FREQUENCY_HZ)
    phase_deg = mod(rad2deg(angle(point.transfer)), 360)
    (
        case_id=case.id,
        periods=1,
        length_mm=case.length_mm,
        gap_mm=case.gap_mm,
        amplitude_mm=amplitude_mm(case.profile),
        sampled_frequency_hz=point.frequency_hz,
        H_real=real(point.transfer),
        H_imag=imag(point.transfer),
        amplitude=point.amplitude,
        phase_deg,
        distance_to_sector_deg=phase_deg < 150 ? 150 - phase_deg : phase_deg > 190 ? phase_deg - 190 : 0.0,
        spectral_valid=point.valid,
    )
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
    figure = scatter(
        getproperty.(rows, :phase_deg),
        getproperty.(rows, :amplitude);
        marker_z=getproperty.(rows, :gap_mm),
        color=:viridis,
        colorbar_title="gap, mm",
        xlabel="Phase at 242 kHz, deg",
        ylabel="|H|",
        title="Targeted search for efficient 150--190 degree elements",
        xlims=(100, 250),
        ylims=(0, max(0.45, 1.15maximum(getproperty.(rows, :amplitude)))),
        markersize=9,
        legend=false,
        gridalpha=0.25,
        size=(900, 620),
    )
    vspan!(figure, [150, 190]; color=:green, alpha=0.08, label=false)
    hline!(figure, [0.4]; color=:gray, linestyle=:dash, label=false)
    for row in rows
        annotate!(figure, row.phase_deg, row.amplitude, text("L=$(row.length_mm)", 7, :left))
    end
    path = joinpath(OUTPUT_ROOT, "target_phase_search.png")
    savefig(figure, path)
    path
end

function run_analysis_stage()
    rows = analyze_case.(search_cases())
    sort!(rows; by=row -> (row.distance_to_sector_deg, -row.amplitude))
    mkpath(OUTPUT_ROOT)
    csv_path = joinpath(OUTPUT_ROOT, "target_phase_search.csv")
    write_csv(csv_path, rows)
    plot_path = save_plot(rows)
    in_sector = filter(row -> row.distance_to_sector_deg == 0 && row.spectral_valid, rows)
    println("[+] $csv_path")
    println("[+] $plot_path")
    if isempty(in_sector)
        println("[-] No geometry reached 150--190 degrees")
    else
        best = in_sector[argmax(getproperty.(in_sector, :amplitude))]
        println("[+] Best in sector: $(best.case_id), phase=$(round(best.phase_deg; digits=1)) deg, |H|=$(round(best.amplitude; digits=3))")
    end
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_conversion_stage()
    elseif stage == "solve"
        run_solver_stage()
    elseif startswith(stage, "solve-one-")
        run_single_solver_stage(parse(Int, split(stage, '-')[3]))
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown target-phase stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Stage: $stage ===")
            run(child_command(stage; threads=1))
        end
    elseif args == ["--list"]
        foreach(case -> println("$(case.id): L=$(case.length_mm), gap=$(case.gap_mm)"), search_cases())
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_target_phase_search.jl [--list | --stage=mesh|convert|solve|analyze]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
