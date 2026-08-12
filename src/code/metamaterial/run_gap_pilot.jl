module GapPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PILOT_ROOT = joinpath(PROJECT_ROOT, "tmp", "gap_pilot_220khz")
const PERIOD_LENGTH_ROOT = joinpath(PROJECT_ROOT, "tmp", "period_length_pilot_220khz")
const PILOT_FREQUENCY_HZ = 220.0e3
const TARGET_FREQUENCIES_HZ = collect(225.0e3:2.0e3:250.0e3)
const PILOT_FINAL_TIME_S = 120.0e-6
const MATERIAL_HEIGHT_MM = 7.0
const ORIGINAL_GAP_VALUES_MM = [1.5, 2.0, 3.0, 4.0]
const TARGET_GAP_VALUES_MM = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
const GAP_VALUES_MM = sort(unique(vcat(ORIGINAL_GAP_VALUES_MM, TARGET_GAP_VALUES_MM)))
const ANCHORS = [
    (periods=1, length_mm=12.0),
    (periods=1, length_mm=22.0),
    (periods=2, length_mm=14.5),
    (periods=2, length_mm=17.0),
    (periods=1, length_mm=33.0),
    (periods=1, length_mm=34.0),
    (periods=1, length_mm=35.0),
]

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

struct GapCase
    id::String
    periods::Int
    length_mm::Float64
    gap_mm::Float64
    profile::SinusoidalProfile
    geometry::GeometryConfig
end

number_label(value::Real) = replace(string(Float64(value)), "." => "p")

function gap_case(periods::Integer, length_mm::Real, gap_mm::Real)
    length_value = Float64(length_mm)
    gap_value = Float64(gap_mm)
    # With staggered walls, A is the vertical channel offset and the constant
    # profiled thickness is H - A (rather than H - 2A).
    amplitude_mm = MATERIAL_HEIGHT_MM - gap_value
    profile = SinusoidalProfile(amplitude_mm; periods)
    geometry = GeometryConfig(length_mm=length_value, height_mm=MATERIAL_HEIGHT_MM)
    GapCase(
        "sin_N_$(periods)_L_$(number_label(length_value))_G_$(number_label(gap_value))",
        Int(periods),
        length_value,
        gap_value,
        profile,
        geometry,
    )
end

function gap_cases()
    [
        gap_case(anchor.periods, anchor.length_mm, gap_mm)
        for anchor in ANCHORS
        for gap_mm in (
            anchor.length_mm >= 33.0 ? TARGET_GAP_VALUES_MM : ORIGINAL_GAP_VALUES_MM
        )
    ]
end

new_cases() = gap_cases()

case_dir(kind::AbstractString, case::GapCase) = joinpath(PILOT_ROOT, kind, case.id)

function mesh_path(case::GapCase)
    joinpath(case_dir("meshes", case), "mesh_$(profile_slug(case.profile)).msh")
end

function model_path(case::GapCase)
    joinpath(case_dir("models", case), "model_$(profile_slug(case.profile)).json")
end

function signal_path(case::GapCase)
    joinpath(
        case_dir("signals", case),
        "data_$(profile_slug(case.profile))_F_220.0.jld2",
    )
end

function reference_signal_path(case::GapCase)
    joinpath(
        PERIOD_LENGTH_ROOT,
        "signals",
        "reference_L_$(number_label(case.length_mm))",
        "data_sin_A_0.0_N_2_F_220.0.jld2",
    )
end

function anchor_cases(periods::Integer, length_mm::Real)
    sort(
        filter(
            case -> case.periods == periods && case.length_mm == length_mm,
            gap_cases(),
        );
        by=case -> case.gap_mm,
    )
end

function run_mesh_stage()
    for case in new_cases()
        if isfile(mesh_path(case))
            println("  [=] Existing mesh: $(mesh_path(case))")
            continue
        end
        mesh = MeshConfig(output_dir=case_dir("meshes", case))
        generate_meshes([case.profile]; geometry=case.geometry, mesh)
    end
end

function run_conversion_stage()
    for case in new_cases()
        if isfile(model_path(case))
            println("  [=] Existing model: $(model_path(case))")
            continue
        end
        convert_mesh(mesh_path(case); output_dir=case_dir("models", case))
    end
end

function simulation(case::GapCase)
    SimulationConfig(
        frequencies_hz=[PILOT_FREQUENCY_HZ],
        final_time_s=PILOT_FINAL_TIME_S,
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
    cases = new_cases()
    1 <= case_index <= length(cases) || error("invalid gap case index: $case_index")
    case = cases[case_index]
    run_acoustic_simulation(
        case.profile,
        PILOT_FREQUENCY_HZ;
        simulation=simulation(case),
    )
end

function process_count()
    requested = parse(Int, get(ENV, "METAMATERIALS_GEOMETRY_JOBS", "4"))
    requested > 0 || error("METAMATERIALS_GEOMETRY_JOBS must be positive")
    min(requested, length(new_cases()), Sys.CPU_THREADS)
end

function child_command(stage; threads::Integer=1)
    julia = Base.julia_cmd()
    script = @__FILE__
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=$threads $script --stage=$stage`
end

function run_solver_stage()
    cases = new_cases()
    pending = findall(case -> !isfile(signal_path(case)), cases)
    if isempty(pending)
        println("[=] All gap FEM results already exist")
        return
    end

    jobs = min(process_count(), length(pending))
    println("=== Gap FEM: $(length(pending)) pending cases, $jobs processes ===")
    semaphore = Base.Semaphore(jobs)
    @sync for case_index in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command("solve-one-$(case_index)"; threads=1))
            finally
                Base.release(semaphore)
            end
        end
    end
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

function result_row(case::GapCase, point, excess_envelope_delay_s)
    transmission_class = if point.amplitude >= 0.1
        "working"
    elseif point.amplitude >= 0.05
        "lossy"
    else
        "deep_minimum"
    end
    (
        case_id=case.id,
        periods=case.periods,
        length_mm=case.length_mm,
        gap_mm=case.gap_mm,
        amplitude_mm=amplitude_mm(case.profile),
        requested_frequency_hz=PILOT_FREQUENCY_HZ,
        sampled_frequency_hz=point.frequency_hz,
        H_real=real(point.transfer),
        H_imag=imag(point.transfer),
        amplitude=point.amplitude,
        magnitude_squared=point.magnitude_squared,
        phase_rad=angle(point.transfer),
        phase_mod_2pi_rad=mod(angle(point.transfer), 2.0 * pi),
        group_delay_s=point.group_delay_s,
        excess_envelope_delay_s,
        valid=point.valid,
        transmission_class,
    )
end

function analyze_case(case::GapCase)
    sample_path = signal_path(case)
    reference_path = reference_signal_path(case)
    isfile(sample_path) || error("sample result is missing: $sample_path")
    isfile(reference_path) || error("reference result is missing: $reference_path")
    output_path = save_analysis(
        sample_path;
        reference_path,
        output_dir=case_dir("characteristics", case),
    )
    saved = JLD2.load(output_path)
    point = SpectralAnalysis.value_at_frequency(
        saved["calibrated_transmission"],
        PILOT_FREQUENCY_HZ,
    )
    result_row(case, point, saved["excess_envelope_delay_s"])
end

function analyze_target_frequency_map(case::GapCase)
    output_path = save_analysis(
        signal_path(case);
        reference_path=reference_signal_path(case),
        output_dir=case_dir("characteristics", case),
    )
    transfer = JLD2.load(output_path)["calibrated_transmission"]
    [
        let point = SpectralAnalysis.value_at_frequency(transfer, frequency_hz)
            (
                case_id=case.id,
                length_mm=case.length_mm,
                gap_mm=case.gap_mm,
                requested_frequency_hz=frequency_hz,
                sampled_frequency_hz=point.frequency_hz,
                amplitude=point.amplitude,
                phase_rad=angle(point.transfer),
                valid=point.valid,
            )
        end
        for frequency_hz in TARGET_FREQUENCIES_HZ
    ]
end

function row_for_case(rows, case::GapCase)
    only(filter(row -> row.case_id == case.id, rows))
end

function unwrap_parameter_phase(rows)
    isempty(rows) && return Float64[]
    phases = Float64[first(rows).phase_rad]
    previous_raw = first(rows).phase_rad
    for row in Iterators.drop(rows, 1)
        raw = row.phase_rad
        push!(phases, last(phases) + mod(raw - previous_raw + pi, 2.0 * pi) - pi)
        previous_raw = raw
    end
    phases
end

anchor_label(periods, length_mm) = "N=$periods, L=$(length_mm) mm"

function save_gap_plot(rows)
    common = (
        linewidth=2,
        markersize=7,
        gridalpha=0.25,
        xlabel="Minimum gap, mm",
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
        top_margin=5Plots.mm,
    )
    amplitude_plot = plot(;
        ylabel="|H|",
        title="Transmission vs minimum gap",
        legend=:outertopright,
        xticks=GAP_VALUES_MM,
        common...,
    )
    phase_plot = plot(;
        ylabel="unwrapped phase, rad",
        title="Phase vs minimum gap",
        legend=:outertopright,
        xticks=GAP_VALUES_MM,
        common...,
    )
    markers = (:circle, :diamond, :square, :utriangle, :dtriangle, :hexagon, :star5)
    for (anchor, marker) in zip(ANCHORS, markers)
        current = [
            row_for_case(rows, case)
            for case in anchor_cases(anchor.periods, anchor.length_mm)
        ]
        label = anchor_label(anchor.periods, anchor.length_mm)
        plot!(
            amplitude_plot,
            getproperty.(current, :gap_mm),
            getproperty.(current, :amplitude);
            marker,
            label,
        )
        plot!(
            phase_plot,
            getproperty.(current, :gap_mm),
            unwrap_parameter_phase(current);
            marker,
            label,
        )
    end
    hline!(amplitude_plot, [0.1]; color=:gray, linestyle=:dash, label=false)
    figure = plot(
        amplitude_plot,
        phase_plot;
        layout=(2, 1),
        size=(1200, 980),
        margin=4Plots.mm,
    )
    figure_dir = joinpath(PILOT_ROOT, "figures")
    mkpath(figure_dir)
    path = joinpath(figure_dir, "gap_sensitivity.png")
    savefig(figure, path)
    println("[+] Gap sensitivity plot: $path")
    path
end

function save_complex_plot(rows)
    plot_radius = max(0.1, 1.25 * maximum(getproperty.(rows, :amplitude)))
    figure = plot(;
        aspect_ratio=:equal,
        xlabel="Re H",
        ylabel="Im H",
        title="Gap trajectories at 220 kHz",
        gridalpha=0.25,
        xlims=(-plot_radius, plot_radius),
        ylims=(-plot_radius, plot_radius),
        size=(850, 740),
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
        legend=:outertopright,
    )
    hline!(figure, [0.0]; color=:gray, linestyle=:dot, label=false)
    vline!(figure, [0.0]; color=:gray, linestyle=:dot, label=false)
    markers = (:circle, :diamond, :square, :utriangle, :dtriangle, :hexagon, :star5)
    for (anchor, marker) in zip(ANCHORS, markers)
        current = [
            row_for_case(rows, case)
            for case in anchor_cases(anchor.periods, anchor.length_mm)
        ]
        plot!(
            figure,
            getproperty.(current, :H_real),
            getproperty.(current, :H_imag);
            marker,
            markersize=7,
            linewidth=2,
            label=anchor_label(anchor.periods, anchor.length_mm),
        )
    end
    figure_dir = joinpath(PILOT_ROOT, "figures")
    mkpath(figure_dir)
    path = joinpath(figure_dir, "gap_complex_H.png")
    savefig(figure, path)
    println("[+] Gap complex-H plot: $path")
    path
end

function save_target_frequency_plots(rows)
    target_lengths = [33.0, 34.0, 35.0]
    panels = Plots.Plot[]
    for length_mm in target_lengths
        panel = plot(
            xlabel="Frequency, kHz",
            ylabel="|H|",
            title="L=$(length_mm) mm",
            gridalpha=0.25,
            legend=length_mm == first(target_lengths) ? :topleft : false,
        )
        for gap_mm in TARGET_GAP_VALUES_MM
            current = sort(filter(rows) do row
                row.length_mm == length_mm && row.gap_mm == gap_mm && row.valid
            end; by=row -> row.sampled_frequency_hz)
            plot!(
                panel,
                getproperty.(current, :sampled_frequency_hz) ./ 1.0e3,
                getproperty.(current, :amplitude);
                linewidth=2,
                label="g=$(gap_mm) mm",
            )
        end
        push!(panels, panel)
    end
    curves = plot(
        panels...;
        layout=(1, 3),
        size=(1450, 520),
        plot_title="Targeted long-branch frequency responses",
        margin=5Plots.mm,
    )
    curves_path = joinpath(PILOT_ROOT, "figures", "target_gap_frequency_curves.png")
    savefig(curves, curves_path)

    best_rows = [
        let current = filter(rows) do row
                row.length_mm == length_mm && row.gap_mm == gap_mm && row.valid
            end
            current[argmax(getproperty.(current, :amplitude))]
        end
        for length_mm in target_lengths, gap_mm in TARGET_GAP_VALUES_MM
    ]
    maxima = plot(
        xlabel="Minimum gap, mm",
        ylabel="maximum |H| in 225--250 kHz",
        title="Peak transmission near the long-branch optimum",
        gridalpha=0.25,
        size=(900, 580),
    )
    for (index, length_mm) in enumerate(target_lengths)
        plot!(
            maxima,
            TARGET_GAP_VALUES_MM,
            getproperty.(best_rows[index, :], :amplitude);
            marker=:circle,
            linewidth=2,
            label="L=$(length_mm) mm",
        )
    end
    maxima_path = joinpath(PILOT_ROOT, "figures", "target_gap_peak_transmission.png")
    savefig(maxima, maxima_path)
    println("[+] Target frequency curves: $curves_path")
    println("[+] Target peak transmission: $maxima_path")
end

function run_analysis_stage()
    rows = [analyze_case(case) for case in gap_cases()]
    target_cases = filter(case -> case.length_mm >= 33.0, gap_cases())
    target_frequency_rows = reduce(vcat, analyze_target_frequency_map.(target_cases))
    mkpath(PILOT_ROOT)
    csv_path = joinpath(PILOT_ROOT, "gap_pilot.csv")
    target_csv_path = joinpath(PILOT_ROOT, "target_gap_frequency_map.csv")
    jld_path = joinpath(PILOT_ROOT, "gap_pilot.jld2")
    write_csv(csv_path, rows)
    write_csv(target_csv_path, target_frequency_rows)
    JLD2.jldsave(
        jld_path;
        format_version=1,
        frequency_hz=PILOT_FREQUENCY_HZ,
        gaps_mm=GAP_VALUES_MM,
        anchors=ANCHORS,
        rows,
    )
    save_gap_plot(rows)
    save_complex_plot(rows)
    save_target_frequency_plots(target_frequency_rows)
    println("[+] Gap pilot: $csv_path")
    println("[+] Target gap/frequency data: $target_csv_path")
    println("[+] Machine-readable gap pilot: $jld_path")
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_conversion_stage()
    elseif stage == "solve"
        run_solver_stage()
    elseif startswith(stage, "solve-one-")
        parts = split(stage, '-')
        length(parts) == 3 || error("invalid single-solver stage: $stage")
        run_single_solver_stage(parse(Int, parts[3]))
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown gap-pilot stage: $stage")
    end
end

function list_pilot()
    println("Pilot root: $PILOT_ROOT")
    println("frequency_hz: $PILOT_FREQUENCY_HZ")
    for case in gap_cases()
        println(
            "  $(case.id): N=$(case.periods), L=$(case.length_mm) mm, ",
            "gap=$(case.gap_mm) mm, A=$(amplitude_mm(case.profile)) mm, staggered",
        )
    end
end

function run_pipeline()
    for stage in ("mesh", "convert", "solve", "analyze")
        println("\n=== Stage: $stage ===")
        run(child_command(stage; threads=1))
    end
end

function main(args=ARGS)
    if isempty(args)
        run_pipeline()
    elseif args == ["--list"]
        list_pilot()
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_gap_pilot.jl [--list | --stage=mesh|convert|solve|analyze]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
