module PeriodLengthPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PILOT_ROOT = joinpath(PROJECT_ROOT, "tmp", "period_length_pilot_220khz")
const PILOT_FREQUENCY_HZ = 220.0e3
const FREQUENCY_MAP_HZ = collect(180.0e3:2.0e3:280.0e3)
const PILOT_FINAL_TIME_S = 120.0e-6
const PILOT_AMPLITUDE_MM = 2.5
const BASELINE_LENGTH_MM = 17.0
const PERIOD_VALUES = [1, 2, 3, 4]
const ORIGINAL_LENGTH_VALUES_MM = [12.0, 14.5, 17.0, 19.5, 22.0]
const EXTENDED_LENGTH_VALUES_MM = [25.0, 28.0, 31.0, 33.0, 34.0, 35.0, 37.0, 40.0, 43.0, 46.0]
const LENGTH_VALUES_MM = vcat(ORIGINAL_LENGTH_VALUES_MM, EXTENDED_LENGTH_VALUES_MM)

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

struct PilotCase
    id::String
    role::Symbol
    profile::WallProfile
    geometry::GeometryConfig
end

length_label(length_mm::Real) = replace(string(Float64(length_mm)), "." => "p")

function sample_case(periods::Integer, length_mm::Real)
    length_value = Float64(length_mm)
    PilotCase(
        "sin_A_2p5_N_$(periods)_L_$(length_label(length_value))",
        :sample,
        SinusoidalProfile(PILOT_AMPLITUDE_MM; periods),
        GeometryConfig(length_mm=length_value),
    )
end

function reference_case(length_mm::Real)
    length_value = Float64(length_mm)
    PilotCase(
        "reference_L_$(length_label(length_value))",
        :reference,
        SinusoidalProfile(0.0; periods=2),
        GeometryConfig(length_mm=length_value),
    )
end

function sample_cases()
    period_cases = [sample_case(periods, BASELINE_LENGTH_MM) for periods in PERIOD_VALUES]
    length_cases = [
        sample_case(periods, length_mm)
        for periods in (1, 2)
        for length_mm in LENGTH_VALUES_MM
        if length_mm != BASELINE_LENGTH_MM &&
           (periods == 1 || length_mm in ORIGINAL_LENGTH_VALUES_MM)
    ]
    vcat(period_cases, length_cases)
end

reference_cases() = [reference_case(length_mm) for length_mm in LENGTH_VALUES_MM]
all_cases() = vcat(reference_cases(), sample_cases())

case_dir(kind::AbstractString, case::PilotCase) = joinpath(PILOT_ROOT, kind, case.id)

function mesh_path(case::PilotCase)
    joinpath(case_dir("meshes", case), "mesh_$(profile_slug(case.profile)).msh")
end

function model_path(case::PilotCase)
    joinpath(case_dir("models", case), "model_$(profile_slug(case.profile)).json")
end

function signal_path(case::PilotCase)
    frequency_khz = PILOT_FREQUENCY_HZ / 1000.0
    joinpath(
        case_dir("signals", case),
        "data_$(profile_slug(case.profile))_F_$(frequency_khz).jld2",
    )
end

function matching_reference(sample::PilotCase)
    only(filter(reference_cases()) do reference
        reference.geometry.length_mm == sample.geometry.length_mm
    end)
end

period_sweep_cases() = sort(
    filter(case -> case.geometry.length_mm == BASELINE_LENGTH_MM, sample_cases());
    by=case -> case.profile.periods,
)

length_sweep_cases(periods::Integer=2) = sort(
    filter(case -> case.profile.periods == periods, sample_cases());
    by=case -> case.geometry.length_mm,
)

function run_mesh_stage()
    for case in all_cases()
        if isfile(mesh_path(case))
            println("  [=] Existing mesh: $(mesh_path(case))")
            continue
        end
        mesh = MeshConfig(output_dir=case_dir("meshes", case))
        generate_meshes([case.profile]; geometry=case.geometry, mesh)
    end
end

function run_conversion_stage()
    for case in all_cases()
        if isfile(model_path(case))
            println("  [=] Existing model: $(model_path(case))")
            continue
        end
        convert_mesh(mesh_path(case); output_dir=case_dir("models", case))
    end
end

function simulation(case::PilotCase)
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
    cases = all_cases()
    1 <= case_index <= length(cases) || error("invalid case index: $case_index")
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
    min(requested, length(all_cases()), Sys.CPU_THREADS)
end

function child_command(stage; threads::Integer=1)
    julia = Base.julia_cmd()
    script = @__FILE__
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=$threads $script --stage=$stage`
end

function run_solver_stage()
    cases = all_cases()
    pending = findall(case -> !isfile(signal_path(case)), cases)
    if isempty(pending)
        println("[=] All period/length FEM results already exist")
        return
    end

    jobs = min(process_count(), length(pending))
    println("=== Period/length FEM: $(length(pending)) pending cases, $jobs processes ===")
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

function result_row(case::PilotCase, point, excess_envelope_delay_s)
    active_length_mm = case.geometry.length_mm - 2.0 * case.geometry.end_margin_mm
    transmission_class = if point.amplitude >= 0.1
        "working"
    elseif point.amplitude >= 0.05
        "lossy"
    else
        "deep_minimum"
    end
    (
        case_id=case.id,
        amplitude_mm=amplitude_mm(case.profile),
        periods=case.profile.periods,
        length_mm=case.geometry.length_mm,
        active_length_mm,
        period_pitch_mm=active_length_mm / case.profile.periods,
        minimum_gap_mm=minimum_gap_mm(case.profile, case.geometry),
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

function analyze_case(case::PilotCase)
    reference = matching_reference(case)
    output_dir = case_dir("characteristics", case)
    output_path = save_analysis(
        signal_path(case);
        reference_path=signal_path(reference),
        output_dir,
    )
    saved = JLD2.load(output_path)
    point = SpectralAnalysis.value_at_frequency(
        saved["calibrated_transmission"],
        PILOT_FREQUENCY_HZ,
    )
    result_row(case, point, saved["excess_envelope_delay_s"])
end

function analyze_frequency_map(case::PilotCase)
    reference = matching_reference(case)
    output_path = save_analysis(
        signal_path(case);
        reference_path=signal_path(reference),
        output_dir=case_dir("characteristics", case),
    )
    transfer = JLD2.load(output_path)["calibrated_transmission"]
    [
        let point = SpectralAnalysis.value_at_frequency(transfer, frequency_hz)
            (
                case_id=case.id,
                periods=case.profile.periods,
                length_mm=case.geometry.length_mm,
                requested_frequency_hz=frequency_hz,
                sampled_frequency_hz=point.frequency_hz,
                amplitude=point.amplitude,
                phase_rad=angle(point.transfer),
                valid=point.valid,
            )
        end
        for frequency_hz in FREQUENCY_MAP_HZ
    ]
end

function row_for_case(rows, case::PilotCase)
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

function save_parameter_plot(rows)
    period_rows = [row_for_case(rows, case) for case in period_sweep_cases()]
    length_rows = Dict(
        periods => [row_for_case(rows, case) for case in length_sweep_cases(periods)]
        for periods in (1, 2)
    )

    common = (
        linewidth=2,
        markersize=7,
        gridalpha=0.25,
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
        top_margin=5Plots.mm,
    )
    p1 = plot(
        getproperty.(period_rows, :periods),
        getproperty.(period_rows, :amplitude);
        marker=:circle,
        xlabel="Number of periods N",
        ylabel="|H|",
        title="Transmission vs periods",
        xticks=PERIOD_VALUES,
        legend=false,
        common...,
    )
    hline!(p1, [0.1]; linestyle=:dash, color=:gray)
    p2 = plot(
        getproperty.(period_rows, :periods),
        getproperty.(period_rows, :phase_mod_2pi_rad);
        marker=:circle,
        xlabel="Number of periods N",
        ylabel="phase mod 2pi, rad",
        title="Phase vs periods",
        xticks=PERIOD_VALUES,
        legend=false,
        common...,
    )
    p3 = plot(;
        xlabel="Corrugated block length, mm",
        ylabel="|H|",
        title="Transmission vs length",
        xticks=LENGTH_VALUES_MM,
        legend=:topright,
        common...,
    )
    p4 = plot(;
        xlabel="Corrugated block length, mm",
        ylabel="unwrapped phase, rad",
        title="Phase vs length",
        xticks=LENGTH_VALUES_MM,
        legend=:topright,
        common...,
    )
    for (periods, marker) in zip((1, 2), (:circle, :diamond))
        current = length_rows[periods]
        plot!(
            p3,
            getproperty.(current, :length_mm),
            getproperty.(current, :amplitude);
            marker,
            label="N=$periods",
        )
        plot!(
            p4,
            getproperty.(current, :length_mm),
            unwrap_parameter_phase(current);
            marker,
            label="N=$periods",
        )
    end
    hline!(p3, [0.1]; linestyle=:dash, color=:gray, label=false)
    figure = plot(
        p1,
        p2,
        p3,
        p4;
        layout=(2, 2),
        size=(1250, 950),
        margin=4Plots.mm,
    )

    figure_dir = joinpath(PILOT_ROOT, "figures")
    mkpath(figure_dir)
    path = joinpath(figure_dir, "period_length_sensitivity.png")
    savefig(figure, path)
    println("[+] Parameter sensitivity plot: $path")
    path
end

function save_complex_plot(rows)
    plot_radius = max(0.1, 1.25 * maximum(getproperty.(rows, :amplitude)))
    figure = plot(
        aspect_ratio=:equal,
        xlabel="Re H",
        ylabel="Im H",
        title="Calibrated transmission at 220 kHz",
        gridalpha=0.25,
        xlims=(-plot_radius, plot_radius),
        ylims=(-plot_radius, plot_radius),
        size=(820, 720),
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    hline!(figure, [0.0]; color=:gray, linestyle=:dot, label=false)
    vline!(figure, [0.0]; color=:gray, linestyle=:dot, label=false)
    for (periods, marker) in zip((1, 2), (:circle, :diamond))
        current = [row_for_case(rows, case) for case in length_sweep_cases(periods)]
        plot!(
            figure,
            getproperty.(current, :H_real),
            getproperty.(current, :H_imag);
            marker,
            markersize=7,
            linewidth=2,
            label="N=$periods, length sweep",
        )
    end
    for periods in (3, 4)
        row = only(filter(row -> row.periods == periods, rows))
        scatter!(
            figure,
            [row.H_real],
            [row.H_imag];
            markersize=8,
            label="N=$periods, L=17 mm",
        )
    end
    figure_dir = joinpath(PILOT_ROOT, "figures")
    mkpath(figure_dir)
    path = joinpath(figure_dir, "period_length_complex_H.png")
    savefig(figure, path)
    println("[+] Complex-H plot: $path")
    path
end

function save_frequency_map(rows)
    long_lengths = filter(>=(minimum(EXTENDED_LENGTH_VALUES_MM)), LENGTH_VALUES_MM)
    amplitude_grid = [
        only(filter(rows) do row
            row.length_mm == length_mm && row.requested_frequency_hz == frequency_hz
        end).amplitude
        for frequency_hz in FREQUENCY_MAP_HZ, length_mm in long_lengths
    ]
    valid_grid = [
        only(filter(rows) do row
            row.length_mm == length_mm && row.requested_frequency_hz == frequency_hz
        end).valid
        for frequency_hz in FREQUENCY_MAP_HZ, length_mm in long_lengths
    ]
    amplitude_grid[.!valid_grid] .= NaN

    heat = heatmap(
        long_lengths,
        FREQUENCY_MAP_HZ ./ 1.0e3,
        amplitude_grid;
        xlabel="Corrugated block length, mm",
        ylabel="Frequency, kHz",
        title="Long N=1 branch: calibrated |H(L,f)|",
        color=:viridis,
        colorbar_title="|H|",
        size=(900, 650),
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    heat_path = joinpath(PILOT_ROOT, "figures", "long_branch_frequency_map.png")
    savefig(heat, heat_path)

    curves = plot(
        xlabel="Frequency, kHz",
        ylabel="|H|",
        title="Frequency response of the long N=1 branch",
        gridalpha=0.25,
        size=(950, 650),
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    for length_mm in long_lengths
        current = sort(filter(row -> row.length_mm == length_mm && row.valid, rows); by=row -> row.sampled_frequency_hz)
        plot!(
            curves,
            getproperty.(current, :sampled_frequency_hz) ./ 1.0e3,
            getproperty.(current, :amplitude);
            linewidth=2,
            label="L=$(length_mm) mm",
        )
    end
    curves_path = joinpath(PILOT_ROOT, "figures", "long_branch_frequency_curves.png")
    savefig(curves, curves_path)
    println("[+] Long-branch frequency map: $heat_path")
    println("[+] Long-branch frequency curves: $curves_path")
end

function run_analysis_stage()
    rows = [analyze_case(case) for case in sample_cases()]
    frequency_cases = filter(
        case -> case.profile.periods == 1 &&
                case.geometry.length_mm >= minimum(EXTENDED_LENGTH_VALUES_MM),
        sample_cases(),
    )
    frequency_rows = reduce(vcat, analyze_frequency_map.(frequency_cases))
    csv_path = joinpath(PILOT_ROOT, "period_length_pilot.csv")
    frequency_csv_path = joinpath(PILOT_ROOT, "long_branch_frequency_map.csv")
    jld_path = joinpath(PILOT_ROOT, "period_length_pilot.jld2")
    mkpath(PILOT_ROOT)
    write_csv(csv_path, rows)
    write_csv(frequency_csv_path, frequency_rows)
    JLD2.jldsave(
        jld_path;
        format_version=1,
        frequency_hz=PILOT_FREQUENCY_HZ,
        amplitude_mm=PILOT_AMPLITUDE_MM,
        period_values=PERIOD_VALUES,
        length_values_mm=LENGTH_VALUES_MM,
        rows,
    )
    save_parameter_plot(rows)
    save_complex_plot(rows)
    save_frequency_map(frequency_rows)
    println("[+] Period/length pilot: $csv_path")
    println("[+] Long-branch frequency data: $frequency_csv_path")
    println("[+] Machine-readable pilot: $jld_path")
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
        error("unknown stage: $stage")
    end
end

function list_pilot()
    println("Pilot root: $PILOT_ROOT")
    println("frequency_hz: $PILOT_FREQUENCY_HZ")
    for case in all_cases()
        println(
            "  $(case.id): role=$(case.role), N=$(case.profile.periods), ",
            "L=$(case.geometry.length_mm) mm, ",
            "gap=$(minimum_gap_mm(case.profile, case.geometry)) mm",
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
        error("usage: julia run_period_length_pilot.jl [--list | --stage=mesh|convert|solve|analyze]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
