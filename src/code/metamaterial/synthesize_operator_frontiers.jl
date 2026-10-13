using Plots

const SYNTHESIS_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const MOTION_ROOT = joinpath(
    SYNTHESIS_ROOT, "results", "operator_space", "motion_frontier_v0",
)
const STATIC_ROOT = joinpath(
    SYNTHESIS_ROOT, "results", "operator_space", "static_frontier_v0",
)
const COMBINED_ROOT = joinpath(
    SYNTHESIS_ROOT, "results", "operator_space", "combined_frontier_v0",
)

function read_simple_csv(path)
    lines = readlines(path)
    isempty(lines) && error("empty CSV: $path")
    headers = Symbol.(split(first(lines), ','))
    [
        Dict(header => value for (header, value) in zip(headers, split(line, ','; keepempty=true)))
        for line in Iterators.drop(lines, 1)
        if !isempty(strip(line))
    ]
end

parse_float(value) = value == "NaN" ? NaN : parse(Float64, value)
parse_int(value) = parse(Int, value)
parse_bool(value) = value == "true"

function csv_value(value)
    value isa AbstractFloat && return isfinite(value) ? string(value) : "NaN"
    value isa Symbol && return String(value)
    value isa Bool && return value ? "true" : "false"
    string(value)
end

function write_csv(path, columns, rows)
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function static_rows()
    raw = read_simple_csv(joinpath(STATIC_ROOT, "frontier_summary.csv"))
    [(
        architecture_id=Symbol(row[:architecture_id]),
        branch=row[:role] in ("baseline", "oracle") ? Symbol(row[:role]) : :static,
        family=Symbol(row[:family]),
        electronic_channels=parse_int(row[:electronic_channels]),
        transmit_events=parse_int(row[:transmit_events]),
        passive_states=parse_int(row[:physical_states]),
        state_switches=parse_int(row[:state_switches]),
        mechanical_positions=1,
        mechanical_moves=0,
        aperture_span_lambda=0.0,
        calibration_complex_coefficients=parse_int(row[:calibration_complex_coefficients]),
        pfa_upper_95=parse_float(row[:pfa_upper_95]),
        pfa_controlled=parse_float(row[:pfa_upper_95]) <= 0.0125,
        pd_reference=parse_float(row[:pd_reference]),
        pdl_reference=parse_float(row[:pdl_reference]),
        alpha90=parse_float(row[:alpha90]),
        reference_gate_passed=parse_bool(row[:intrinsic_gate_passed]),
        source_campaign=:static_frontier_v0,
    ) for row in raw]
end

function motion_rows()
    raw = read_simple_csv(joinpath(MOTION_ROOT, "motion_summary.csv"))
    comparator_ids = Set((
        "a_scan_uniform",
        "four_channel_two_static_states",
        "paut_fmc_oracle",
    ))
    [(
        architecture_id=Symbol(row[:architecture_id]),
        branch=:motion,
        family=Symbol(row[:strategy]),
        electronic_channels=parse_int(row[:electronic_channels]),
        transmit_events=parse_int(row[:transmit_events]),
        passive_states=parse_int(row[:physical_states]),
        state_switches=max(parse_int(row[:physical_states]) - 1, 0),
        mechanical_positions=parse_int(row[:mechanical_positions]),
        mechanical_moves=parse_int(row[:mechanical_moves]),
        aperture_span_lambda=parse_float(row[:aperture_span_lambda]),
        calibration_complex_coefficients=parse_int(row[:complex_features]),
        pfa_upper_95=parse_float(row[:pfa_upper_95]),
        pfa_controlled=parse_bool(row[:pfa_controlled]),
        pd_reference=parse_float(row[:pd_reference]),
        pdl_reference=parse_float(row[:pdl_reference]),
        alpha90=parse_float(row[:alpha90]),
        reference_gate_passed=parse_bool(row[:pfa_controlled]) &&
            parse_float(row[:pd_reference]) >= 0.90 &&
            parse_float(row[:pdl_reference]) >= 0.80,
        source_campaign=:motion_frontier_v0,
    ) for row in raw if !(row[:architecture_id] in comparator_ids)]
end

function time_value(row, scenario)
    scenario == :events_only && return Float64(row.transmit_events)
    scenario == :balanced_overhead && return Float64(
        row.transmit_events + row.state_switches + row.mechanical_moves +
        row.aperture_span_lambda,
    )
    scenario == :motion_expensive && return Float64(
        row.transmit_events + row.state_switches + 4 * row.mechanical_moves +
        2 * row.aperture_span_lambda,
    )
    scenario == :state_expensive && return Float64(
        row.transmit_events + 4 * row.state_switches + row.mechanical_moves +
        row.aperture_span_lambda,
    )
    throw(ArgumentError("unknown scenario: $scenario"))
end

function dominated(row, other, row_time, other_time)
    objectives = (
        row.electronic_channels,
        row_time,
        row.calibration_complex_coefficients,
        row.alpha90,
    )
    other_objectives = (
        other.electronic_channels,
        other_time,
        other.calibration_complex_coefficients,
        other.alpha90,
    )
    all(other_objectives[index] <= objectives[index] for index in eachindex(objectives)) &&
        any(other_objectives[index] < objectives[index] for index in eachindex(objectives))
end

function pareto_flag(row, rows, scenario; require_reference_gate)
    eligible = filter(other ->
        other.pfa_controlled && isfinite(other.alpha90) &&
        (!require_reference_gate || other.reference_gate_passed),
        rows,
    )
    row in eligible || return false
    row_time = time_value(row, scenario)
    !any(other -> other !== row && dominated(
        row, other, row_time, time_value(other, scenario),
    ), eligible)
end

function resource_frontiers(rows)
    scenarios = (:events_only, :balanced_overhead, :motion_expensive, :state_expensive)
    [(
        scenario,
        architecture_id=row.architecture_id,
        branch=row.branch,
        family=row.family,
        electronic_channels=row.electronic_channels,
        normalized_time=time_value(row, scenario),
        transmit_events=row.transmit_events,
        passive_states=row.passive_states,
        mechanical_positions=row.mechanical_positions,
        calibration_complex_coefficients=row.calibration_complex_coefficients,
        alpha90=row.alpha90,
        pfa_controlled=row.pfa_controlled,
        reference_gate_passed=row.reference_gate_passed,
        pareto_pfa=pareto_flag(row, rows, scenario; require_reference_gate=false),
        pareto_reference_gate=pareto_flag(row, rows, scenario; require_reference_gate=true),
    ) for scenario in scenarios for row in rows]
end

function best_at_budget(rows, branch, budget; require_gate)
    eligible = filter(row ->
        row.branch == branch && row.transmit_events <= budget && row.pfa_controlled &&
        isfinite(row.alpha90) && (!require_gate || row.reference_gate_passed),
        rows,
    )
    isempty(eligible) && return nothing
    eligible[argmin(getproperty.(eligible, :alpha90))]
end

function event_budget_envelope(rows)
    output = NamedTuple[]
    for budget in (1, 2, 4, 8, 9, 12, 17), branch in (:static, :motion)
        for require_gate in (false, true)
            best = best_at_budget(rows, branch, budget; require_gate)
            push!(output, (
                event_budget=budget,
                branch,
                eligibility=require_gate ? :reference_gate : :pfa_only,
                architecture_id=isnothing(best) ? :none : best.architecture_id,
                electronic_channels=isnothing(best) ? 0 : best.electronic_channels,
                alpha90=isnothing(best) ? NaN : best.alpha90,
            ))
        end
    end
    output
end

function plot_combined(rows, frontier_rows)
    colors = Dict(:static => :steelblue, :motion => :darkorange, :baseline => :black,
                  :oracle => :firebrick)
    quality = plot(
        xlabel="transmit events",
        ylabel="alpha90 (lower is better)",
        title="Fresh combined operator frontier",
        legend=:topright,
        size=(950, 620),
    )
    for branch in (:static, :motion, :baseline, :oracle)
        selected = filter(row -> row.branch == branch && isfinite(row.alpha90), rows)
        isempty(selected) && continue
        scatter!(
            quality,
            getproperty.(selected, :transmit_events),
            getproperty.(selected, :alpha90);
            marker=[row.reference_gate_passed ? :circle : :xcross for row in selected],
            markersize=[row.reference_gate_passed ? 7 : 5 for row in selected],
            color=colors[branch],
            alpha=0.75,
            label=String(branch),
        )
    end
    savefig(quality, joinpath(COMBINED_ROOT, "combined_quality_frontier.png"))

    figures = Plots.Plot[]
    for scenario in (:events_only, :balanced_overhead, :motion_expensive)
        selected = filter(row ->
            row.scenario == scenario && row.pfa_controlled && isfinite(row.alpha90),
            frontier_rows,
        )
        figure = plot(
            xlabel="normalized time",
            ylabel="alpha90",
            title=replace(String(scenario), '_' => ' '),
            legend=false,
        )
        for branch in (:static, :motion)
            branch_rows = filter(row -> row.branch == branch, selected)
            scatter!(
                figure,
                getproperty.(branch_rows, :normalized_time),
                getproperty.(branch_rows, :alpha90);
                color=colors[branch],
                alpha=0.55,
                markersize=4,
            )
            pareto = filter(row -> row.branch == branch && row.pareto_pfa, branch_rows)
            scatter!(
                figure,
                getproperty.(pareto, :normalized_time),
                getproperty.(pareto, :alpha90);
                marker=:star5,
                color=colors[branch],
                markersize=9,
            )
        end
        push!(figures, figure)
    end
    savefig(
        plot(figures...; layout=(1, 3), size=(1550, 500), margin=6 * Plots.mm),
        joinpath(COMBINED_ROOT, "combined_cost_sensitivity.png"),
    )
end

function main()
    mkpath(COMBINED_ROOT)
    rows = vcat(static_rows(), motion_rows())
    frontiers = resource_frontiers(rows)
    envelopes = event_budget_envelope(rows)
    write_csv(
        joinpath(COMBINED_ROOT, "combined_architectures.csv"),
        propertynames(first(rows)), rows,
    )
    write_csv(
        joinpath(COMBINED_ROOT, "combined_resource_frontiers.csv"),
        propertynames(first(frontiers)), frontiers,
    )
    write_csv(
        joinpath(COMBINED_ROOT, "event_budget_envelope.csv"),
        propertynames(first(envelopes)), envelopes,
    )
    plot_combined(rows, frontiers)
    println("[+] combined $(length(rows)) architecture points")
    println("[+] wrote synthesis to $COMBINED_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
