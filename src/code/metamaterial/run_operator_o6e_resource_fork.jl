module RunOperatorO6EResourceFork

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o6d_isolated_bus_bound.jl"))
include(joinpath(@__DIR__, "run_operator_f1d_expanded_calibration.jl"))

using .RunOperatorO6DIsolatedBusBound
using .RunOperatorF1DExpandedCalibration

export run_o6e_resource_fork

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const O6D_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6d_isolated_bus_bound_v0",
    "o6d_isolated_bus_bound.jld2",
)
const F1D_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "fullwave_f1d_expanded_v0",
    "f1d_expanded.jld2",
)
const MOTION_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "motion_frontier_v0",
    "motion_frontier.jld2",
)
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6e_resource_fork_v0",
)
const TOLERANCE_SEED = UInt64(0x06e50117)
const TOLERANCE_DRAWS = 50_000
const TOLERANCE_SCALES = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5)
const LOG_AMPLITUDE_RMS_REFERENCE = 0.03
const DELAY_RMS_REFERENCE_S = 82.6e-9
const PASS_PROBABILITY_LOWER_GATE = 0.99
const ELEMENT_COUNT = 17
const PORT_COUNT = 3

function wilson_lower_bound(successes, trials; z=1.6448536269514722)
    trials > 0 || return 0.0
    proportion = successes / trials
    denominator = 1 + z^2 / trials
    center = proportion + z^2 / (2trials)
    spread = z * sqrt(proportion * (1 - proportion) / trials + z^2 / (4trials^2))
    max(0.0, (center - spread) / denominator)
end

function delay_rms(delay, target)
    active = abs.(target) .>= 0.10maximum(abs, target)
    weights = abs2.(target[active])
    values = delay[active]
    common = sum(weights .* values) / sum(weights)
    sqrt(sum(weights .* (values .- common) .^ 2) / sum(weights))
end

function summarize_draws(
    scale, base_transfer, target, frequencies_hz, amplitude_standard, delay_standard,
)
    draws = size(amplitude_standard, 3)
    passed = falses(draws)
    minimum_fidelity = zeros(draws)
    maximum_coherence = zeros(draws)
    maximum_imbalance_db = zeros(draws)
    maximum_amplitude_rms = zeros(draws)
    maximum_delay_rms_ns = zeros(draws)
    minimum_useful_amplitude = zeros(draws)
    maximum_uncontracted_sigma = zeros(draws)
    contractive_scale = zeros(draws)
    for draw in 1:draws
        log_amplitude_error = scale * LOG_AMPLITUDE_RMS_REFERENCE .* @view(
            amplitude_standard[:, :, draw],
        )
        tap_delay = scale * DELAY_RMS_REFERENCE_S .* @view(
            delay_standard[:, :, draw],
        )
        magnitude = base_transfer .* exp.(log_amplitude_error)
        raw = [
            complex.(magnitude) .* exp.(-2pi * im * frequency_hz .* tap_delay)
            for frequency_hz in frequencies_hz
        ]
        maximum_sigma = maximum(opnorm.(raw))
        contraction = inv(max(1.0, maximum_sigma))
        metrics = [
            RunOperatorO6DIsolatedBusBound.transfer_metrics(
                contraction .* transfer, target,
            ) for transfer in raw
        ]
        delay_value = delay_rms(tap_delay, target)
        minimum_fidelity[draw] = minimum(getproperty.(metrics, :fidelity))
        maximum_coherence[draw] = maximum(getproperty.(metrics, :coherence))
        maximum_imbalance_db[draw] = maximum(getproperty.(metrics, :imbalance_db))
        maximum_amplitude_rms[draw] = maximum(getproperty.(metrics, :amplitude_rms))
        maximum_delay_rms_ns[draw] = 1e9delay_value
        minimum_useful_amplitude[draw] = minimum(
            getproperty.(metrics, :useful_amplitude),
        )
        maximum_uncontracted_sigma[draw] = maximum_sigma
        contractive_scale[draw] = contraction
        passed[draw] =
            minimum_fidelity[draw] >= RunOperatorO6DIsolatedBusBound.FIDELITY_GATE &&
            maximum_coherence[draw] <= RunOperatorO6DIsolatedBusBound.COHERENCE_GATE &&
            maximum_imbalance_db[draw] <= RunOperatorO6DIsolatedBusBound.IMBALANCE_GATE_DB &&
            maximum_amplitude_rms[draw] <= RunOperatorO6DIsolatedBusBound.AMPLITUDE_RMS_GATE &&
            delay_value <= RunOperatorO6DIsolatedBusBound.DELAY_RMS_GATE_S &&
            minimum_useful_amplitude[draw] >= RunOperatorO6DIsolatedBusBound.USEFUL_AMPLITUDE_GATE
    end
    successes = count(passed)
    (
        scale,
        log_amplitude_rms=scale * LOG_AMPLITUDE_RMS_REFERENCE,
        delay_rms_ns=1e9scale * DELAY_RMS_REFERENCE_S,
        draws,
        successes,
        observed_joint_pass_probability=successes / draws,
        wilson_lower_95=wilson_lower_bound(successes, draws),
        q01_minimum_fidelity=quantile(minimum_fidelity, 0.01),
        q99_maximum_coherence=quantile(maximum_coherence, 0.99),
        q99_maximum_imbalance_db=quantile(maximum_imbalance_db, 0.99),
        q99_maximum_amplitude_rms=quantile(maximum_amplitude_rms, 0.99),
        q99_maximum_delay_rms_ns=quantile(maximum_delay_rms_ns, 0.99),
        q01_minimum_useful_amplitude=quantile(minimum_useful_amplitude, 0.01),
        q99_uncontracted_sigma=quantile(maximum_uncontracted_sigma, 0.99),
        q01_contractive_scale=quantile(contractive_scale, 0.01),
        probability_gate_passed=wilson_lower_bound(successes, draws) >=
                                PASS_PROBABILITY_LOWER_GATE,
    )
end

function tolerance_threshold(
    base_transfer, target, frequencies_hz, amplitude_standard, delay_standard,
    frontier_rows,
)
    passing = filter(row -> row.probability_gate_passed, frontier_rows)
    isempty(passing) && error("zero-error beamformer does not pass the probability gate")
    lower = maximum(getproperty.(passing, :scale))
    failing = filter(row -> !row.probability_gate_passed && row.scale > lower, frontier_rows)
    upper = isempty(failing) ? 2maximum(getproperty.(frontier_rows, :scale)) :
            minimum(getproperty.(failing, :scale))
    if isempty(failing)
        while summarize_draws(
            upper, base_transfer, target, frequencies_hz,
            amplitude_standard, delay_standard,
        ).probability_gate_passed
            lower = upper
            upper *= 2
        end
    end
    for _ in 1:12
        middle = (lower + upper) / 2
        result = summarize_draws(
            middle, base_transfer, target, frequencies_hz,
            amplitude_standard, delay_standard,
        )
        if result.probability_gate_passed
            lower = middle
        else
            upper = middle
        end
    end
    summarize_draws(
        lower, base_transfer, target, frequencies_hz,
        amplitude_standard, delay_standard,
    )
end

function open_tap_audit(base_transfer, target)
    rows = NamedTuple[]
    for linear_index in eachindex(base_transfer)
        abs(base_transfer[linear_index]) > 1e-12 || continue
        opened = copy(base_transfer)
        opened[linear_index] = 0.0
        calibrated = RunOperatorO6DIsolatedBusBound.balance_and_contract(opened)
        metrics = RunOperatorO6DIsolatedBusBound.transfer_metrics(calibrated, target)
        gates = RunOperatorO6DIsolatedBusBound.electrical_gates(metrics)
        cartesian = CartesianIndices(base_transfer)[linear_index]
        push!(rows, (
            element_index=cartesian[1],
            port_index=cartesian[2],
            nominal_tap=base_transfer[linear_index],
            normalized_transfer_fidelity=metrics.fidelity,
            maximum_port_coherence=metrics.coherence,
            differential_amplitude_rms=metrics.amplitude_rms,
            useful_amplitude=metrics.useful_amplitude,
            sigma_max=metrics.sigma_max,
            passed=all(values(gates)),
        ))
    end
    rows
end

function replay_q8_f1d(state, frozen_f1d)
    runner = RunOperatorF1DExpandedCalibration
    stochastic = runner.RunOperatorO3F1CStochastic
    f1c = load(runner.F1C_PATH)
    electronics_sigma = Float64(f1c["electronics_sigma"])
    definition = (
        id=:o5a_dct_k3_q1,
        kind=:static,
        electronic_channels=3,
        physical_states=1,
        transmit_events=3,
        indices=Int[],
        states=[state],
    )
    validation_h0_specs = stochastic.nuisance_specs(
        runner.SEEDS.validation_h0,
        runner.VALIDATION_CONTEXTS,
        runner.VALIDATION_H0_REPLICATES;
        with_targets=false,
    )
    validation_signal_specs = stochastic.nuisance_specs(
        runner.SEEDS.validation_signal,
        runner.VALIDATION_CONTEXTS,
        runner.VALIDATION_SIGNAL_REPLICATES;
        with_targets=true,
    )
    test_h0_specs = stochastic.nuisance_specs(
        runner.SEEDS.test_h0,
        runner.TEST_CONTEXTS,
        runner.TEST_H0_REPLICATES;
        with_targets=false,
    )
    test_signal_specs = stochastic.nuisance_specs(
        runner.SEEDS.test_signal,
        runner.TEST_CONTEXTS,
        runner.TEST_SIGNAL_REPLICATES;
        with_targets=true,
    )
    evaluation = runner.evaluate_definition(
        definition,
        electronics_sigma,
        validation_h0_specs,
        validation_signal_specs,
        test_h0_specs,
        test_signal_specs,
    )
    reference_index = only(findall(==(runner.REFERENCE_ALPHA), stochastic.AMPLITUDES))
    reference = evaluation.test_metrics[reference_index]
    old_rows = frozen_f1d["summary_rows"]
    a_scan = only(filter(row -> row.architecture_id == :uniform_a_scan, old_rows))
    fmc = only(filter(row -> row.architecture_id == :paut_fmc_oracle, old_rows))
    ideal = only(filter(row -> row.architecture_id == :o5a_dct_k3_q1, old_rows))
    gap_closure = (a_scan.alpha90 - evaluation.alpha90) /
                  (a_scan.alpha90 - fmc.alpha90)
    gate_passed = Bool(frozen_f1d["positive_control_passed"]) &&
        evaluation.pfa_upper_95 <= 0.0125 && reference.pd >= 0.90 &&
        reference.pdl >= 0.80 && reference.minimum_case_pd >= 0.75 &&
        isfinite(gap_closure) && gap_closure >= 0.50
    (
        architecture_id="beamformer_q8_frozen_f1d_replay",
        evidence_level="fullwave_f1d_frozen_replay",
        pfa_upper_95=evaluation.pfa_upper_95,
        pd_reference=reference.pd,
        pdl_reference=reference.pdl,
        minimum_case_pd_reference=reference.minimum_case_pd,
        alpha90=evaluation.alpha90,
        alpha90_localized=evaluation.alpha90_localized,
        gap_closure,
        ideal_dct_alpha90=ideal.alpha90,
        alpha90_ratio_to_ideal_dct=evaluation.alpha90 / ideal.alpha90,
        passed=gate_passed,
    )
end

function write_namedtuple_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function resource_row_beamformer(replay)
    (
        architecture_id="beamformer_q8",
        branch="passive_beamformer",
        evidence_level=replay.evidence_level,
        electronic_channels=3,
        transmit_events=3,
        mechanical_positions=1,
        mechanical_moves=0,
        aperture_span_lambda=0.0,
        active_transducer_columns=17,
        passive_bus_taps=50,
        weight_levels_per_port=8,
        distinct_weight_ratios=17,
        calibration_sites_or_taps=50,
        time_electronics_only=3.0,
        time_one_event_per_settle=3.0,
        time_settle_plus_travel=3.0,
        time_slow_motion=3.0,
        pfa_upper_95=replay.pfa_upper_95,
        pd_reference=replay.pd_reference,
        pdl_reference=replay.pdl_reference,
        alpha90_within_evidence=replay.alpha90,
        normalized_gap_closure_within_evidence=replay.gap_closure,
        pfa_controlled=replay.pfa_upper_95 <= 0.0125,
        unresolved_burden="network parasitics; isolated-column fabrication",
    )
end

function resource_row_motion(row)
    (
        architecture_id=String(row.architecture_id),
        branch="mechanical_motion",
        evidence_level="reduced_born_o3_fresh_test",
        electronic_channels=row.electronic_channels,
        transmit_events=row.transmit_events,
        mechanical_positions=row.mechanical_positions,
        mechanical_moves=row.mechanical_moves,
        aperture_span_lambda=row.aperture_span_lambda,
        active_transducer_columns=1,
        passive_bus_taps=0,
        weight_levels_per_port=0,
        distinct_weight_ratios=0,
        calibration_sites_or_taps=row.mechanical_positions,
        time_electronics_only=row.time_electronics_only,
        time_one_event_per_settle=row.time_one_event_per_settle,
        time_settle_plus_travel=row.time_settle_plus_travel,
        time_slow_motion=row.time_slow_motion,
        pfa_upper_95=row.pfa_upper_95,
        pd_reference=row.pd_reference,
        pdl_reference=row.pdl_reference,
        alpha90_within_evidence=row.alpha90,
        normalized_gap_closure_within_evidence=row.aggregate_paut_gap_closure,
        pfa_controlled=row.pfa_controlled,
        unresolved_burden="encoder error; couplant reset; probe footprint",
    )
end

function exchange_rows(beamformer, motions)
    [(
        reference=beamformer.architecture_id,
        alternative=motion.architecture_id,
        extra_electronic_channels_beamformer=
            beamformer.electronic_channels - motion.electronic_channels,
        events_saved_by_beamformer=
            motion.transmit_events - beamformer.transmit_events,
        moves_saved_by_beamformer=
            motion.mechanical_moves - beamformer.mechanical_moves,
        extra_active_columns_beamformer=
            beamformer.active_transducer_columns - motion.active_transducer_columns,
        extra_passive_taps_beamformer=
            beamformer.passive_bus_taps - motion.passive_bus_taps,
        electronics_time_saved_by_beamformer=
            motion.time_electronics_only - beamformer.time_electronics_only,
        settle_time_saved_by_beamformer=
            motion.time_one_event_per_settle - beamformer.time_one_event_per_settle,
        travel_time_saved_by_beamformer=
            motion.time_settle_plus_travel - beamformer.time_settle_plus_travel,
        slow_motion_time_saved_by_beamformer=
            motion.time_slow_motion - beamformer.time_slow_motion,
        detection_metrics_commensurate=false,
    ) for motion in motions]
end

function run_o6e_resource_fork()
    for path in (O6D_PATH, F1D_PATH, MOTION_PATH)
        isfile(path) || error("required frozen result not found: $path")
    end
    mkpath(OUTPUT_ROOT)
    o6d = load(O6D_PATH)
    frozen_f1d = load(F1D_PATH)
    motion = load(MOTION_PATH)
    target = Float64.(o6d["target"])
    base_transfer = Float64.(o6d["quantized_transfers"][8])
    frequencies_hz = Float64.(o6d["frequencies_hz"])
    active_taps = count(abs.(base_transfer) .> 1e-12)
    active_taps == 50 || error("O6e beamformer no longer has 50 active taps")
    nominal_metrics = RunOperatorO6DIsolatedBusBound.transfer_metrics(
        base_transfer, target,
    )
    nominal_loss_amplitude_floor =
        RunOperatorO6DIsolatedBusBound.USEFUL_AMPLITUDE_GATE /
        nominal_metrics.useful_amplitude
    nominal_loss_ceiling_db = -20log10(nominal_loss_amplitude_floor)

    rng = Xoshiro(TOLERANCE_SEED)
    amplitude_standard = randn(rng, ELEMENT_COUNT, PORT_COUNT, TOLERANCE_DRAWS)
    delay_standard = randn(rng, ELEMENT_COUNT, PORT_COUNT, TOLERANCE_DRAWS)
    println("[>] O6e beamformer tolerance frontier")
    tolerance_rows = [
        summarize_draws(
            scale, base_transfer, target, frequencies_hz,
            amplitude_standard, delay_standard,
        ) for scale in TOLERANCE_SCALES
    ]
    threshold = tolerance_threshold(
        base_transfer, target, frequencies_hz,
        amplitude_standard, delay_standard, tolerance_rows,
    )
    robust_loss_amplitude_floor =
        RunOperatorO6DIsolatedBusBound.USEFUL_AMPLITUDE_GATE /
        threshold.q01_minimum_useful_amplitude
    robust_loss_ceiling_db = robust_loss_amplitude_floor <= 1 ?
        -20log10(robust_loss_amplitude_floor) : -Inf

    println("[>] O6e open-tap audit")
    open_tap_rows = open_tap_audit(base_transfer, target)
    println("[>] O6e q8 frozen F1d replay")
    replay = replay_q8_f1d(base_transfer, frozen_f1d)

    motion_rows = motion["summary_rows"]
    p12 = only(filter(
        row -> row.architecture_id == :motion_compact_center_p12_native,
        motion_rows,
    ))
    p17 = only(filter(
        row -> row.architecture_id == :motion_uniform_full_p17_a8,
        motion_rows,
    ))
    beamformer_resource = resource_row_beamformer(replay)
    motion_resources = resource_row_motion.((p12, p17))
    resource_rows = [beamformer_resource, motion_resources...]
    exchanges = exchange_rows(beamformer_resource, motion_resources)
    branch_ready = replay.passed && nominal_loss_ceiling_db > 0 &&
                   threshold.wilson_lower_95 >= PASS_PROBABILITY_LOWER_GATE

    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "tolerance_frontier.csv"), tolerance_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "tolerance_threshold.csv"), [threshold])
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "open_tap_audit.csv"), open_tap_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "f1d_q8_replay.csv"), [replay])
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "resource_frontier.csv"), resource_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "pairwise_resource_exchange.csv"), exchanges)
    jldsave(
        joinpath(OUTPUT_ROOT, "o6e_resource_fork.jld2");
        format_version=1,
        source_o6d_path=abspath(O6D_PATH),
        source_f1d_path=abspath(F1D_PATH),
        source_motion_path=abspath(MOTION_PATH),
        tolerance_seed=TOLERANCE_SEED,
        tolerance_draws=TOLERANCE_DRAWS,
        target,
        base_transfer,
        frequencies_hz,
        nominal_metrics,
        nominal_loss_amplitude_floor,
        nominal_loss_ceiling_db,
        threshold,
        robust_loss_amplitude_floor,
        robust_loss_ceiling_db,
        tolerance_rows,
        open_tap_rows,
        replay,
        resource_rows,
        exchanges,
        beamformer_branch_ready=branch_ready,
    )

    tap_passes = count(getproperty.(open_tap_rows, :passed))
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6e passive-beamformer versus motion resource fork")
        println(io)
        println(io, branch_ready ?
            "**The passive-beamformer branch clears its algebraic pre-geometry gates; neither branch dominates the full resource vector.**" :
            "**The passive-beamformer branch fails its pre-geometry gate; motion remains the surviving branch.**")
        println(io)
        println(io, "## Passive 50-tap branch")
        println(io)
        @printf(io,
            "The nominal q8 map has useful amplitude `%.6f`, leaving at most `%.3f dB` one-way common insertion loss before the O4 floor. Its frozen F1d replay %s: `Pfa95=%.3f%%`, `Pd=%.3f`, `Pdl=%.3f`, `alpha90=%.4f`, gap closure `%.3f`; alpha90 changes by `%.3f%%` from ideal DCT3.\n\n",
            nominal_metrics.useful_amplitude, nominal_loss_ceiling_db,
            replay.passed ? "passes" : "fails", 100replay.pfa_upper_95,
            replay.pd_reference, replay.pdl_reference, replay.alpha90,
            replay.gap_closure, 100(replay.alpha90_ratio_to_ideal_dct - 1),
        )
        @printf(io,
            "The largest common manufacturing-error scale with Wilson lower-95 pass probability >=99%% is `%.4f`: log-amplitude RMS `%.3f%%`, delay RMS `%.3f ns`. At that point observed pass is `%.3f%%`, lower bound `%.3f%%`, q01 useful amplitude `%.6f`, leaving `%.3f dB` additional common-loss headroom.\n\n",
            threshold.scale, 100threshold.log_amplitude_rms,
            threshold.delay_rms_ns, 100threshold.observed_joint_pass_probability,
            100threshold.wilson_lower_95, threshold.q01_minimum_useful_amplitude,
            robust_loss_ceiling_db,
        )
        reference_tolerance = only(filter(row -> row.scale == 1.0, tolerance_rows))
        @printf(io,
            "At the full O4 reference perturbation (`3%%`, `82.6 ns`) only `%.3f%%` of draws pass jointly, so the q8 network does not inherit the entire coefficient tolerance budget: fabrication and network errors together must stay below the reported scale.\n\n",
            100reference_tolerance.observed_joint_pass_probability,
        )
        println(io, "Single-open-tap audit: `$tap_passes/$(length(open_tap_rows))` faults still pass O4 after optimistic fixed port recalibration. This is diagnostic only; no fault tolerance is credited to the architecture.")
        println(io)
        println(io, "## Unsummed resource frontier")
        println(io)
        println(io, "| Branch | Ch | Events | Positions / moves | T electronics / settle / travel / slow | Columns | Passive taps | Gap closure* |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in resource_rows
            @printf(io,
                "| `%s` | %d | %d | %d / %d | %.1f / %.1f / %.1f / %.1f | %d | %d | %.3f |\n",
                row.architecture_id, row.electronic_channels, row.transmit_events,
                row.mechanical_positions, row.mechanical_moves,
                row.time_electronics_only, row.time_one_event_per_settle,
                row.time_settle_plus_travel, row.time_slow_motion,
                row.active_transducer_columns, row.passive_bus_taps,
                row.normalized_gap_closure_within_evidence,
            )
        end
        println(io)
        println(io, "`*` Gap closures are normalized inside different evidence campaigns: full-wave frozen F1d for beamformer q8, reduced fresh Born O3 for motion. They are not a shared ROC number.")
        println(io)
        println(io, "Relative to P12/P17, the static beamformer spends two extra electronic channels and 50 calibrated passive taps, but saves `9/14` transmit events and `11/16` moves. Under the existing time models it saves `9--64` normalized units versus P12 and `14--94` versus P17. P12/P17 avoid the distributed network but retain unresolved encoder, couplant-reset and probe-footprint risks.")
        println(io)
        println(io, "## Decision")
        println(io)
        println(io, "Both branches remain open. The next paired physical-risk stage is: one realizable passive-network topology with measured/estimated parasitics for q8, and one positioning/contact perturbation campaign for P12/P17. A scalar winner remains forbidden until those missing coordinates exist.")
    end

    tolerance_plot = plot(
        getproperty.(tolerance_rows, :scale),
        100 .* getproperty.(tolerance_rows, :wilson_lower_95);
        marker=:circle, linewidth=2, label="Wilson lower 95%",
        xlabel="O4 tolerance scale", ylabel="joint pass probability, %",
        title="O6e passive-beamformer tolerance",
    )
    hline!(tolerance_plot, [100PASS_PROBABILITY_LOWER_GATE];
        linestyle=:dash, color=:black, label="99% gate")
    vline!(tolerance_plot, [threshold.scale];
        linestyle=:dot, color=:red, label="threshold")
    savefig(tolerance_plot, joinpath(OUTPUT_ROOT, "beamformer_tolerance_frontier.png"))

    resource_plot = scatter(
        getproperty.(resource_rows, :time_settle_plus_travel),
        getproperty.(resource_rows, :normalized_gap_closure_within_evidence);
        markersize=8, xlabel="normalized settle+travel time",
        ylabel="within-evidence gap closure", label="architectures",
        title="O6e resource fork (evidence levels differ)",
    )
    for row in resource_rows
        annotate!(resource_plot,
            row.time_settle_plus_travel,
            row.normalized_gap_closure_within_evidence,
            text(row.architecture_id, 7, :left),
        )
    end
    savefig(resource_plot, joinpath(OUTPUT_ROOT, "resource_fork.png"))
    println("[+] O6e resource fork complete")
    println("[+] $(joinpath(OUTPUT_ROOT, "README.md"))")
    (; branch_ready, threshold, replay, resource_rows, exchanges)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6e_resource_fork()

end
