module RunOperatorO5AFreshTest

using JLD2
using LinearAlgebra
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o3_f1c_stochastic.jl"))

using .RunOperatorO3F1CStochastic

export run_o5a_fresh_test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o5a_smooth_projection_v0",
    "fresh_test_v0",
)
const F1C_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "fullwave_f1c_stochastic_v0",
    "f1c_stochastic.jld2",
)
const TEST_CONTEXT = :o5a_test_1
const TEST_H0_SEED = UInt64(0x05a59001)
const TEST_SIGNAL_SEED = UInt64(0x05a59002)
const TEST_NOISE_SEED = UInt64(0x05a59003)
const REFERENCE_ALPHA = 3.0

function dct_basis(count)
    basis = zeros(Float64, count, count)
    basis[:, 1] .= inv(sqrt(count))
    for mode in 2:count
        basis[:, mode] .= sqrt(2 / count) .* cos.(
            (mode - 1) .* ((1:count) .- 0.5) .* pi / count,
        )
    end
    basis
end

function smooth_definition(channel_count)
    basis = dct_basis(17)
    (
        id=Symbol("o5a_dct_k$(channel_count)_q1"),
        kind=:static,
        electronic_channels=channel_count,
        physical_states=1,
        transmit_events=channel_count,
        indices=Int[],
        states=[basis[:, 1:channel_count]],
    )
end

function evaluate_definition(
    definition, electronics_sigma, validation_h0_specs, validation_signal_specs,
    test_h0_specs, test_signal_specs,
)
    model = RunOperatorO3F1CStochastic.build_architecture_model(definition)
    factor = RunOperatorO3F1CStochastic.noise_factor(model)
    validation_base = RunOperatorO3F1CStochastic.load_architecture_scenes(
        definition, [:validation_1],
    )
    validation_h0 = RunOperatorO3F1CStochastic.score_components(
        model, validation_base, validation_h0_specs, electronics_sigma,
    )
    validation_signal = RunOperatorO3F1CStochastic.score_components(
        model, validation_base, validation_signal_specs, electronics_sigma,
    )
    validation_statistics = RunOperatorO3F1CStochastic.h0_statistics(
        Xoshiro(
            RunOperatorO3F1CStochastic.SEEDS.validation_noise +
            RunOperatorO3F1CStochastic.architecture_offset(definition.id),
        ),
        validation_h0, factor, RunOperatorO3F1CStochastic.H0_VALIDATION_TRIALS,
    )
    threshold = quantile(validation_statistics, 0.99)
    validation_metrics = [
        RunOperatorO3F1CStochastic.signal_metrics(
            Xoshiro(
                RunOperatorO3F1CStochastic.SEEDS.validation_noise +
                RunOperatorO3F1CStochastic.architecture_offset(definition.id) + UInt64(index),
            ),
            validation_signal, factor, threshold, amplitude,
        )
        for (index, amplitude) in enumerate(RunOperatorO3F1CStochastic.AMPLITUDES)
    ]

    test_base = RunOperatorO3F1CStochastic.load_architecture_scenes(
        definition, [TEST_CONTEXT],
    )
    test_h0 = RunOperatorO3F1CStochastic.score_components(
        model, test_base, test_h0_specs, electronics_sigma,
    )
    test_signal = RunOperatorO3F1CStochastic.score_components(
        model, test_base, test_signal_specs, electronics_sigma,
    )
    test_statistics = RunOperatorO3F1CStochastic.h0_statistics(
        Xoshiro(TEST_NOISE_SEED + RunOperatorO3F1CStochastic.architecture_offset(definition.id)),
        test_h0, factor, RunOperatorO3F1CStochastic.H0_TEST_TRIALS,
    )
    false_alarms = count(>(threshold), test_statistics)
    test_metrics = [
        RunOperatorO3F1CStochastic.signal_metrics(
            Xoshiro(
                TEST_NOISE_SEED + RunOperatorO3F1CStochastic.architecture_offset(definition.id) +
                UInt64(100 + index),
            ),
            test_signal, factor, threshold, amplitude,
        )
        for (index, amplitude) in enumerate(RunOperatorO3F1CStochastic.AMPLITUDES)
    ]
    (
        model,
        threshold,
        observed_pfa=false_alarms / RunOperatorO3F1CStochastic.H0_TEST_TRIALS,
        pfa_upper_95=RunOperatorO3F1CStochastic.wilson_upper_bound(
            false_alarms, RunOperatorO3F1CStochastic.H0_TEST_TRIALS,
        ),
        validation_metrics,
        test_metrics,
        alpha90=RunOperatorO3F1CStochastic.crossing_amplitude(
            RunOperatorO3F1CStochastic.AMPLITUDES, getproperty.(test_metrics, :pd),
        ),
        alpha90_localized=RunOperatorO3F1CStochastic.crossing_amplitude(
            RunOperatorO3F1CStochastic.AMPLITUDES, getproperty.(test_metrics, :pdl),
        ),
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

function run_o5a_fresh_test()
    mkpath(OUTPUT_ROOT)
    f1c = JLD2.load(F1C_PATH)
    electronics_sigma = Float64(f1c["electronics_sigma"])
    frozen = RunOperatorO3F1CStochastic.RunOperatorO3FullWaveCarrier.frozen_definitions()
    a_scan = only(filter(definition -> definition.id == :uniform_a_scan, frozen))
    fmc = only(filter(definition -> definition.id == :paut_fmc_oracle, frozen))
    definitions = [a_scan, smooth_definition(3), smooth_definition(4), fmc]

    validation_h0_specs = RunOperatorO3F1CStochastic.nuisance_specs(
        RunOperatorO3F1CStochastic.SEEDS.validation_h0,
        [:validation_1], RunOperatorO3F1CStochastic.VALIDATION_H0_REPLICATES;
        with_targets=false,
    )
    validation_signal_specs = RunOperatorO3F1CStochastic.nuisance_specs(
        RunOperatorO3F1CStochastic.SEEDS.validation_signal,
        [:validation_1], RunOperatorO3F1CStochastic.VALIDATION_SIGNAL_REPLICATES;
        with_targets=true,
    )
    test_h0_specs = RunOperatorO3F1CStochastic.nuisance_specs(
        TEST_H0_SEED, [TEST_CONTEXT], RunOperatorO3F1CStochastic.TEST_H0_REPLICATES;
        with_targets=false,
    )
    test_signal_specs = RunOperatorO3F1CStochastic.nuisance_specs(
        TEST_SIGNAL_SEED, [TEST_CONTEXT], RunOperatorO3F1CStochastic.TEST_SIGNAL_REPLICATES;
        with_targets=true,
    )
    evaluations = Dict{Symbol,Any}()
    for definition in definitions
        println("[+] O5a fresh test $(definition.id)")
        evaluations[definition.id] = evaluate_definition(
            definition, electronics_sigma,
            validation_h0_specs, validation_signal_specs,
            test_h0_specs, test_signal_specs,
        )
    end
    reference_index = only(findall(
        ==(REFERENCE_ALPHA), RunOperatorO3F1CStochastic.AMPLITUDES,
    ))
    a_alpha90 = evaluations[:uniform_a_scan].alpha90
    fmc_evaluation = evaluations[:paut_fmc_oracle]
    fmc_alpha90 = fmc_evaluation.alpha90
    fmc_reference = fmc_evaluation.test_metrics[reference_index]
    fmc_model = fmc_evaluation.model
    positive_control_passed =
        fmc_model.stable_rank >= 3.0 &&
        fmc_model.maximum_template_coherence <= 0.80 &&
        isfinite(a_alpha90) && fmc_alpha90 <= 0.80a_alpha90 &&
        fmc_reference.pdl >= 0.90 && fmc_evaluation.pfa_upper_95 <= 0.0125
    valid_gap = isfinite(a_alpha90) && isfinite(fmc_alpha90) && a_alpha90 > fmc_alpha90

    summary_rows = NamedTuple[]
    target_rows = NamedTuple[]
    for definition in definitions
        evaluation = evaluations[definition.id]
        reference = evaluation.test_metrics[reference_index]
        gap_closure = valid_gap ?
            (a_alpha90 - evaluation.alpha90) / (a_alpha90 - fmc_alpha90) : NaN
        candidate = definition.id in (:o5a_dct_k3_q1, :o5a_dct_k4_q1)
        gate_passed = positive_control_passed && candidate &&
            evaluation.pfa_upper_95 <= 0.0125 && reference.pd >= 0.90 &&
            reference.pdl >= 0.80 && reference.minimum_case_pd >= 0.75 &&
            isfinite(gap_closure) && gap_closure >= 0.50
        push!(summary_rows, (
            architecture_id=definition.id,
            kind=definition.kind,
            electronic_channels=definition.electronic_channels,
            physical_states=definition.physical_states,
            transmit_events=definition.transmit_events,
            nuisance_rank=evaluation.model.nuisance_rank,
            template_stable_rank=evaluation.model.stable_rank,
            maximum_template_coherence=evaluation.model.maximum_template_coherence,
            observed_pfa=evaluation.observed_pfa,
            pfa_upper_95=evaluation.pfa_upper_95,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            minimum_case_pd_reference=reference.minimum_case_pd,
            minimum_case_pdl_reference=reference.minimum_case_pdl,
            alpha90=evaluation.alpha90,
            alpha90_localized=evaluation.alpha90_localized,
            a_scan_to_fmc_gap_closure=gap_closure,
            positive_control_passed,
            o5a_gate_passed=gate_passed,
        ))
        for target_index in eachindex(RunOperatorO3F1CStochastic.TARGET_IDS)
            push!(target_rows, (
                architecture_id=definition.id,
                target_id=RunOperatorO3F1CStochastic.TARGET_IDS[target_index],
                target_y_mm=RunOperatorO3F1CStochastic.TARGET_Y_MM[target_index],
                pd_reference=reference.case_pd[target_index],
                pdl_reference=reference.case_pdl[target_index],
            ))
        end
    end
    passing = filter(row -> row.o5a_gate_passed, summary_rows)
    priority = Dict(:o5a_dct_k3_q1 => 1, :o5a_dct_k4_q1 => 2)
    selected_id = isempty(passing) ? :none :
        first(sort(passing; by=row -> priority[row.architecture_id])).architecture_id
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "summary.csv"), summary_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "target_metrics.csv"), target_rows)
    jldsave(
        joinpath(OUTPUT_ROOT, "o5a_fresh_test.jld2");
        format_version=1,
        test_context=TEST_CONTEXT,
        test_h0_seed=TEST_H0_SEED,
        test_signal_seed=TEST_SIGNAL_SEED,
        test_noise_seed=TEST_NOISE_SEED,
        reference_alpha=REFERENCE_ALPHA,
        electronics_sigma,
        positive_control_passed,
        selected_id,
        summary_rows,
        target_rows,
    )
    readme_path = joinpath(OUTPUT_ROOT, "README.md")
    open(readme_path, "w") do io
        println(io, "# O5a fresh full-wave test")
        println(io)
        println(io, "> Candidate set and priority were frozen on train/validation before `o5a_test_1` was solved.")
        println(io)
        println(io, "Positive control: **$(positive_control_passed ? "PASS" : "FAIL")**. Selected geometry-aware operator: **`$(selected_id)`**.")
        println(io)
        println(io, "| Architecture | Channels | States | Events | Pfa upper | Pd@3 | Pdl@3 | Min Pd | alpha90 | Gap closure | Gate |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
        for row in summary_rows
            closure = isfinite(row.a_scan_to_fmc_gap_closure) ? @sprintf("%.3f", row.a_scan_to_fmc_gap_closure) : "n/a"
            alpha = isfinite(row.alpha90) ? @sprintf("%.3f", row.alpha90) : "inf"
            println(io, @sprintf(
                "| `%s` | %d | %d | %d | %.3f%% | %.3f | %.3f | %.3f | %s | %s | %s |",
                row.architecture_id, row.electronic_channels, row.physical_states,
                row.transmit_events, 100row.pfa_upper_95, row.pd_reference,
                row.pdl_reference, row.minimum_case_pd_reference, alpha, closure,
                row.o5a_gate_passed ? "pass" : "fail",
            ))
        end
    end
    println("[+] O5a fresh-test report: $readme_path")
    summary_rows
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o5a_fresh_test()

end
