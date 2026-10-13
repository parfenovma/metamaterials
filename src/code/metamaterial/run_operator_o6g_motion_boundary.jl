module RunOperatorO6GMotionBoundary

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o6f_motion_risk.jl"))
using .RunOperatorO6FMotionRisk

export run_o6g_motion_boundary

const M = RunOperatorO6FMotionRisk
const O6G_MOTION_SOURCE = joinpath(
    M.PROJECT_ROOT, "results", "operator_space", "o6f_motion_risk_v0",
    "o6f_motion_risk.jld2",
)
const O6G_MOTION_OUTPUT_ROOT = joinpath(
    M.PROJECT_ROOT, "results", "operator_space", "o6g_motion_boundary_v0",
)
const COUPLED_COARSE_SCALES = (0.0, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0)
const BISECTION_STEPS = 12
const MAXIMUM_SEARCH_SCALE = 64.0
const SCENARIOS = (
    (id=:coupled, position=1.0, amplitude=1.0, delay=1.0),
    (id=:position_only, position=1.0, amplitude=0.0, delay=0.0),
    (id=:amplitude_only, position=0.0, amplitude=1.0, delay=0.0),
    (id=:delay_only, position=0.0, amplitude=0.0, delay=1.0),
)

function setup_frozen_replay()
    isfile(M.O6F_MOTION_SOURCE) || error("motion frontier source is missing")
    isfile(O6G_MOTION_SOURCE) || error("O6f motion result is missing")
    saved_motion = load(M.O6F_MOTION_SOURCE)
    saved_o6f = load(O6G_MOTION_SOURCE)
    config = M.SurrogateConfig()
    hypotheses = M.localization_hypotheses()
    architectures = M.scan_architectures(config)
    split_libraries = Dict(
        split => M.fresh_nuisance_libraries(config, split)
        for split in (:train, :test)
    )
    oracle = only(filter(architecture -> architecture.kind == :fmc,
        M.fixed_architectures(config)))
    scales, scale_rows = M.fresh_nuisance_scales(
        config, split_libraries, oracle,
    )
    models = M.build_frozen_models(
        architectures,
        config,
        hypotheses,
        split_libraries,
        scales,
        saved_motion,
    )
    amplitudes = Float64.(saved_motion["amplitudes"])
    reference_alpha = Float64(saved_motion["reference_alpha"])
    contexts, replay_rows = M.nominal_replay(
        architectures,
        models,
        hypotheses,
        config,
        amplitudes,
        reference_alpha,
        saved_motion,
    )
    baseline_alpha = Dict(row.architecture_id => row.nominal_replay_alpha90
        for row in saved_o6f["summary"] if row.scale == 0.0)
    (; saved_motion, saved_o6f, config, hypotheses, architectures,
       split_libraries, scales, scale_rows, models, amplitudes,
       reference_alpha, contexts, replay_rows, baseline_alpha)
end

function evaluate_scale(setup, scenario, scale)
    config = setup.config
    hypotheses = setup.hypotheses
    frequency_count = length(config.frequency_ratios)
    physical = setup.saved_o6f
    position_standard = Float64.(physical["position_standard"])
    amplitude_standard = Float64.(physical["amplitude_standard"])
    delay_standard = Float64.(physical["delay_standard"])
    signal_noise_rng = Xoshiro(M.O6F_SIGNAL_NOISE_SEED)
    h0_noise_rng = Xoshiro(M.O6F_H0_NOISE_SEED)
    test_hypotheses = (
        material=M.material_clutter_hypotheses(
            M.MOTION_SPLIT_SEEDS[:nuisance_test],
        ),
        backwall=M.backwall_hypotheses(
            M.MOTION_SPLIT_SEEDS[:nuisance_test],
        ),
    )
    ringing = setup.split_libraries[:test][:ringing]
    zero_standard = zeros(config.element_count)
    counts = Dict{Tuple{Symbol,Float64,Float64},NamedTuple}()
    false_alarms = Dict(architecture.id => 0 for architecture in setup.architectures)
    trials = Dict(architecture.id => 0 for architecture in setup.architectures)
    for (context_index, context) in enumerate(setup.contexts)
        signal_latent = M.complex_noise(
            signal_noise_rng, length(hypotheses), length(hypotheses),
        )
        h0_latent = M.complex_noise(
            h0_noise_rng, length(hypotheses), M.H0_TRIALS_PER_CONTEXT,
        )
        pstandard = scenario.position > 0 ?
            view(position_standard, :, context_index) : zero_standard
        astandard = scenario.amplitude > 0 ?
            view(amplitude_standard, :, context_index) : zero_standard
        dstandard = scenario.delay > 0 ?
            view(delay_standard, :, context_index) : zero_standard
        full = M.physical_realization_features(
            context,
            config,
            test_hypotheses,
            ringing,
            scenario.position .* pstandard,
            scenario.amplitude .* astandard,
            scenario.delay .* dstandard,
            scale,
        )
        for architecture in setup.architectures
            id = architecture.id
            model = setup.models[id]
            indices = architecture.mechanical_indices
            signal_dictionary = M.select_stops(
                full.signal, indices, frequency_count,
            )
            clutter = hcat(
                setup.scales[:backwall] .* M.select_stops(
                    full.backwall, indices, frequency_count,
                ),
                setup.scales[:ringing] .* M.select_stops(
                    full.ringing, indices, frequency_count,
                ),
                setup.scales[:material] .* M.select_stops(
                    full.material, indices, frequency_count,
                ),
            )
            factor = M.covariance_factor(M.matched_output_covariance(
                model.precision_applied, clutter,
            ))
            batch = M.evaluate_with_common_noise(
                factor * signal_latent,
                M.matched_mean(model.precision_applied, signal_dictionary),
                model.normalization,
                model.threshold,
                setup.amplitudes,
                hypotheses,
                context.hypotheses,
            )
            M.add_batch!(counts, id, scale, batch)
            scores, _ = M.normalized_scores(
                factor * h0_latent, model.normalization,
            )
            false_alarms[id] += count(>(model.threshold), scores)
            trials[id] += length(scores)
        end
    end
    curves = M.aggregate_curve_rows(
        counts, setup.architectures, (scale,), setup.amplitudes,
    )
    rows = NamedTuple[]
    for architecture in setup.architectures
        id = architecture.id
        alpha90 = M.curve_alpha90(curves, id, scale)
        reference = M.curve_reference(
            curves, id, scale, setup.reference_alpha,
        )
        pfa_upper = M.wilson_upper_bound(false_alarms[id], trials[id])
        ratio = alpha90 / setup.baseline_alpha[id]
        gates = (
            pfa=pfa_upper <= M.PFA_GATE,
            alpha90=isfinite(ratio) && ratio <= M.ALPHA90_RATIO_GATE,
            pdl=reference.pdl >= M.PDL_REFERENCE_GATE,
        )
        push!(rows, (
            scenario=scenario.id,
            architecture_id=id,
            scale,
            position_rms_lambda=scenario.position * scale * M.POSITION_RMS_LAMBDA,
            contact_log_amplitude_rms=
                scenario.amplitude * scale * M.CONTACT_LOG_AMPLITUDE_RMS,
            contact_delay_rms_cycles=
                scenario.delay * scale * M.CONTACT_DELAY_RMS_CYCLES,
            h0_trials=trials[id],
            false_alarms=false_alarms[id],
            observed_pfa=false_alarms[id] / trials[id],
            pfa_upper_95=pfa_upper,
            alpha90,
            nominal_replay_alpha90=setup.baseline_alpha[id],
            alpha90_ratio_to_nominal=ratio,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            pfa_gate=gates.pfa,
            alpha90_gate=gates.alpha90,
            pdl_gate=gates.pdl,
            joint_pass=all(values(gates)),
        ))
    end
    rows, [merge((scenario=scenario.id,), row) for row in curves]
end

function cached_evaluation!(cache, curve_cache, setup, scenario, scale)
    key = (scenario.id, Float64(scale))
    if !haskey(cache, key)
        println("[>] O6g-B $(scenario.id) scale=$(round(scale; digits=6))")
        rows, curves = evaluate_scale(setup, scenario, Float64(scale))
        cache[key] = rows
        curve_cache[key] = curves
    end
    cache[key]
end

function architecture_row(rows, architecture_id)
    only(filter(row -> row.architecture_id == architecture_id, rows))
end

function first_failure_bracket!(
    cache, curve_cache, setup, scenario, architecture_id, initial_scales,
)
    scales = sort(unique(Float64.(collect(initial_scales))))
    first(scales) == 0.0 || pushfirst!(scales, 0.0)
    rows = [architecture_row(cached_evaluation!(
        cache, curve_cache, setup, scenario, scale,
    ), architecture_id) for scale in scales]
    rows[1].joint_pass || error("zero physical error must pass")
    failure_index = findfirst(row -> !row.joint_pass, rows)
    while isnothing(failure_index)
        next_scale = 2last(scales)
        next_scale <= MAXIMUM_SEARCH_SCALE || error(
            "no failure below maximum search scale for $architecture_id",
        )
        push!(scales, next_scale)
        push!(rows, architecture_row(cached_evaluation!(
            cache, curve_cache, setup, scenario, next_scale,
        ), architecture_id))
        failure_index = rows[end].joint_pass ? nothing : length(rows)
    end
    upper_index = something(failure_index)
    lower_index = upper_index - 1
    while lower_index > 1 && !rows[lower_index].joint_pass
        lower_index -= 1
    end
    rows[lower_index].joint_pass || error("passing lower bracket not found")
    (scales[lower_index], scales[upper_index])
end

function bisect_threshold!(
    cache, curve_cache, setup, scenario, architecture_id, lower, upper,
)
    for _ in 1:BISECTION_STEPS
        middle = (lower + upper) / 2
        row = architecture_row(cached_evaluation!(
            cache, curve_cache, setup, scenario, middle,
        ), architecture_id)
        if row.joint_pass
            lower = middle
        else
            upper = middle
        end
    end
    passing = architecture_row(cache[(scenario.id, Float64(lower))], architecture_id)
    failing = architecture_row(cache[(scenario.id, Float64(upper))], architecture_id)
    failed_gates = join((String(name) for (name, passed) in pairs((
        pfa=failing.pfa_gate,
        alpha90=failing.alpha90_gate,
        pdl=failing.pdl_gate,
    )) if !passed), ";")
    (
        scenario=scenario.id,
        architecture_id,
        maximum_passing_scale=lower,
        first_failing_scale=upper,
        threshold_width=upper - lower,
        position_rms_lambda=passing.position_rms_lambda,
        contact_log_amplitude_rms=passing.contact_log_amplitude_rms,
        contact_delay_rms_cycles=passing.contact_delay_rms_cycles,
        passing_pfa_upper_95=passing.pfa_upper_95,
        passing_alpha90_ratio=passing.alpha90_ratio_to_nominal,
        passing_pdl_reference=passing.pdl_reference,
        first_failed_gates=failed_gates,
    )
end

function write_rows(path, rows)
    M.write_rows(path, rows)
end

function plot_thresholds(threshold_rows)
    plot_object = plot(
        xlabel="error scenario",
        ylabel="maximum passing scale",
        title="O6g-B motion tolerance attribution",
        legend=:topright,
    )
    scenario_order = getproperty.(SCENARIOS, :id)
    for (architecture_id, marker, color) in zip(
        M.ARCHITECTURE_IDS, (:circle, :diamond), (:steelblue, :darkorange),
    )
        selected = [only(filter(row ->
            row.architecture_id == architecture_id && row.scenario == scenario,
            threshold_rows,
        )) for scenario in scenario_order]
        plot!(
            plot_object,
            1:length(scenario_order),
            getproperty.(selected, :maximum_passing_scale);
            marker,
            color,
            linewidth=2,
            label=String(architecture_id),
        )
    end
    plot!(plot_object; xticks=(1:length(scenario_order), String.(scenario_order)))
    savefig(plot_object, joinpath(
        O6G_MOTION_OUTPUT_ROOT, "motion_tolerance_attribution.png",
    ))
end

function run_o6g_motion_boundary()
    mkpath(O6G_MOTION_OUTPUT_ROOT)
    setup = setup_frozen_replay()
    cache = Dict{Tuple{Symbol,Float64},Vector{NamedTuple}}()
    curve_cache = Dict{Tuple{Symbol,Float64},Vector{NamedTuple}}()
    threshold_rows = NamedTuple[]
    coupled = first(SCENARIOS)
    println("[>] O6g-B coupled failure brackets")
    coupled_thresholds = Dict{Symbol,Float64}()
    for architecture_id in M.ARCHITECTURE_IDS
        lower, upper = first_failure_bracket!(
            cache,
            curve_cache,
            setup,
            coupled,
            architecture_id,
            COUPLED_COARSE_SCALES,
        )
        threshold = bisect_threshold!(
            cache,
            curve_cache,
            setup,
            coupled,
            architecture_id,
            lower,
            upper,
        )
        coupled_thresholds[architecture_id] = threshold.maximum_passing_scale
        push!(threshold_rows, threshold)
    end
    println("[>] O6g-B one-at-a-time attribution")
    for scenario in Base.tail(SCENARIOS), architecture_id in M.ARCHITECTURE_IDS
        initial = (0.0, coupled_thresholds[architecture_id])
        lower, upper = first_failure_bracket!(
            cache,
            curve_cache,
            setup,
            scenario,
            architecture_id,
            initial,
        )
        push!(threshold_rows, bisect_threshold!(
            cache,
            curve_cache,
            setup,
            scenario,
            architecture_id,
            lower,
            upper,
        ))
    end
    all_rows = reduce(vcat, values(cache))
    all_curves = reduce(vcat, values(curve_cache))
    sort!(all_rows; by=row -> (String(row.scenario), String(row.architecture_id), row.scale))
    sort!(all_curves; by=row ->
        (String(row.scenario), String(row.architecture_id), row.scale, row.amplitude))
    sort!(threshold_rows; by=row -> (String(row.architecture_id), String(row.scenario)))
    write_rows(joinpath(O6G_MOTION_OUTPUT_ROOT, "boundary_evaluations.csv"), all_rows)
    write_rows(joinpath(O6G_MOTION_OUTPUT_ROOT, "boundary_curves.csv"), all_curves)
    write_rows(joinpath(O6G_MOTION_OUTPUT_ROOT, "motion_tolerance_thresholds.csv"), threshold_rows)
    jldsave(
        joinpath(O6G_MOTION_OUTPUT_ROOT, "o6g_motion_boundary.jld2");
        format_version=1,
        source_o6f_path=abspath(O6G_MOTION_SOURCE),
        architecture_ids=M.ARCHITECTURE_IDS,
        scenarios=SCENARIOS,
        coarse_scales=COUPLED_COARSE_SCALES,
        bisection_steps=BISECTION_STEPS,
        replay_rows=setup.replay_rows,
        all_rows,
        all_curves,
        threshold_rows,
    )
    open(joinpath(O6G_MOTION_OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6g-B frozen motion failure boundary")
        println(io)
        println(io, "The exact O6f detector reconstruction, 60 signal contexts, physical standard draws and fresh-noise streams are reused. Every evaluated point contains 12,000 aggregate H0 trials; no layout or threshold was reoptimized.")
        println(io)
        println(io, "| Layout | Scenario | Max passing scale | Position RMS | Contact log-amp RMS | Delay RMS | First failed gate |")
        println(io, "|---|---|---:|---:|---:|---:|---|")
        for row in threshold_rows
            @printf(io,
                "| `%s` | `%s` | `%.5f` | `%.4f lambda` | `%.2f%%` | `%.4f cycles` | `%s` |\n",
                String(row.architecture_id), String(row.scenario),
                row.maximum_passing_scale, row.position_rms_lambda,
                100row.contact_log_amplitude_rms,
                row.contact_delay_rms_cycles, row.first_failed_gates,
            )
        end
        println(io)
        println(io, "These are finite boundaries of the reduced multiplicative/delay contact model, not real couplant specifications. O6g exhausts this surrogate; surviving motion risk is promoted to full-wave contact/lift-off in O6h.")
    end
    plot_thresholds(threshold_rows)
    for row in threshold_rows
        println("[+] $(row.architecture_id) $(row.scenario) threshold=$(round(row.maximum_passing_scale; digits=5))")
    end
    println("[+] wrote O6g-B to $O6G_MOTION_OUTPUT_ROOT")
    (; threshold_rows, all_rows)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6g_motion_boundary()

end
