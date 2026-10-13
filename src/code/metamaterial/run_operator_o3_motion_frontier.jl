using LinearAlgebra
using Random
using Statistics
using JLD2
using Plots

include(joinpath(@__DIR__, "run_operator_o3_paired_clutter.jl"))

const MOTION_OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "motion_frontier_v0",
)
const MOTION_BASE_SEED = UInt64(0x04f0a117)
const MOTION_SPLIT_SEEDS = Dict(
    :nuisance_train => MOTION_BASE_SEED + 0x101,
    :nuisance_validation => MOTION_BASE_SEED + 0x202,
    :nuisance_test => MOTION_BASE_SEED + 0x303,
    :signal_validation => MOTION_BASE_SEED + 0x404,
    :signal_test => MOTION_BASE_SEED + 0x505,
    :h0_validation => MOTION_BASE_SEED + 0x606,
    :h0_test => MOTION_BASE_SEED + 0x707,
)
const MOTION_VALIDATION_H0_TRIALS = 10_000
const MOTION_TEST_H0_TRIALS = 10_000
const MOTION_VALIDATION_CONTEXTS = 24
const MOTION_TEST_CONTEXTS = 60
const MOTION_POSITION_COUNTS = [1, 2, 3, 4, 6, 8, 12, 17]
const MOTION_AMPLITUDE_MULTIPLIERS = O3_AMPLITUDE_MULTIPLIERS
const MOTION_TIME_SCENARIOS = [
    (id=:electronics_only, step_cost=0.0, span_cost=0.0),
    (id=:one_event_per_settle, step_cost=1.0, span_cost=0.0),
    (id=:settle_plus_travel, step_cost=1.0, span_cost=1.0),
    (id=:slow_motion, step_cost=4.0, span_cost=2.0),
]

function fixed_architectures(config)
    isfile(O3_CODEBOOK_PATH) ||
        error("run optimize_operator_multistate_codebooks.jl first")
    definitions = JLD2.load(O3_CODEBOOK_PATH)["definitions"]
    state_ports = definitions[O3_ENTRY_ID].state_ports
    [
        (
            id=:a_scan_uniform,
            kind=:a_scan,
            electronic_channels=1,
            transmit_events=1,
            physical_states=1,
            mechanical_positions=1,
            state_blocks=1,
            block_size=1,
            state_ports=(),
            mechanical_indices=Int[],
            strategy=:stationary_uniform,
            role=:baseline,
        ),
        (
            id=:four_channel_two_static_states,
            kind=:two_static_states,
            electronic_channels=4,
            transmit_events=8,
            physical_states=2,
            mechanical_positions=1,
            state_blocks=2,
            block_size=16,
            state_ports,
            mechanical_indices=Int[],
            strategy=:frozen_o2c,
            role=:static_comparator,
        ),
        (
            id=:paut_fmc_oracle,
            kind=:fmc,
            electronic_channels=config.element_count,
            transmit_events=config.element_count,
            physical_states=1,
            mechanical_positions=1,
            state_blocks=1,
            block_size=config.element_count^2,
            state_ports=(),
            mechanical_indices=Int[],
            strategy=:full_matrix_capture,
            role=:oracle,
        ),
    ]
end

function quantized_indices(targets, positions)
    unused = collect(eachindex(positions))
    selected = Int[]
    for target in targets
        distances = abs.(positions[unused] .- target)
        local_index = argmin(distances)
        push!(selected, unused[local_index])
        deleteat!(unused, local_index)
    end
    sort!(selected)
end

function span_uniform_indices(position_count, span, positions)
    position_count == 1 && return [argmin(abs.(positions))]
    eligible = findall(abs.(positions) .<= span / 2 + 1.0e-12)
    position_count <= length(eligible) || return Int[]
    quantized_indices(
        collect(range(-span / 2, span / 2; length=position_count)),
        positions,
    )
end

function compact_indices(position_count, positions)
    sort(sortperm(abs.(positions))[1:position_count])
end

function chebyshev_indices(position_count, positions)
    position_count == 1 && return [argmin(abs.(positions))]
    targets = reverse(4 .* cos.((2 .* (1:position_count) .- 1) .* pi ./ (2position_count)))
    quantized_indices(targets, positions)
end

function motion_architecture(id, indices, strategy)
    (
        id,
        kind=:mechanical_scan,
        electronic_channels=1,
        transmit_events=length(indices),
        physical_states=1,
        mechanical_positions=length(indices),
        state_blocks=length(indices),
        block_size=1,
        state_ports=(),
        mechanical_indices=collect(indices),
        strategy,
        role=:motion_candidate,
    )
end

function predefined_motion_architectures(config)
    positions = element_positions(config)
    architectures = NamedTuple[]
    seen = Set{Tuple{Vararg{Int}}}()
    function register(indices, strategy, suffix)
        isempty(indices) && return
        key = Tuple(indices)
        key in seen && return
        push!(seen, key)
        id = Symbol("motion_$(strategy)_p$(length(indices))_$(suffix)")
        push!(architectures, motion_architecture(id, indices, strategy))
    end
    for count in MOTION_POSITION_COUNTS
        register(
            span_uniform_indices(count, config.aperture_wavelengths, positions),
            :uniform_full,
            "a8",
        )
        register(compact_indices(count, positions), :compact_center, "native")
        register(chebyshev_indices(count, positions), :chebyshev_full, "a8")
    end
    for count in (4, 8), span in (2.0, 4.0, 6.0, 8.0)
        indices = span_uniform_indices(count, span, positions)
        register(indices, :aperture_sweep, "a$(Int(span))")
    end
    architectures
end

function fresh_nuisance_libraries(config, split)
    seed = MOTION_SPLIT_SEEDS[Symbol("nuisance_$(split)")]
    Dict(
        :material => response_library(material_clutter_hypotheses(seed), config),
        :backwall => response_library(backwall_hypotheses(seed), config),
        :ringing => ringing_library(config, seed),
    )
end

function fresh_nuisance_scales(config, split_libraries, oracle)
    scales = Dict{Symbol,Float64}()
    rows = NamedTuple[]
    for category in (:backwall, :ringing, :material)
        features = electronics_whiten(
            compress_pulse_library(split_libraries[:train][category], oracle, config),
            oracle,
            config,
        )
        target = O3_CATEGORY_POWER[category]
        scale = sqrt(target * size(features, 1) / sum(abs2, features))
        scales[category] = scale
        push!(rows, (
            category,
            target_fmc_power_ratio=target,
            physical_scale=scale,
            achieved_fmc_power_ratio=sum(abs2, scale .* features) / size(features, 1),
        ))
    end
    scales, rows
end

function architecture_dictionary(library, architecture, config; state_errors=nothing)
    electronics_whiten(
        compress_pulse_library(
            library,
            architecture,
            config;
            state_errors,
        ),
        architecture,
        config,
    )
end

function architecture_clutter(libraries, architecture, config, scales)
    hcat((
        scales[category] .* architecture_dictionary(
            libraries[category], architecture, config,
        )
        for category in (:backwall, :ringing, :material)
    )...)
end

function d_optimal_score(dictionary, clutter, gamma)
    precision_applied = precision_dictionary(dictionary, clutter)
    singular_values = svdvals(precision_applied)
    sum(log1p.(gamma .* abs2.(singular_values)))
end

function greedy_motion_architectures(
    config,
    nominal_library,
    train_libraries,
    scales,
    gamma,
)
    selected = Int[]
    available = collect(1:config.element_count)
    architectures = NamedTuple[]
    trace = NamedTuple[]
    for step in 1:config.element_count
        candidates = NamedTuple[]
        for index in available
            indices = sort(vcat(selected, index))
            architecture = motion_architecture(:temporary, indices, :train_d_optimal)
            dictionary = architecture_dictionary(nominal_library, architecture, config)
            clutter = architecture_clutter(train_libraries, architecture, config, scales)
            push!(candidates, (
                added_index=index,
                score=d_optimal_score(dictionary, clutter, gamma),
            ))
        end
        sort!(candidates; by=row -> (-row.score, row.added_index))
        winner = first(candidates)
        push!(selected, winner.added_index)
        sort!(selected)
        deleteat!(available, findfirst(==(winner.added_index), available))
        push!(trace, (
            step,
            added_index=winner.added_index,
            selected_indices=join(selected, ";"),
            train_d_optimal_score=winner.score,
        ))
        if step in MOTION_POSITION_COUNTS
            push!(architectures, motion_architecture(
                Symbol("motion_train_d_optimal_p$(step)"),
                selected,
                :train_d_optimal,
            ))
        end
    end
    architectures, trace
end

function unique_architectures(architectures)
    result = NamedTuple[]
    seen = Set{Tuple{Vararg{Int}}}()
    for architecture in architectures
        architecture.kind == :mechanical_scan || begin
            push!(result, architecture)
            continue
        end
        key = Tuple(architecture.mechanical_indices)
        key in seen && continue
        push!(seen, key)
        push!(result, architecture)
    end
    result
end

function build_fresh_feature_cache(architectures, config, split_libraries, scales)
    cache = Dict{Tuple{Symbol,Symbol,Symbol},Matrix{ComplexF64}}()
    for architecture in architectures, split in (:train, :validation, :test)
        for category in (:backwall, :ringing, :material)
            cache[(architecture.id, split, category)] = scales[category] .*
                architecture_dictionary(
                    split_libraries[split][category], architecture, config,
                )
        end
    end
    cache
end

function cached_clutter(cache, architecture, split)
    hcat((cache[(architecture.id, split, category)]
          for category in (:backwall, :ringing, :material))...)
end

function build_fresh_models(architectures, nominal_dictionaries, cache)
    models = Dict{Symbol,NamedTuple}()
    rows = NamedTuple[]
    for (architecture_index, architecture) in enumerate(architectures)
        dictionary = nominal_dictionaries[architecture.id]
        precision_applied = precision_dictionary(
            dictionary, cached_clutter(cache, architecture, :train),
        )
        gram = matched_mean(precision_applied, dictionary)
        normalization = max.(real.(diag(gram)), eps(Float64))
        validation_factor = covariance_factor(matched_output_covariance(
            precision_applied, cached_clutter(cache, architecture, :validation),
        ))
        test_factor = covariance_factor(matched_output_covariance(
            precision_applied, cached_clutter(cache, architecture, :test),
        ))
        validation_rng = Xoshiro(
            MOTION_SPLIT_SEEDS[:h0_validation] + UInt64(10_000 * architecture_index),
        )
        validation_scores = no_defect_scores(
            validation_rng,
            validation_factor,
            normalization,
            MOTION_VALIDATION_H0_TRIALS,
        )
        threshold = quantile(validation_scores, 0.99)
        test_rng = Xoshiro(
            MOTION_SPLIT_SEEDS[:h0_test] + UInt64(10_000 * architecture_index),
        )
        test_scores = no_defect_scores(
            test_rng,
            test_factor,
            normalization,
            MOTION_TEST_H0_TRIALS,
        )
        false_alarms = count(>(threshold), test_scores)
        models[architecture.id] = (;
            precision_applied,
            normalization,
            validation_factor,
            test_factor,
            threshold,
        )
        push!(rows, (
            architecture_id=architecture.id,
            threshold,
            validation_trials=MOTION_VALIDATION_H0_TRIALS,
            test_trials=MOTION_TEST_H0_TRIALS,
            false_alarms,
            observed_pfa=false_alarms / MOTION_TEST_H0_TRIALS,
            pfa_upper_95=wilson_upper_bound(false_alarms, MOTION_TEST_H0_TRIALS),
        ))
    end
    models, rows
end

function evaluate_with_common_noise(
    noise,
    mean_matrix,
    normalization,
    threshold,
    amplitudes,
    templates,
    true_hypotheses,
)
    rows = NamedTuple[]
    for amplitude in amplitudes
        scores, predicted = normalized_scores(
            noise .+ amplitude .* mean_matrix,
            normalization,
        )
        detected = scores .> threshold
        localized = 0
        errors = Float64[]
        for true_index in eachindex(true_hypotheses)
            detected[true_index] || continue
            error = localization_error(
                predicted[true_index], true_index, templates, true_hypotheses,
            )
            push!(errors, error)
            localized += error <= 0.5 + 1.0e-12
        end
        push!(rows, (
            amplitude,
            trials=length(true_hypotheses),
            detections=count(detected),
            localized,
            pd=count(detected) / length(true_hypotheses),
            pdl=localized / length(true_hypotheses),
            median_localization_error=isempty(errors) ? NaN : median(errors),
            q90_localization_error=isempty(errors) ? NaN : quantile(errors, 0.90),
        ))
    end
    rows
end

function fresh_signal_context(rng, base_hypotheses, config, context_index)
    signal_context(rng, base_hypotheses, config, context_index)
end

function validation_amplitudes(oracle, model, nominal_dictionary, hypotheses, config)
    nominal_mean = matched_mean(model.precision_applied, nominal_dictionary)
    requirements = [
        sqrt(model.threshold * model.normalization[index]) /
            max(abs(nominal_mean[index, index]), eps(Float64))
        for index in eachindex(model.normalization)
    ]
    anchor = median(requirements)
    amplitudes = anchor .* MOTION_AMPLITUDE_MULTIPLIERS
    detections = zeros(Int, length(amplitudes))
    trials = zeros(Int, length(amplitudes))
    rng = Xoshiro(MOTION_SPLIT_SEEDS[:signal_validation])
    for context_index in 1:MOTION_VALIDATION_CONTEXTS
        context = fresh_signal_context(rng, hypotheses, config, context_index)
        dictionary = architecture_dictionary(
            response_library(context.hypotheses, context.config),
            oracle,
            context.config,
        )
        mean_matrix = matched_mean(model.precision_applied, dictionary)
        latent = complex_noise(rng, length(hypotheses), length(hypotheses))
        batch = evaluate_with_common_noise(
            model.validation_factor * latent,
            mean_matrix,
            model.normalization,
            model.threshold,
            amplitudes,
            hypotheses,
            context.hypotheses,
        )
        for (index, row) in enumerate(batch)
            detections[index] += row.detections
            trials[index] += row.trials
        end
    end
    probabilities = detections ./ trials
    index = findfirst(>=(0.95), accumulate(max, probabilities))
    isnothing(index) && error("fresh FMC validation Pd did not reach 0.95")
    amplitudes, probabilities, amplitudes[index], anchor
end

function fresh_test_campaign(architectures, models, hypotheses, config, amplitudes)
    counts = Dict{Tuple{Symbol,Float64},NamedTuple}()
    errors = Dict{Tuple{Symbol,Float64},Vector{Float64}}()
    for architecture in architectures, amplitude in amplitudes
        counts[(architecture.id, amplitude)] = (trials=0, detections=0, localized=0)
        errors[(architecture.id, amplitude)] = Float64[]
    end
    context_rows = NamedTuple[]
    rng = Xoshiro(MOTION_SPLIT_SEEDS[:signal_test])
    for context_index in 1:MOTION_TEST_CONTEXTS
        context = fresh_signal_context(rng, hypotheses, config, 100_000 + context_index)
        library = response_library(context.hypotheses, context.config)
        latent = complex_noise(rng, length(hypotheses), length(hypotheses))
        for architecture in architectures
            model = models[architecture.id]
            dictionary = architecture_dictionary(
                library,
                architecture,
                context.config;
                state_errors=architecture_state_errors(architecture, context),
            )
            mean_matrix = matched_mean(model.precision_applied, dictionary)
            batch = evaluate_with_common_noise(
                model.test_factor * latent,
                mean_matrix,
                model.normalization,
                model.threshold,
                amplitudes,
                hypotheses,
                context.hypotheses,
            )
            for row in batch
                key = (architecture.id, row.amplitude)
                previous = counts[key]
                counts[key] = (
                    trials=previous.trials + row.trials,
                    detections=previous.detections + row.detections,
                    localized=previous.localized + row.localized,
                )
                # Context-level quantiles are already sufficient for the audit.
                push!(context_rows, merge((;
                    context_index,
                    architecture_id=architecture.id,
                ), row))
            end
        end
        context_index % 10 == 0 && println(
            "[+] fresh motion test contexts: $context_index / $MOTION_TEST_CONTEXTS",
        )
    end
    aggregate = NamedTuple[]
    for architecture in architectures, amplitude in amplitudes
        count_row = counts[(architecture.id, amplitude)]
        selected = filter(row ->
            row.architecture_id == architecture.id &&
            isapprox(row.amplitude, amplitude; rtol=1.0e-12),
            context_rows,
        )
        medians = filter(isfinite, getproperty.(selected, :median_localization_error))
        q90s = filter(isfinite, getproperty.(selected, :q90_localization_error))
        push!(aggregate, (
            architecture_id=architecture.id,
            amplitude,
            trials=count_row.trials,
            detections=count_row.detections,
            localized=count_row.localized,
            pd=count_row.detections / count_row.trials,
            pdl=count_row.localized / count_row.trials,
            median_localization_error=isempty(medians) ? NaN : median(medians),
            q90_localization_error=isempty(q90s) ? NaN : quantile(q90s, 0.90),
        ))
    end
    aggregate, context_rows
end

function reference_row(rows, architecture_id, reference_alpha)
    only(filter(row ->
        row.architecture_id == architecture_id &&
        isapprox(row.amplitude, reference_alpha; rtol=1.0e-12),
        rows,
    ))
end

function alpha90_for(rows, architecture_id)
    selected = sort(filter(row -> row.architecture_id == architecture_id, rows);
                    by=row -> row.amplitude)
    crossing_amplitude(
        getproperty.(selected, :amplitude), getproperty.(selected, :pd), 0.90,
    )
end

function context_alpha_rows(context_rows, architectures, reference_alpha)
    rows = NamedTuple[]
    for context_index in 1:MOTION_TEST_CONTEXTS, architecture in architectures
        selected = sort(filter(row ->
            row.context_index == context_index &&
            row.architecture_id == architecture.id,
            context_rows,
        ); by=row -> row.amplitude)
        reference = only(filter(row ->
            isapprox(row.amplitude, reference_alpha; rtol=1.0e-12), selected,
        ))
        push!(rows, (
            context_index,
            architecture_id=architecture.id,
            alpha90=crossing_amplitude(
                getproperty.(selected, :amplitude), getproperty.(selected, :pd), 0.90,
            ),
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
        ))
    end
    rows
end

function motion_resource_metadata(architectures, config)
    positions = element_positions(config)
    rows = NamedTuple[]
    for architecture in architectures
        indices = architecture.mechanical_indices
        span = isempty(indices) ? 0.0 : maximum(positions[indices]) - minimum(positions[indices])
        move_count = architecture.kind == :mechanical_scan ? max(length(indices) - 1, 0) : 0
        base = (
            architecture_id=architecture.id,
            role=architecture.role,
            strategy=architecture.strategy,
            electronic_channels=architecture.electronic_channels,
            transmit_events=architecture.transmit_events,
            physical_states=architecture.physical_states,
            mechanical_positions=architecture.mechanical_positions,
            mechanical_moves=move_count,
            aperture_span_lambda=span,
            position_indices=join(indices, ";"),
            total_incident_energy=1.0,
            complex_features=architecture.state_blocks * architecture.block_size *
                length(config.frequency_ratios),
        )
        times = NamedTuple{Tuple(Symbol("time_$(scenario.id)") for scenario in MOTION_TIME_SCENARIOS)}(
            Tuple(
                architecture.transmit_events + scenario.step_cost * move_count +
                    scenario.span_cost * span
                for scenario in MOTION_TIME_SCENARIOS
            ),
        )
        push!(rows, merge(base, times))
    end
    rows
end

function summarize_frontier(
    architectures,
    metadata,
    thresholds,
    curves,
    context_alpha,
    reference_alpha,
)
    alpha_values = Dict(architecture.id => alpha90_for(curves, architecture.id)
                        for architecture in architectures)
    a_alpha = alpha_values[:a_scan_uniform]
    fmc_alpha = alpha_values[:paut_fmc_oracle]
    rows = NamedTuple[]
    for architecture in architectures
        id = architecture.id
        meta = only(filter(row -> row.architecture_id == id, metadata))
        threshold = only(filter(row -> row.architecture_id == id, thresholds))
        reference = reference_row(curves, id, reference_alpha)
        selected_contexts = filter(row -> row.architecture_id == id, context_alpha)
        finite_contexts = filter(row -> isfinite(row.alpha90), selected_contexts)
        gap = a_alpha > fmc_alpha ?
            (a_alpha - alpha_values[id]) / (a_alpha - fmc_alpha) : NaN
        context_gaps = Float64[]
        for row in finite_contexts
            a_context = only(filter(other ->
                other.context_index == row.context_index &&
                other.architecture_id == :a_scan_uniform,
                context_alpha,
            )).alpha90
            fmc_context = only(filter(other ->
                other.context_index == row.context_index &&
                other.architecture_id == :paut_fmc_oracle,
                context_alpha,
            )).alpha90
            if all(isfinite, (a_context, row.alpha90, fmc_context)) &&
               a_context > fmc_context + eps(Float64)
                push!(context_gaps,
                    (a_context - row.alpha90) / (a_context - fmc_context),
                )
            end
        end
        push!(rows, merge(meta, (
            observed_pfa=threshold.observed_pfa,
            pfa_upper_95=threshold.pfa_upper_95,
            pfa_controlled=threshold.pfa_upper_95 <= 0.0125,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            alpha90=alpha_values[id],
            aggregate_paut_gap_closure=gap,
            q10_context_gap_closure=isempty(context_gaps) ? NaN : quantile(context_gaps, 0.10),
            finite_context_fraction=length(finite_contexts) / MOTION_TEST_CONTEXTS,
        )))
    end
    rows
end

function dominates(first, second, time_field)
    first_time = getproperty(first, time_field)
    second_time = getproperty(second, time_field)
    first.alpha90 <= second.alpha90 && first_time <= second_time &&
        (first.alpha90 < second.alpha90 || first_time < second_time)
end

function dominates_with_channels(first, second, time_field)
    first_time = getproperty(first, time_field)
    second_time = getproperty(second, time_field)
    first.alpha90 <= second.alpha90 && first_time <= second_time &&
        first.electronic_channels <= second.electronic_channels &&
        (first.alpha90 < second.alpha90 || first_time < second_time ||
         first.electronic_channels < second.electronic_channels)
end

function resource_frontiers(summary_rows)
    rows = NamedTuple[]
    for scenario in MOTION_TIME_SCENARIOS
        time_field = Symbol("time_$(scenario.id)")
        for row in summary_rows
            pareto = row.pfa_controlled && !any(other ->
                other.architecture_id != row.architecture_id &&
                other.pfa_controlled &&
                dominates(other, row, time_field),
                summary_rows,
            )
            pareto_time_channels = row.pfa_controlled && !any(other ->
                other.architecture_id != row.architecture_id &&
                other.pfa_controlled &&
                dominates_with_channels(other, row, time_field),
                summary_rows,
            )
            efficiency = row.aggregate_paut_gap_closure / getproperty(row, time_field)
            channel_time = row.electronic_channels * getproperty(row, time_field)
            push!(rows, (
                scenario=scenario.id,
                step_cost=scenario.step_cost,
                span_cost=scenario.span_cost,
                architecture_id=row.architecture_id,
                normalized_time=getproperty(row, time_field),
                alpha90=row.alpha90,
                paut_gap_closure=row.aggregate_paut_gap_closure,
                gap_closure_per_time=efficiency,
                channel_time_product=channel_time,
                gap_closure_per_channel_time=row.aggregate_paut_gap_closure / channel_time,
                pfa_controlled=row.pfa_controlled,
                pareto_time_only=pareto,
                pareto_time_channels,
                within_static_time_budget=getproperty(row, time_field) <= 8.0 + 1.0e-12,
            ))
        end
    end
    rows
end

function pairwise_static_comparisons(summary_rows, context_summary)
    static = only(filter(row ->
        row.architecture_id == :four_channel_two_static_states, summary_rows,
    ))
    static_contexts = sort(filter(row ->
        row.architecture_id == :four_channel_two_static_states, context_summary,
    ); by=row -> row.context_index)
    rows = NamedTuple[]
    for motion in filter(row -> row.role == :motion_candidate, summary_rows)
        motion_contexts = sort(filter(row ->
            row.architecture_id == motion.architecture_id, context_summary,
        ); by=row -> row.context_index)
        ratios = [
            motion_contexts[index].alpha90 / static_contexts[index].alpha90
            for index in eachindex(static_contexts)
            if isfinite(motion_contexts[index].alpha90) &&
               isfinite(static_contexts[index].alpha90)
        ]
        push!(rows, (
            architecture_id=motion.architecture_id,
            mechanical_positions=motion.mechanical_positions,
            aperture_span_lambda=motion.aperture_span_lambda,
            pfa_controlled=motion.pfa_controlled,
            aggregate_alpha90_ratio_to_static=motion.alpha90 / static.alpha90,
            median_context_alpha90_ratio=median(ratios),
            q10_context_alpha90_ratio=quantile(ratios, 0.10),
            q90_context_alpha90_ratio=quantile(ratios, 0.90),
            context_win_fraction=mean(ratios .< 1.0),
            valid_contexts=length(ratios),
        ))
    end
    rows
end

function split_contract_rows()
    [(
        split=key,
        seed=string("0x", string(value; base=16)),
        base_seed=string("0x", string(MOTION_BASE_SEED; base=16)),
        reused_old_o3_test=false,
    ) for (key, value) in sort(collect(MOTION_SPLIT_SEEDS); by=first)]
end

function self_check(architectures, metadata, split_rows)
    @assert length(unique(getproperty.(architectures, :id))) == length(architectures)
    @assert length(unique(getproperty.(split_rows, :seed))) == length(split_rows)
    @assert all(row -> isapprox(row.total_incident_energy, 1.0), metadata)
    @assert all(architecture ->
        architecture.kind != :mechanical_scan ||
        length(unique(architecture.mechanical_indices)) == architecture.mechanical_positions,
        architectures,
    )
    @assert MOTION_SPLIT_SEEDS[:nuisance_test] != UInt64(0x3303)
    true
end

function plot_motion_frontier(summary_rows, frontier_rows)
    mechanical = filter(row -> row.role == :motion_candidate, summary_rows)
    positions_plot = plot(
        xlabel="mechanical positions P",
        ylabel="alpha90 (lower is better)",
        title="Fresh held-out mechanical operator frontier",
        legend=:topright,
        xscale=:log10,
    )
    colors = Dict(
        :uniform_full => :steelblue,
        :compact_center => :darkorange,
        :chebyshev_full => :seagreen,
        :aperture_sweep => :purple,
        :train_d_optimal => :firebrick,
    )
    for strategy in keys(colors)
        selected = sort(filter(row -> row.strategy == strategy, mechanical);
                        by=row -> row.mechanical_positions)
        isempty(selected) && continue
        scatter!(
            positions_plot,
            getproperty.(selected, :mechanical_positions),
            getproperty.(selected, :alpha90);
            marker=:circle,
            color=colors[strategy],
            label=String(strategy),
        )
    end
    static_alpha = only(filter(row ->
        row.architecture_id == :four_channel_two_static_states, summary_rows,
    )).alpha90
    fmc_alpha = only(filter(row ->
        row.architecture_id == :paut_fmc_oracle, summary_rows,
    )).alpha90
    hline!(positions_plot, [static_alpha]; color=:black, linestyle=:dash, label="frozen 4ch x 2state")
    hline!(positions_plot, [fmc_alpha]; color=:gray, linestyle=:dot, label="FMC")
    savefig(positions_plot, joinpath(MOTION_OUTPUT_ROOT, "motion_alpha90_vs_positions.png"))

    frontier_plot = plot(
        xlabel="normalized acquisition time",
        ylabel="alpha90 (lower is better)",
        title="Parametric time/performance fronts",
        legend=:topright,
        xscale=:log10,
    )
    scenario_colors = [:steelblue, :darkorange, :seagreen, :firebrick]
    for (scenario, color) in zip(MOTION_TIME_SCENARIOS, scenario_colors)
        selected = sort(filter(row -> row.scenario == scenario.id && row.pareto_time_only, frontier_rows);
                        by=row -> row.normalized_time)
        plot!(
            frontier_plot,
            getproperty.(selected, :normalized_time),
            getproperty.(selected, :alpha90);
            marker=:circle,
            linewidth=2,
            color,
            label=String(scenario.id),
        )
    end
    savefig(frontier_plot, joinpath(MOTION_OUTPUT_ROOT, "motion_resource_frontiers.png"))
end

function main_motion_frontier()
    mkpath(MOTION_OUTPUT_ROOT)
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    nominal_library = response_library(hypotheses, config)
    fixed = fixed_architectures(config)
    oracle = only(filter(architecture -> architecture.kind == :fmc, fixed))
    split_libraries = Dict(
        split => fresh_nuisance_libraries(config, split)
        for split in (:train, :validation, :test)
    )
    scales, scale_rows = fresh_nuisance_scales(config, split_libraries, oracle)

    a_scan = only(filter(architecture -> architecture.kind == :a_scan, fixed))
    a_dictionary = architecture_dictionary(nominal_library, a_scan, config)
    a_clutter = architecture_clutter(split_libraries[:train], a_scan, config, scales)
    a_precision = precision_dictionary(a_dictionary, a_clutter)
    gamma = inv(median(vec(sum(abs2, a_precision; dims=1))))
    greedy, greedy_trace = greedy_motion_architectures(
        config,
        nominal_library,
        split_libraries[:train],
        scales,
        gamma,
    )
    architectures = unique_architectures(vcat(
        fixed, predefined_motion_architectures(config), greedy,
    ))
    println("[+] frozen $(length(architectures)) architectures; no old O3 test data loaded")

    metadata = motion_resource_metadata(architectures, config)
    splits = split_contract_rows()
    self_check(architectures, metadata, splits)
    nominal_dictionaries = Dict(
        architecture.id => architecture_dictionary(
            nominal_library, architecture, config,
        )
        for architecture in architectures
    )
    cache = build_fresh_feature_cache(
        architectures, config, split_libraries, scales,
    )
    models, thresholds = build_fresh_models(
        architectures, nominal_dictionaries, cache,
    )
    println("[+] fresh C3 GLRT models and held-out Pfa estimates assembled")

    amplitudes, validation_pd, reference_alpha, amplitude_anchor =
        validation_amplitudes(
            oracle,
            models[oracle.id],
            nominal_dictionaries[oracle.id],
            hypotheses,
            config,
        )
    println("[+] validation-only FMC reference alpha=$(round(reference_alpha; digits=4))")
    curves, context_curves = fresh_test_campaign(
        architectures, models, hypotheses, config, amplitudes,
    )
    context_summary = context_alpha_rows(
        context_curves, architectures, reference_alpha,
    )
    summary_rows = summarize_frontier(
        architectures,
        metadata,
        thresholds,
        curves,
        context_summary,
        reference_alpha,
    )
    frontier_rows = resource_frontiers(summary_rows)
    pairwise_rows = pairwise_static_comparisons(summary_rows, context_summary)

    write_csv(joinpath(MOTION_OUTPUT_ROOT, "split_contract.csv"),
              propertynames(first(splits)), splits)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "nuisance_scales.csv"),
              propertynames(first(scale_rows)), scale_rows)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "greedy_trace.csv"),
              propertynames(first(greedy_trace)), greedy_trace)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "architecture_metadata.csv"),
              propertynames(first(metadata)), metadata)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "thresholds.csv"),
              propertynames(first(thresholds)), thresholds)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "detection_curves.csv"),
              propertynames(first(curves)), curves)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "context_summary.csv"),
              propertynames(first(context_summary)), context_summary)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "motion_summary.csv"),
              propertynames(first(summary_rows)), summary_rows)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "resource_frontiers.csv"),
              propertynames(first(frontier_rows)), frontier_rows)
    write_csv(joinpath(MOTION_OUTPUT_ROOT, "pairwise_vs_static.csv"),
              propertynames(first(pairwise_rows)), pairwise_rows)
    jldsave(
        joinpath(MOTION_OUTPUT_ROOT, "motion_frontier.jld2");
        config,
        hypotheses,
        splits,
        scale_rows,
        greedy_trace,
        metadata,
        thresholds,
        amplitudes,
        validation_pd,
        reference_alpha,
        amplitude_anchor,
        curves,
        context_summary,
        summary_rows,
        frontier_rows,
        pairwise_rows,
        base_seed=MOTION_BASE_SEED,
        validation_h0_trials=MOTION_VALIDATION_H0_TRIALS,
        test_h0_trials=MOTION_TEST_H0_TRIALS,
        validation_contexts=MOTION_VALIDATION_CONTEXTS,
        test_contexts=MOTION_TEST_CONTEXTS,
    )
    plot_motion_frontier(summary_rows, frontier_rows)

    static = only(filter(row ->
        row.architecture_id == :four_channel_two_static_states, summary_rows,
    ))
    motion = sort(filter(row -> row.role == :motion_candidate, summary_rows);
                  by=row -> row.alpha90)
    println("[+] frozen static alpha90=$(round(static.alpha90; digits=3))")
    for row in first(motion, min(5, length(motion)))
        println(
            "[+] $(row.architecture_id): P=$(row.mechanical_positions), " *
            "span=$(row.aperture_span_lambda), alpha90=$(round(row.alpha90; digits=3)), " *
            "Pfa95=$(round(row.pfa_upper_95; digits=4))",
        )
    end
    println("[+] wrote fresh motion frontier to $MOTION_OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_motion_frontier()
end
