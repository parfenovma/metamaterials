using LinearAlgebra
using Random
using Statistics
using JLD2
using Plots

include(joinpath(@__DIR__, "optimize_operator_multistate_codebooks.jl"))
include(joinpath(@__DIR__, "operator_o3_statistics.jl"))

using .OperatorO3Statistics

const O3_OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "paired_clutter_born_v0",
)
const O3_CODEBOOK_PATH = joinpath(OUTPUT_ROOT, "o2c_multistate_codebooks.jld2")
const O3_ENTRY_ID = :real_static_two_state_6
const O3_VALIDATION_H0_TRIALS = 20_000
const O3_TEST_H0_TRIALS = 20_000
const O3_VALIDATION_CONTEXTS = 40
const O3_TEST_CONTEXTS = 100
const O3_NOMINAL_REPLICATES = 20
const O3_RANDOM_SEED = 0x03c1a77e
const O3_AMPLITUDE_MULTIPLIERS = [
    0.4, 0.6, 0.8, 1.0, 1.3, 1.7, 2.2, 3.0, 4.0, 6.0, 8.0, 12.0,
    16.0, 24.0, 32.0,
]
const O3_CATEGORY_POWER = Dict(
    :backwall => 2.0,
    :ringing => 2.0,
    :material => 3.0,
)

function o3_architectures(config)
    isfile(O3_CODEBOOK_PATH) ||
        error("run optimize_operator_multistate_codebooks.jl first")
    bundle = JLD2.load(O3_CODEBOOK_PATH)
    definitions = bundle["definitions"]
    haskey(definitions, O3_ENTRY_ID) || error("frozen O3 entrant is missing")
    state_ports = definitions[O3_ENTRY_ID].state_ports
    mechanical_indices = unique(round.(Int, range(
        1, config.element_count; length=8,
    )))
    length(mechanical_indices) == 8 || error("mechanical scan positions collapsed")
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
            role=:baseline,
        ),
        (
            id=:mechanical_scan_equal_time,
            kind=:mechanical_scan,
            electronic_channels=1,
            transmit_events=8,
            physical_states=1,
            mechanical_positions=8,
            state_blocks=8,
            block_size=1,
            state_ports=(),
            mechanical_indices,
            role=:baseline,
        ),
        (
            id=:four_channel_one_static_state,
            kind=:one_static_state,
            electronic_channels=4,
            transmit_events=4,
            physical_states=1,
            mechanical_positions=1,
            state_blocks=1,
            block_size=16,
            state_ports=(state_ports[1],),
            mechanical_indices=Int[],
            role=:negative_control,
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
            role=:candidate,
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
            role=:oracle,
        ),
    ]
end

function pulse_features(raw, config, block_size)
    features, _, _ = apply_pulse_and_select(
        raw,
        config.frequency_ratios,
        collect(eachindex(config.frequency_ratios)),
        block_size,
    )
    features
end

function compress_pulse_library(
    library,
    architecture,
    config;
    state_errors=nothing,
)
    if architecture.kind == :a_scan
        port = ComplexF64.(uniform_port(config))
        raw = one_port_dictionary(library, _ -> port, config)
        return pulse_features(raw, config, 1)
    elseif architecture.kind == :mechanical_scan
        state_features_list = Matrix{ComplexF64}[]
        for element_index in architecture.mechanical_indices
            port = zeros(ComplexF64, config.element_count)
            port[element_index] = 1
            raw = one_port_dictionary(library, _ -> port, config)
            push!(state_features_list, pulse_features(raw, config, 1))
        end
        return reduce(vcat, state_features_list) ./ sqrt(length(state_features_list))
    elseif architecture.kind in (:one_static_state, :two_static_states)
        errors = isnothing(state_errors) ? [(
            amplitude=ones(size(ports)),
            delay=zeros(size(ports)),
        ) for ports in architecture.state_ports] : state_errors
        length(errors) == length(architecture.state_ports) ||
            throw(DimensionMismatch("one port-error realization per state is required"))
        state_features_list = Matrix{ComplexF64}[]
        for state_index in eachindex(architecture.state_ports)
            ports = architecture.state_ports[state_index]
            raw, _ = causal_multiport_dictionary(
                library,
                ports,
                config;
                model=:constant_contractive,
                amplitude_errors=errors[state_index].amplitude,
                delay_errors=errors[state_index].delay,
            )
            push!(state_features_list, pulse_features(raw, config, 16))
        end
        return reduce(vcat, state_features_list) ./ sqrt(length(state_features_list))
    elseif architecture.kind == :fmc
        raw = fmc_dictionary(library, config)
        return pulse_features(raw, config, config.element_count^2)
    end
    throw(ArgumentError("unknown O3 architecture: $(architecture.kind)"))
end

function electronics_frequency_covariance(frequencies)
    frequency_covariance(frequencies, 1; model=:correlated)
end

function electronics_whiten(features, architecture, config)
    frequency_count = length(config.frequency_ratios)
    rows_per_state = frequency_count * architecture.block_size
    expected_rows = architecture.state_blocks * rows_per_state
    size(features, 1) == expected_rows ||
        throw(DimensionMismatch("feature rows do not match architecture metadata"))
    covariance = electronics_frequency_covariance(config.frequency_ratios)
    lower = cholesky(Hermitian(covariance)).L
    result = similar(complex.(features))
    for state in 1:architecture.state_blocks
        rows = ((state - 1) * rows_per_state + 1):(state * rows_per_state)
        source = reshape(
            view(features, rows, :),
            architecture.block_size,
            frequency_count,
            size(features, 2),
        )
        destination = reshape(
            view(result, rows, :),
            architecture.block_size,
            frequency_count,
            size(features, 2),
        )
        for channel in 1:architecture.block_size
            destination[channel, :, :] .= lower \ source[channel, :, :]
        end
    end
    result
end

function material_clutter_hypotheses(seed; count=48)
    rng = Xoshiro(seed)
    hypotheses = DefectHypothesis[]
    for index in 1:count
        depth = 1.0 + 6.2 * rand(rng)
        lateral = -4.0 + 8.0 * rand(rng)
        angle = 2pi * rand(rng)
        secondary = -0.35 + 1.0 * rand(rng)
        rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)]
        tensor = rotation * Diagonal([1.0, secondary]) * transpose(rotation)
        push!(hypotheses, DefectHypothesis(
            Symbol("material_$(seed)_$(index)"),
            depth,
            lateral,
            Matrix{Float64}(tensor),
        ))
    end
    hypotheses
end

function backwall_hypotheses(seed; count=17)
    rng = Xoshiro(seed)
    laterals = range(-4.25, 4.25; length=count)
    [
        DefectHypothesis(
            Symbol("backwall_$(seed)_$(index)"),
            7.25 + 0.08 * randn(rng),
            lateral + 0.08 * randn(rng),
            [1.0 0.0; 0.0 0.08 + 0.04 * rand(rng)],
        )
        for (index, lateral) in enumerate(laterals)
    ]
end

function ringing_library(config, seed; atom_count=6)
    rng = Xoshiro(seed)
    count = config.element_count
    positions = element_positions(config)
    normalized_positions = positions ./ maximum(abs, positions)
    spatial_atoms = Matrix{ComplexF64}[]
    push!(spatial_atoms, Matrix{ComplexF64}(I, count, count))
    push!(spatial_atoms, ones(ComplexF64, count, count))
    decay = [exp(-abs(first - second) / 1.5) for first in positions, second in positions]
    push!(spatial_atoms, ComplexF64.(decay))
    for mode in 1:(atom_count - length(spatial_atoms))
        vector = @. cos(pi * mode * (normalized_positions + 1) / 2) +
            0.25 * sin(pi * (mode + 1) * normalized_positions)
        push!(spatial_atoms, ComplexF64.(vector * transpose(vector)))
    end
    for index in eachindex(spatial_atoms)
        spatial_atoms[index] ./= norm(spatial_atoms[index])
    end

    library = Array{ComplexF64}(undef,
        count, count, length(config.frequency_ratios), atom_count,
    )
    for atom in 1:atom_count
        delay = 0.12 + 1.5 * rand(rng)
        resonance = 0.85 + 0.30 * rand(rng)
        bandwidth = 0.10 + 0.12 * rand(rng)
        for (frequency_index, frequency) in enumerate(config.frequency_ratios)
            transfer = exp(-2pi * im * frequency * delay) /
                (1 + im * (frequency - resonance) / bandwidth)
            library[:, :, frequency_index, atom] = transfer .* spatial_atoms[atom]
        end
    end
    library
end

function nuisance_libraries(config, split)
    seed = Dict(:train => 0x1101, :validation => 0x2202, :test => 0x3303)[split]
    Dict(
        :material => response_library(material_clutter_hypotheses(seed), config),
        :backwall => response_library(backwall_hypotheses(seed), config),
        :ringing => ringing_library(config, seed),
    )
end

function nuisance_categories(level)
    level == :C0 && return Symbol[]
    level == :C1 && return [:backwall, :ringing]
    level in (:C2, :C3) && return [:backwall, :ringing, :material]
    throw(ArgumentError("unknown nuisance level: $level"))
end

function build_nuisance_features(architectures, config)
    split_libraries = Dict(
        split => nuisance_libraries(config, split)
        for split in (:train, :validation, :test)
    )
    full = only(filter(architecture -> architecture.kind == :fmc, architectures))
    category_scales = Dict{Symbol,Float64}()
    scale_rows = NamedTuple[]
    for category in keys(O3_CATEGORY_POWER)
        full_features = electronics_whiten(
            compress_pulse_library(split_libraries[:train][category], full, config),
            full,
            config,
        )
        target_power = O3_CATEGORY_POWER[category]
        scale = sqrt(target_power * size(full_features, 1) / sum(abs2, full_features))
        category_scales[category] = scale
        push!(scale_rows, (
            category,
            target_fmc_power_ratio=target_power,
            physical_scale=scale,
            achieved_fmc_power_ratio=sum(abs2, scale .* full_features) /
                size(full_features, 1),
        ))
    end

    features = Dict{Tuple{Symbol,Symbol,Symbol},Matrix{ComplexF64}}()
    for architecture in architectures, split in (:train, :validation, :test)
        for category in keys(O3_CATEGORY_POWER)
            compressed = compress_pulse_library(
                split_libraries[split][category], architecture, config,
            )
            features[(architecture.id, split, category)] =
                category_scales[category] .* electronics_whiten(
                    compressed, architecture, config,
                )
        end
    end
    features, scale_rows
end

function clutter_matrix(feature_cache, architecture, split, level, measurement_count)
    categories = nuisance_categories(level)
    isempty(categories) && return zeros(ComplexF64, measurement_count, 0)
    hcat((feature_cache[(architecture.id, split, category)] for category in categories)...)
end

complex_noise(rng, rows, columns) =
    (randn(rng, rows, columns) .+ im .* randn(rng, rows, columns)) ./ sqrt(2)

function no_defect_scores(rng, factor, normalization, trial_count; batch_size=2000)
    scores = Float64[]
    remaining = trial_count
    while remaining > 0
        count = min(batch_size, remaining)
        outputs = factor * complex_noise(rng, size(factor, 1), count)
        maxima, _ = normalized_scores(outputs, normalization)
        append!(scores, maxima)
        remaining -= count
    end
    scores
end

function build_detection_models(
    architectures,
    nominal_dictionaries,
    nuisance_features,
)
    levels = (:C0, :C1, :C2, :C3)
    models = Dict{Tuple{Symbol,Symbol},NamedTuple}()
    threshold_rows = NamedTuple[]
    for (architecture_index, architecture) in enumerate(architectures)
        dictionary = nominal_dictionaries[architecture.id]
        for (level_index, level) in enumerate(levels)
            train_clutter = clutter_matrix(
                nuisance_features, architecture, :train, level, size(dictionary, 1),
            )
            validation_clutter = clutter_matrix(
                nuisance_features, architecture, :validation, level, size(dictionary, 1),
            )
            test_clutter = clutter_matrix(
                nuisance_features, architecture, :test, level, size(dictionary, 1),
            )
            precision_applied = precision_dictionary(dictionary, train_clutter)
            gram = matched_mean(precision_applied, dictionary)
            normalization = max.(real.(diag(gram)), eps(Float64))
            validation_covariance = matched_output_covariance(
                precision_applied, validation_clutter,
            )
            test_covariance = matched_output_covariance(
                precision_applied, test_clutter,
            )
            validation_factor = OperatorO3Statistics.covariance_factor(validation_covariance)
            test_factor = OperatorO3Statistics.covariance_factor(test_covariance)
            validation_rng = Xoshiro(
                O3_RANDOM_SEED + UInt(10_000 * architecture_index + 100 * level_index + 1),
            )
            test_rng = Xoshiro(
                O3_RANDOM_SEED + UInt(10_000 * architecture_index + 100 * level_index + 2),
            )
            validation_scores = no_defect_scores(
                validation_rng,
                validation_factor,
                normalization,
                O3_VALIDATION_H0_TRIALS,
            )
            threshold = quantile(validation_scores, 0.99)
            test_scores = no_defect_scores(
                test_rng,
                test_factor,
                normalization,
                O3_TEST_H0_TRIALS,
            )
            false_alarms = count(>(threshold), test_scores)
            observed_pfa = false_alarms / O3_TEST_H0_TRIALS
            pfa_upper = wilson_upper_bound(false_alarms, O3_TEST_H0_TRIALS)
            models[(architecture.id, level)] = (;
                precision_applied,
                normalization,
                validation_factor,
                test_factor,
                threshold,
            )
            push!(threshold_rows, (
                architecture_id=architecture.id,
                nuisance_level=level,
                validation_trials=O3_VALIDATION_H0_TRIALS,
                test_trials=O3_TEST_H0_TRIALS,
                threshold,
                false_alarms,
                observed_pfa,
                pfa_upper_95=pfa_upper,
            ))
        end
    end
    models, threshold_rows
end

function scattering_tensors()
    rotation = [cos(pi / 4) -sin(pi / 4); sin(pi / 4) cos(pi / 4)]
    [
        Matrix{Float64}(I, 2, 2),
        [1.0 0.0; 0.0 0.15],
        Matrix{Float64}(rotation * Diagonal([1.0, 0.15]) * transpose(rotation)),
    ]
end

function signal_context(rng, base_hypotheses, nominal_config, context_index)
    perturbed = perturbed_config(rng, nominal_config)
    lateral_shift = 0.10 * (2 * rand(rng) - 1)
    tensors = scattering_tensors()
    hypotheses = DefectHypothesis[]
    for (index, hypothesis) in enumerate(base_hypotheses)
        tensor = tensors[mod1(index + context_index - 1, length(tensors))]
        push!(hypotheses, DefectHypothesis(
            hypothesis.id,
            hypothesis.x_wavelengths,
            hypothesis.y_wavelengths - lateral_shift,
            copy(tensor),
        ))
    end
    state_errors = [(
        amplitude=exp.(0.03 .* randn(rng, nominal_config.element_count, 4)),
        delay=0.02 .* randn(rng, nominal_config.element_count, 4),
    ) for _ in 1:2]
    (
        config=perturbed,
        hypotheses,
        state_errors,
        lateral_shift,
    )
end

function architecture_state_errors(architecture, context)
    architecture.kind == :one_static_state && return context.state_errors[1:1]
    architecture.kind == :two_static_states && return context.state_errors
    nothing
end

function localization_error(predicted_index, true_index, templates, true_hypotheses)
    hypot(
        templates[predicted_index].x_wavelengths -
            true_hypotheses[true_index].x_wavelengths,
        templates[predicted_index].y_wavelengths -
            true_hypotheses[true_index].y_wavelengths,
    )
end

function evaluate_signal_batch(
    rng,
    factor,
    mean_matrix,
    normalization,
    threshold,
    amplitudes,
    templates,
    true_hypotheses;
    replicates=1,
)
    true_count = size(mean_matrix, 2)
    true_count == length(true_hypotheses) ||
        throw(DimensionMismatch("one true hypothesis per signal column is required"))
    true_indices = repeat(collect(1:true_count); outer=replicates)
    noise = factor * complex_noise(rng, size(factor, 1), length(true_indices))
    means = mean_matrix[:, true_indices]
    rows = NamedTuple[]
    for amplitude in amplitudes
        scores, predicted = normalized_scores(
            noise .+ amplitude .* means,
            normalization,
        )
        detected = scores .> threshold
        errors = Float64[]
        localized = 0
        for trial in eachindex(true_indices)
            detected[trial] || continue
            error = localization_error(
                predicted[trial],
                true_indices[trial],
                templates,
                true_hypotheses,
            )
            push!(errors, error)
            localized += error <= 0.5 + 1.0e-12
        end
        push!(rows, (
            amplitude,
            trials=length(true_indices),
            detections=count(detected),
            localized,
            errors,
        ))
    end
    rows
end

function summarize_batch_rows(rows)
    [(
        amplitude=row.amplitude,
        trials=row.trials,
        detections=row.detections,
        localized=row.localized,
        pd=row.detections / row.trials,
        pdl=row.localized / row.trials,
        median_localization_error=isempty(row.errors) ? NaN : median(row.errors),
        q90_localization_error=isempty(row.errors) ? NaN : quantile(row.errors, 0.90),
    ) for row in rows]
end

function validation_reference_amplitudes(
    architectures,
    models,
    nominal_dictionaries,
    hypotheses,
    config,
)
    oracle = only(filter(architecture -> architecture.kind == :fmc, architectures))
    model = models[(oracle.id, :C3)]
    nominal_mean = matched_mean(
        model.precision_applied, nominal_dictionaries[oracle.id],
    )
    single_template_requirements = [
        sqrt(model.threshold * model.normalization[index]) /
            max(abs(nominal_mean[index, index]), eps(Float64))
        for index in eachindex(model.normalization)
    ]
    anchor = median(single_template_requirements)
    amplitudes = anchor .* O3_AMPLITUDE_MULTIPLIERS
    detections = zeros(Int, length(amplitudes))
    trials = zeros(Int, length(amplitudes))
    rng = Xoshiro(O3_RANDOM_SEED + UInt(0x5100))
    for context_index in 1:O3_VALIDATION_CONTEXTS
        context = signal_context(rng, hypotheses, config, context_index)
        library = response_library(context.hypotheses, context.config)
        dictionary = electronics_whiten(
            compress_pulse_library(library, oracle, context.config),
            oracle,
            context.config,
        )
        mean_matrix = matched_mean(model.precision_applied, dictionary)
        batch = evaluate_signal_batch(
            rng,
            model.validation_factor,
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
    reference_index = findfirst(>=(0.95), accumulate(max, probabilities))
    isnothing(reference_index) && error(
        "FMC validation Pd did not reach 0.95; extend the preregistered amplitude grid",
    )
    amplitudes, probabilities, amplitudes[reference_index], anchor
end

function nominal_detection_curves(
    architectures,
    models,
    nominal_dictionaries,
    hypotheses,
    amplitudes,
)
    rows = NamedTuple[]
    for (architecture_index, architecture) in enumerate(architectures)
        for (level_index, level) in enumerate((:C0, :C1, :C2))
            model = models[(architecture.id, level)]
            mean_matrix = matched_mean(
                model.precision_applied, nominal_dictionaries[architecture.id],
            )
            rng = Xoshiro(
                O3_RANDOM_SEED + UInt(0x60_000 + 1000 * architecture_index + level_index),
            )
            batch = summarize_batch_rows(evaluate_signal_batch(
                rng,
                model.test_factor,
                mean_matrix,
                model.normalization,
                model.threshold,
                amplitudes,
                hypotheses,
                hypotheses;
                replicates=O3_NOMINAL_REPLICATES,
            ))
            for row in batch
                push!(rows, merge((;
                    architecture_id=architecture.id,
                    nuisance_level=level,
                ), row))
            end
        end
    end
    rows
end

function c3_detection_campaign(
    architectures,
    models,
    hypotheses,
    config,
    amplitudes,
)
    counts = Dict{Tuple{Symbol,Float64},NamedTuple}()
    error_cache = Dict{Tuple{Symbol,Float64},Vector{Float64}}()
    for architecture in architectures, amplitude in amplitudes
        counts[(architecture.id, amplitude)] = (trials=0, detections=0, localized=0)
        error_cache[(architecture.id, amplitude)] = Float64[]
    end
    context_rows = NamedTuple[]
    rng = Xoshiro(O3_RANDOM_SEED + UInt(0x7100))
    for context_index in 1:O3_TEST_CONTEXTS
        context = signal_context(rng, hypotheses, config, 10_000 + context_index)
        library = response_library(context.hypotheses, context.config)
        for architecture in architectures
            model = models[(architecture.id, :C3)]
            dictionary = electronics_whiten(
                compress_pulse_library(
                    library,
                    architecture,
                    context.config;
                    state_errors=architecture_state_errors(architecture, context),
                ),
                architecture,
                context.config,
            )
            mean_matrix = matched_mean(model.precision_applied, dictionary)
            raw_batch = evaluate_signal_batch(
                rng,
                model.test_factor,
                mean_matrix,
                model.normalization,
                model.threshold,
                amplitudes,
                hypotheses,
                context.hypotheses,
            )
            summarized = summarize_batch_rows(raw_batch)
            for (raw, summary) in zip(raw_batch, summarized)
                key = (architecture.id, raw.amplitude)
                previous = counts[key]
                counts[key] = (
                    trials=previous.trials + raw.trials,
                    detections=previous.detections + raw.detections,
                    localized=previous.localized + raw.localized,
                )
                append!(error_cache[key], raw.errors)
                push!(context_rows, (
                    context_index,
                    architecture_id=architecture.id,
                    amplitude=raw.amplitude,
                    pd=summary.pd,
                    pdl=summary.pdl,
                    median_localization_error=summary.median_localization_error,
                    q90_localization_error=summary.q90_localization_error,
                ))
            end
        end
        context_index % 20 == 0 &&
            println("[+] O3 C3 contexts completed: $context_index / $O3_TEST_CONTEXTS")
    end

    aggregate_rows = NamedTuple[]
    for architecture in architectures, amplitude in amplitudes
        key = (architecture.id, amplitude)
        count_row = counts[key]
        errors = error_cache[key]
        push!(aggregate_rows, (
            architecture_id=architecture.id,
            nuisance_level=:C3,
            amplitude,
            trials=count_row.trials,
            detections=count_row.detections,
            localized=count_row.localized,
            pd=count_row.detections / count_row.trials,
            pdl=count_row.localized / count_row.trials,
            median_localization_error=isempty(errors) ? NaN : median(errors),
            q90_localization_error=isempty(errors) ? NaN : quantile(errors, 0.90),
        ))
    end
    aggregate_rows, context_rows
end

function context_alpha90_rows(context_curve_rows, architectures, amplitudes, reference_alpha)
    rows = NamedTuple[]
    for context_index in 1:O3_TEST_CONTEXTS, architecture in architectures
        selected = filter(row ->
            row.context_index == context_index &&
            row.architecture_id == architecture.id,
            context_curve_rows,
        )
        sort!(selected; by=row -> row.amplitude)
        alpha90 = crossing_amplitude(
            getproperty.(selected, :amplitude), getproperty.(selected, :pd), 0.90,
        )
        reference = only(filter(row -> isapprox(
            row.amplitude, reference_alpha; rtol=1.0e-12,
        ), selected))
        push!(rows, (
            context_index,
            architecture_id=architecture.id,
            alpha90,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            median_localization_error_reference=reference.median_localization_error,
            q90_localization_error_reference=reference.q90_localization_error,
        ))
    end
    rows
end

function aggregate_alpha90(curve_rows, architecture_id, level)
    selected = sort(
        filter(row ->
            row.architecture_id == architecture_id && row.nuisance_level == level,
            curve_rows,
        );
        by=row -> row.amplitude,
    )
    crossing_amplitude(
        getproperty.(selected, :amplitude), getproperty.(selected, :pd), 0.90,
    )
end

function o3_summary(
    architectures,
    threshold_rows,
    c3_rows,
    context_summaries,
    amplitudes,
    reference_alpha,
)
    alpha90_values = Dict(
        architecture.id => aggregate_alpha90(c3_rows, architecture.id, :C3)
        for architecture in architectures
    )
    a_scan_alpha = alpha90_values[:a_scan_uniform]
    oracle_alpha = alpha90_values[:paut_fmc_oracle]
    candidate_alpha = alpha90_values[:four_channel_two_static_states]
    mechanical_alpha = alpha90_values[:mechanical_scan_equal_time]
    aggregate_gap = isfinite(a_scan_alpha) && isfinite(oracle_alpha) &&
        a_scan_alpha > oracle_alpha ?
        (a_scan_alpha - candidate_alpha) / (a_scan_alpha - oracle_alpha) : NaN

    context_gaps = Float64[]
    for context_index in 1:O3_TEST_CONTEXTS
        lookup(id) = only(filter(row ->
            row.context_index == context_index && row.architecture_id == id,
            context_summaries,
        )).alpha90
        a_value = lookup(:a_scan_uniform)
        candidate_value = lookup(:four_channel_two_static_states)
        oracle_value = lookup(:paut_fmc_oracle)
        if all(isfinite, (a_value, candidate_value, oracle_value)) &&
           a_value > oracle_value + eps(Float64)
            push!(context_gaps,
                (a_value - candidate_value) / (a_value - oracle_value),
            )
        end
    end
    q10_gap = isempty(context_gaps) ? NaN : quantile(context_gaps, 0.10)

    rows = NamedTuple[]
    for architecture in architectures
        threshold = only(filter(row ->
            row.architecture_id == architecture.id && row.nuisance_level == :C3,
            threshold_rows,
        ))
        reference = only(filter(row ->
            row.architecture_id == architecture.id &&
            row.nuisance_level == :C3 &&
            isapprox(row.amplitude, reference_alpha; rtol=1.0e-12),
            c3_rows,
        ))
        is_candidate = architecture.id == :four_channel_two_static_states
        gate_pfa = threshold.pfa_upper_95 <= 0.0125
        gate_pd = reference.pd >= 0.90
        gate_pdl = reference.pdl >= 0.80
        gap_context_coverage = length(context_gaps) / O3_TEST_CONTEXTS
        gate_gap = is_candidate && gap_context_coverage >= 0.90 &&
            isfinite(q10_gap) && q10_gap >= 0.50
        gate_localization = reference.q90_localization_error <= 1.0
        gate_mechanical = is_candidate && isfinite(candidate_alpha) &&
            isfinite(mechanical_alpha) && candidate_alpha <= mechanical_alpha
        passed = is_candidate && gate_pfa && gate_pd && gate_pdl && gate_gap &&
            gate_localization && gate_mechanical
        push!(rows, (
            architecture_id=architecture.id,
            electronic_channels=architecture.electronic_channels,
            transmit_events=architecture.transmit_events,
            physical_states=architecture.physical_states,
            mechanical_positions=architecture.mechanical_positions,
            reference_alpha,
            observed_pfa=threshold.observed_pfa,
            pfa_upper_95=threshold.pfa_upper_95,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            median_localization_error_reference=reference.median_localization_error,
            q90_localization_error_reference=reference.q90_localization_error,
            alpha90=alpha90_values[architecture.id],
            aggregate_paut_gap_closure=is_candidate ? aggregate_gap : NaN,
            q10_context_paut_gap_closure=is_candidate ? q10_gap : NaN,
            valid_context_gap_count=is_candidate ? length(context_gaps) : 0,
            context_gap_coverage=is_candidate ? gap_context_coverage : NaN,
            gate_pfa,
            gate_pd,
            gate_pdl,
            gate_gap,
            gate_localization,
            gate_mechanical_noninferiority=gate_mechanical,
            passed,
        ))
    end
    rows, context_gaps
end

function architecture_metadata(architectures)
    [(
        architecture_id=architecture.id,
        role=architecture.role,
        electronic_channels=architecture.electronic_channels,
        transmit_events=architecture.transmit_events,
        physical_states=architecture.physical_states,
        mechanical_positions=architecture.mechanical_positions,
        complex_features=architecture.state_blocks * architecture.block_size * 7,
        energy_normalization=architecture.kind == :two_static_states ?
            "1/sqrt(2) across states" :
            architecture.kind == :mechanical_scan ? "1/sqrt(8) across positions" :
            "unit total transmit energy",
        lift_status=architecture.kind in (:one_static_state, :two_static_states) ?
            "group_sector+translation_kz" : "plane_strain_surrogate",
    ) for architecture in architectures]
end

function plot_o3_results(curve_rows, summary_rows, amplitudes, reference_alpha)
    colors = Dict(
        :a_scan_uniform => :black,
        :mechanical_scan_equal_time => :darkorange,
        :four_channel_one_static_state => :gray,
        :four_channel_two_static_states => :steelblue,
        :paut_fmc_oracle => :firebrick,
    )
    detection = plot(
        xscale=:log10,
        xlabel="scattering amplitude alpha",
        ylabel="Pd at frozen Pfa threshold",
        title="O3 C3 detection",
        legend=:bottomright,
        ylim=(0, 1.02),
    )
    localization = plot(
        xscale=:log10,
        xlabel="scattering amplitude alpha",
        ylabel="P(detect and localize <= 0.5 lambda_P)",
        title="O3 C3 joint decision",
        legend=false,
        ylim=(0, 1.02),
    )
    for architecture_id in keys(colors)
        selected = sort(filter(row ->
            row.architecture_id == architecture_id && row.nuisance_level == :C3,
            curve_rows,
        ); by=row -> row.amplitude)
        plot!(
            detection,
            getproperty.(selected, :amplitude),
            getproperty.(selected, :pd);
            marker=:circle,
            linewidth=2,
            color=colors[architecture_id],
            label=String(architecture_id),
        )
        plot!(
            localization,
            getproperty.(selected, :amplitude),
            getproperty.(selected, :pdl);
            marker=:circle,
            linewidth=2,
            color=colors[architecture_id],
        )
    end
    hline!(detection, [0.90]; linestyle=:dash, color=:black, label="Pd gate")
    vline!(detection, [reference_alpha]; linestyle=:dot, color=:purple, label="FMC validation alpha")
    hline!(localization, [0.80]; linestyle=:dash, color=:black)
    vline!(localization, [reference_alpha]; linestyle=:dot, color=:purple)
    savefig(
        plot(detection, localization; layout=(1, 2), size=(1450, 560), margin=7 * Plots.mm),
        joinpath(O3_OUTPUT_ROOT, "o3_detection_curves.png"),
    )

    finite_rows = filter(row -> isfinite(row.alpha90), summary_rows)
    alpha_plot = bar(
        String.(getproperty.(finite_rows, :architecture_id)),
        getproperty.(finite_rows, :alpha90);
        ylabel="alpha90 (lower is better)",
        title="O3 C3 detection threshold",
        label=false,
        xrotation=24,
        size=(1050, 520),
    )
    savefig(alpha_plot, joinpath(O3_OUTPUT_ROOT, "o3_alpha90.png"))
end

function main_o3()
    mkpath(O3_OUTPUT_ROOT)
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    architectures = o3_architectures(config)
    nominal_library = response_library(hypotheses, config)
    nominal_dictionaries = Dict(
        architecture.id => electronics_whiten(
            compress_pulse_library(nominal_library, architecture, config),
            architecture,
            config,
        )
        for architecture in architectures
    )
    println("[+] frozen O3 architectures and nominal dictionaries assembled")

    nuisance_features, nuisance_scale_rows = build_nuisance_features(
        architectures, config,
    )
    println("[+] independent train/validation/test nuisance libraries assembled")
    models, threshold_rows = build_detection_models(
        architectures, nominal_dictionaries, nuisance_features,
    )
    println("[+] frozen GLRT thresholds calibrated on $(O3_VALIDATION_H0_TRIALS) H0 trials")

    amplitudes, oracle_validation_pd, reference_alpha, amplitude_anchor =
        validation_reference_amplitudes(
            architectures, models, nominal_dictionaries, hypotheses, config,
        )
    println(
        "[+] FMC validation reference alpha=$(round(reference_alpha; digits=4)); " *
        "Pd=$(round(oracle_validation_pd[findfirst(==(reference_alpha), amplitudes)]; digits=3))",
    )

    nominal_rows = nominal_detection_curves(
        architectures, models, nominal_dictionaries, hypotheses, amplitudes,
    )
    c3_rows, context_curve_rows = c3_detection_campaign(
        architectures, models, hypotheses, config, amplitudes,
    )
    all_curve_rows = vcat(nominal_rows, c3_rows)
    context_summaries = context_alpha90_rows(
        context_curve_rows, architectures, amplitudes, reference_alpha,
    )
    summary_rows, context_gaps = o3_summary(
        architectures,
        threshold_rows,
        c3_rows,
        context_summaries,
        amplitudes,
        reference_alpha,
    )

    metadata_rows = architecture_metadata(architectures)
    write_csv(
        joinpath(O3_OUTPUT_ROOT, "architecture_metadata.csv"),
        propertynames(first(metadata_rows)),
        metadata_rows,
    )
    write_csv(
        joinpath(O3_OUTPUT_ROOT, "nuisance_scales.csv"),
        propertynames(first(nuisance_scale_rows)),
        nuisance_scale_rows,
    )
    write_csv(
        joinpath(O3_OUTPUT_ROOT, "thresholds.csv"),
        propertynames(first(threshold_rows)),
        threshold_rows,
    )
    write_csv(
        joinpath(O3_OUTPUT_ROOT, "detection_curves.csv"),
        propertynames(first(all_curve_rows)),
        all_curve_rows,
    )
    write_csv(
        joinpath(O3_OUTPUT_ROOT, "c3_context_curves.csv"),
        propertynames(first(context_curve_rows)),
        context_curve_rows,
    )
    write_csv(
        joinpath(O3_OUTPUT_ROOT, "c3_context_summary.csv"),
        propertynames(first(context_summaries)),
        context_summaries,
    )
    write_csv(
        joinpath(O3_OUTPUT_ROOT, "o3_summary.csv"),
        propertynames(first(summary_rows)),
        summary_rows,
    )
    jldsave(
        joinpath(O3_OUTPUT_ROOT, "o3_paired_clutter.jld2");
        config,
        hypotheses,
        metadata_rows,
        nuisance_scale_rows,
        threshold_rows,
        amplitudes,
        amplitude_anchor,
        oracle_validation_pd,
        reference_alpha,
        all_curve_rows,
        context_summaries,
        summary_rows,
        context_gaps,
        random_seed=O3_RANDOM_SEED,
        validation_h0_trials=O3_VALIDATION_H0_TRIALS,
        test_h0_trials=O3_TEST_H0_TRIALS,
        validation_contexts=O3_VALIDATION_CONTEXTS,
        test_contexts=O3_TEST_CONTEXTS,
    )
    plot_o3_results(all_curve_rows, summary_rows, amplitudes, reference_alpha)

    for row in summary_rows
        println(
            "[+] $(row.architecture_id): Pfa95=$(round(row.pfa_upper_95; digits=4)), " *
            "Pd_ref=$(round(row.pd_reference; digits=3)), " *
            "Pdl_ref=$(round(row.pdl_reference; digits=3)), " *
            "alpha90=$(round(row.alpha90; digits=3)), passed=$(row.passed)",
        )
    end
    candidate = only(filter(row ->
        row.architecture_id == :four_channel_two_static_states,
        summary_rows,
    ))
    println(
        "[+] candidate Q10 context PAUT-gap closure=" *
        "$(round(candidate.q10_context_paut_gap_closure; digits=3))",
    )
    println("[+] wrote reduced O3 results to $O3_OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_o3()
end
