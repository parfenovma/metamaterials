using LinearAlgebra
using Random
using Statistics
using JLD2
using Plots

include(joinpath(@__DIR__, "run_operator_o3_static_frontier.jl"))

const O7_OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o7_global_state_v0",
)
const O7_STATIC_BUNDLE = joinpath(
    PROJECT_ROOT, "results", "operator_space", "static_frontier_v0",
    "static_frontier.jld2",
)
const O7_TARGET_ID = :trained_k3_q4_reflection_balanced
const O7_BASE_SEED = UInt(0x0707a117)
const O7_SEEDS = (
    nuisance_train=O7_BASE_SEED + UInt(0x101),
    nuisance_validation=O7_BASE_SEED + UInt(0x202),
    nuisance_test=O7_BASE_SEED + UInt(0x303),
    signal_validation=O7_BASE_SEED + UInt(0x404),
    signal_test=O7_BASE_SEED + UInt(0x505),
    h0_validation=O7_BASE_SEED + UInt(0x606),
    h0_test=O7_BASE_SEED + UInt(0x707),
    candidate_search=O7_BASE_SEED + UInt(0x808),
    proxy_h0=O7_BASE_SEED + UInt(0x909),
)
const O7_VALIDATION_H0_TRIALS = 10_000
const O7_TEST_H0_TRIALS = 10_000
const O7_VALIDATION_CONTEXTS = 24
const O7_TEST_CONTEXTS = 60
const O7_PROXY_CONTEXTS = 4
const O7_SCREEN_SAMPLES = 800
const O7_OPERATOR_FINALISTS = 8
const O7_CARRIER_HZ = 242_000.0
const O7_PAIR_CAPACITANCE_F = 4.0e-9
const O7_BITS = ((0, 0), (1, 0), (0, 1), (1, 1))
const O7_PERMUTATIONS = let values = NTuple{4,Int}[]
    for a in 1:4, b in 1:4, c in 1:4, d in 1:4
        length(Set((a, b, c, d))) == 4 && push!(values, (a, b, c, d))
    end
    Tuple(values)
end

o7_complex_noise(rng, rows, columns) =
    (randn(rng, rows, columns) .+ im .* randn(rng, rows, columns)) ./ sqrt(2)

function o7_architecture(id, family, frequency_states, role; mechanism=:none)
    Q = length(frequency_states)
    K = size(first(first(frequency_states)), 2)
    (
        id=Symbol(id), kind=:frequency_static, electronic_channels=K,
        physical_states=Q, transmit_events=K * Q, block_size=K^2,
        state_blocks=Q, family, states=(), frequency_states, role, mechanism,
    )
end

function o7_reflection_matrix(count)
    Matrix{Float64}(I, count, count)[reverse(1:count), :]
end

function o7_zone_operator(position, width, basis, config)
    x = element_positions(config) ./ maximum(abs, element_positions(config))
    profile = @. exp(-((abs(x) - position) / width)^2)
    profile ./= maximum(profile)
    Matrix(Symmetric(basis' * Diagonal(profile) * basis))
end

function o7_input_coupling(mode_count)
    raw = zeros(Float64, mode_count, 3)
    raw[[1, 3, 5], 1] .= [1.0, 0.30, 0.08]
    raw[[3, 5, 7], 2] .= [1.0, -0.32, 0.12]
    raw[[2, 4, 6], 3] .= [1.0, 0.35, -0.10]
    result = zeros(Float64, mode_count, 3)
    result[:, 1:2] .= Matrix(qr(raw[:, 1:2]).Q)[:, 1:2]
    result[:, 3] .= raw[:, 3] ./ norm(raw[:, 3])
    result
end

function o7_feedback(kind, bit, frequency, strength, resonance, quality, loss)
    bit == 0 && return 0.0 + 0.0im
    if kind == :pzt
        denominator = resonance^2 - frequency^2 -
            im * resonance * frequency / quality
        strength * resonance^2 / denominator
    elseif kind == :mechanical
        strength * (1 + im * loss * frequency)
    else
        throw(ArgumentError("unknown O7 mechanism: $kind"))
    end
end

function o7_transfer_states(kind, parameters, config; zero_control=false)
    mode_count = 11
    full_basis = cosine_ports(config.element_count, config)
    basis = full_basis[:, 1:mode_count]
    coupling = o7_input_coupling(mode_count)
    zone_a = o7_zone_operator(
        parameters.position_a, parameters.width_a, basis, config,
    )
    zone_b = o7_zone_operator(
        parameters.position_b, parameters.width_b, basis, config,
    )
    mode_frequencies = collect(range(0.70, 1.30; length=mode_count)) .+
        parameters.host_shift .* sin.(range(0, 2pi; length=mode_count))
    raw_states = Vector{Vector{Matrix{ComplexF64}}}()
    for bits in O7_BITS
        state = Matrix{ComplexF64}[]
        for frequency in config.frequency_ratios
            dynamic = Diagonal(complex.(
                mode_frequencies .^ 2 .- frequency^2,
                2 .* parameters.host_damping .* mode_frequencies .* frequency,
            ))
            strength_a = zero_control ? 0.0 : parameters.strength_a
            strength_b = zero_control ? 0.0 : parameters.strength_b
            feedback_a = o7_feedback(
                kind, bits[1], frequency, strength_a, parameters.resonance_a,
                parameters.quality_a, parameters.loss_a,
            )
            feedback_b = o7_feedback(
                kind, bits[2], frequency, strength_b, parameters.resonance_b,
                parameters.quality_b, parameters.loss_b,
            )
            resolvent = dynamic + feedback_a .* zone_a + feedback_b .* zone_b
            push!(state, ComplexF64.(basis * (resolvent \ coupling)))
        end
        push!(raw_states, state)
    end
    maximum_sigma = maximum(opnorm(ports) for state in raw_states for ports in state)
    common_scale = min(1.0, 0.98 / maximum_sigma)
    states = [[common_scale .* ports for ports in state] for state in raw_states]
    states, common_scale, maximum_sigma
end

function o7_projector(ports)
    decomposition = svd(ports; full=false)
    rank = count(>(maximum(decomposition.S) * 1.0e-10), decomposition.S)
    rank == size(ports, 2) || return zeros(ComplexF64, size(ports, 1), size(ports, 1))
    vectors = decomposition.U[:, 1:rank]
    vectors * vectors'
end

function o7_target_fit(states, target_states)
    target_projectors = o7_projector.(target_states)
    costs = zeros(Float64, 4, 4)
    for physical in 1:4, target in 1:4
        costs[physical, target] = mean(
            norm(o7_projector(ports) - target_projectors[target])^2 / 6
            for ports in states[physical]
        )
    end
    best = minimum(
        mean(costs[index, permutation[index]] for index in 1:4)
        for permutation in O7_PERMUTATIONS
    )
    center = cld(length(first(states)), 2)
    projectors = [o7_projector(state[center]) for state in states]
    diversity = minimum(
        norm(projectors[first] - projectors[second]) / sqrt(6)
        for first in 1:4 for second in (first + 1):4
    )
    best, diversity
end

function o7_random_parameters(rng, kind)
    position_a = 0.15 + 0.30 * rand(rng)
    position_b = 0.55 + 0.38 * rand(rng)
    width_a = 0.08 + 0.18 * rand(rng)
    width_b = 0.08 + 0.18 * rand(rng)
    host_shift = -0.035 + 0.070 * rand(rng)
    host_damping = 0.025 + 0.055 * rand(rng)
    if kind == :pzt
        strength_a = 10.0^(-2.0 + 1.45 * rand(rng))
        strength_b = 10.0^(-2.0 + 1.45 * rand(rng))
        resonance_a = 0.84 + 0.32 * rand(rng)
        resonance_b = 0.84 + 0.32 * rand(rng)
        quality_a = 1.4 + 5.6 * rand(rng)
        quality_b = 1.4 + 5.6 * rand(rng)
        loss_a = 0.0
        loss_b = 0.0
    else
        strength_a = 10.0^(-1.5 + 1.20 * rand(rng))
        strength_b = 10.0^(-1.5 + 1.20 * rand(rng))
        resonance_a = 1.0
        resonance_b = 1.0
        quality_a = 1.0
        quality_b = 1.0
        loss_a = 0.02 + 0.13 * rand(rng)
        loss_b = 0.02 + 0.13 * rand(rng)
    end
    (;
        position_a, position_b, width_a, width_b, host_shift, host_damping,
        strength_a, strength_b, resonance_a, resonance_b, quality_a,
        quality_b, loss_a, loss_b,
    )
end

function o7_apply_errors(ports, frequency, errors)
    raw = ports .* errors.amplitude .* exp.(-2pi * im * frequency .* errors.delay)
    raw ./ max(1.0, opnorm(raw))
end

function o7_compress_library(library, architecture, config; state_errors=nothing)
    architecture.kind == :frequency_static ||
        return sf_compress_library(
            library, architecture, config; state_errors=state_errors,
        )
    Q = architecture.physical_states
    K = architecture.electronic_channels
    errors = isnothing(state_errors) ? [(
        amplitude=ones(config.element_count, K),
        delay=zeros(config.element_count, K),
    ) for _ in 1:Q] : state_errors
    length(errors) == Q || throw(DimensionMismatch("one error draw per state"))
    state_features = Matrix{ComplexF64}[]
    for state_index in 1:Q
        raw = Matrix{ComplexF64}(undef,
            length(config.frequency_ratios) * K^2, size(library, 4),
        )
        for (frequency_index, frequency) in enumerate(config.frequency_ratios)
            ports = o7_apply_errors(
                architecture.frequency_states[state_index][frequency_index],
                frequency, errors[state_index],
            )
            opnorm(ports) <= 1 + 1.0e-10 || error("O7 contraction failed")
            row = (frequency_index - 1) * K^2 + 1
            for hypothesis_index in axes(library, 4)
                response = view(library, :, :, frequency_index, hypothesis_index)
                compressed = transpose(ports) * response * ports / sqrt(K)
                destination = row
                for transmit in 1:K, receive in 1:K
                    raw[destination, hypothesis_index] = compressed[receive, transmit]
                    destination += 1
                end
            end
        end
        push!(state_features, sf_pulse_features(raw, config, K^2))
    end
    features = reduce(vcat, state_features) ./ sqrt(Q)
    sf_electronics_whiten(features, K^2, Q, config.frequency_ratios)
end

function o7_context_errors(architecture, context)
    architecture.kind in (:static, :frequency_static) || return nothing
    K = architecture.electronic_channels
    [(
        amplitude=context.amplitude_errors[state][:, 1:K],
        delay=context.delay_errors[state][:, 1:K],
    ) for state in 1:architecture.physical_states]
end

function o7_signal_contexts(config, hypotheses, split, count)
    seed = split == :validation ? O7_SEEDS.signal_validation : O7_SEEDS.signal_test
    rng = Xoshiro(seed)
    rotation = [cos(pi / 4) -sin(pi / 4); sin(pi / 4) cos(pi / 4)]
    tensors = [
        Matrix{Float64}(I, 2, 2), [1.0 0.0; 0.0 0.15],
        Matrix{Float64}(rotation * Diagonal([1.0, 0.15]) * transpose(rotation)),
    ]
    contexts = NamedTuple[]
    for context_index in 1:count
        perturbed = perturbed_config(rng, config)
        lateral_shift = 0.10 * (2 * rand(rng) - 1)
        true_hypotheses = [
            DefectHypothesis(
                hypothesis.id, hypothesis.x_wavelengths,
                hypothesis.y_wavelengths - lateral_shift,
                copy(tensors[mod1(index + context_index - 1, length(tensors))]),
            )
            for (index, hypothesis) in enumerate(hypotheses)
        ]
        amplitude_errors = [
            exp.(0.03 .* randn(rng, config.element_count, 4)) for _ in 1:4
        ]
        delay_errors = [
            0.02 .* randn(rng, config.element_count, 4) for _ in 1:4
        ]
        push!(contexts, (
            context_index, config=perturbed, hypotheses=true_hypotheses,
            library=response_library(true_hypotheses, perturbed),
            amplitude_errors, delay_errors, lateral_shift,
        ))
    end
    contexts
end

function o7_nuisance_contract(config)
    libraries = Dict(
        :train => sf_nuisance_libraries(config, O7_SEEDS.nuisance_train),
        :validation => sf_nuisance_libraries(config, O7_SEEDS.nuisance_validation),
        :test => sf_nuisance_libraries(config, O7_SEEDS.nuisance_test),
    )
    fmc = sf_architecture(:o7_scaling_fmc, :fmc, 17, 1, :full_aperture, (), :scaling)
    scales = Dict{Symbol,Float64}()
    rows = NamedTuple[]
    for category in (:backwall, :ringing, :material)
        features = sf_compress_library(libraries[:train][category], fmc, config)
        target = SF_CATEGORY_POWER[category]
        scale = sqrt(target * size(features, 1) / sum(abs2, features))
        scales[category] = scale
        push!(rows, (
            category, target_fmc_power_ratio=target, physical_scale=scale,
            achieved_fmc_power_ratio=sum(abs2, scale .* features) / size(features, 1),
        ))
    end
    libraries, scales, rows
end

function o7_clutter_features(architecture, split_libraries, scales, config)
    hcat((
        scales[category] .* o7_compress_library(
            split_libraries[category], architecture, config,
        ) for category in (:backwall, :ringing, :material)
    )...)
end

function o7_proxy(
    architecture, config, hypotheses, nuisance_libraries, nuisance_scales,
    contexts, proxy_index,
)
    nominal = o7_compress_library(response_library(hypotheses, config), architecture, config)
    train_clutter = o7_clutter_features(
        architecture, nuisance_libraries[:train], nuisance_scales, config,
    )
    validation_clutter = o7_clutter_features(
        architecture, nuisance_libraries[:validation], nuisance_scales, config,
    )
    precision_applied = precision_dictionary(nominal, train_clutter)
    gram = matched_mean(precision_applied, nominal)
    normalization = max.(real.(diag(gram)), eps(Float64))
    covariance = matched_output_covariance(precision_applied, validation_clutter)
    factor = OperatorO3Statistics.covariance_factor(covariance)
    rng = Xoshiro(O7_SEEDS.proxy_h0 + UInt(10_000 * proxy_index))
    scores, _ = normalized_scores(
        factor * o7_complex_noise(rng, size(factor, 1), 2000), normalization,
    )
    threshold = quantile(scores, 0.99)
    requirements = Float64[]
    fidelities = Float64[]
    for context in contexts
        dictionary = o7_compress_library(
            context.library, architecture, context.config;
            state_errors=o7_context_errors(architecture, context),
        )
        mean_matrix = matched_mean(precision_applied, dictionary)
        append!(requirements, [
            sqrt(threshold * normalization[index]) /
                max(abs(mean_matrix[index, index]), eps(Float64))
            for index in eachindex(normalization)
        ])
        push!(fidelities, calibration_fidelity(
            nominal, dictionary,
            Matrix{Float64}(I, size(dictionary, 1), size(dictionary, 1)),
        ).q10)
    end
    summary = sf_cheap_summary(precision_applied, sf_nonlocal_mask(hypotheses))
    (
        proxy_alpha90=quantile(requirements, 0.90),
        q10_fidelity=quantile(fidelities, 0.10),
        residual_minimum_detectability=summary.minimum_detectability,
        residual_nonlocal_coherence=summary.nonlocal_coherence,
        residual_effective_rank=summary.effective_rank,
        validation_objective=quantile(requirements, 0.90) /
            max(quantile(fidelities, 0.10), 0.50),
    )
end

function o7_screen_branch(kind, target_states, config, nominal_library, hypotheses)
    rng = Xoshiro(O7_SEEDS.candidate_search + (kind == :pzt ? UInt(0) : UInt(1_000_000)))
    screen_rows = NamedTuple[]
    parameter_cache = Dict{Int,NamedTuple}()
    state_cache = Dict{Int,Vector{Vector{Matrix{ComplexF64}}}}()
    for sample in 1:O7_SCREEN_SAMPLES
        parameters = o7_random_parameters(rng, kind)
        states, common_scale, raw_maximum_sigma = o7_transfer_states(kind, parameters, config)
        target_fit, state_diversity = o7_target_fit(states, target_states)
        objective = target_fit + 0.10 / max(state_diversity, 1.0e-4)
        row = (;
            mechanism=kind, sample, target_fit, state_diversity,
            screen_objective=objective, common_scale, raw_maximum_sigma,
        )
        push!(screen_rows, row)
        parameter_cache[sample] = parameters
        state_cache[sample] = states
    end
    ranked = sort(screen_rows; by=row -> row.screen_objective)
    operator_rows = NamedTuple[]
    architectures = NamedTuple[]
    for row in first(ranked, O7_OPERATOR_FINALISTS)
        architecture = o7_architecture(
            Symbol("$(kind)_screen_$(row.sample)"), Symbol("$(kind)_2bit"),
            state_cache[row.sample], :search; mechanism=kind,
        )
        dictionary = o7_compress_library(nominal_library, architecture, config)
        summary = sf_cheap_summary(dictionary, sf_nonlocal_mask(hypotheses))
        push!(operator_rows, merge(row, (
            minimum_detectability=summary.minimum_detectability,
            nonlocal_coherence=summary.nonlocal_coherence,
            effective_rank=summary.effective_rank,
            nominal_score=summary.nominal_score,
        )))
        push!(architectures, architecture)
    end
    screen_rows, operator_rows, architectures, parameter_cache
end

function o7_freeze_branch(
    kind, target_states, config, nominal_library, hypotheses,
    nuisance_libraries, nuisance_scales, validation_contexts,
)
    screen_rows, operator_rows, architectures, parameter_cache = o7_screen_branch(
        kind, target_states, config, nominal_library, hypotheses,
    )
    proxy_rows = NamedTuple[]
    for (index, architecture) in enumerate(architectures)
        proxy = o7_proxy(
            architecture, config, hypotheses, nuisance_libraries,
            nuisance_scales, first(validation_contexts, O7_PROXY_CONTEXTS), index,
        )
        sample = parse(Int, last(split(String(architecture.id), '_')))
        operator = only(filter(row -> row.sample == sample, operator_rows))
        push!(proxy_rows, merge(operator, proxy))
        println(
            "[+] O7 $(kind) finalist $(index)/$(length(architectures)): " *
            "sample=$(sample), proxy=$(round(proxy.validation_objective; digits=2))",
        )
    end
    winner_row = first(sort(proxy_rows; by=row -> row.validation_objective))
    winner_parameters = parameter_cache[winner_row.sample]
    winner_states, common_scale, raw_maximum_sigma = o7_transfer_states(
        kind, winner_parameters, config,
    )
    winner = o7_architecture(
        Symbol("$(kind)_k3q4_frozen"), Symbol("$(kind)_2bit"), winner_states,
        :frozen_physical; mechanism=kind,
    )
    winner, winner_parameters, winner_row, screen_rows, operator_rows, proxy_rows,
        common_scale, raw_maximum_sigma
end

function o7_detection_model(
    architecture, dictionary, train_clutter, validation_clutter, test_clutter,
    architecture_index,
)
    precision_applied = precision_dictionary(dictionary, train_clutter)
    gram = matched_mean(precision_applied, dictionary)
    normalization = max.(real.(diag(gram)), eps(Float64))
    validation_covariance = matched_output_covariance(precision_applied, validation_clutter)
    test_covariance = matched_output_covariance(precision_applied, test_clutter)
    validation_factor = OperatorO3Statistics.covariance_factor(validation_covariance)
    test_factor = OperatorO3Statistics.covariance_factor(test_covariance)
    validation_rng = Xoshiro(O7_SEEDS.h0_validation + UInt(10_000 * architecture_index))
    test_rng = Xoshiro(O7_SEEDS.h0_test + UInt(10_000 * architecture_index))
    validation_scores, _ = normalized_scores(
        validation_factor * o7_complex_noise(
            validation_rng, size(validation_factor, 1), O7_VALIDATION_H0_TRIALS,
        ), normalization,
    )
    threshold = quantile(validation_scores, 0.99)
    test_scores, _ = normalized_scores(
        test_factor * o7_complex_noise(
            test_rng, size(test_factor, 1), O7_TEST_H0_TRIALS,
        ), normalization,
    )
    false_alarms = count(>(threshold), test_scores)
    (
        precision_applied, normalization, validation_factor, test_factor,
        threshold, false_alarms,
        observed_pfa=false_alarms / O7_TEST_H0_TRIALS,
        pfa_upper_95=wilson_upper_bound(false_alarms, O7_TEST_H0_TRIALS),
    )
end

function o7_build_models(architectures, config, hypotheses, libraries, scales)
    nominal_library = response_library(hypotheses, config)
    dictionaries = Dict{Symbol,Matrix{ComplexF64}}()
    models = Dict{Symbol,NamedTuple}()
    rows = NamedTuple[]
    for (index, architecture) in enumerate(architectures)
        dictionary = o7_compress_library(nominal_library, architecture, config)
        dictionaries[architecture.id] = dictionary
        clutters = Dict(split => o7_clutter_features(
            architecture, libraries[split], scales, config,
        ) for split in (:train, :validation, :test))
        model = o7_detection_model(
            architecture, dictionary, clutters[:train], clutters[:validation],
            clutters[:test], index,
        )
        models[architecture.id] = model
        push!(rows, (
            architecture_id=architecture.id,
            validation_trials=O7_VALIDATION_H0_TRIALS,
            test_trials=O7_TEST_H0_TRIALS, threshold=model.threshold,
            false_alarms=model.false_alarms, observed_pfa=model.observed_pfa,
            pfa_upper_95=model.pfa_upper_95,
        ))
        println(
            "[+] O7 model $(index)/$(length(architectures)): $(architecture.id), " *
            "Pfa95=$(round(model.pfa_upper_95; digits=4))",
        )
    end
    dictionaries, models, rows
end

function o7_reference_amplitudes(oracle, model, dictionary, contexts, hypotheses)
    nominal_mean = matched_mean(model.precision_applied, dictionary)
    requirements = [
        sqrt(model.threshold * model.normalization[index]) /
            max(abs(nominal_mean[index, index]), eps(Float64))
        for index in eachindex(model.normalization)
    ]
    anchor = median(requirements)
    amplitudes = anchor .* SF_AMPLITUDE_MULTIPLIERS
    detections = zeros(Int, length(amplitudes))
    trials = zeros(Int, length(amplitudes))
    for context in contexts
        dictionary_true = o7_compress_library(context.library, oracle, context.config)
        mean_matrix = matched_mean(model.precision_applied, dictionary_true)
        rng = Xoshiro(O7_SEEDS.signal_validation + UInt(context.context_index))
        batch = sf_signal_batch(
            rng, model.validation_factor, mean_matrix, model.normalization,
            model.threshold, amplitudes, hypotheses, context.hypotheses,
        )
        for (index, row) in enumerate(batch)
            detections[index] += row.detections
            trials[index] += row.trials
        end
    end
    probabilities = detections ./ trials
    reference_index = findfirst(>=(0.95), accumulate(max, probabilities))
    isnothing(reference_index) && error("O7 FMC validation Pd did not reach 0.95")
    amplitudes, probabilities, amplitudes[reference_index], anchor
end

function o7_test_campaign(architectures, models, contexts, hypotheses, amplitudes)
    counts = Dict{Tuple{Symbol,Float64},NamedTuple}()
    errors = Dict{Tuple{Symbol,Float64},Vector{Float64}}()
    for architecture in architectures, amplitude in amplitudes
        counts[(architecture.id, amplitude)] = (trials=0, detections=0, localized=0)
        errors[(architecture.id, amplitude)] = Float64[]
    end
    context_rows = NamedTuple[]
    for context in contexts
        for (architecture_index, architecture) in enumerate(architectures)
            model = models[architecture.id]
            dictionary = o7_compress_library(
                context.library, architecture, context.config;
                state_errors=o7_context_errors(architecture, context),
            )
            mean_matrix = matched_mean(model.precision_applied, dictionary)
            rng = Xoshiro(
                O7_SEEDS.signal_test + UInt(10_000 * context.context_index + architecture_index),
            )
            batch = sf_signal_batch(
                rng, model.test_factor, mean_matrix, model.normalization,
                model.threshold, amplitudes, hypotheses, context.hypotheses,
            )
            for row in batch
                key = (architecture.id, row.amplitude)
                previous = counts[key]
                counts[key] = (
                    trials=previous.trials + row.trials,
                    detections=previous.detections + row.detections,
                    localized=previous.localized + row.localized,
                )
                append!(errors[key], row.errors)
                push!(context_rows, (
                    context_index=context.context_index,
                    architecture_id=architecture.id, amplitude=row.amplitude,
                    pd=row.pd, pdl=row.pdl,
                    median_localization_error=row.median_localization_error,
                    q90_localization_error=row.q90_localization_error,
                ))
            end
        end
        context.context_index % 10 == 0 && println(
            "[+] O7 final-test contexts $(context.context_index)/$(length(contexts))",
        )
    end
    aggregate_rows = NamedTuple[]
    for architecture in architectures, amplitude in amplitudes
        key = (architecture.id, amplitude)
        value = counts[key]
        error = errors[key]
        push!(aggregate_rows, (
            architecture_id=architecture.id, nuisance_level=:C3, amplitude,
            trials=value.trials, detections=value.detections,
            localized=value.localized, pd=value.detections / value.trials,
            pdl=value.localized / value.trials,
            median_localization_error=isempty(error) ? NaN : median(error),
            q90_localization_error=isempty(error) ? NaN : quantile(error, 0.90),
        ))
    end
    aggregate_rows, context_rows
end

function o7_physical_checks(architecture, parameters, config)
    reflection = o7_reflection_matrix(config.element_count)
    maximum_sigma = maximum(
        opnorm(ports) for state in architecture.frequency_states for ports in state
    )
    reflection_error = maximum(
        norm(reflection * (ports * ports') - (ports * ports') * reflection) /
            max(norm(ports * ports'), eps(Float64))
        for state in architecture.frequency_states for ports in state
    )
    reciprocity_error = maximum(
        norm(transpose(ports) * ports - transpose(transpose(ports) * ports)) /
            max(norm(transpose(ports) * ports), eps(Float64))
        for state in architecture.frequency_states for ports in state
    )
    center = cld(length(config.frequency_ratios), 2)
    projectors = [o7_projector(state[center]) for state in architecture.frequency_states]
    minimum_state_distance = minimum(
        norm(projectors[first] - projectors[second]) / sqrt(6)
        for first in 1:4 for second in (first + 1):4
    )
    zero_states, _, _ = o7_transfer_states(
        architecture.mechanism, parameters, config; zero_control=true,
    )
    zero_control_collapse = maximum(
        norm(zero_states[index][frequency] - zero_states[1][frequency])
        for index in 2:4 for frequency in eachindex(config.frequency_ratios)
    )
    (
        maximum_sigma, reflection_error, reciprocity_error,
        minimum_state_distance, zero_control_collapse,
        nonnegative_loss=parameters.host_damping >= 0 && parameters.quality_a > 0 &&
            parameters.quality_b > 0 && parameters.loss_a >= 0 && parameters.loss_b >= 0,
        passed=maximum_sigma <= 1 + 1.0e-10 && reflection_error <= 1.0e-10 &&
            reciprocity_error <= 1.0e-10 && minimum_state_distance > 1.0e-4 &&
            zero_control_collapse <= 1.0e-12,
    )
end

function o7_summary(
    architectures, threshold_rows, aggregate_rows, reference_alpha,
    pzt_comparator, mechanical_comparator,
)
    thresholds = Dict(row.architecture_id => row for row in threshold_rows)
    alpha = Dict{Symbol,Float64}()
    base_rows = NamedTuple[]
    for architecture in architectures
        selected = filter(row -> row.architecture_id == architecture.id, aggregate_rows)
        alpha[architecture.id] = sf_alpha90(selected)
        reference = only(filter(row -> isapprox(
            row.amplitude, reference_alpha; rtol=1.0e-12,
        ), selected))
        intrinsic = thresholds[architecture.id].pfa_upper_95 <= 0.012 &&
            reference.pd >= 0.90 && reference.pdl >= 0.85 &&
            isfinite(alpha[architecture.id])
        push!(base_rows, (
            architecture_id=architecture.id, role=architecture.role,
            family=architecture.family,
            electronic_channels=architecture.electronic_channels,
            physical_states=architecture.physical_states,
            transmit_events=architecture.transmit_events,
            pfa_upper_95=thresholds[architecture.id].pfa_upper_95,
            alpha90=alpha[architecture.id], pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            q90_localization_error_reference=reference.q90_localization_error,
            intrinsic_gate_passed=intrinsic,
        ))
    end
    [begin
        comparator = row.architecture_id == :pzt_k3q4_frozen ? pzt_comparator :
            row.architecture_id == :mechanical_k3q4_frozen ? mechanical_comparator : nothing
        switching_gain = isnothing(comparator) ? NaN :
            1 - row.alpha90 / alpha[comparator]
        final_pass = row.intrinsic_gate_passed &&
            (isnothing(comparator) || switching_gain >= 0.05)
        merge(row, (; switching_gain_vs_matched_q1=switching_gain, final_gate_passed=final_pass))
    end for row in base_rows]
end

function o7_component_rows(parameters, mechanism)
    rows = NamedTuple[]
    for zone in (:a, :b)
        resonance_ratio = getproperty(parameters, Symbol("resonance_$(zone)"))
        quality = getproperty(parameters, Symbol("quality_$(zone)"))
        frequency_hz = resonance_ratio * O7_CARRIER_HZ
        inductance_h = 1 / ((2pi * frequency_hz)^2 * O7_PAIR_CAPACITANCE_F)
        resistance_ohm = 2pi * frequency_hz * inductance_h / quality
        push!(rows, (
            mechanism, zone, resonance_ratio, frequency_hz,
            assumed_pair_capacitance_nf=O7_PAIR_CAPACITANCE_F * 1.0e9,
            equivalent_inductance_uh=mechanism == :pzt ? inductance_h * 1.0e6 : NaN,
            equivalent_series_resistance_ohm=mechanism == :pzt ? resistance_ohm : NaN,
            quality, position=getproperty(parameters, Symbol("position_$(zone)")),
            width=getproperty(parameters, Symbol("width_$(zone)")),
            feedback_strength=getproperty(parameters, Symbol("strength_$(zone)")),
            loss=getproperty(parameters, Symbol("loss_$(zone)")),
        ))
    end
    rows
end

function o7_write_readme(summary_rows, pzt_check, mechanical_check, component_rows)
    pzt = only(filter(row -> row.architecture_id == :pzt_k3q4_frozen, summary_rows))
    mechanical = only(filter(
        row -> row.architecture_id == :mechanical_k3q4_frozen, summary_rows,
    ))
    open(joinpath(O7_OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O7 stationary global-state gate")
        println(io)
        println(io, "Evidence level: coupled-modal surrogate; not coupled piezo FEM.")
        println(io)
        println(io, "| Frozen candidate | alpha90 | Pfa95 | Pd@ref | Pdl@ref | Q4 gain over matched Q1 | Final gate |")
        println(io, "|---|---:|---:|---:|---:|---:|---|")
        for row in (pzt, mechanical)
            println(io, "| $(row.architecture_id) | $(round(row.alpha90; digits=3)) | " *
                "$(round(100row.pfa_upper_95; digits=3))% | $(round(row.pd_reference; digits=3)) | " *
                "$(round(row.pdl_reference; digits=3)) | " *
                "$(round(100row.switching_gain_vs_matched_q1; digits=2))% | $(row.final_gate_passed) |")
        end
        println(io)
        println(io, "Positive controls: PZT `sigma=$(pzt_check.maximum_sigma)`, reflection error " *
            "`$(pzt_check.reflection_error)`, minimum state distance `$(pzt_check.minimum_state_distance)`; " *
            "mechanical `sigma=$(mechanical_check.maximum_sigma)`, reflection error " *
            "`$(mechanical_check.reflection_error)`, minimum state distance `$(mechanical_check.minimum_state_distance)`." )
        println(io)
        println(io, "The RLC values in `component_translation.csv` assume two 2 nF mirrored patches in parallel per bit. " *
            "They translate the reduced resonance only; they are not a matched BVD terminal design.")
        println(io)
        if pzt.final_gate_passed || mechanical.final_gate_passed
            println(io, "Verdict: at least one stationary two-bit mechanism passes the reduced O7 gate. " *
                "The next authorized calculation is a frozen coupled plane-strain physical model.")
        else
            println(io, "Verdict: neither stationary two-bit mechanism passes. Stop before FEM geometry synthesis.")
        end
    end
end

function o7_plot(summary_rows)
    rows = filter(row -> row.role != :oracle, summary_rows)
    colors = [row.final_gate_passed ? :seagreen : :firebrick for row in rows]
    figure = scatter(
        getproperty.(rows, :transmit_events), getproperty.(rows, :alpha90);
        markercolor=colors, markersize=9, xlabel="Tx/state events",
        ylabel="alpha90 (lower is better)", title="O7 stationary global-state gate",
        label=false, size=(980, 650),
    )
    for row in rows
        annotate!(figure, row.transmit_events, row.alpha90, text(String(row.architecture_id), 7, :left))
    end
    savefig(figure, joinpath(O7_OUTPUT_ROOT, "stationary_state_gate.png"))
end

function o7_plot_state_profiles(pzt, mechanical, config)
    x = element_positions(config)
    center = cld(length(config.frequency_ratios), 2)
    panels = Any[]
    for (architecture, label) in ((pzt, "PZT"), (mechanical, "mechanical"))
        for (state_index, bits) in enumerate(O7_BITS)
            ports = architecture.frequency_states[state_index][center]
            panel = plot(
                xlabel=label == "mechanical" ? "aperture / lambdaP" : "",
                ylabel=state_index == 1 ? "normalized real field" : "",
                title="$label q$(bits[1])$(bits[2])",
                legend=label == "PZT" && state_index == 1 ? :bottomleft : false,
                ylim=(-1.1, 1.1), titlefontsize=11, guidefontsize=9,
            )
            for port_index in 1:size(ports, 2)
                vector = copy(ports[:, port_index])
                anchor = argmax(abs.(vector))
                vector .*= exp(-im * angle(vector[anchor]))
                vector ./= maximum(abs, vector)
                port_label = port_index == 1 ? "even 1" :
                    port_index == 2 ? "even 2" : "odd"
                plot!(panel, x, real.(vector); linewidth=2, label=port_label)
            end
            push!(panels, panel)
        end
    end
    figure = plot(panels...; layout=(2, 4), size=(1400, 760), margin=4Plots.mm)
    savefig(figure, joinpath(O7_OUTPUT_ROOT, "frozen_state_aperture_profiles.png"))
end

function main_o7_global_states()
    mkpath(O7_OUTPUT_ROOT)
    isfile(O7_STATIC_BUNDLE) || error("missing frozen static bundle: $O7_STATIC_BUNDLE")
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    target_states = JLD2.load(O7_STATIC_BUNDLE, "architecture_states")[O7_TARGET_ID]
    nominal_library = response_library(hypotheses, config)
    nuisance_libraries, nuisance_scales, scale_rows = o7_nuisance_contract(config)
    validation_contexts = o7_signal_contexts(
        config, hypotheses, :validation, O7_VALIDATION_CONTEXTS,
    )
    test_contexts = o7_signal_contexts(config, hypotheses, :test, O7_TEST_CONTEXTS)
    println("[+] O7 fresh splits assembled; final test remains unopened")

    pzt, pzt_parameters, pzt_winner, pzt_screen, pzt_operator, pzt_proxy,
        pzt_scale, pzt_raw_sigma = o7_freeze_branch(
            :pzt, target_states, config, nominal_library, hypotheses,
            nuisance_libraries, nuisance_scales, validation_contexts,
        )
    mechanical, mechanical_parameters, mechanical_winner, mechanical_screen,
        mechanical_operator, mechanical_proxy, mechanical_scale,
        mechanical_raw_sigma = o7_freeze_branch(
            :mechanical, target_states, config, nominal_library, hypotheses,
            nuisance_libraries, nuisance_scales, validation_contexts,
        )
    println("[+] O7 froze one PZT and one mechanical candidate before test")

    pzt_q1 = o7_architecture(
        :pzt_k3q1_matched, :frequency_only,
        [pzt.frequency_states[1]], :matched_comparator; mechanism=:pzt,
    )
    mechanical_q1 = o7_architecture(
        :mechanical_k3q1_matched, :frequency_only,
        [mechanical.frequency_states[1]], :matched_comparator; mechanism=:mechanical,
    )
    ideal = sf_architecture(
        :ideal_k3q4_frozen, :static, 3, 4, :reflection_balanced,
        target_states, :operator_target,
    )
    a_scan = sf_architecture(:a_scan_fresh_o7, :a_scan, 1, 1, :uniform, (), :baseline)
    fmc = sf_architecture(:fmc_fresh_o7, :fmc, 17, 1, :full_aperture, (), :oracle)
    architectures = [a_scan, pzt_q1, mechanical_q1, pzt, mechanical, ideal, fmc]

    pzt_check = o7_physical_checks(pzt, pzt_parameters, config)
    mechanical_check = o7_physical_checks(mechanical, mechanical_parameters, config)
    pzt_check.passed || error("PZT positive controls failed: $pzt_check")
    mechanical_check.passed || error("mechanical positive controls failed: $mechanical_check")

    dictionaries, models, threshold_rows = o7_build_models(
        architectures, config, hypotheses, nuisance_libraries, nuisance_scales,
    )
    amplitudes, oracle_validation_pd, reference_alpha, amplitude_anchor =
        o7_reference_amplitudes(
            fmc, models[fmc.id], dictionaries[fmc.id], validation_contexts, hypotheses,
        )
    println("[+] O7 opened frozen final test; FMC reference alpha=$(reference_alpha)")
    aggregate_rows, context_rows = o7_test_campaign(
        architectures, models, test_contexts, hypotheses, amplitudes,
    )
    summary_rows = o7_summary(
        architectures, threshold_rows, aggregate_rows, reference_alpha,
        pzt_q1.id, mechanical_q1.id,
    )

    component_rows = vcat(
        o7_component_rows(pzt_parameters, :pzt),
        o7_component_rows(mechanical_parameters, :mechanical),
    )
    seed_rows = [(purpose=key, seed=string(value)) for (key, value) in pairs(O7_SEEDS)]
    contract_rows = [
        (setting=:evidence_level, value="coupled_modal_surrogate"),
        (setting=:screen_samples_per_branch, value=string(O7_SCREEN_SAMPLES)),
        (setting=:operator_finalists_per_branch, value=string(O7_OPERATOR_FINALISTS)),
        (setting=:validation_h0_trials, value=string(O7_VALIDATION_H0_TRIALS)),
        (setting=:test_h0_trials, value=string(O7_TEST_H0_TRIALS)),
        (setting=:validation_contexts, value=string(O7_VALIDATION_CONTEXTS)),
        (setting=:test_contexts, value=string(O7_TEST_CONTEXTS)),
        (setting=:state_energy, value="1/sqrt(Q)"),
        (setting=:port_energy, value="1/sqrt(K)"),
        (setting=:motion, value="forbidden"),
        (setting=:pzt_resources, value="4 patches; 2 symmetric pairs; 2 passive RLC branches; 2 switches"),
    ]
    physical_rows = [
        merge((architecture_id=:pzt_k3q4_frozen,), pzt_check),
        merge((architecture_id=:mechanical_k3q4_frozen,), mechanical_check),
    ]

    write_csv(joinpath(O7_OUTPUT_ROOT, "seed_contract.csv"), propertynames(first(seed_rows)), seed_rows)
    write_csv(joinpath(O7_OUTPUT_ROOT, "campaign_contract.csv"), propertynames(first(contract_rows)), contract_rows)
    write_csv(joinpath(O7_OUTPUT_ROOT, "nuisance_scales.csv"), propertynames(first(scale_rows)), scale_rows)
    write_csv(joinpath(O7_OUTPUT_ROOT, "pzt_screen.csv"), propertynames(first(pzt_screen)), pzt_screen)
    write_csv(joinpath(O7_OUTPUT_ROOT, "mechanical_screen.csv"), propertynames(first(mechanical_screen)), mechanical_screen)
    write_csv(joinpath(O7_OUTPUT_ROOT, "pzt_operator_finalists.csv"), propertynames(first(pzt_operator)), pzt_operator)
    write_csv(joinpath(O7_OUTPUT_ROOT, "mechanical_operator_finalists.csv"), propertynames(first(mechanical_operator)), mechanical_operator)
    write_csv(joinpath(O7_OUTPUT_ROOT, "pzt_validation_proxy.csv"), propertynames(first(pzt_proxy)), pzt_proxy)
    write_csv(joinpath(O7_OUTPUT_ROOT, "mechanical_validation_proxy.csv"), propertynames(first(mechanical_proxy)), mechanical_proxy)
    write_csv(joinpath(O7_OUTPUT_ROOT, "thresholds.csv"), propertynames(first(threshold_rows)), threshold_rows)
    write_csv(joinpath(O7_OUTPUT_ROOT, "detection_curves.csv"), propertynames(first(aggregate_rows)), aggregate_rows)
    write_csv(joinpath(O7_OUTPUT_ROOT, "test_context_curves.csv"), propertynames(first(context_rows)), context_rows)
    write_csv(joinpath(O7_OUTPUT_ROOT, "frontier_summary.csv"), propertynames(first(summary_rows)), summary_rows)
    write_csv(joinpath(O7_OUTPUT_ROOT, "physical_checks.csv"), propertynames(first(physical_rows)), physical_rows)
    write_csv(joinpath(O7_OUTPUT_ROOT, "component_translation.csv"), propertynames(first(component_rows)), component_rows)
    jldsave(
        joinpath(O7_OUTPUT_ROOT, "o7_global_states.jld2"); config, hypotheses,
        base_seed=O7_BASE_SEED, seed_contract=O7_SEEDS, pzt_parameters,
        mechanical_parameters, pzt_winner, mechanical_winner, pzt_scale,
        pzt_raw_sigma, mechanical_scale, mechanical_raw_sigma,
        pzt_frequency_states=pzt.frequency_states,
        mechanical_frequency_states=mechanical.frequency_states,
        physical_rows, threshold_rows, amplitudes, amplitude_anchor,
        oracle_validation_pd, reference_alpha, aggregate_rows, summary_rows,
    )
    o7_write_readme(summary_rows, pzt_check, mechanical_check, component_rows)
    o7_plot(summary_rows)
    o7_plot_state_profiles(pzt, mechanical, config)
    for row in summary_rows
        println(
            "[+] O7 $(row.architecture_id): alpha90=$(round(row.alpha90; digits=2)), " *
            "Pfa95=$(round(row.pfa_upper_95; digits=4)), pass=$(row.final_gate_passed)",
        )
    end
    println("[+] wrote O7 stationary-state gate to $O7_OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_o7_global_states()
end
