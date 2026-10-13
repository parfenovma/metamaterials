using LinearAlgebra
using Random
using Statistics
using JLD2
using Plots

include(joinpath(@__DIR__, "run_operator_causal_robustness.jl"))
include(joinpath(@__DIR__, "operator_o3_statistics.jl"))

using .OperatorO3Statistics

const SF_OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "static_frontier_v0",
)
const SF_FROZEN_CODEBOOK_PATH = joinpath(
    OUTPUT_ROOT, "o2c_multistate_codebooks.jld2",
)

# Shared fresh-split contract with the parallel motion-allowed campaign.
const SF_BASE_SEED = UInt(0x04f0a117)
const SF_SEEDS = (
    nuisance_train=SF_BASE_SEED + UInt(0x101),
    nuisance_validation=SF_BASE_SEED + UInt(0x202),
    nuisance_test=SF_BASE_SEED + UInt(0x303),
    signal_validation=SF_BASE_SEED + UInt(0x404),
    signal_test=SF_BASE_SEED + UInt(0x505),
    h0_validation=SF_BASE_SEED + UInt(0x606),
    h0_test=SF_BASE_SEED + UInt(0x707),
    operator_search=SF_BASE_SEED + UInt(0x808),
    proxy_h0=SF_BASE_SEED + UInt(0x909),
)
const SF_VALIDATION_H0_TRIALS = 10_000
const SF_TEST_H0_TRIALS = 10_000
const SF_VALIDATION_CONTEXTS = 24
const SF_TEST_CONTEXTS = 60
const SF_SEARCH_SAMPLES = 36
const SF_PROXY_FINALISTS = 3
const SF_PROXY_CONTEXTS = 4
const SF_AMPLITUDE_MULTIPLIERS = [
    0.4, 0.6, 0.8, 1.0, 1.3, 1.7, 2.2, 3.0, 4.0, 6.0, 8.0, 12.0,
    16.0, 24.0, 32.0,
]
const SF_CATEGORY_POWER = Dict(
    :backwall => 2.0,
    :ringing => 2.0,
    :material => 3.0,
)
const SF_FAMILIES = (
    :reflection_balanced,
    :reflection_even,
    :reflection_odd,
    :unconstrained,
)

sf_complex_noise(rng, rows, columns) =
    (randn(rng, rows, columns) .+ im .* randn(rng, rows, columns)) ./ sqrt(2)

function sf_orthonormal_columns(raw, count)
    count == 0 && return zeros(eltype(raw), size(raw, 1), 0)
    Matrix(qr(raw).Q)[:, 1:count]
end

function sf_subspace_ports(rng, basis, count)
    count == 0 && return zeros(Float64, size(basis, 1), 0)
    coefficients = randn(rng, size(basis, 2), count)
    basis * sf_orthonormal_columns(coefficients, count)
end

"""A real, reciprocal, contractive K-port state in a prescribed C_s sector."""
function sf_random_state(rng, config, electronic_channels, family)
    1 <= electronic_channels <= 4 || throw(ArgumentError("K must be in 1:4"))
    full_basis = cosine_ports(config.element_count, config)
    even_basis = full_basis[:, 1:2:end]
    odd_basis = full_basis[:, 2:2:end]
    ports = if family == :reflection_balanced
        even_count = cld(electronic_channels, 2)
        odd_count = electronic_channels - even_count
        hcat(
            sf_subspace_ports(rng, even_basis, even_count),
            sf_subspace_ports(rng, odd_basis, odd_count),
        )
    elseif family == :reflection_even
        sf_subspace_ports(rng, even_basis, electronic_channels)
    elseif family == :reflection_odd
        sf_subspace_ports(rng, odd_basis, electronic_channels)
    elseif family == :unconstrained
        sf_orthonormal_columns(
            randn(rng, config.element_count, electronic_channels),
            electronic_channels,
        )
    else
        throw(ArgumentError("unknown static family: $family"))
    end
    isapprox(ports' * ports, Matrix{Float64}(I, electronic_channels, electronic_channels);
        atol=1.0e-10) || error("state ports lost orthonormality")
    ports
end

function sf_architecture(id, kind, K, Q, family, states, role)
    block_size = kind == :fmc ? 17^2 : kind == :a_scan ? 1 : K^2
    state_blocks = kind == :static ? Q : 1
    (
        id=Symbol(id), kind, electronic_channels=K, physical_states=Q,
        transmit_events=kind == :static ? K * Q : kind == :fmc ? 17 : 1,
        block_size, state_blocks, family, states, role,
    )
end

function sf_pulse_features(raw, config, block_size)
    features, _, _ = apply_pulse_and_select(
        raw,
        config.frequency_ratios,
        collect(eachindex(config.frequency_ratios)),
        block_size,
    )
    features
end

function sf_electronics_whiten(features, block_size, state_blocks, frequencies)
    frequency_count = length(frequencies)
    rows_per_state = frequency_count * block_size
    size(features, 1) == rows_per_state * state_blocks ||
        throw(DimensionMismatch("feature rows do not match state metadata"))
    covariance = frequency_covariance(frequencies, 1; model=:correlated)
    lower = cholesky(Hermitian(covariance)).L
    result = similar(complex.(features))
    for state in 1:state_blocks
        rows = ((state - 1) * rows_per_state + 1):(state * rows_per_state)
        source = reshape(
            view(features, rows, :), block_size, frequency_count, size(features, 2),
        )
        destination = reshape(
            view(result, rows, :), block_size, frequency_count, size(features, 2),
        )
        for channel in 1:block_size
            destination[channel, :, :] .= lower \ source[channel, :, :]
        end
    end
    result
end

function sf_compress_library(library, architecture, config; state_errors=nothing)
    if architecture.kind == :a_scan
        port = ComplexF64.(uniform_port(config))
        raw = one_port_dictionary(library, _ -> port, config)
        features = sf_pulse_features(raw, config, 1)
    elseif architecture.kind == :fmc
        raw = fmc_dictionary(library, config)
        features = sf_pulse_features(raw, config, config.element_count^2)
    elseif architecture.kind == :static
        errors = isnothing(state_errors) ? [(
            amplitude=ones(config.element_count, architecture.electronic_channels),
            delay=zeros(config.element_count, architecture.electronic_channels),
        ) for _ in 1:architecture.physical_states] : state_errors
        length(errors) == architecture.physical_states ||
            throw(DimensionMismatch("one error realization per state is required"))
        state_features = Matrix{ComplexF64}[]
        for state_index in 1:architecture.physical_states
            raw, sigma_rows = causal_multiport_dictionary(
                library,
                architecture.states[state_index],
                config;
                model=:constant_contractive,
                amplitude_errors=errors[state_index].amplitude,
                delay_errors=errors[state_index].delay,
            )
            maximum(getproperty.(sigma_rows, :sigma_after)) <= 1 + 1.0e-10 ||
                error("passivity contraction failed")
            push!(state_features, sf_pulse_features(
                raw, config, architecture.electronic_channels^2,
            ))
        end
        # Fixed total incident energy, independent of the number of states.
        features = reduce(vcat, state_features) ./ sqrt(architecture.physical_states)
    else
        throw(ArgumentError("unknown architecture kind: $(architecture.kind)"))
    end
    sf_electronics_whiten(
        features, architecture.block_size, architecture.state_blocks,
        config.frequency_ratios,
    )
end

function sf_material_hypotheses(seed; count=48)
    rng = Xoshiro(seed)
    hypotheses = DefectHypothesis[]
    for index in 1:count
        depth = 1.0 + 6.2 * rand(rng)
        lateral = -4.0 + 8.0 * rand(rng)
        angle = 2pi * rand(rng)
        secondary = -0.35 + rand(rng)
        rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)]
        tensor = rotation * Diagonal([1.0, secondary]) * transpose(rotation)
        push!(hypotheses, DefectHypothesis(
            Symbol("sf_material_$(seed)_$(index)"), depth, lateral,
            Matrix{Float64}(tensor),
        ))
    end
    hypotheses
end

function sf_backwall_hypotheses(seed; count=17)
    rng = Xoshiro(seed)
    laterals = range(-4.25, 4.25; length=count)
    [
        DefectHypothesis(
            Symbol("sf_backwall_$(seed)_$(index)"),
            7.25 + 0.08 * randn(rng),
            lateral + 0.08 * randn(rng),
            [1.0 0.0; 0.0 0.08 + 0.04 * rand(rng)],
        )
        for (index, lateral) in enumerate(laterals)
    ]
end

function sf_ringing_library(config, seed; atom_count=6)
    rng = Xoshiro(seed)
    count = config.element_count
    positions = element_positions(config)
    normalized_positions = positions ./ maximum(abs, positions)
    spatial_atoms = Matrix{ComplexF64}[
        Matrix{ComplexF64}(I, count, count),
        ones(ComplexF64, count, count),
        ComplexF64.([exp(-abs(first - second) / 1.5)
            for first in positions, second in positions]),
    ]
    for mode in 1:(atom_count - length(spatial_atoms))
        vector = @. cos(pi * mode * (normalized_positions + 1) / 2) +
            0.25 * sin(pi * (mode + 1) * normalized_positions)
        push!(spatial_atoms, ComplexF64.(vector * transpose(vector)))
    end
    spatial_atoms ./= norm.(spatial_atoms)
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

function sf_nuisance_libraries(config, seed)
    Dict(
        :material => response_library(sf_material_hypotheses(seed), config),
        :backwall => response_library(sf_backwall_hypotheses(seed), config),
        :ringing => sf_ringing_library(config, seed),
    )
end

function sf_nuisance_contract(config)
    libraries = Dict(
        :train => sf_nuisance_libraries(config, SF_SEEDS.nuisance_train),
        :validation => sf_nuisance_libraries(config, SF_SEEDS.nuisance_validation),
        :test => sf_nuisance_libraries(config, SF_SEEDS.nuisance_test),
    )
    full = sf_architecture(:scaling_fmc, :fmc, 17, 1, :full_aperture, (), :scaling)
    scales = Dict{Symbol,Float64}()
    scale_rows = NamedTuple[]
    for category in (:backwall, :ringing, :material)
        features = sf_compress_library(libraries[:train][category], full, config)
        target = SF_CATEGORY_POWER[category]
        scale = sqrt(target * size(features, 1) / sum(abs2, features))
        scales[category] = scale
        push!(scale_rows, (
            category,
            target_fmc_power_ratio=target,
            physical_scale=scale,
            achieved_fmc_power_ratio=sum(abs2, scale .* features) / size(features, 1),
        ))
    end
    libraries, scales, scale_rows
end

function sf_clutter_features(architecture, split_libraries, scales, config)
    hcat((
        scales[category] .* sf_compress_library(
            split_libraries[category], architecture, config,
        )
        for category in (:backwall, :ringing, :material)
    )...)
end

function sf_nonlocal_mask(hypotheses; minimum_separation=1.0)
    hypothesis_separation_mask(hypotheses; minimum_separation)
end

function sf_cheap_summary(dictionary, mask)
    norms = vec(sqrt.(sum(abs2, dictionary; dims=1)))
    normalized = dictionary ./ reshape(max.(norms, eps(Float64)), 1, :)
    gram = abs.(normalized' * normalized)
    gram_energy = sum(abs2, gram)
    effective_rank = length(norms)^2 / gram_energy
    (
        minimum_detectability=minimum(norms),
        nonlocal_coherence=maximum(gram[mask]),
        effective_rank,
        nominal_score=minimum(norms) * sqrt(effective_rank) /
            (1 + 2 * maximum(gram[mask])),
    )
end

function sf_signal_contexts(config, hypotheses, split, count)
    seed = split == :validation ? SF_SEEDS.signal_validation : SF_SEEDS.signal_test
    rng = Xoshiro(seed)
    tensors = begin
        rotation = [cos(pi / 4) -sin(pi / 4); sin(pi / 4) cos(pi / 4)]
        [
            Matrix{Float64}(I, 2, 2),
            [1.0 0.0; 0.0 0.15],
            Matrix{Float64}(rotation * Diagonal([1.0, 0.15]) * transpose(rotation)),
        ]
    end
    contexts = NamedTuple[]
    for context_index in 1:count
        perturbed = perturbed_config(rng, config)
        lateral_shift = 0.10 * (2 * rand(rng) - 1)
        true_hypotheses = DefectHypothesis[]
        for (index, hypothesis) in enumerate(hypotheses)
            tensor = tensors[mod1(index + context_index - 1, length(tensors))]
            push!(true_hypotheses, DefectHypothesis(
                hypothesis.id,
                hypothesis.x_wavelengths,
                hypothesis.y_wavelengths - lateral_shift,
                copy(tensor),
            ))
        end
        amplitude_errors = [
            exp.(0.03 .* randn(rng, config.element_count, 4)) for _ in 1:4
        ]
        delay_errors = [
            0.02 .* randn(rng, config.element_count, 4) for _ in 1:4
        ]
        push!(contexts, (
            context_index,
            config=perturbed,
            hypotheses=true_hypotheses,
            library=response_library(true_hypotheses, perturbed),
            amplitude_errors,
            delay_errors,
            lateral_shift,
        ))
    end
    contexts
end

function sf_context_errors(architecture, context)
    architecture.kind == :static || return nothing
    K = architecture.electronic_channels
    [(
        amplitude=context.amplitude_errors[state][:, 1:K],
        delay=context.delay_errors[state][:, 1:K],
    ) for state in 1:architecture.physical_states]
end

function sf_proxy_threshold(precision_dictionary_value, validation_clutter, normalization)
    covariance = matched_output_covariance(
        precision_dictionary_value, validation_clutter,
    )
    factor = OperatorO3Statistics.covariance_factor(covariance)
    rng = Xoshiro(SF_SEEDS.proxy_h0)
    outputs = factor * sf_complex_noise(rng, size(factor, 1), 2000)
    scores, _ = normalized_scores(outputs, normalization)
    quantile(scores, 0.99), covariance
end

function sf_validation_proxy(
    architecture,
    nominal_dictionary,
    train_clutter,
    validation_clutter,
    proxy_contexts,
    config,
    mask,
)
    precision_applied = precision_dictionary(nominal_dictionary, train_clutter)
    gram = matched_mean(precision_applied, nominal_dictionary)
    normalization = max.(real.(diag(gram)), eps(Float64))
    threshold, validation_covariance = sf_proxy_threshold(
        precision_applied, validation_clutter, normalization,
    )
    requirements = Float64[]
    fidelities = Float64[]
    for context in proxy_contexts
        dictionary = sf_compress_library(
            context.library,
            architecture,
            context.config;
            state_errors=sf_context_errors(architecture, context),
        )
        mean_matrix = matched_mean(precision_applied, dictionary)
        append!(requirements, [
            sqrt(threshold * normalization[index]) /
                max(abs(mean_matrix[index, index]), eps(Float64))
            for index in eachindex(normalization)
        ])
        push!(fidelities, calibration_fidelity(
            nominal_dictionary, dictionary,
            Matrix{Float64}(I, size(dictionary, 1), size(dictionary, 1)),
        ).q10)
    end
    summary = sf_cheap_summary(precision_applied, mask)
    proxy_alpha90 = quantile(requirements, 0.90)
    q10_fidelity = quantile(fidelities, 0.10)
    validation_objective = proxy_alpha90 *
        (1 + max(0.0, summary.nonlocal_coherence - 0.60)) /
        max(q10_fidelity, 0.50)
    (
        proxy_alpha90,
        q10_fidelity,
        residual_minimum_detectability=summary.minimum_detectability,
        residual_nonlocal_coherence=summary.nonlocal_coherence,
        residual_effective_rank=summary.effective_rank,
        validation_objective,
        validation_output_variance=median(real.(diag(validation_covariance))),
    )
end

function sf_candidate_search(
    config,
    hypotheses,
    nuisance_libraries,
    nuisance_scales,
    proxy_contexts,
)
    rng = Xoshiro(SF_SEEDS.operator_search)
    nominal_library = response_library(hypotheses, config)
    mask = sf_nonlocal_mask(hypotheses)
    stage1_rows = NamedTuple[]
    state_cache = Dict{Tuple{Int,Int,Symbol,Int},Vector{Matrix{Float64}}}()
    stage1_selected = NamedTuple[]

    for K in 1:4, Q in 1:4, family in SF_FAMILIES
        family_rows = NamedTuple[]
        for sample in 1:SF_SEARCH_SAMPLES
            states = [sf_random_state(rng, config, K, family) for _ in 1:Q]
            id = Symbol("search_k$(K)_q$(Q)_$(family)_$(sample)")
            architecture = sf_architecture(id, :static, K, Q, family, states, :search)
            dictionary = sf_compress_library(nominal_library, architecture, config)
            summary = sf_cheap_summary(dictionary, mask)
            row = (
                electronic_channels=K,
                physical_states=Q,
                family,
                sample,
                minimum_detectability=summary.minimum_detectability,
                nonlocal_coherence=summary.nonlocal_coherence,
                effective_rank=summary.effective_rank,
                nominal_score=summary.nominal_score,
            )
            push!(stage1_rows, row)
            push!(family_rows, row)
            state_cache[(K, Q, family, sample)] = states
        end
        sort!(family_rows; by=row -> row.nominal_score, rev=true)
        append!(stage1_selected, first(family_rows, SF_PROXY_FINALISTS))
    end

    proxy_rows = NamedTuple[]
    proxy_architectures = Dict{Tuple{Int,Int,Symbol,Int},NamedTuple}()
    for (index, row) in enumerate(stage1_selected)
        key = (row.electronic_channels, row.physical_states, row.family, row.sample)
        id = Symbol("proxy_k$(key[1])_q$(key[2])_$(key[3])_$(key[4])")
        architecture = sf_architecture(
            id, :static, key[1], key[2], key[3], state_cache[key], :proxy,
        )
        dictionary = sf_compress_library(nominal_library, architecture, config)
        train_clutter = sf_clutter_features(
            architecture, nuisance_libraries[:train], nuisance_scales, config,
        )
        validation_clutter = sf_clutter_features(
            architecture, nuisance_libraries[:validation], nuisance_scales, config,
        )
        proxy = sf_validation_proxy(
            architecture,
            dictionary,
            train_clutter,
            validation_clutter,
            proxy_contexts,
            config,
            mask,
        )
        push!(proxy_rows, merge(row, proxy))
        proxy_architectures[key] = architecture
        index % 24 == 0 && println(
            "[+] static proxy finalists evaluated: $index / $(length(stage1_selected))",
        )
    end

    family_winners = NamedTuple[]
    for K in 1:4, Q in 1:4, family in SF_FAMILIES
        rows = filter(row ->
            row.electronic_channels == K && row.physical_states == Q &&
            row.family == family,
        proxy_rows)
        push!(family_winners, rows[argmin(getproperty.(rows, :validation_objective))])
    end
    cell_winners = NamedTuple[]
    architectures = NamedTuple[]
    for K in 1:4, Q in 1:4
        rows = filter(row ->
            row.electronic_channels == K && row.physical_states == Q,
        family_winners)
        winner = rows[argmin(getproperty.(rows, :validation_objective))]
        key = (K, Q, winner.family, winner.sample)
        source = proxy_architectures[key]
        id = Symbol("trained_k$(K)_q$(Q)_$(winner.family)")
        architecture = merge(source, (id=id, role=:trained_cell_winner))
        push!(cell_winners, merge((architecture_id=id,), winner))
        push!(architectures, architecture)
    end
    stage1_rows, proxy_rows, family_winners, cell_winners, architectures
end

function sf_frozen_architecture(config)
    isfile(SF_FROZEN_CODEBOOK_PATH) ||
        error("frozen O2c codebook is missing: $SF_FROZEN_CODEBOOK_PATH")
    bundle = JLD2.load(SF_FROZEN_CODEBOOK_PATH)
    definition = bundle["definitions"][:real_static_two_state_6]
    sf_architecture(
        :frozen_o2c_k4_q2,
        :static,
        4,
        2,
        :reflection_balanced,
        collect(definition.state_ports),
        :frozen_entrant,
    )
end

function sf_detection_model(
    architecture,
    dictionary,
    train_clutter,
    validation_clutter,
    test_clutter,
    architecture_index,
)
    precision_applied = precision_dictionary(dictionary, train_clutter)
    gram = matched_mean(precision_applied, dictionary)
    normalization = max.(real.(diag(gram)), eps(Float64))
    validation_covariance = matched_output_covariance(
        precision_applied, validation_clutter,
    )
    test_covariance = matched_output_covariance(precision_applied, test_clutter)
    validation_factor = OperatorO3Statistics.covariance_factor(validation_covariance)
    test_factor = OperatorO3Statistics.covariance_factor(test_covariance)
    validation_rng = Xoshiro(SF_SEEDS.h0_validation + UInt(10_000 * architecture_index))
    test_rng = Xoshiro(SF_SEEDS.h0_test + UInt(10_000 * architecture_index))
    validation_outputs = validation_factor * sf_complex_noise(
        validation_rng, size(validation_factor, 1), SF_VALIDATION_H0_TRIALS,
    )
    validation_scores, _ = normalized_scores(validation_outputs, normalization)
    threshold = quantile(validation_scores, 0.99)
    test_outputs = test_factor * sf_complex_noise(
        test_rng, size(test_factor, 1), SF_TEST_H0_TRIALS,
    )
    test_scores, _ = normalized_scores(test_outputs, normalization)
    false_alarms = count(>(threshold), test_scores)
    (
        precision_applied,
        normalization,
        validation_factor,
        test_factor,
        threshold,
        false_alarms,
        observed_pfa=false_alarms / SF_TEST_H0_TRIALS,
        pfa_upper_95=wilson_upper_bound(false_alarms, SF_TEST_H0_TRIALS),
    )
end

function sf_build_models(architectures, config, hypotheses, libraries, scales)
    nominal_library = response_library(hypotheses, config)
    dictionaries = Dict{Symbol,Matrix{ComplexF64}}()
    models = Dict{Symbol,NamedTuple}()
    threshold_rows = NamedTuple[]
    for (index, architecture) in enumerate(architectures)
        dictionary = sf_compress_library(nominal_library, architecture, config)
        dictionaries[architecture.id] = dictionary
        train_clutter = sf_clutter_features(
            architecture, libraries[:train], scales, config,
        )
        validation_clutter = sf_clutter_features(
            architecture, libraries[:validation], scales, config,
        )
        test_clutter = sf_clutter_features(
            architecture, libraries[:test], scales, config,
        )
        model = sf_detection_model(
            architecture,
            dictionary,
            train_clutter,
            validation_clutter,
            test_clutter,
            index,
        )
        models[architecture.id] = model
        push!(threshold_rows, (
            architecture_id=architecture.id,
            validation_trials=SF_VALIDATION_H0_TRIALS,
            test_trials=SF_TEST_H0_TRIALS,
            threshold=model.threshold,
            false_alarms=model.false_alarms,
            observed_pfa=model.observed_pfa,
            pfa_upper_95=model.pfa_upper_95,
        ))
        println(
            "[+] static model $(index)/$(length(architectures)): $(architecture.id), " *
            "Pfa95=$(round(model.pfa_upper_95; digits=4))",
        )
    end
    dictionaries, models, threshold_rows
end

function sf_localization_error(predicted, truth, templates, true_hypotheses)
    hypot(
        templates[predicted].x_wavelengths - true_hypotheses[truth].x_wavelengths,
        templates[predicted].y_wavelengths - true_hypotheses[truth].y_wavelengths,
    )
end

function sf_signal_batch(
    rng,
    factor,
    mean_matrix,
    normalization,
    threshold,
    amplitudes,
    templates,
    true_hypotheses,
)
    hypothesis_count = size(mean_matrix, 2)
    noise = factor * sf_complex_noise(rng, size(factor, 1), hypothesis_count)
    rows = NamedTuple[]
    for amplitude in amplitudes
        scores, predicted = normalized_scores(
            noise .+ amplitude .* mean_matrix, normalization,
        )
        detected = scores .> threshold
        errors = Float64[]
        localized = 0
        for truth in 1:hypothesis_count
            detected[truth] || continue
            error = sf_localization_error(
                predicted[truth], truth, templates, true_hypotheses,
            )
            push!(errors, error)
            localized += error <= 0.5 + 1.0e-12
        end
        push!(rows, (
            amplitude,
            trials=hypothesis_count,
            detections=count(detected),
            localized,
            pd=count(detected) / hypothesis_count,
            pdl=localized / hypothesis_count,
            median_localization_error=isempty(errors) ? NaN : median(errors),
            q90_localization_error=isempty(errors) ? NaN : quantile(errors, 0.90),
            errors,
        ))
    end
    rows
end

function sf_reference_amplitudes(
    oracle,
    model,
    dictionary,
    validation_contexts,
    hypotheses,
)
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
    for context in validation_contexts
        true_dictionary = sf_compress_library(
            context.library, oracle, context.config,
        )
        mean_matrix = matched_mean(model.precision_applied, true_dictionary)
        rng = Xoshiro(SF_SEEDS.signal_validation + UInt(context.context_index))
        batch = sf_signal_batch(
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
    isnothing(reference_index) && error("FMC validation Pd did not reach 0.95")
    amplitudes, probabilities, amplitudes[reference_index], anchor
end

function sf_test_campaign(
    architectures,
    models,
    contexts,
    hypotheses,
    amplitudes,
)
    counts = Dict{Tuple{Symbol,Float64},NamedTuple}()
    error_cache = Dict{Tuple{Symbol,Float64},Vector{Float64}}()
    for architecture in architectures, amplitude in amplitudes
        counts[(architecture.id, amplitude)] = (trials=0, detections=0, localized=0)
        error_cache[(architecture.id, amplitude)] = Float64[]
    end
    context_rows = NamedTuple[]
    for context in contexts
        for (architecture_index, architecture) in enumerate(architectures)
            model = models[architecture.id]
            true_dictionary = sf_compress_library(
                context.library,
                architecture,
                context.config;
                state_errors=sf_context_errors(architecture, context),
            )
            mean_matrix = matched_mean(model.precision_applied, true_dictionary)
            rng = Xoshiro(
                SF_SEEDS.signal_test + UInt(10_000 * context.context_index + architecture_index),
            )
            batch = sf_signal_batch(
                rng,
                model.test_factor,
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
                append!(error_cache[key], row.errors)
                push!(context_rows, (
                    context_index=context.context_index,
                    architecture_id=architecture.id,
                    amplitude=row.amplitude,
                    pd=row.pd,
                    pdl=row.pdl,
                    median_localization_error=row.median_localization_error,
                    q90_localization_error=row.q90_localization_error,
                ))
            end
        end
        context.context_index % 10 == 0 && println(
            "[+] static final-test contexts: $(context.context_index) / $(length(contexts))",
        )
    end
    aggregate_rows = NamedTuple[]
    for architecture in architectures, amplitude in amplitudes
        key = (architecture.id, amplitude)
        value = counts[key]
        errors = error_cache[key]
        push!(aggregate_rows, (
            architecture_id=architecture.id,
            nuisance_level=:C3,
            amplitude,
            trials=value.trials,
            detections=value.detections,
            localized=value.localized,
            pd=value.detections / value.trials,
            pdl=value.localized / value.trials,
            median_localization_error=isempty(errors) ? NaN : median(errors),
            q90_localization_error=isempty(errors) ? NaN : quantile(errors, 0.90),
        ))
    end
    aggregate_rows, context_rows
end

function sf_alpha90(rows)
    selected = sort(rows; by=row -> row.amplitude)
    crossing_amplitude(
        getproperty.(selected, :amplitude), getproperty.(selected, :pd), 0.90,
    )
end

function sf_context_summaries(context_rows, architectures, reference_alpha)
    rows = NamedTuple[]
    for context_index in 1:SF_TEST_CONTEXTS, architecture in architectures
        selected = filter(row ->
            row.context_index == context_index &&
            row.architecture_id == architecture.id,
        context_rows)
        reference = only(filter(row -> isapprox(
            row.amplitude, reference_alpha; rtol=1.0e-12,
        ), selected))
        push!(rows, (
            context_index,
            architecture_id=architecture.id,
            alpha90=sf_alpha90(selected),
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            q90_localization_error_reference=reference.q90_localization_error,
        ))
    end
    rows
end

function sf_resource_metadata(architectures)
    [(
        architecture_id=architecture.id,
        role=architecture.role,
        family=architecture.family,
        electronic_channels=architecture.electronic_channels,
        physical_states=architecture.physical_states,
        transmit_events=architecture.transmit_events,
        state_switches=max(architecture.physical_states - 1, 0),
        calibration_complex_coefficients=architecture.kind == :static ?
            length(SurrogateConfig().frequency_ratios) *
                architecture.electronic_channels^2 * architecture.physical_states :
            architecture.kind == :fmc ?
                length(SurrogateConfig().frequency_ratios) * 17^2 : 7,
        complex_features=architecture.block_size * architecture.state_blocks * 7,
        total_energy_rule=architecture.kind == :static ?
            "1/sqrt(Q) across states and 1/sqrt(K) per Tx basis" :
            "unit total incident energy",
        mechanical_positions=1,
        motion_forbidden=true,
        lift_status=architecture.kind == :static ?
            "C_s group sectors + translation_kz" : "plane_strain_surrogate",
    ) for architecture in architectures]
end

function sf_dominates(first, second)
    costs_first = (
        first.alpha90,
        first.electronic_channels,
        first.physical_states,
        first.calibration_complex_coefficients,
        first.transmit_events,
    )
    costs_second = (
        second.alpha90,
        second.electronic_channels,
        second.physical_states,
        second.calibration_complex_coefficients,
        second.transmit_events,
    )
    all(first_value <= second_value + 1.0e-12
        for (first_value, second_value) in zip(costs_first, costs_second)) &&
        any(first_value < second_value - 1.0e-12
            for (first_value, second_value) in zip(costs_first, costs_second))
end

function sf_frontier_summary(
    architectures,
    metadata,
    thresholds,
    aggregate_rows,
    context_rows,
    reference_alpha,
)
    alpha = Dict(
        architecture.id => sf_alpha90(filter(row ->
            row.architecture_id == architecture.id,
        aggregate_rows)) for architecture in architectures
    )
    a_value = alpha[:a_scan_uniform_fresh]
    fmc_value = alpha[:paut_fmc_oracle_fresh]
    rows = NamedTuple[]
    for architecture in architectures
        meta = only(filter(row -> row.architecture_id == architecture.id, metadata))
        threshold = only(filter(row -> row.architecture_id == architecture.id, thresholds))
        reference = only(filter(row ->
            row.architecture_id == architecture.id &&
            isapprox(row.amplitude, reference_alpha; rtol=1.0e-12),
        aggregate_rows))
        context_values = filter(row -> row.architecture_id == architecture.id, context_rows)
        context_gaps = Float64[]
        for context in 1:SF_TEST_CONTEXTS
            lookup(id) = only(filter(row ->
                row.context_index == context && row.architecture_id == id,
            context_rows)).alpha90
            av = lookup(:a_scan_uniform_fresh)
            cv = lookup(architecture.id)
            fv = lookup(:paut_fmc_oracle_fresh)
            if all(isfinite, (av, cv, fv)) && av > fv + eps(Float64)
                push!(context_gaps, (av - cv) / (av - fv))
            end
        end
        gap = isfinite(a_value) && isfinite(fmc_value) && a_value > fmc_value ?
            (a_value - alpha[architecture.id]) / (a_value - fmc_value) : NaN
        q10_gap = isempty(context_gaps) ? NaN : quantile(context_gaps, 0.10)
        intrinsic_pass = architecture.kind == :static &&
            threshold.pfa_upper_95 <= 0.0125 &&
            reference.pd >= 0.90 &&
            reference.pdl >= 0.80 &&
            reference.q90_localization_error <= 1.0 &&
            length(context_gaps) >= ceil(Int, 0.9 * SF_TEST_CONTEXTS) &&
            q10_gap >= 0.50
        push!(rows, (
            architecture_id=architecture.id,
            role=architecture.role,
            family=architecture.family,
            electronic_channels=architecture.electronic_channels,
            physical_states=architecture.physical_states,
            transmit_events=architecture.transmit_events,
            state_switches=meta.state_switches,
            calibration_complex_coefficients=meta.calibration_complex_coefficients,
            pfa_upper_95=threshold.pfa_upper_95,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            q90_localization_error_reference=reference.q90_localization_error,
            alpha90=alpha[architecture.id],
            aggregate_fmc_gap_closure=gap,
            q10_context_fmc_gap_closure=q10_gap,
            valid_context_gap_count=length(context_gaps),
            intrinsic_gate_passed=intrinsic_pass,
            pareto=false,
        ))
    end
    static_indices = findall(row -> row.role in (:trained_cell_winner, :frozen_entrant), rows)
    for index in static_indices
        row = rows[index]
        feasible = isfinite(row.alpha90) && row.pfa_upper_95 <= 0.0125
        pareto = feasible && !any(
            sf_dominates(rows[other], row)
            for other in static_indices
            if other != index && isfinite(rows[other].alpha90) &&
                rows[other].pfa_upper_95 <= 0.0125
        )
        rows[index] = merge(row, (pareto=pareto,))
    end
    rows
end

function sf_plot_results(summary_rows, family_winners)
    candidates = filter(row -> row.role in (:trained_cell_winner, :frozen_entrant), summary_rows)
    figure = scatter(
        getproperty.(candidates, :transmit_events),
        getproperty.(candidates, :alpha90);
        marker_z=getproperty.(candidates, :electronic_channels),
        markersize=5 .+ 2 .* getproperty.(candidates, :physical_states),
        xlabel="Tx/state events KQ",
        ylabel="alpha90 (lower is better)",
        title="Motion-forbidden static resource frontier",
        colorbar_title="K",
        label="frozen cell winners",
        size=(980, 650),
    )
    pareto = filter(row -> row.pareto, candidates)
    isempty(pareto) || scatter!(
        figure,
        getproperty.(pareto, :transmit_events),
        getproperty.(pareto, :alpha90);
        marker=:star5,
        markersize=14,
        color=:black,
        label="feasible Pareto",
    )
    savefig(figure, joinpath(SF_OUTPUT_ROOT, "static_resource_frontier.png"))

    family_plot = scatter(
        getproperty.(family_winners, :electronic_channels) .*
            getproperty.(family_winners, :physical_states),
        getproperty.(family_winners, :validation_objective);
        group=String.(getproperty.(family_winners, :family)),
        yscale=:log10,
        xlabel="Tx/state events KQ",
        ylabel="validation proxy alpha (lower is better)",
        title="Train+validation symmetry-family comparison",
        markersize=6,
        size=(980, 650),
    )
    savefig(family_plot, joinpath(SF_OUTPUT_ROOT, "symmetry_family_validation.png"))
end

function sf_seed_rows()
    [(purpose=key, seed=string(value)) for (key, value) in pairs(SF_SEEDS)]
end

function sf_campaign_contract_rows()
    [
        (setting=:base_seed, value=string(SF_BASE_SEED)),
        (setting=:validation_h0_trials, value=string(SF_VALIDATION_H0_TRIALS)),
        (setting=:test_h0_trials, value=string(SF_TEST_H0_TRIALS)),
        (setting=:validation_signal_contexts, value=string(SF_VALIDATION_CONTEXTS)),
        (setting=:test_signal_contexts, value=string(SF_TEST_CONTEXTS)),
        (setting=:amplitude_multipliers,
            value=join(string.(SF_AMPLITUDE_MULTIPLIERS), '|')),
        (setting=:reference_rule, value="first FMC validation point with Pd>=0.95"),
        (setting=:background_levels, value="C3=backwall+ringing+material"),
        (setting=:backwall_atoms, value="17"),
        (setting=:ringing_atoms, value="6"),
        (setting=:material_atoms, value="48"),
        (setting=:fmc_power_backwall, value="2.0"),
        (setting=:fmc_power_ringing, value="2.0"),
        (setting=:fmc_power_material, value="3.0"),
        (setting=:search_samples_per_family_cell, value=string(SF_SEARCH_SAMPLES)),
        (setting=:proxy_finalists_per_family_cell, value=string(SF_PROXY_FINALISTS)),
        (setting=:proxy_signal_contexts, value=string(SF_PROXY_CONTEXTS)),
        (setting=:selection_rule,
            value="one KQ winner frozen on train+validation before test"),
        (setting=:energy_rule,
            value="1/sqrt(Q) across states and 1/sqrt(K) per Tx basis"),
        (setting=:mechanical_positions, value="1"),
    ]
end

function sf_self_check(architectures, config)
    for architecture in architectures
        architecture.kind == :static || continue
        length(architecture.states) == architecture.physical_states ||
            error("state-count mismatch")
        for ports in architecture.states
            size(ports) == (config.element_count, architecture.electronic_channels) ||
                error("port shape mismatch")
            isapprox(ports' * ports,
                Matrix{Float64}(I, architecture.electronic_channels,
                    architecture.electronic_channels); atol=1.0e-9) ||
                error("non-orthonormal ports")
            opnorm(ports) <= 1 + 1.0e-9 || error("non-contractive state")
        end
    end
    true
end

function main_static_frontier()
    mkpath(SF_OUTPUT_ROOT)
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    nuisance_libraries, nuisance_scales, scale_rows = sf_nuisance_contract(config)
    validation_contexts = sf_signal_contexts(
        config, hypotheses, :validation, SF_VALIDATION_CONTEXTS,
    )
    test_contexts = sf_signal_contexts(
        config, hypotheses, :test, SF_TEST_CONTEXTS,
    )
    println("[+] fresh static-frontier nuisance and signal splits assembled")

    stage1_rows, proxy_rows, family_winners, cell_winners, trained_architectures =
        sf_candidate_search(
            config,
            hypotheses,
            nuisance_libraries,
            nuisance_scales,
            first(validation_contexts, SF_PROXY_CONTEXTS),
        )
    println("[+] froze 16 train+validation-selected (K,Q) cell winners")

    a_scan = sf_architecture(
        :a_scan_uniform_fresh, :a_scan, 1, 1, :uniform, (), :baseline,
    )
    fmc = sf_architecture(
        :paut_fmc_oracle_fresh, :fmc, 17, 1, :full_aperture, (), :oracle,
    )
    frozen = sf_frozen_architecture(config)
    architectures = vcat([a_scan], trained_architectures, [frozen, fmc])
    sf_self_check(architectures, config)

    dictionaries, models, threshold_rows = sf_build_models(
        architectures,
        config,
        hypotheses,
        nuisance_libraries,
        nuisance_scales,
    )
    amplitudes, oracle_validation_pd, reference_alpha, amplitude_anchor =
        sf_reference_amplitudes(
            fmc,
            models[fmc.id],
            dictionaries[fmc.id],
            validation_contexts,
            hypotheses,
        )
    println(
        "[+] fresh FMC reference alpha=$(round(reference_alpha; digits=4)), " *
        "validation Pd=$(round(oracle_validation_pd[findfirst(==(reference_alpha), amplitudes)]; digits=3))",
    )

    aggregate_rows, context_curve_rows = sf_test_campaign(
        architectures, models, test_contexts, hypotheses, amplitudes,
    )
    context_summary_rows = sf_context_summaries(
        context_curve_rows, architectures, reference_alpha,
    )
    metadata_rows = sf_resource_metadata(architectures)
    summary_rows = sf_frontier_summary(
        architectures,
        metadata_rows,
        threshold_rows,
        aggregate_rows,
        context_summary_rows,
        reference_alpha,
    )

    write_csv(joinpath(SF_OUTPUT_ROOT, "seed_contract.csv"),
        propertynames(first(sf_seed_rows())), sf_seed_rows())
    write_csv(joinpath(SF_OUTPUT_ROOT, "campaign_contract.csv"),
        propertynames(first(sf_campaign_contract_rows())), sf_campaign_contract_rows())
    write_csv(joinpath(SF_OUTPUT_ROOT, "nuisance_scales.csv"),
        propertynames(first(scale_rows)), scale_rows)
    write_csv(joinpath(SF_OUTPUT_ROOT, "search_stage1.csv"),
        propertynames(first(stage1_rows)), stage1_rows)
    write_csv(joinpath(SF_OUTPUT_ROOT, "validation_proxy.csv"),
        propertynames(first(proxy_rows)), proxy_rows)
    write_csv(joinpath(SF_OUTPUT_ROOT, "family_winners.csv"),
        propertynames(first(family_winners)), family_winners)
    write_csv(joinpath(SF_OUTPUT_ROOT, "cell_winners.csv"),
        propertynames(first(cell_winners)), cell_winners)
    write_csv(joinpath(SF_OUTPUT_ROOT, "architecture_metadata.csv"),
        propertynames(first(metadata_rows)), metadata_rows)
    write_csv(joinpath(SF_OUTPUT_ROOT, "thresholds.csv"),
        propertynames(first(threshold_rows)), threshold_rows)
    write_csv(joinpath(SF_OUTPUT_ROOT, "detection_curves.csv"),
        propertynames(first(aggregate_rows)), aggregate_rows)
    write_csv(joinpath(SF_OUTPUT_ROOT, "test_context_curves.csv"),
        propertynames(first(context_curve_rows)), context_curve_rows)
    write_csv(joinpath(SF_OUTPUT_ROOT, "test_context_summary.csv"),
        propertynames(first(context_summary_rows)), context_summary_rows)
    write_csv(joinpath(SF_OUTPUT_ROOT, "frontier_summary.csv"),
        propertynames(first(summary_rows)), summary_rows)
    jldsave(
        joinpath(SF_OUTPUT_ROOT, "static_frontier.jld2");
        config,
        hypotheses,
        base_seed=SF_BASE_SEED,
        seed_contract=SF_SEEDS,
        category_power=SF_CATEGORY_POWER,
        stage1_rows,
        proxy_rows,
        family_winners,
        cell_winners,
        architecture_states=Dict(
            architecture.id => architecture.states
            for architecture in architectures if architecture.kind == :static
        ),
        metadata_rows,
        threshold_rows,
        amplitudes,
        amplitude_anchor,
        oracle_validation_pd,
        reference_alpha,
        aggregate_rows,
        context_summary_rows,
        summary_rows,
    )
    sf_plot_results(summary_rows, family_winners)

    for row in filter(row -> row.pareto, summary_rows)
        println(
            "[+] Pareto $(row.architecture_id): K=$(row.electronic_channels), " *
            "Q=$(row.physical_states), alpha90=$(round(row.alpha90; digits=2)), " *
            "Pfa95=$(round(row.pfa_upper_95; digits=4)), " *
            "gap=$(round(row.aggregate_fmc_gap_closure; digits=3))",
        )
    end
    println("[+] wrote fresh static O3 frontier to $SF_OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_static_frontier()
end
