using LinearAlgebra
using Random
using Statistics
using JLD2
using Plots

include(joinpath(@__DIR__, "optimize_operator_causal_codebooks.jl"))

const MULTISTATE_SEARCH_SEED = 0x51a7e002
const MULTISTATE_SAMPLES_PER_FAMILY = 300
const MULTISTATE_POOL_PER_FAMILY = 16
const MULTISTATE_FINALISTS_PER_TRANSFER = 6

function unconstrained_ports(rng, config; complex_weights)
    coefficient_type = complex_weights ? ComplexF64 : Float64
    raw = randn(rng, coefficient_type, config.element_count, 4)
    causal_orthonormal_columns(raw, 4)
end

function state_spec(family)
    if family == :real_sector
        return (transfer=:constant_contractive, complex_weights=false, sector=true)
    elseif family == :real_unconstrained
        return (transfer=:constant_contractive, complex_weights=false, sector=false)
    elseif family == :complex_sector
        return (transfer=:delay_contractive, complex_weights=true, sector=true)
    elseif family == :complex_unconstrained
        return (transfer=:delay_contractive, complex_weights=true, sector=false)
    end
    throw(ArgumentError("unknown multistate family: $family"))
end

function random_state_ports(rng, config, family)
    spec = state_spec(family)
    spec.sector ?
        causal_sector_ports(rng, config; complex_weights=spec.complex_weights) :
        unconstrained_ports(rng, config; complex_weights=spec.complex_weights)
end

function state_features(ports, family, library, config; amplitude_errors=ones(size(ports)), delay_errors=zeros(size(ports)))
    raw, sigma_rows = causal_multiport_dictionary(
        library,
        ports,
        config;
        model=state_spec(family).transfer,
        amplitude_errors,
        delay_errors,
    )
    features, _, _ = apply_pulse_and_select(
        raw,
        config.frequency_ratios,
        collect(eachindex(config.frequency_ratios)),
        16,
    )
    features, maximum(getproperty.(sigma_rows, :sigma_after))
end

function combined_state_operator(feature_matrices, base_covariance, hypotheses)
    state_count = length(feature_matrices)
    dictionary = reduce(vcat, feature_matrices) ./ sqrt(state_count)
    covariance = kron(Matrix{Float64}(I, state_count, state_count), base_covariance)
    operator_summary(dictionary, hypotheses, covariance), dictionary, covariance
end

function search_multistate_basis(library, hypotheses, config, baseline_minimum)
    rng = Xoshiro(MULTISTATE_SEARCH_SEED)
    covariance = frequency_covariance(config.frequency_ratios, 16; model=:correlated)
    legacy = JLD2.load(CODEBOOK_PATH)["best_ports"]
    legacy_ports = Dict(
        :real_sector => legacy[(4, :real_sector)],
        :complex_sector => legacy[(4, :complex_sector)],
    )
    families = (
        :real_sector,
        :real_unconstrained,
        :complex_sector,
        :complex_unconstrained,
    )
    rows = NamedTuple[]
    ports_by_key = Dict{Tuple{Symbol,Int},Matrix}()
    features_by_key = Dict{Tuple{Symbol,Int},Matrix{ComplexF64}}()

    for family in families
        first_sample = haskey(legacy_ports, family) ? 0 : 1
        last_sample = MULTISTATE_SAMPLES_PER_FAMILY
        for sample in first_sample:last_sample
            ports = sample == 0 ? legacy_ports[family] : random_state_ports(rng, config, family)
            features, maximum_sigma = state_features(ports, family, library, config)
            summary = operator_summary(features, hypotheses, covariance)
            score = causal_gate_score(summary, baseline_minimum)
            push!(rows, (
                family,
                sample,
                minimum_detectability=summary.minimum_detectability,
                nonlocal_coherence=summary.nonlocal_coherence,
                stable_rank=summary.stable_rank,
                maximum_sigma,
                gate_score=score,
                nominal_gate_passed=score >= 1.0,
            ))
            key = (family, sample)
            ports_by_key[key] = ports
            features_by_key[key] = features
        end
        best = argmax(row -> row.gate_score, filter(row -> row.family == family, rows))
        println(
            "[+] one-state $(family): score=$(round(best.gate_score; digits=3)), " *
            "min=$(round(best.minimum_detectability; digits=3)), " *
            "mu=$(round(best.nonlocal_coherence; digits=3)), " *
            "rank=$(round(best.stable_rank; digits=3))",
        )
    end
    rows, ports_by_key, features_by_key, covariance
end

function multistate_pairs(
    search_rows,
    ports_by_key,
    features_by_key,
    covariance,
    hypotheses,
    baseline_minimum,
)
    transfer_groups = Dict(
        :real_static => (:real_sector, :real_unconstrained),
        :complex_delay_contractive => (:complex_sector, :complex_unconstrained),
    )
    pair_rows = NamedTuple[]
    pair_definitions = Dict{Tuple{Symbol,Int},NamedTuple}()

    for (transfer_group, families) in transfer_groups
        pool = NamedTuple[]
        for family in families
            selected = sort(
                filter(row -> row.family == family, search_rows);
                by=row -> row.gate_score,
                rev=true,
            )
            append!(pool, first(selected, MULTISTATE_POOL_PER_FAMILY))
        end
        pair_index = 0
        for first_index in 1:(length(pool) - 1), second_index in (first_index + 1):length(pool)
            first_state = pool[first_index]
            second_state = pool[second_index]
            first_key = (first_state.family, first_state.sample)
            second_key = (second_state.family, second_state.sample)
            summary, _, _ = combined_state_operator(
                [features_by_key[first_key], features_by_key[second_key]],
                covariance,
                hypotheses,
            )
            pair_index += 1
            score = causal_gate_score(summary, baseline_minimum)
            push!(pair_rows, (
                transfer_group,
                pair_index,
                first_family=first_state.family,
                first_sample=first_state.sample,
                second_family=second_state.family,
                second_sample=second_state.sample,
                sector_preserving=first_state.family == second_state.family &&
                    first_state.family in (:real_sector, :complex_sector),
                minimum_detectability=summary.minimum_detectability,
                nonlocal_coherence=summary.nonlocal_coherence,
                stable_rank=summary.stable_rank,
                gate_score=score,
                nominal_gate_passed=score >= 1.0,
            ))
            pair_definitions[(transfer_group, pair_index)] = (
                state_keys=(first_key, second_key),
                state_ports=(ports_by_key[first_key], ports_by_key[second_key]),
            )
        end
        group_rows = filter(row -> row.transfer_group == transfer_group, pair_rows)
        best = group_rows[argmax(getproperty.(group_rows, :gate_score))]
        println(
            "[+] two-state $(transfer_group): score=$(round(best.gate_score; digits=3)), " *
            "min=$(round(best.minimum_detectability; digits=3)), " *
            "mu=$(round(best.nonlocal_coherence; digits=3)), " *
            "rank=$(round(best.stable_rank; digits=3)), " *
            "nominal_pass=$(best.nominal_gate_passed)",
        )
    end
    pair_rows, pair_definitions
end

function select_multistate_finalists(pair_rows, pair_definitions)
    finalists = NamedTuple[]
    definitions = Dict{Symbol,NamedTuple}()
    for transfer_group in (:real_static, :complex_delay_contractive)
        selected = sort(
            filter(row -> row.transfer_group == transfer_group, pair_rows);
            by=row -> row.gate_score,
            rev=true,
        )
        for (rank_index, row) in enumerate(first(selected, MULTISTATE_FINALISTS_PER_TRANSFER))
            id = Symbol(transfer_group, "_two_state_", rank_index)
            push!(finalists, merge((; architecture_id=id), row))
            definitions[id] = pair_definitions[(transfer_group, row.pair_index)]
        end
    end
    finalists, definitions
end

function evaluate_multistate_finalists(
    finalists,
    definitions,
    nominal_library,
    hypotheses,
    config,
)
    base_covariance = frequency_covariance(config.frequency_ratios, 16; model=:correlated)
    nominal_dictionaries = Dict{Symbol,Matrix{ComplexF64}}()
    combined_covariance = kron(Matrix{Float64}(I, 2, 2), base_covariance)
    for finalist in finalists
        definition = definitions[finalist.architecture_id]
        feature_matrices = Matrix{ComplexF64}[]
        for (state_index, state_key) in enumerate(definition.state_keys)
            features, _ = state_features(
                definition.state_ports[state_index],
                state_key[1],
                nominal_library,
                config,
            )
            push!(feature_matrices, features)
        end
        nominal_dictionaries[finalist.architecture_id] = reduce(vcat, feature_matrices) ./ sqrt(2)
    end

    rng = Xoshiro(RANDOM_SEED)
    rows = NamedTuple[]
    for sample in 1:ENSEMBLE_SIZE
        perturbed = perturbed_config(rng, config)
        lateral_shift = 0.10 * (2 * rand(rng) - 1)
        perturbed_hypotheses = shifted_hypotheses(hypotheses, lateral_shift)
        perturbed_library = response_library(perturbed_hypotheses, perturbed)
        common_errors = [(
            amplitude=exp.(0.03 .* randn(rng, config.element_count, 4)),
            delay=0.02 .* randn(rng, config.element_count, 4),
        ) for _ in 1:2]
        for finalist in finalists
            definition = definitions[finalist.architecture_id]
            feature_matrices = Matrix{ComplexF64}[]
            maximum_sigma = 0.0
            for (state_index, state_key) in enumerate(definition.state_keys)
                features, sigma = state_features(
                    definition.state_ports[state_index],
                    state_key[1],
                    perturbed_library,
                    perturbed;
                    amplitude_errors=common_errors[state_index].amplitude,
                    delay_errors=common_errors[state_index].delay,
                )
                push!(feature_matrices, features)
                maximum_sigma = max(maximum_sigma, sigma)
            end
            dictionary = reduce(vcat, feature_matrices) ./ sqrt(2)
            summary = operator_summary(dictionary, hypotheses, combined_covariance)
            fidelity = calibration_fidelity(
                nominal_dictionaries[finalist.architecture_id],
                dictionary,
                combined_covariance,
            )
            push!(rows, (
                architecture_id=finalist.architecture_id,
                transfer_group=finalist.transfer_group,
                sample,
                minimum_detectability=summary.minimum_detectability,
                nonlocal_coherence=summary.nonlocal_coherence,
                stable_rank=summary.stable_rank,
                q10_calibration_fidelity=fidelity.q10,
                max_sigma_after=maximum_sigma,
            ))
        end
    end
    rows, nominal_dictionaries
end

function summarize_multistate(finalists, ensemble_rows, baseline_minimum)
    summaries = NamedTuple[]
    for finalist in finalists
        rows = filter(row -> row.architecture_id == finalist.architecture_id, ensemble_rows)
        q10_minimum = quantile(getproperty.(rows, :minimum_detectability), 0.10)
        q90_coherence = quantile(getproperty.(rows, :nonlocal_coherence), 0.90)
        q10_rank = quantile(getproperty.(rows, :stable_rank), 0.10)
        q10_fidelity = quantile(getproperty.(rows, :q10_calibration_fidelity), 0.10)
        passivity = maximum(getproperty.(rows, :max_sigma_after)) <= 1 + 1.0e-12
        passed = finalist.minimum_detectability >= baseline_minimum &&
            q10_minimum >= 0.8 * baseline_minimum &&
            finalist.nonlocal_coherence <= 0.50 &&
            q90_coherence <= 0.60 &&
            finalist.stable_rank >= 10 &&
            q10_rank >= 8 &&
            q10_fidelity >= 0.90 &&
            passivity
        push!(summaries, (
            architecture_id=finalist.architecture_id,
            transfer_group=finalist.transfer_group,
            first_family=finalist.first_family,
            second_family=finalist.second_family,
            sector_preserving=finalist.sector_preserving,
            nominal_minimum_detectability=finalist.minimum_detectability,
            q10_minimum_detectability=q10_minimum,
            q10_minimum_over_ascan=q10_minimum / baseline_minimum,
            nominal_nonlocal_coherence=finalist.nonlocal_coherence,
            q90_nonlocal_coherence=q90_coherence,
            nominal_stable_rank=finalist.stable_rank,
            q10_stable_rank=q10_rank,
            q10_calibration_fidelity=q10_fidelity,
            maximum_sigma_after=maximum(getproperty.(rows, :max_sigma_after)),
            passed,
        ))
    end
    sort(summaries; by=row -> (row.passed, row.q10_stable_rank), rev=true)
end

function plot_multistate(single_rows, pair_rows, summaries)
    figure = plot(
        xlabel="minimum detectability",
        ylabel="nonlocal coherence",
        title="One versus two passive code states (four channels)",
        legend=:topright,
        size=(950, 650),
    )
    scatter!(
        figure,
        getproperty.(single_rows, :minimum_detectability),
        getproperty.(single_rows, :nonlocal_coherence);
        markersize=2.5,
        alpha=0.20,
        color=:gray,
        label="one state",
    )
    scatter!(
        figure,
        getproperty.(pair_rows, :minimum_detectability),
        getproperty.(pair_rows, :nonlocal_coherence);
        marker_z=getproperty.(pair_rows, :stable_rank),
        markersize=4,
        alpha=0.45,
        colorbar_title="stable rank",
        label="two states",
    )
    passed = filter(row -> row.passed, summaries)
    if !isempty(passed)
        scatter!(
            figure,
            getproperty.(passed, :nominal_minimum_detectability),
            getproperty.(passed, :nominal_nonlocal_coherence);
            marker=:star5,
            markersize=13,
            color=:black,
            label="robust pass",
        )
    end
    savefig(figure, joinpath(OUTPUT_ROOT, "o2c_multistate_search.png"))
end

function main_multistate_optimization()
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    library = response_library(hypotheses, config)
    baseline_rows = baseline_cases(library, hypotheses, config)
    baseline = only(filter(row ->
        row.frequency_grid == :seven && row.covariance_model == :correlated,
        baseline_rows,
    ))

    single_rows, ports_by_key, features_by_key, covariance = search_multistate_basis(
        library, hypotheses, config, baseline.minimum_detectability,
    )
    pair_rows, pair_definitions = multistate_pairs(
        single_rows,
        ports_by_key,
        features_by_key,
        covariance,
        hypotheses,
        baseline.minimum_detectability,
    )
    finalists, definitions = select_multistate_finalists(pair_rows, pair_definitions)
    ensemble_rows, nominal_dictionaries = evaluate_multistate_finalists(
        finalists, definitions, library, hypotheses, config,
    )
    summaries = summarize_multistate(
        finalists, ensemble_rows, baseline.minimum_detectability,
    )

    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_multistate_single_search.csv"),
        propertynames(first(single_rows)),
        single_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_multistate_pair_search.csv"),
        propertynames(first(pair_rows)),
        pair_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_multistate_ensemble.csv"),
        propertynames(first(ensemble_rows)),
        ensemble_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_multistate_summary.csv"),
        propertynames(first(summaries)),
        summaries,
    )
    jldsave(
        joinpath(OUTPUT_ROOT, "o2c_multistate_codebooks.jld2");
        config,
        hypotheses,
        definitions,
        finalists,
        nominal_dictionaries,
        summaries,
        random_seed=MULTISTATE_SEARCH_SEED,
        samples_per_family=MULTISTATE_SAMPLES_PER_FAMILY,
        ensemble_size=ENSEMBLE_SIZE,
    )
    plot_multistate(single_rows, pair_rows, summaries)

    for summary in summaries
        println(
            "[+] $(summary.architecture_id): Q10 min/A=" *
            "$(round(summary.q10_minimum_over_ascan; digits=3)), Q90 mu=" *
            "$(round(summary.q90_nonlocal_coherence; digits=3)), Q10 rank=" *
            "$(round(summary.q10_stable_rank; digits=3)), fidelity=" *
            "$(round(summary.q10_calibration_fidelity; digits=3)), " *
            "passed=$(summary.passed)",
        )
    end
    println("[+] wrote two-state O2c search to $OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_multistate_optimization()
end
