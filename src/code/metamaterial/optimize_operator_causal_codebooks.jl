using LinearAlgebra
using Random
using Statistics
using JLD2
using Plots

include(joinpath(@__DIR__, "run_operator_causal_robustness.jl"))

const CAUSAL_SEARCH_SEED = 0x73c24a10
const CAUSAL_SAMPLES_PER_FAMILY = 500
const CAUSAL_FINALISTS_PER_FAMILY = 4

function causal_orthonormal_columns(matrix, count)
    factor = qr(matrix)
    Matrix(factor.Q)[:, 1:count]
end

function causal_sector_ports(rng, config; complex_weights)
    port_count = 4
    full_basis = cosine_ports(config.element_count, config)
    even_basis = full_basis[:, 1:2:end]
    odd_basis = full_basis[:, 2:2:end]
    coefficient_type = complex_weights ? ComplexF64 : Float64
    even_coefficients = randn(rng, coefficient_type, size(even_basis, 2), 2)
    odd_coefficients = randn(rng, coefficient_type, size(odd_basis, 2), 2)
    even_ports = even_basis * causal_orthonormal_columns(even_coefficients, 2)
    odd_ports = odd_basis * causal_orthonormal_columns(odd_coefficients, 2)
    hcat(even_ports, odd_ports)
end

function causal_gate_score(summary, baseline_minimum)
    detectability_margin = summary.minimum_detectability / baseline_minimum
    coherence_margin = 0.50 / summary.nonlocal_coherence
    rank_margin = summary.stable_rank / 10.0
    minimum((detectability_margin, coherence_margin, rank_margin))
end

function nominal_causal_candidate(ports, family, library, hypotheses, config, covariance)
    transfer = family == :real_static ? :constant_contractive : :delay_contractive
    raw, sigma_rows = causal_multiport_dictionary(
        library, ports, config; model=transfer,
    )
    features, _, _ = apply_pulse_and_select(
        raw, config.frequency_ratios, collect(eachindex(config.frequency_ratios)), 16,
    )
    summary = operator_summary(features, hypotheses, covariance)
    summary, features, maximum(getproperty.(sigma_rows, :sigma_after))
end

function search_causal_codebooks(library, hypotheses, config, baseline_minimum)
    rng = Xoshiro(CAUSAL_SEARCH_SEED)
    covariance = frequency_covariance(config.frequency_ratios, 16; model=:correlated)
    families = (:real_static, :complex_delay_contractive)
    rows = NamedTuple[]
    ports_by_key = Dict{Tuple{Symbol,Int},Matrix}()

    for family in families
        complex_weights = family == :complex_delay_contractive
        for sample in 1:CAUSAL_SAMPLES_PER_FAMILY
            ports = causal_sector_ports(rng, config; complex_weights)
            summary, _, maximum_sigma = nominal_causal_candidate(
                ports, family, library, hypotheses, config, covariance,
            )
            score = causal_gate_score(summary, baseline_minimum)
            push!(rows, (
                family,
                sample,
                minimum_detectability=summary.minimum_detectability,
                nonlocal_coherence=summary.nonlocal_coherence,
                stable_rank=summary.stable_rank,
                numerical_rank=summary.numerical_rank,
                maximum_sigma,
                gate_score=score,
                nominal_gate_passed=score >= 1.0,
            ))
            ports_by_key[(family, sample)] = ports
        end
        family_rows = filter(row -> row.family == family, rows)
        best = family_rows[argmax(getproperty.(family_rows, :gate_score))]
        println(
            "[+] nominal causal search $(family): best gate score=" *
            "$(round(best.gate_score; digits=3)), min=" *
            "$(round(best.minimum_detectability; digits=3)), mu=" *
            "$(round(best.nonlocal_coherence; digits=3)), rank=" *
            "$(round(best.stable_rank; digits=3))",
        )
    end
    rows, ports_by_key
end

function select_causal_finalists(search_rows, ports_by_key)
    finalists = NamedTuple[]
    finalist_ports = Dict{Symbol,Matrix}()
    for family in (:real_static, :complex_delay_contractive)
        family_rows = sort(
            filter(row -> row.family == family, search_rows);
            by=row -> row.gate_score,
            rev=true,
        )
        for (rank_index, row) in enumerate(first(family_rows, CAUSAL_FINALISTS_PER_FAMILY))
            id = Symbol(family, "_", rank_index)
            push!(finalists, merge((; architecture_id=id), row))
            finalist_ports[id] = ports_by_key[(family, row.sample)]
        end
    end
    finalists, finalist_ports
end

function causal_perturbation_cases(config, hypotheses)
    rng = Xoshiro(RANDOM_SEED)
    cases = NamedTuple[]
    for sample in 1:ENSEMBLE_SIZE
        perturbed = perturbed_config(rng, config)
        lateral_shift = 0.10 * (2 * rand(rng) - 1)
        amplitude_errors = exp.(0.03 .* randn(rng, config.element_count, 4))
        delay_errors = 0.02 .* randn(rng, config.element_count, 4)
        push!(cases, (
            sample,
            config=perturbed,
            hypotheses=shifted_hypotheses(hypotheses, lateral_shift),
            lateral_shift,
            amplitude_errors,
            delay_errors,
        ))
    end
    cases
end

function evaluate_causal_finalists(
    finalists,
    finalist_ports,
    cases,
    nominal_library,
    hypotheses,
    config,
)
    covariance = frequency_covariance(config.frequency_ratios, 16; model=:correlated)
    nominal_dictionaries = Dict{Symbol,Matrix{ComplexF64}}()
    for finalist in finalists
        _, features, _ = nominal_causal_candidate(
            finalist_ports[finalist.architecture_id],
            finalist.family,
            nominal_library,
            hypotheses,
            config,
            covariance,
        )
        nominal_dictionaries[finalist.architecture_id] = features
    end

    rows = NamedTuple[]
    for perturbation in cases
        perturbed_library = response_library(
            perturbation.hypotheses, perturbation.config,
        )
        for finalist in finalists
            transfer = finalist.family == :real_static ?
                :constant_contractive : :delay_contractive
            raw, sigma_rows = causal_multiport_dictionary(
                perturbed_library,
                finalist_ports[finalist.architecture_id],
                perturbation.config;
                model=transfer,
                amplitude_errors=perturbation.amplitude_errors,
                delay_errors=perturbation.delay_errors,
            )
            features, _, _ = apply_pulse_and_select(
                raw,
                config.frequency_ratios,
                collect(eachindex(config.frequency_ratios)),
                16,
            )
            summary = operator_summary(features, hypotheses, covariance)
            fidelity = calibration_fidelity(
                nominal_dictionaries[finalist.architecture_id], features, covariance,
            )
            push!(rows, (
                architecture_id=finalist.architecture_id,
                family=finalist.family,
                sample=perturbation.sample,
                minimum_detectability=summary.minimum_detectability,
                nonlocal_coherence=summary.nonlocal_coherence,
                stable_rank=summary.stable_rank,
                q10_calibration_fidelity=fidelity.q10,
                max_sigma_after=maximum(getproperty.(sigma_rows, :sigma_after)),
            ))
        end
    end
    rows, nominal_dictionaries
end

function summarize_causal_finalists(finalists, ensemble_rows, baseline_minimum)
    summaries = NamedTuple[]
    for finalist in finalists
        rows = filter(
            row -> row.architecture_id == finalist.architecture_id,
            ensemble_rows,
        )
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
            family=finalist.family,
            source_sample=finalist.sample,
            nominal_gate_score=finalist.gate_score,
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
    sort(summaries; by=row -> (
        row.passed,
        minimum((
            row.q10_minimum_over_ascan / 0.8,
            0.60 / row.q90_nonlocal_coherence,
            row.q10_stable_rank / 8,
            row.q10_calibration_fidelity / 0.90,
        )),
    ), rev=true)
end

function plot_causal_search(search_rows, summaries)
    figure = plot(
        xlabel="minimum detectability",
        ylabel="nonlocal coherence",
        title="O2c causal codebook search",
        legend=:topright,
        size=(850, 620),
    )
    colors = Dict(:real_static => :steelblue, :complex_delay_contractive => :firebrick)
    for family in keys(colors)
        rows = filter(row -> row.family == family, search_rows)
        scatter!(
            figure,
            getproperty.(rows, :minimum_detectability),
            getproperty.(rows, :nonlocal_coherence);
            marker_z=getproperty.(rows, :stable_rank),
            color=colors[family],
            colorbar=false,
            alpha=0.35,
            markersize=3,
            label=String(family),
        )
    end
    passed = filter(row -> row.passed, summaries)
    if !isempty(passed)
        scatter!(
            figure,
            getproperty.(passed, :nominal_minimum_detectability),
            getproperty.(passed, :nominal_nonlocal_coherence);
            marker=:star5,
            markersize=12,
            color=:black,
            label="robust pass",
        )
    end
    savefig(figure, joinpath(OUTPUT_ROOT, "o2c_codebook_search.png"))
end

function main_causal_optimization()
    mkpath(OUTPUT_ROOT)
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    library = response_library(hypotheses, config)
    baseline_rows = baseline_cases(library, hypotheses, config)
    baseline = only(filter(row ->
        row.frequency_grid == :seven && row.covariance_model == :correlated,
        baseline_rows,
    ))

    search_rows, ports_by_key = search_causal_codebooks(
        library, hypotheses, config, baseline.minimum_detectability,
    )
    finalists, finalist_ports = select_causal_finalists(search_rows, ports_by_key)
    perturbation_cases = causal_perturbation_cases(config, hypotheses)
    ensemble_rows, nominal_dictionaries = evaluate_causal_finalists(
        finalists,
        finalist_ports,
        perturbation_cases,
        library,
        hypotheses,
        config,
    )
    summaries = summarize_causal_finalists(
        finalists, ensemble_rows, baseline.minimum_detectability,
    )

    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_codebook_search.csv"),
        propertynames(first(search_rows)),
        search_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_codebook_finalists.csv"),
        propertynames(first(summaries)),
        summaries,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_codebook_ensemble.csv"),
        propertynames(first(ensemble_rows)),
        ensemble_rows,
    )
    jldsave(
        joinpath(OUTPUT_ROOT, "o2c_optimized_codebooks.jld2");
        config,
        hypotheses,
        finalist_ports,
        finalists,
        nominal_dictionaries,
        summaries,
        random_seed=CAUSAL_SEARCH_SEED,
        samples_per_family=CAUSAL_SAMPLES_PER_FAMILY,
        ensemble_size=ENSEMBLE_SIZE,
    )
    plot_causal_search(search_rows, summaries)

    for summary in summaries
        println(
            "[+] $(summary.architecture_id): nominal score=" *
            "$(round(summary.nominal_gate_score; digits=3)), Q10 min/A=" *
            "$(round(summary.q10_minimum_over_ascan; digits=3)), Q90 mu=" *
            "$(round(summary.q90_nonlocal_coherence; digits=3)), Q10 rank=" *
            "$(round(summary.q10_stable_rank; digits=3)), fidelity=" *
            "$(round(summary.q10_calibration_fidelity; digits=3)), " *
            "passed=$(summary.passed)",
        )
    end
    println("[+] wrote robust causal codebook optimization to $OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_causal_optimization()
end
