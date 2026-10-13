using LinearAlgebra
using Random
using Statistics
using JLD2
using Plots

include(joinpath(@__DIR__, "measurement_operator_space.jl"))
include(joinpath(@__DIR__, "elastic_operator_surrogate.jl"))

using .MeasurementOperatorSpace
using .ElasticOperatorSurrogate

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "dimensionless_elastic_born_v0",
)
const RANDOM_SEED = 0x0f3a51c2
const SAMPLES_PER_FAMILY = 300

function csv_value(value)
    value isa AbstractFloat && return isfinite(value) ? string(value) : "NaN"
    value isa Symbol && return String(value)
    string(value)
end

function write_csv(path::AbstractString, columns, rows)
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function separation_mask(hypotheses; minimum_separation=1.0)
    count = length(hypotheses)
    mask = falses(count, count)
    for first_index in 1:count, second_index in (first_index + 1):count
        separation = hypot(
            hypotheses[first_index].x_wavelengths - hypotheses[second_index].x_wavelengths,
            hypotheses[first_index].y_wavelengths - hypotheses[second_index].y_wavelengths,
        )
        if separation + 1.0e-12 >= minimum_separation
            mask[first_index, second_index] = true
            mask[second_index, first_index] = true
        end
    end
    any(mask) || error("no hypothesis pairs satisfy the separation gate")
    mask
end

function cheap_metrics(dictionary, nonlocal_mask)
    norms = vec(sqrt.(sum(abs2, dictionary; dims=1)))
    all(>(0), norms) || return (
        minimum_detectability=0.0,
        median_detectability=0.0,
        nonlocal_coherence=1.0,
    )
    normalized = dictionary ./ reshape(norms, 1, :)
    gram = abs.(normalized' * normalized)
    (
        minimum_detectability=Float64(minimum(norms)),
        median_detectability=Float64(median(norms)),
        nonlocal_coherence=Float64(maximum(gram[nonlocal_mask])),
    )
end

function orthonormal_columns(matrix, count)
    factor = qr(matrix)
    Matrix(factor.Q)[:, 1:count]
end

function random_unconstrained_ports(rng, config, port_count; complex_weights)
    raw = complex_weights ?
        randn(rng, ComplexF64, config.element_count, port_count) :
        randn(rng, config.element_count, port_count)
    orthonormal_columns(raw, port_count)
end

function random_sector_ports(rng, config, port_count; complex_weights)
    full_basis = cosine_ports(config.element_count, config)
    even_basis = full_basis[:, 1:2:end]
    odd_basis = full_basis[:, 2:2:end]
    even_count = cld(port_count, 2)
    odd_count = fld(port_count, 2)
    coefficient_type = complex_weights ? ComplexF64 : Float64
    even_coefficients = randn(rng, coefficient_type, size(even_basis, 2), even_count)
    even_ports = even_basis * orthonormal_columns(even_coefficients, even_count)
    if odd_count == 0
        return even_ports
    end
    odd_coefficients = randn(rng, coefficient_type, size(odd_basis, 2), odd_count)
    odd_ports = odd_basis * orthonormal_columns(odd_coefficients, odd_count)
    hcat(even_ports, odd_ports)
end

function candidate_ports(rng, family, config, port_count)
    if family == :real_unconstrained
        random_unconstrained_ports(rng, config, port_count; complex_weights=false)
    elseif family == :complex_unconstrained
        random_unconstrained_ports(rng, config, port_count; complex_weights=true)
    elseif family == :real_sector
        random_sector_ports(rng, config, port_count; complex_weights=false)
    elseif family == :complex_sector
        random_sector_ports(rng, config, port_count; complex_weights=true)
    else
        error("unknown codebook family: $family")
    end
end

function balanced_score(metrics, baseline_minimum, oracle_nonlocal_coherence)
    detectability_fraction = min(metrics.minimum_detectability / baseline_minimum, 1.0)
    ambiguity_closure = clamp(
        (1 - metrics.nonlocal_coherence) / (1 - oracle_nonlocal_coherence), 0.0, 1.0,
    )
    sqrt(detectability_fraction * ambiguity_closure)
end

function search_codebooks(library, hypotheses, config, baseline_minimum, oracle_coherence)
    rng = Xoshiro(RANDOM_SEED)
    nonlocal_mask = separation_mask(hypotheses)
    families = (:real_sector, :complex_sector, :real_unconstrained, :complex_unconstrained)
    rows = NamedTuple[]
    best_ports = Dict{Tuple{Int,Symbol},Matrix}()
    best_scores = Dict{Tuple{Int,Symbol},Float64}()

    for port_count in (2, 4), family in families
        best_scores[(port_count, family)] = -Inf
        for sample in 1:SAMPLES_PER_FAMILY
            ports = candidate_ports(rng, family, config, port_count)
            dictionary = multiport_dictionary(library, ports, config)
            metrics = cheap_metrics(dictionary, nonlocal_mask)
            score = balanced_score(metrics, baseline_minimum, oracle_coherence)
            push!(rows, (
                port_count,
                family,
                sample,
                minimum_detectability=metrics.minimum_detectability,
                median_detectability=metrics.median_detectability,
                nonlocal_coherence=metrics.nonlocal_coherence,
                balanced_score=score,
            ))
            if score > best_scores[(port_count, family)]
                best_scores[(port_count, family)] = score
                best_ports[(port_count, family)] = ports
            end
        end
    end
    rows, best_ports
end

function selected_rows(best_ports, library, hypotheses, config, baseline_minimum, oracle_coherence)
    nonlocal_mask = separation_mask(hypotheses)
    rows = NamedTuple[]
    dictionaries = Dict{Tuple{Int,Symbol},Matrix{ComplexF64}}()
    for key in sort(collect(keys(best_ports)); by=string)
        port_count, family = key
        dictionary = multiport_dictionary(library, best_ports[key], config)
        dictionaries[key] = dictionary
        cheap = cheap_metrics(dictionary, nonlocal_mask)
        full = operator_metrics(dictionary)
        push!(rows, (
            port_count,
            family,
            minimum_detectability=cheap.minimum_detectability,
            median_detectability=cheap.median_detectability,
            nonlocal_coherence=cheap.nonlocal_coherence,
            stable_rank=full.stable_rank,
            numerical_rank=full.numerical_rank,
            information_logdet_nats=full.information_logdet_nats,
            balanced_score=balanced_score(cheap, baseline_minimum, oracle_coherence),
        ))
    end
    rows, dictionaries
end

function cosine_baseline_rows(library, hypotheses, config, baseline_minimum, oracle_coherence)
    mask = separation_mask(hypotheses)
    rows = NamedTuple[]
    for port_count in (2, 4)
        dictionary = multiport_dictionary(library, cosine_ports(port_count, config), config)
        cheap = cheap_metrics(dictionary, mask)
        full = operator_metrics(dictionary)
        push!(rows, (
            port_count,
            family=:cosine_baseline,
            minimum_detectability=cheap.minimum_detectability,
            median_detectability=cheap.median_detectability,
            nonlocal_coherence=cheap.nonlocal_coherence,
            stable_rank=full.stable_rank,
            numerical_rank=full.numerical_rank,
            information_logdet_nats=full.information_logdet_nats,
            balanced_score=balanced_score(cheap, baseline_minimum, oracle_coherence),
        ))
    end
    rows
end

function plot_search(search_rows, selected, baselines)
    family_colors = Dict(
        :real_sector => :steelblue,
        :complex_sector => :navy,
        :real_unconstrained => :orange,
        :complex_unconstrained => :firebrick,
    )
    figures = Plots.Plot[]
    for port_count in (2, 4)
        figure = plot(
            xlabel="minimum detectability",
            ylabel="max coherence, separation >= 1 lambda_P",
            title="$port_count reciprocal ports",
            legend=:topright,
            size=(650, 520),
        )
        for family in keys(family_colors)
            rows = filter(row -> row.port_count == port_count && row.family == family, search_rows)
            scatter!(
                figure,
                getproperty.(rows, :minimum_detectability),
                getproperty.(rows, :nonlocal_coherence);
                markersize=2.5,
                alpha=0.35,
                color=family_colors[family],
                label=String(family),
            )
        end
        baseline = only(filter(row -> row.port_count == port_count, baselines))
        scatter!(
            figure,
            [baseline.minimum_detectability],
            [baseline.nonlocal_coherence];
            marker=:star5,
            markersize=10,
            color=:black,
            label="cosine baseline",
        )
        chosen = filter(row -> row.port_count == port_count, selected)
        scatter!(
            figure,
            getproperty.(chosen, :minimum_detectability),
            getproperty.(chosen, :nonlocal_coherence);
            marker=:diamond,
            markersize=7,
            markerstrokewidth=1.5,
            color=[family_colors[row.family] for row in chosen],
            label="selected per family",
        )
        push!(figures, figure)
    end
    savefig(
        plot(figures...; layout=(1, 2), size=(1300, 540), margin=7 * Plots.mm),
        joinpath(OUTPUT_ROOT, "codebook_search.png"),
    )
end

function main()
    mkpath(OUTPUT_ROOT)
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    library = response_library(hypotheses, config)
    uniform_dictionary = one_port_dictionary(
        library, _ -> ComplexF64.(uniform_port(config)), config,
    )
    oracle_dictionary = fmc_dictionary(library, config)
    nonlocal_mask = separation_mask(hypotheses)
    baseline_metrics = cheap_metrics(uniform_dictionary, nonlocal_mask)
    oracle_metrics = cheap_metrics(oracle_dictionary, nonlocal_mask)

    search_rows, best_ports = search_codebooks(
        library,
        hypotheses,
        config,
        baseline_metrics.minimum_detectability,
        oracle_metrics.nonlocal_coherence,
    )
    selected, selected_dictionaries = selected_rows(
        best_ports,
        library,
        hypotheses,
        config,
        baseline_metrics.minimum_detectability,
        oracle_metrics.nonlocal_coherence,
    )
    baselines = cosine_baseline_rows(
        library,
        hypotheses,
        config,
        baseline_metrics.minimum_detectability,
        oracle_metrics.nonlocal_coherence,
    )

    write_csv(
        joinpath(OUTPUT_ROOT, "codebook_search.csv"),
        propertynames(first(search_rows)),
        search_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "codebook_selected.csv"),
        propertynames(first(selected)),
        selected,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "codebook_baselines.csv"),
        propertynames(first(baselines)),
        baselines,
    )
    jldsave(
        joinpath(OUTPUT_ROOT, "optimized_codebooks.jld2");
        config,
        hypotheses,
        best_ports,
        selected_dictionaries,
        selected,
        baselines,
        random_seed=RANDOM_SEED,
        samples_per_family=SAMPLES_PER_FAMILY,
    )
    plot_search(search_rows, selected, baselines)

    for port_count in (2, 4)
        baseline = only(filter(row -> row.port_count == port_count, baselines))
        println("[+] $port_count-port cosine score=$(round(baseline.balanced_score; digits=4))")
        for row in filter(row -> row.port_count == port_count, selected)
            println(
                "[+] $port_count-port $(row.family): min_d=$(round(row.minimum_detectability; digits=4)), " *
                "mu_nonlocal=$(round(row.nonlocal_coherence; digits=4)), " *
                "rank=$(round(row.stable_rank; digits=3)), " *
                "score=$(round(row.balanced_score; digits=4))",
            )
        end
    end
    println("[+] wrote O2b codebook search to $OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
