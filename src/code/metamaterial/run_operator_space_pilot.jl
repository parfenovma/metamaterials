using LinearAlgebra
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

function csv_value(value)
    if value isa AbstractFloat
        return isfinite(value) ? string(value) : "NaN"
    elseif value isa Symbol
        return String(value)
    elseif value isa Tuple || value isa AbstractVector
        return join(string.(value), "|")
    end
    string(value)
end

function write_csv(path::AbstractString, columns, rows)
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, column)) for column in columns), ','))
        end
    end
    path
end

function worst_pair(dictionary, hypotheses; minimum_separation=0.0)
    whitened = whiten_dictionary(dictionary)
    gram = normalized_gram(whitened)
    best_value = -Inf
    best_pair = (0, 0)
    for first_index in axes(gram, 1), second_index in (first_index + 1):size(gram, 2)
        separation = hypot(
            hypotheses[first_index].x_wavelengths - hypotheses[second_index].x_wavelengths,
            hypotheses[first_index].y_wavelengths - hypotheses[second_index].y_wavelengths,
        )
        separation + 1.0e-12 >= minimum_separation || continue
        if gram[first_index, second_index] > best_value
            best_value = gram[first_index, second_index]
            best_pair = (first_index, second_index)
        end
    end
    best_pair..., Float64(best_value), gram
end

function mirror_coherence(dictionary, mirror_indices)
    gram = normalized_gram(whiten_dictionary(dictionary))
    values = Float64[
        gram[index, mirror_indices[index]]
        for index in eachindex(mirror_indices)
        if index != mirror_indices[index]
    ]
    isempty(values) ? (0.0, 0.0) : (maximum(values), median(values))
end

function center_hypothesis_index(hypotheses)
    distances = [
        hypot(hypothesis.x_wavelengths - 4.25, hypothesis.y_wavelengths)
        for hypothesis in hypotheses
    ]
    argmin(distances)
end

function architecture_row(
    id,
    dictionary,
    hypotheses,
    mirror_indices;
    electronic_channels,
    acquisitions,
    txrx_constraint,
    code_family,
    frequency_anchors,
    role=:candidate,
    lift_status=:translation_kz,
)
    metrics = operator_metrics(dictionary)
    first_index, second_index, pair_coherence, _ = worst_pair(dictionary, hypotheses)
    nonlocal_first, nonlocal_second, nonlocal_coherence, _ = worst_pair(
        dictionary, hypotheses; minimum_separation=1.0,
    )
    mirror_maximum, mirror_median = mirror_coherence(dictionary, mirror_indices)
    center_index = center_hypothesis_index(hypotheses)
    center_detectability = norm(view(whiten_dictionary(dictionary), :, center_index))
    (
        architecture_id=id,
        role,
        electronic_channels,
        acquisitions,
        frequency_anchors,
        measurement_count=metrics.measurement_count,
        hypothesis_count=metrics.hypothesis_count,
        txrx_constraint,
        code_family,
        lift_status,
        center_detectability=Float64(center_detectability),
        minimum_detectability=metrics.minimum_detectability,
        median_detectability=metrics.median_detectability,
        maximum_pair_coherence=metrics.maximum_pair_coherence,
        maximum_nonlocal_coherence_1lambda=nonlocal_coherence,
        maximum_mirror_coherence=mirror_maximum,
        median_mirror_coherence=mirror_median,
        stable_rank=metrics.stable_rank,
        numerical_rank=metrics.numerical_rank,
        information_logdet_nats=metrics.information_logdet_nats,
        worst_first_index=first_index,
        worst_second_index=second_index,
        worst_pair_coherence=pair_coherence,
        nonlocal_first_index=nonlocal_first,
        nonlocal_second_index=nonlocal_second,
        worst_nonlocal_coherence=nonlocal_coherence,
    )
end

function physical_dictionaries(library, config)
    focus_x = 4.25
    focus_y = 0.0
    Dict{Symbol,Matrix{ComplexF64}}(
        :a_scan_uniform => one_port_dictionary(
            library, _ -> ComplexF64.(uniform_port(config)), config,
        ),
        :fixed_focus_reciprocal => one_port_dictionary(
            library,
            frequency -> fixed_focus_port(frequency, focus_x, focus_y, config),
            config,
        ),
        :coded_even_single_port => one_port_dictionary(
            library,
            frequency -> dispersive_port(frequency, config; symmetry=:even),
            config,
        ),
        :coded_asymmetric_single_port => one_port_dictionary(
            library,
            frequency -> dispersive_port(frequency, config; symmetry=:asymmetric),
            config,
        ),
        :few_port_2 => multiport_dictionary(library, cosine_ports(2, config), config),
        :few_port_4 => multiport_dictionary(library, cosine_ports(4, config), config),
        :paut_fmc_oracle => fmc_dictionary(library, config),
    )
end

function resource_contract(id, config)
    if id in (:a_scan_uniform, :fixed_focus_reciprocal, :coded_even_single_port,
              :coded_asymmetric_single_port)
        return (1, 1, :common_monostatic, id == :a_scan_uniform ? :uniform : :reciprocal_code)
    elseif id == :few_port_2
        return (2, 2, :reciprocal_matched, :cosine_symmetry_ports)
    elseif id == :few_port_4
        return (4, 4, :reciprocal_matched, :cosine_symmetry_ports)
    elseif id == :paut_fmc_oracle
        return (config.element_count, config.element_count, :independent, :full_matrix_capture)
    end
    error("missing resource contract for $id")
end

function plot_frontiers(architecture_rows, ideal_rows)
    labels = [
        "A-scan", "fixed focus", "even code", "asym. code", "2 ports", "4 ports", "FMC",
    ]
    measurements = getproperty.(architecture_rows, :measurement_count)
    colors = [:gray, :orange, :orange, :orange, :steelblue, :steelblue, :firebrick]

    ambiguity_plot = bar(
        labels,
        getproperty.(architecture_rows, :maximum_nonlocal_coherence_1lambda);
        color=colors,
        xlabel="architecture",
        ylabel="max coherence, separation >= 1 lambda_P",
        title="Localization ambiguity (lower is better)",
        label=false,
        xrotation=22,
        ylim=(0, 1.05),
    )
    for (index, count) in enumerate(measurements)
        annotate!(ambiguity_plot, index, architecture_rows[index].maximum_nonlocal_coherence_1lambda + 0.035,
                  text("m=$count", 7))
    end

    rank_plot = bar(
        labels,
        getproperty.(architecture_rows, :stable_rank);
        color=colors,
        xlabel="architecture",
        ylabel="stable rank",
        title="Stable measurement degrees of freedom",
        label=false,
        xrotation=22,
    )
    for (index, count) in enumerate(measurements)
        annotate!(rank_plot, index, architecture_rows[index].stable_rank + 0.8, text("m=$count", 7))
    end
    savefig(
        plot(ambiguity_plot, rank_plot; layout=(2, 1), size=(1000, 850), margin=8 * Plots.mm),
        joinpath(OUTPUT_ROOT, "architecture_frontier.png"),
    )

    compression_plot = plot(
        getproperty.(ideal_rows, :ideal_measurements),
        getproperty.(ideal_rows, :retained_dictionary_energy);
        xscale=:log2,
        marker=:circle,
        linewidth=2,
        xlabel="ideal linear measurements K",
        ylabel="retained dictionary energy",
        ylim=(0, 1.02),
        label="SVD compression ceiling",
        title="O2 ideal-K compression frontier",
        size=(850, 500),
        xticks=(getproperty.(ideal_rows, :ideal_measurements),
                string.(getproperty.(ideal_rows, :ideal_measurements))),
    )
    savefig(compression_plot, joinpath(OUTPUT_ROOT, "ideal_k_compression.png"))
end

function plot_mode_ablation(mode_rows)
    labels = ["A-scan", "fixed focus", "even code", "asym. code", "2 ports", "4 ports", "FMC"]
    full_rank = getproperty.(mode_rows, :full_elastic_stable_rank)
    pressure_rank = getproperty.(mode_rows, :pressure_only_stable_rank)
    positions = collect(eachindex(labels))
    figure = bar(
        positions .- 0.18,
        pressure_rank;
        bar_width=0.34,
        label="P only",
        ylabel="stable rank",
        xlabel="architecture",
        title="Elastic-mode diversity ablation",
        size=(1000, 560),
        margin=8 * Plots.mm,
        xticks=(positions, labels),
        xrotation=22,
    )
    bar!(
        figure,
        positions .+ 0.18,
        full_rank;
        bar_width=0.34,
        label="P + S",
    )
    savefig(figure, joinpath(OUTPUT_ROOT, "mode_ablation.png"))
end

function main()
    mkpath(OUTPUT_ROOT)
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    mirror_indices = mirror_hypothesis_index(hypotheses)
    library = response_library(hypotheses, config)
    dictionaries = physical_dictionaries(library, config)
    pressure_only_config = SurrogateConfig(shear_amplitude=0.0)
    pressure_only_library = response_library(hypotheses, pressure_only_config)
    pressure_only_dictionaries = physical_dictionaries(pressure_only_library, pressure_only_config)

    reciprocity_error = maximum_reciprocity_error(library)
    mirror_error = maximum_mirror_error(library, hypotheses)
    reciprocity_error <= 1.0e-12 || error("reciprocity gate failed: $reciprocity_error")
    mirror_error <= 1.0e-12 || error("mirror equivariance gate failed: $mirror_error")

    architecture_rows = NamedTuple[]
    ambiguity_rows = NamedTuple[]
    for id in (
        :a_scan_uniform,
        :fixed_focus_reciprocal,
        :coded_even_single_port,
        :coded_asymmetric_single_port,
        :few_port_2,
        :few_port_4,
        :paut_fmc_oracle,
    )
        channels, acquisitions, txrx, code_family = resource_contract(id, config)
        role = id == :a_scan_uniform ? :baseline : id == :paut_fmc_oracle ? :bound : :candidate
        row = architecture_row(
            id,
            dictionaries[id],
            hypotheses,
            mirror_indices;
            electronic_channels=channels,
            acquisitions,
            txrx_constraint=txrx,
            code_family,
            frequency_anchors=length(config.frequency_ratios),
            role,
        )
        push!(architecture_rows, row)
        first_hypothesis = hypotheses[row.worst_first_index]
        second_hypothesis = hypotheses[row.worst_second_index]
        push!(ambiguity_rows, (
            architecture_id=id,
            pair_kind=:all,
            first_id=first_hypothesis.id,
            first_x=first_hypothesis.x_wavelengths,
            first_y=first_hypothesis.y_wavelengths,
            second_id=second_hypothesis.id,
            second_x=second_hypothesis.x_wavelengths,
            second_y=second_hypothesis.y_wavelengths,
            coherence=row.worst_pair_coherence,
        ))
        nonlocal_first = hypotheses[row.nonlocal_first_index]
        nonlocal_second = hypotheses[row.nonlocal_second_index]
        push!(ambiguity_rows, (
            architecture_id=id,
            pair_kind=:separation_ge_1lambda,
            first_id=nonlocal_first.id,
            first_x=nonlocal_first.x_wavelengths,
            first_y=nonlocal_first.y_wavelengths,
            second_id=nonlocal_second.id,
            second_x=nonlocal_second.x_wavelengths,
            second_y=nonlocal_second.y_wavelengths,
            coherence=row.worst_nonlocal_coherence,
        ))
    end

    full_dictionary = dictionaries[:paut_fmc_oracle]
    maximum_ideal_count = min(size(full_dictionary)...)
    ideal_counts = filter(<=(maximum_ideal_count), [1, 2, 4, 8, 16, 32, 64, 128])
    ideal_rows = NamedTuple[]
    ideal_dictionaries = Dict{Int,Matrix{ComplexF64}}()
    singular_values = svdvals(full_dictionary)
    for count in ideal_counts
        compressed, retained_energy, _ = ideal_compression_dictionary(full_dictionary, count)
        ideal_dictionaries[count] = compressed
        metrics = operator_metrics(compressed)
        push!(ideal_rows, (
            ideal_measurements=count,
            retained_dictionary_energy=retained_energy,
            minimum_detectability=metrics.minimum_detectability,
            median_detectability=metrics.median_detectability,
            maximum_pair_coherence=metrics.maximum_pair_coherence,
            stable_rank=metrics.stable_rank,
            numerical_rank=metrics.numerical_rank,
            information_logdet_nats=metrics.information_logdet_nats,
        ))
    end


    mode_rows = NamedTuple[]
    for id in (
        :a_scan_uniform,
        :fixed_focus_reciprocal,
        :coded_even_single_port,
        :coded_asymmetric_single_port,
        :few_port_2,
        :few_port_4,
        :paut_fmc_oracle,
    )
        full_metrics = operator_metrics(dictionaries[id])
        pressure_metrics = operator_metrics(pressure_only_dictionaries[id])
        _, _, full_nonlocal, _ = worst_pair(
            dictionaries[id], hypotheses; minimum_separation=1.0,
        )
        _, _, pressure_nonlocal, _ = worst_pair(
            pressure_only_dictionaries[id], hypotheses; minimum_separation=1.0,
        )
        push!(mode_rows, (
            architecture_id=id,
            pressure_only_stable_rank=pressure_metrics.stable_rank,
            full_elastic_stable_rank=full_metrics.stable_rank,
            stable_rank_ratio=full_metrics.stable_rank / pressure_metrics.stable_rank,
            pressure_only_information_logdet=pressure_metrics.information_logdet_nats,
            full_elastic_information_logdet=full_metrics.information_logdet_nats,
            information_delta_nats=full_metrics.information_logdet_nats -
                pressure_metrics.information_logdet_nats,
            pressure_only_nonlocal_coherence=pressure_nonlocal,
            full_elastic_nonlocal_coherence=full_nonlocal,
        ))
    end

    write_csv(
        joinpath(OUTPUT_ROOT, "architecture_metrics.csv"),
        propertynames(first(architecture_rows)),
        architecture_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "ideal_k_frontier.csv"),
        propertynames(first(ideal_rows)),
        ideal_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "ambiguity_pairs.csv"),
        propertynames(first(ambiguity_rows)),
        ambiguity_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "mode_ablation.csv"),
        propertynames(first(mode_rows)),
        mode_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "hypothesis_metadata.csv"),
        (:index, :id, :x_wavelengths, :y_wavelengths, :mirror_index),
        [(
            index=index,
            id=hypothesis.id,
            x_wavelengths=hypothesis.x_wavelengths,
            y_wavelengths=hypothesis.y_wavelengths,
            mirror_index=mirror_indices[index],
        ) for (index, hypothesis) in enumerate(hypotheses)],
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "singular_values.csv"),
        (:index, :singular_value, :cumulative_energy),
        [(
            index=index,
            singular_value=value,
            cumulative_energy=sum(abs2, singular_values[1:index]) / sum(abs2, singular_values),
        ) for (index, value) in enumerate(singular_values)],
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "operator_checks.csv"),
        (:check, :value, :tolerance, :passed),
        [
            (check=:reciprocity_error, value=reciprocity_error, tolerance=1.0e-12,
             passed=reciprocity_error <= 1.0e-12),
            (check=:mirror_equivariance_error, value=mirror_error, tolerance=1.0e-12,
             passed=mirror_error <= 1.0e-12),
        ],
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "surrogate_config.csv"),
        (:parameter, :value),
        [
            (parameter=:pressure_speed, value=config.pressure_speed),
            (parameter=:shear_speed, value=config.shear_speed),
            (parameter=:reference_frequency, value=config.reference_frequency),
            (parameter=:frequency_ratios, value=config.frequency_ratios),
            (parameter=:aperture_wavelengths, value=config.aperture_wavelengths),
            (parameter=:element_count, value=config.element_count),
            (parameter=:pressure_attenuation_per_wavelength,
             value=config.pressure_attenuation_per_wavelength),
            (parameter=:shear_attenuation_per_wavelength,
             value=config.shear_attenuation_per_wavelength),
            (parameter=:shear_amplitude, value=config.shear_amplitude),
        ],
    )

    jldsave(
        joinpath(OUTPUT_ROOT, "operator_dictionary.jld2");
        config,
        hypotheses,
        mirror_indices,
        library,
        dictionaries,
        pressure_only_config,
        pressure_only_dictionaries,
        ideal_dictionaries,
        singular_values,
        architecture_rows,
        ideal_rows,
        mode_rows,
        reciprocity_error,
        mirror_error,
    )
    plot_frontiers(architecture_rows, ideal_rows)
    plot_mode_ablation(mode_rows)

    println("[+] O1 reciprocity error: $reciprocity_error")
    println("[+] O1 mirror equivariance error: $mirror_error")
    for row in architecture_rows
        println(
            "[+] $(row.architecture_id): measurements=$(row.measurement_count), " *
            "coherence=$(round(row.maximum_pair_coherence; digits=5)), " *
            "stable_rank=$(round(row.stable_rank; digits=3))",
        )
    end
    println("[+] wrote operator-space pilot to $OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
