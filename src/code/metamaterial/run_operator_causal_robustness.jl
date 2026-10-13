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
const CODEBOOK_PATH = joinpath(OUTPUT_ROOT, "optimized_codebooks.jld2")
const ENSEMBLE_SIZE = 120
const RANDOM_SEED = 0xc2a51b17

function csv_value(value)
    value isa AbstractFloat && return isfinite(value) ? string(value) : "NaN"
    value isa Symbol && return String(value)
    value isa Bool && return value ? "true" : "false"
    value isa AbstractVector && return join(string.(value), '|')
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

function hypothesis_separation_mask(hypotheses; minimum_separation=1.0)
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
    mask
end

function hann_burst_spectrum(frequency_ratios; cycles=5.0, sample_count=12001)
    times = collect(range(0.0, cycles; length=sample_count))
    step = times[2] - times[1]
    waveform = @. sin(2pi * times) * 0.5 * (1 - cos(2pi * times / cycles))
    amplitudes = Float64[
        abs(sum(waveform .* exp.(-2pi * im * frequency * times)) * step)
        for frequency in frequency_ratios
    ]
    maximum(amplitudes) > 0 || error("pulse spectrum vanished")
    amplitudes
end

function frequency_covariance(frequencies, block_size; model=:white)
    count = length(frequencies)
    frequency_covariance = if model == :white
        Matrix{Float64}(I, count, count)
    elseif model == :correlated
        correlation_length = 0.12
        correlated_fraction = 0.60
        kernel = [
            exp(-0.5 * ((first - second) / correlation_length)^2)
            for first in frequencies, second in frequencies
        ]
        (1 - correlated_fraction) .* Matrix{Float64}(I, count, count) .+
            correlated_fraction .* kernel
    else
        throw(ArgumentError("unknown frequency covariance model: $model"))
    end
    kron(frequency_covariance, Matrix{Float64}(I, block_size, block_size))
end

function apply_pulse_and_select(dictionary, all_frequencies, indices, block_size)
    selected_frequencies = all_frequencies[indices]
    raw_weights = hann_burst_spectrum(selected_frequencies)
    weights = raw_weights / norm(raw_weights)
    selected = Matrix{ComplexF64}(undef,
        length(indices) * block_size,
        size(dictionary, 2),
    )
    destination = 1
    for (weight_index, frequency_index) in enumerate(indices)
        source_rows = ((frequency_index - 1) * block_size + 1):(frequency_index * block_size)
        destination_rows = destination:(destination + block_size - 1)
        selected[destination_rows, :] = weight_index <= length(weights) ?
            weights[weight_index] .* dictionary[source_rows, :] : dictionary[source_rows, :]
        destination += block_size
    end
    selected, selected_frequencies, weights
end

function transfer_parameters(base_ports)
    amplitudes = abs.(base_ports)
    delays = mod.(-angle.(complex.(base_ports)) ./ (2pi), 1.0)
    amplitudes, delays
end

function frequency_ports(
    base_ports,
    frequency;
    model,
    amplitude_errors=ones(size(base_ports)),
    delay_errors=zeros(size(base_ports)),
)
    size(amplitude_errors) == size(base_ports) || throw(DimensionMismatch("amplitude error shape mismatch"))
    size(delay_errors) == size(base_ports) || throw(DimensionMismatch("delay error shape mismatch"))
    if model in (:carrier_oracle, :constant_contractive)
        raw = complex.(base_ports) .* amplitude_errors .* exp.(-2pi * im * frequency .* delay_errors)
    elseif model in (:delay_unscaled, :delay_contractive)
        amplitudes, delays = transfer_parameters(base_ports)
        raw = amplitudes .* amplitude_errors .* exp.(-2pi * im * frequency .* (delays .+ delay_errors))
    else
        throw(ArgumentError("unknown transfer model: $model"))
    end
    sigma_before = opnorm(raw)
    contractive = model in (:constant_contractive, :delay_contractive)
    ports = contractive ? raw ./ max(1.0, sigma_before) : raw
    ports, Float64(sigma_before), Float64(opnorm(ports))
end

function causal_multiport_dictionary(
    library,
    base_ports,
    config;
    model,
    amplitude_errors=ones(size(base_ports)),
    delay_errors=zeros(size(base_ports)),
)
    port_count = size(base_ports, 2)
    hypothesis_count = size(library, 4)
    block_size = port_count^2
    dictionary = Matrix{ComplexF64}(undef,
        length(config.frequency_ratios) * block_size,
        hypothesis_count,
    )
    sigma_rows = NamedTuple[]
    energy_scale = inv(sqrt(port_count))
    for (frequency_index, frequency) in enumerate(config.frequency_ratios)
        ports, sigma_before, sigma_after = frequency_ports(
            base_ports,
            frequency;
            model,
            amplitude_errors,
            delay_errors,
        )
        push!(sigma_rows, (; frequency, sigma_before, sigma_after))
        row = (frequency_index - 1) * block_size + 1
        for hypothesis_index in 1:hypothesis_count
            response = view(library, :, :, frequency_index, hypothesis_index)
            compressed = transpose(ports) * response * ports .* energy_scale
            destination = row
            for transmit_port in 1:port_count, receive_port in 1:port_count
                dictionary[destination, hypothesis_index] = compressed[receive_port, transmit_port]
                destination += 1
            end
        end
    end
    dictionary, sigma_rows
end

function operator_summary(dictionary, hypotheses, covariance)
    whitened = whiten_dictionary(dictionary; noise_covariance=covariance)
    norms = vec(sqrt.(sum(abs2, whitened; dims=1)))
    normalized = whitened ./ reshape(norms, 1, :)
    gram = abs.(normalized' * normalized)
    mask = hypothesis_separation_mask(hypotheses)
    singular_values = svdvals(whitened)
    (
        minimum_detectability=Float64(minimum(norms)),
        median_detectability=Float64(median(norms)),
        nonlocal_coherence=Float64(maximum(gram[mask])),
        stable_rank=Float64(sum(abs2, singular_values) / first(singular_values)^2),
        numerical_rank=count(>(max(size(whitened)...) * eps(Float64) * first(singular_values)), singular_values),
    )
end

function calibration_fidelity(nominal, perturbed, covariance)
    nominal_white = whiten_dictionary(nominal; noise_covariance=covariance)
    perturbed_white = whiten_dictionary(perturbed; noise_covariance=covariance)
    values = Float64[]
    for index in axes(nominal, 2)
        first = view(nominal_white, :, index)
        second = view(perturbed_white, :, index)
        push!(values, Float64(abs(dot(first, second)) / (norm(first) * norm(second))))
    end
    (
        minimum=minimum(values),
        q10=quantile(values, 0.10),
        median=median(values),
    )
end

function shifted_hypotheses(hypotheses, lateral_shift)
    [
        DefectHypothesis(
            hypothesis.id,
            hypothesis.x_wavelengths,
            hypothesis.y_wavelengths - lateral_shift,
            copy(hypothesis.scattering_tensor),
        )
        for hypothesis in hypotheses
    ]
end

function perturbed_config(rng, nominal)
    pressure_scale = 1 + 0.01 * (2 * rand(rng) - 1)
    shear_scale = 1 + 0.01 * (2 * rand(rng) - 1)
    attenuation_scale_p = 1 + 0.20 * (2 * rand(rng) - 1)
    attenuation_scale_s = 1 + 0.20 * (2 * rand(rng) - 1)
    SurrogateConfig(
        pressure_speed=nominal.pressure_speed * pressure_scale,
        shear_speed=nominal.shear_speed * shear_scale,
        reference_frequency=nominal.reference_frequency,
        frequency_ratios=copy(nominal.frequency_ratios),
        aperture_wavelengths=nominal.aperture_wavelengths,
        element_count=nominal.element_count,
        pressure_attenuation_per_wavelength=
            nominal.pressure_attenuation_per_wavelength * attenuation_scale_p,
        shear_attenuation_per_wavelength=
            nominal.shear_attenuation_per_wavelength * attenuation_scale_s,
        shear_amplitude=nominal.shear_amplitude,
    )
end

function nominal_cases(library, hypotheses, config, real_ports, complex_ports)
    models = [
        (id=:real_static, ports=real_ports, transfer=:constant_contractive),
        (id=:real_delay, ports=real_ports, transfer=:delay_contractive),
        (id=:complex_carrier_oracle, ports=complex_ports, transfer=:carrier_oracle),
        (id=:complex_delay_unscaled, ports=complex_ports, transfer=:delay_unscaled),
        (id=:complex_delay_contractive, ports=complex_ports, transfer=:delay_contractive),
    ]
    grids = [
        (id=:center, indices=[4]),
        (id=:three, indices=[1, 4, 7]),
        (id=:seven, indices=collect(1:7)),
    ]
    rows = NamedTuple[]
    dictionaries = Dict{Tuple{Symbol,Symbol,Symbol},Matrix{ComplexF64}}()
    passivity_rows = NamedTuple[]
    for model in models
        raw_dictionary, sigma_rows = causal_multiport_dictionary(
            library, model.ports, config; model=model.transfer,
        )
        for sigma in sigma_rows
            push!(passivity_rows, (
                architecture_id=model.id,
                transfer_model=model.transfer,
                frequency=sigma.frequency,
                sigma_before=sigma.sigma_before,
                sigma_after=sigma.sigma_after,
                passive=sigma.sigma_after <= 1 + 1.0e-12,
            ))
        end
        for grid in grids
            features, selected_frequencies, weights = apply_pulse_and_select(
                raw_dictionary,
                config.frequency_ratios,
                grid.indices,
                16,
            )
            covariance_models = grid.id == :center ? (:white,) : (:white, :correlated)
            for covariance_model in covariance_models
                covariance = frequency_covariance(
                    selected_frequencies, 16; model=covariance_model,
                )
                summary = operator_summary(features, hypotheses, covariance)
                dictionaries[(model.id, grid.id, covariance_model)] = features
                push!(rows, (
                    architecture_id=model.id,
                    transfer_model=model.transfer,
                    frequency_grid=grid.id,
                    covariance_model,
                    frequency_count=length(grid.indices),
                    pulse_weights=weights,
                    minimum_detectability=summary.minimum_detectability,
                    median_detectability=summary.median_detectability,
                    nonlocal_coherence=summary.nonlocal_coherence,
                    stable_rank=summary.stable_rank,
                    numerical_rank=summary.numerical_rank,
                    max_sigma_before=maximum(getproperty.(sigma_rows, :sigma_before)),
                    max_sigma_after=maximum(getproperty.(sigma_rows, :sigma_after)),
                ))
            end
        end
    end
    rows, dictionaries, passivity_rows
end

function baseline_cases(library, hypotheses, config)
    uniform = ComplexF64.(uniform_port(config))
    raw = one_port_dictionary(library, _ -> uniform, config)
    grids = [
        (id=:center, indices=[4]),
        (id=:three, indices=[1, 4, 7]),
        (id=:seven, indices=collect(1:7)),
    ]
    rows = NamedTuple[]
    for grid in grids
        features, selected_frequencies, weights = apply_pulse_and_select(
            raw, config.frequency_ratios, grid.indices, 1,
        )
        covariance_models = grid.id == :center ? (:white,) : (:white, :correlated)
        for covariance_model in covariance_models
            covariance = frequency_covariance(selected_frequencies, 1; model=covariance_model)
            summary = operator_summary(features, hypotheses, covariance)
            push!(rows, (
                architecture_id=:a_scan_uniform,
                frequency_grid=grid.id,
                covariance_model,
                frequency_count=length(grid.indices),
                pulse_weights=weights,
                minimum_detectability=summary.minimum_detectability,
                median_detectability=summary.median_detectability,
                nonlocal_coherence=summary.nonlocal_coherence,
                stable_rank=summary.stable_rank,
                numerical_rank=summary.numerical_rank,
            ))
        end
    end
    rows
end

function robustness_ensemble(
    nominal_library,
    hypotheses,
    config,
    real_ports,
    complex_ports,
)
    rng = Xoshiro(RANDOM_SEED)
    selected = [
        (id=:real_static, ports=real_ports, transfer=:constant_contractive),
        (id=:complex_delay_contractive, ports=complex_ports, transfer=:delay_contractive),
    ]
    nominal = Dict{Symbol,Matrix{ComplexF64}}()
    frequencies = config.frequency_ratios
    covariance = frequency_covariance(frequencies, 16; model=:correlated)
    for candidate in selected
        raw, _ = causal_multiport_dictionary(
            nominal_library, candidate.ports, config; model=candidate.transfer,
        )
        features, _, _ = apply_pulse_and_select(raw, frequencies, collect(1:7), 16)
        nominal[candidate.id] = features
    end

    rows = NamedTuple[]
    for sample in 1:ENSEMBLE_SIZE
        perturbed = perturbed_config(rng, config)
        lateral_shift = 0.10 * (2 * rand(rng) - 1)
        perturbed_hypotheses = shifted_hypotheses(hypotheses, lateral_shift)
        perturbed_library = response_library(perturbed_hypotheses, perturbed)
        for candidate in selected
            amplitude_errors = exp.(0.03 .* randn(rng, size(candidate.ports)))
            delay_errors = 0.02 .* randn(rng, size(candidate.ports))
            raw, sigma_rows = causal_multiport_dictionary(
                perturbed_library,
                candidate.ports,
                perturbed;
                model=candidate.transfer,
                amplitude_errors,
                delay_errors,
            )
            features, _, _ = apply_pulse_and_select(raw, frequencies, collect(1:7), 16)
            summary = operator_summary(features, hypotheses, covariance)
            fidelity = calibration_fidelity(nominal[candidate.id], features, covariance)
            push!(rows, (
                architecture_id=candidate.id,
                sample,
                pressure_speed=perturbed.pressure_speed,
                shear_speed=perturbed.shear_speed,
                pressure_attenuation=perturbed.pressure_attenuation_per_wavelength,
                shear_attenuation=perturbed.shear_attenuation_per_wavelength,
                lateral_shift,
                minimum_detectability=summary.minimum_detectability,
                median_detectability=summary.median_detectability,
                nonlocal_coherence=summary.nonlocal_coherence,
                stable_rank=summary.stable_rank,
                numerical_rank=summary.numerical_rank,
                minimum_calibration_fidelity=fidelity.minimum,
                q10_calibration_fidelity=fidelity.q10,
                median_calibration_fidelity=fidelity.median,
                max_sigma_before=maximum(getproperty.(sigma_rows, :sigma_before)),
                max_sigma_after=maximum(getproperty.(sigma_rows, :sigma_after)),
            ))
        end
    end
    rows
end

function robustness_summary(ensemble_rows, baseline_minimum, nominal_rows)
    rows = NamedTuple[]
    for id in (:real_static, :complex_delay_contractive)
        selected = filter(row -> row.architecture_id == id, ensemble_rows)
        nominal = only(filter(row ->
            row.architecture_id == id &&
            row.frequency_grid == :seven &&
            row.covariance_model == :correlated,
            nominal_rows,
        ))
        minimum_values = getproperty.(selected, :minimum_detectability)
        coherence_values = getproperty.(selected, :nonlocal_coherence)
        rank_values = getproperty.(selected, :stable_rank)
        fidelity_values = getproperty.(selected, :q10_calibration_fidelity)
        q10_minimum = quantile(minimum_values, 0.10)
        q90_coherence = quantile(coherence_values, 0.90)
        q10_rank = quantile(rank_values, 0.10)
        q10_fidelity = quantile(fidelity_values, 0.10)
        passivity = maximum(getproperty.(selected, :max_sigma_after)) <= 1 + 1.0e-12
        push!(rows, (
            architecture_id=id,
            ensemble_size=length(selected),
            nominal_minimum_detectability=nominal.minimum_detectability,
            q10_minimum_detectability=q10_minimum,
            q10_minimum_over_ascan=q10_minimum / baseline_minimum,
            nominal_nonlocal_coherence=nominal.nonlocal_coherence,
            q90_nonlocal_coherence=q90_coherence,
            nominal_stable_rank=nominal.stable_rank,
            q10_stable_rank=q10_rank,
            q10_calibration_fidelity=q10_fidelity,
            maximum_sigma_after=maximum(getproperty.(selected, :max_sigma_after)),
            gate_nominal_minimum=nominal.minimum_detectability >= baseline_minimum,
            gate_q10_minimum=q10_minimum >= 0.8 * baseline_minimum,
            gate_nominal_coherence=nominal.nonlocal_coherence <= 0.50,
            gate_q90_coherence=q90_coherence <= 0.60,
            gate_nominal_rank=nominal.stable_rank >= 10,
            gate_q10_rank=q10_rank >= 8,
            gate_q10_fidelity=q10_fidelity >= 0.90,
            gate_passivity=passivity,
            passed=nominal.minimum_detectability >= baseline_minimum &&
                q10_minimum >= 0.8 * baseline_minimum &&
                nominal.nonlocal_coherence <= 0.50 &&
                q90_coherence <= 0.60 &&
                nominal.stable_rank >= 10 &&
                q10_rank >= 8 &&
                q10_fidelity >= 0.90 &&
                passivity,
        ))
    end
    rows
end

function plot_nominal(nominal_rows, baseline_rows)
    selected = filter(row ->
        row.frequency_grid == :seven && row.covariance_model == :correlated,
        nominal_rows,
    )
    labels = String.(getproperty.(selected, :architecture_id))
    baseline = only(filter(row ->
        row.frequency_grid == :seven && row.covariance_model == :correlated,
        baseline_rows,
    ))
    detection = bar(
        labels,
        getproperty.(selected, :minimum_detectability);
        label=false,
        ylabel="minimum detectability",
        title="Causal/pulse nominal comparison",
        xrotation=25,
        size=(1100, 520),
    )
    hline!(detection, [baseline.minimum_detectability]; linestyle=:dash, linewidth=2,
           label="A-scan minimum")
    coherence = bar(
        labels,
        getproperty.(selected, :nonlocal_coherence);
        label=false,
        ylabel="nonlocal coherence",
        title="Ambiguity after pulse covariance",
        xrotation=25,
    )
    savefig(
        plot(detection, coherence; layout=(1, 2), size=(1450, 540), margin=7 * Plots.mm),
        joinpath(OUTPUT_ROOT, "o2c_nominal.png"),
    )
end

function plot_robustness(summary_rows)
    labels = String.(getproperty.(summary_rows, :architecture_id))
    detection = bar(
        labels,
        getproperty.(summary_rows, :q10_minimum_over_ascan);
        label=false,
        ylabel="Q10 minimum / A-scan minimum",
        title="Robust detectability",
        ylim=(0, max(1.2, 1.1 * maximum(getproperty.(summary_rows, :q10_minimum_over_ascan)))),
    )
    hline!(detection, [0.8]; linestyle=:dash, linewidth=2, label="gate")
    fidelity = bar(
        labels,
        getproperty.(summary_rows, :q10_calibration_fidelity);
        label=false,
        ylabel="Q10 calibration fidelity",
        title="Nominal-to-perturbed dictionary fidelity",
        ylim=(0, 1.02),
    )
    hline!(fidelity, [0.9]; linestyle=:dash, linewidth=2, label="gate")
    savefig(
        plot(detection, fidelity; layout=(1, 2), size=(1150, 470), margin=7 * Plots.mm),
        joinpath(OUTPUT_ROOT, "o2c_robustness.png"),
    )
end

function main()
    isfile(CODEBOOK_PATH) || error("run optimize_operator_codebooks.jl first")
    codebooks = JLD2.load(CODEBOOK_PATH)
    best_ports = codebooks["best_ports"]
    real_ports = best_ports[(4, :real_sector)]
    complex_ports = best_ports[(4, :complex_sector)]

    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    library = response_library(hypotheses, config)
    nominal_rows, nominal_dictionaries, passivity_rows = nominal_cases(
        library, hypotheses, config, real_ports, complex_ports,
    )
    baseline_rows = baseline_cases(library, hypotheses, config)
    ensemble_rows = robustness_ensemble(
        library, hypotheses, config, real_ports, complex_ports,
    )
    baseline = only(filter(row ->
        row.frequency_grid == :seven && row.covariance_model == :correlated,
        baseline_rows,
    ))
    summary_rows = robustness_summary(
        ensemble_rows, baseline.minimum_detectability, nominal_rows,
    )

    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_nominal.csv"),
        propertynames(first(nominal_rows)),
        nominal_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_ascan_baseline.csv"),
        propertynames(first(baseline_rows)),
        baseline_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_passivity.csv"),
        propertynames(first(passivity_rows)),
        passivity_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_ensemble.csv"),
        propertynames(first(ensemble_rows)),
        ensemble_rows,
    )
    write_csv(
        joinpath(OUTPUT_ROOT, "o2c_summary.csv"),
        propertynames(first(summary_rows)),
        summary_rows,
    )
    jldsave(
        joinpath(OUTPUT_ROOT, "o2c_causal_robustness.jld2");
        config,
        hypotheses,
        real_ports,
        complex_ports,
        nominal_rows,
        nominal_dictionaries,
        passivity_rows,
        baseline_rows,
        ensemble_rows,
        summary_rows,
        ensemble_size=ENSEMBLE_SIZE,
        random_seed=RANDOM_SEED,
    )
    plot_nominal(nominal_rows, baseline_rows)
    plot_robustness(summary_rows)

    println("[+] O2c A-scan correlated-pulse minimum=$(baseline.minimum_detectability)")
    for row in summary_rows
        println(
            "[+] $(row.architecture_id): Q10 min/A=$(round(row.q10_minimum_over_ascan; digits=3)), " *
            "Q90 mu=$(round(row.q90_nonlocal_coherence; digits=3)), " *
            "Q10 rank=$(round(row.q10_stable_rank; digits=3)), " *
            "Q10 fidelity=$(round(row.q10_calibration_fidelity; digits=3)), " *
            "passed=$(row.passed)",
        )
    end
    println("[+] wrote O2c results to $OUTPUT_ROOT")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
