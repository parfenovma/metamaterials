module RunOperatorO8DetectorReplay

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o7_global_states.jl"))

export run_o8_detector_replay

const O8_OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o8_plane_strain_v0",
)
const O8_PZT_BUNDLE = joinpath(O8_OUTPUT_ROOT, "pzt_physical_bundle.jld2")
const O8_BRANCH_SUMMARY = joinpath(O8_OUTPUT_ROOT, "branch_summary.csv")
const O8_BASE_SEED = UInt(0x0808a117)
const O8_SEEDS = (
    nuisance_train=O8_BASE_SEED + UInt(0x101),
    nuisance_validation=O8_BASE_SEED + UInt(0x202),
    nuisance_test=O8_BASE_SEED + UInt(0x303),
    signal_validation=O8_BASE_SEED + UInt(0x404),
    signal_test=O8_BASE_SEED + UInt(0x505),
    h0_validation=O8_BASE_SEED + UInt(0x606),
    h0_test=O8_BASE_SEED + UInt(0x707),
)
const O8_VALIDATION_H0_TRIALS = 10_000
const O8_TEST_H0_TRIALS = 10_000
const O8_VALIDATION_CONTEXTS = 24
const O8_TEST_CONTEXTS = 60

o8_noise(rng, rows, columns) =
    (randn(rng, rows, columns) .+ im .* randn(rng, rows, columns)) ./ sqrt(2)

function o8_signal_contexts(config, hypotheses, split, count)
    seed = split == :validation ? O8_SEEDS.signal_validation : O8_SEEDS.signal_test
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

function o8_nuisance_contract(config)
    libraries = Dict(
        :train => sf_nuisance_libraries(config, O8_SEEDS.nuisance_train),
        :validation => sf_nuisance_libraries(config, O8_SEEDS.nuisance_validation),
        :test => sf_nuisance_libraries(config, O8_SEEDS.nuisance_test),
    )
    fmc = sf_architecture(:o8_scale_fmc, :fmc, 17, 1, :full_aperture, (), :scale)
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

function o8_clutter(architecture, libraries, scales, config)
    hcat((
        scales[category] .* o7_compress_library(
            libraries[category], architecture, config,
        ) for category in (:backwall, :ringing, :material)
    )...)
end

function o8_detection_model(
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
    validation_rng = Xoshiro(O8_SEEDS.h0_validation + UInt(10_000 * architecture_index))
    test_rng = Xoshiro(O8_SEEDS.h0_test + UInt(10_000 * architecture_index))
    validation_scores, _ = normalized_scores(
        validation_factor * o8_noise(
            validation_rng, size(validation_factor, 1), O8_VALIDATION_H0_TRIALS,
        ), normalization,
    )
    threshold = quantile(validation_scores, 0.99)
    test_scores, _ = normalized_scores(
        test_factor * o8_noise(
            test_rng, size(test_factor, 1), O8_TEST_H0_TRIALS,
        ), normalization,
    )
    false_alarms = count(>(threshold), test_scores)
    (
        precision_applied, normalization, validation_factor, test_factor,
        threshold, false_alarms,
        observed_pfa=false_alarms / O8_TEST_H0_TRIALS,
        pfa_upper_95=wilson_upper_bound(false_alarms, O8_TEST_H0_TRIALS),
    )
end

function o8_build_models(architectures, config, hypotheses, libraries, scales)
    nominal_library = response_library(hypotheses, config)
    dictionaries = Dict{Symbol,Matrix{ComplexF64}}()
    models = Dict{Symbol,NamedTuple}()
    rows = NamedTuple[]
    for (index, architecture) in enumerate(architectures)
        dictionary = o7_compress_library(nominal_library, architecture, config)
        dictionaries[architecture.id] = dictionary
        clutters = Dict(split => o8_clutter(
            architecture, libraries[split], scales, config,
        ) for split in (:train, :validation, :test))
        model = o8_detection_model(
            architecture, dictionary, clutters[:train], clutters[:validation],
            clutters[:test], index,
        )
        models[architecture.id] = model
        push!(rows, (
            architecture_id=architecture.id,
            validation_trials=O8_VALIDATION_H0_TRIALS,
            test_trials=O8_TEST_H0_TRIALS, threshold=model.threshold,
            false_alarms=model.false_alarms, observed_pfa=model.observed_pfa,
            pfa_upper_95=model.pfa_upper_95,
        ))
        println(
            "[+] O8 detector $(index)/$(length(architectures)): " *
            "$(architecture.id), Pfa95=$(round(model.pfa_upper_95; digits=4))",
        )
    end
    dictionaries, models, rows
end

function o8_reference_amplitudes(oracle, model, dictionary, contexts, hypotheses)
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
        true_dictionary = o7_compress_library(context.library, oracle, context.config)
        mean_matrix = matched_mean(model.precision_applied, true_dictionary)
        rng = Xoshiro(O8_SEEDS.signal_validation + UInt(context.context_index))
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
    isnothing(reference_index) && error("O8 FMC validation Pd did not reach 0.95")
    amplitudes, probabilities, amplitudes[reference_index], anchor
end

function o8_test_campaign(architectures, models, contexts, hypotheses, amplitudes)
    counts = Dict{Tuple{Symbol,Float64},NamedTuple}()
    errors = Dict{Tuple{Symbol,Float64},Vector{Float64}}()
    for architecture in architectures, amplitude in amplitudes
        counts[(architecture.id, amplitude)] = (trials=0, detections=0, localized=0)
        errors[(architecture.id, amplitude)] = Float64[]
    end
    for context in contexts
        for (architecture_index, architecture) in enumerate(architectures)
            model = models[architecture.id]
            dictionary = o7_compress_library(
                context.library, architecture, context.config;
                state_errors=o7_context_errors(architecture, context),
            )
            mean_matrix = matched_mean(model.precision_applied, dictionary)
            rng = Xoshiro(
                O8_SEEDS.signal_test + UInt(10_000 * context.context_index + architecture_index),
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
            end
        end
        context.context_index % 10 == 0 && println(
            "[+] O8 detector contexts $(context.context_index)/$(length(contexts))",
        )
    end
    rows = NamedTuple[]
    for architecture in architectures, amplitude in amplitudes
        key = (architecture.id, amplitude)
        value = counts[key]
        error = errors[key]
        push!(rows, (
            architecture_id=architecture.id, nuisance_level=:C3, amplitude,
            trials=value.trials, detections=value.detections,
            localized=value.localized, pd=value.detections / value.trials,
            pdl=value.localized / value.trials,
            median_localization_error=isempty(error) ? NaN : median(error),
            q90_localization_error=isempty(error) ? NaN : quantile(error, 0.90),
        ))
    end
    rows
end

function o8_summary(architectures, threshold_rows, curve_rows, reference_alpha)
    thresholds = Dict(row.architecture_id => row for row in threshold_rows)
    alpha = Dict{Symbol,Float64}()
    base = NamedTuple[]
    for architecture in architectures
        selected = filter(row -> row.architecture_id == architecture.id, curve_rows)
        alpha[architecture.id] = sf_alpha90(selected)
        reference = only(filter(row -> isapprox(
            row.amplitude, reference_alpha; rtol=1.0e-12,
        ), selected))
        intrinsic = thresholds[architecture.id].pfa_upper_95 <= 0.012 &&
            reference.pd >= 0.90 && reference.pdl >= 0.85 &&
            isfinite(alpha[architecture.id])
        push!(base, (
            architecture_id=architecture.id, role=architecture.role,
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
        gain = row.architecture_id == :o8_pzt_k3q4 ?
            1 - row.alpha90 / alpha[:o8_pzt_k3q1] : NaN
        final = row.intrinsic_gate_passed &&
            (row.architecture_id != :o8_pzt_k3q4 || gain >= 0.05)
        merge(row, (; switching_gain_vs_q1=gain, final_gate_passed=final))
    end for row in base]
end

function append_readme(summary_rows)
    pzt = only(filter(row -> row.architecture_id == :o8_pzt_k3q4, summary_rows))
    comparator = only(filter(row -> row.architecture_id == :o8_pzt_k3q1, summary_rows))
    open(joinpath(O8_OUTPUT_ROOT, "README.md"), "a") do io
        println(io)
        println(io, "## Fresh detector replay")
        println(io)
        println(io, "| Architecture | alpha90 | Pfa95 | Pd@ref | Pdl@ref | final gate |")
        println(io, "|---|---:|---:|---:|---:|---|")
        for row in summary_rows
            println(io, "| $(row.architecture_id) | $(round(row.alpha90; digits=3)) | " *
                "$(round(100 * row.pfa_upper_95; digits=3))% | $(round(row.pd_reference; digits=3)) | " *
                "$(round(row.pdl_reference; digits=3)) | $(row.final_gate_passed) |")
        end
        println(io)
        println(io, "PZT Q4 alpha90 gain over the same FEM head in q00 only: " *
            "$(round(100 * (1 - pzt.alpha90 / comparator.alpha90); digits=2))%.")
        println(io)
        println(io, pzt.final_gate_passed ?
            "Verdict: the frozen PZT boundary passes O8 and may enter an explicit patch/adhesive piezoelectric FEM." :
            "Verdict: the frozen PZT boundary fails O8; stop before domain-resolved piezoelectric FEM.")
    end
end

function plot_curves(curves, summary)
    figure = plot(
        xscale=:log10, xlabel="signal amplitude", ylabel="Pd",
        title="O8 fresh detector replay", legend=:bottomright, size=(980, 650),
    )
    for row in summary
        selected = sort(filter(curve -> curve.architecture_id == row.architecture_id, curves);
            by=curve -> curve.amplitude)
        plot!(figure, getproperty.(selected, :amplitude), getproperty.(selected, :pd);
            marker=:circle, label=String(row.architecture_id))
    end
    hline!(figure, [0.90]; color=:black, linestyle=:dash, label="Pd=0.90")
    savefig(figure, joinpath(O8_OUTPUT_ROOT, "detector_replay.png"))
end

function run_o8_detector_replay()
    isfile(O8_PZT_BUNDLE) || error("missing O8 PZT physical bundle")
    physical = JLD2.load(O8_PZT_BUNDLE)
    physical["passed"] || error("O8 PZT failed physical gates; detector replay forbidden")
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    frequency_states = physical["frequency_states"]
    target_states = JLD2.load(O7_STATIC_BUNDLE, "architecture_states")[O7_TARGET_ID]
    pzt_q1 = o7_architecture(
        :o8_pzt_k3q1, :frequency_only, [frequency_states[1]],
        :matched_comparator; mechanism=:pzt,
    )
    pzt_q4 = o7_architecture(
        :o8_pzt_k3q4, :pzt_plane_strain, frequency_states,
        :frozen_physical; mechanism=:pzt,
    )
    ideal = sf_architecture(
        :o8_ideal_k3q4, :static, 3, 4, :reflection_balanced,
        target_states, :operator_target,
    )
    a_scan = sf_architecture(:o8_a_scan, :a_scan, 1, 1, :uniform, (), :baseline)
    fmc = sf_architecture(:o8_fmc, :fmc, 17, 1, :full_aperture, (), :oracle)
    architectures = [a_scan, pzt_q1, pzt_q4, ideal, fmc]

    libraries, scales, scale_rows = o8_nuisance_contract(config)
    validation_contexts = o8_signal_contexts(
        config, hypotheses, :validation, O8_VALIDATION_CONTEXTS,
    )
    test_contexts = o8_signal_contexts(config, hypotheses, :test, O8_TEST_CONTEXTS)
    println("[+] O8 fresh detector splits assembled after physical freeze")
    dictionaries, models, threshold_rows = o8_build_models(
        architectures, config, hypotheses, libraries, scales,
    )
    amplitudes, oracle_validation_pd, reference_alpha, amplitude_anchor =
        o8_reference_amplitudes(
            fmc, models[fmc.id], dictionaries[fmc.id], validation_contexts, hypotheses,
        )
    println("[+] O8 final detector test opened; FMC reference alpha=$reference_alpha")
    curves = o8_test_campaign(
        architectures, models, test_contexts, hypotheses, amplitudes,
    )
    summary = o8_summary(architectures, threshold_rows, curves, reference_alpha)
    seed_rows = [(purpose=key, seed=string(value)) for (key, value) in pairs(O8_SEEDS)]
    write_csv(joinpath(O8_OUTPUT_ROOT, "detector_seed_contract.csv"),
        propertynames(first(seed_rows)), seed_rows)
    write_csv(joinpath(O8_OUTPUT_ROOT, "detector_nuisance_scales.csv"),
        propertynames(first(scale_rows)), scale_rows)
    write_csv(joinpath(O8_OUTPUT_ROOT, "detector_thresholds.csv"),
        propertynames(first(threshold_rows)), threshold_rows)
    write_csv(joinpath(O8_OUTPUT_ROOT, "detector_curves.csv"),
        propertynames(first(curves)), curves)
    write_csv(joinpath(O8_OUTPUT_ROOT, "detector_summary.csv"),
        propertynames(first(summary)), summary)
    jldsave(
        joinpath(O8_OUTPUT_ROOT, "o8_detector_replay.jld2"); config, hypotheses,
        base_seed=O8_BASE_SEED, seed_contract=O8_SEEDS, threshold_rows,
        amplitudes, amplitude_anchor, oracle_validation_pd, reference_alpha,
        curves, summary,
    )
    append_readme(summary)
    plot_curves(curves, summary)
    for row in summary
        println(
            "[+] O8 $(row.architecture_id): alpha90=$(round(row.alpha90; digits=2)), " *
            "Pfa95=$(round(row.pfa_upper_95; digits=4)), pass=$(row.final_gate_passed)",
        )
    end
    summary
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_o8_detector_replay()
end

end

