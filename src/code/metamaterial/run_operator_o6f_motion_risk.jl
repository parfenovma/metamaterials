module RunOperatorO6FMotionRisk

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o3_motion_frontier.jl"))

export scan_features_at_positions,
       run_o6f_motion_risk

const O6F_MOTION_SOURCE = joinpath(
    PROJECT_ROOT, "results", "operator_space", "motion_frontier_v0",
    "motion_frontier.jld2",
)
const O6F_MOTION_OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6f_motion_risk_v0",
)
const O6F_MOTION_BASE_SEED = UInt64(0x06f0b117)
const O6F_PHYSICAL_SEED = O6F_MOTION_BASE_SEED + 0x101
const O6F_SIGNAL_NOISE_SEED = O6F_MOTION_BASE_SEED + 0x202
const O6F_H0_NOISE_SEED = O6F_MOTION_BASE_SEED + 0x303
const O6F_SCALES = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
const POSITION_RMS_LAMBDA = 0.05
const CONTACT_LOG_AMPLITUDE_RMS = 0.05
const CONTACT_DELAY_RMS_CYCLES = 0.02
const H0_TRIALS_PER_CONTEXT = 200
const PFA_GATE = 0.0125
const ALPHA90_RATIO_GATE = 1.10
const PDL_REFERENCE_GATE = 0.80
const ARCHITECTURE_IDS = (
    :motion_compact_center_p12_native,
    :motion_uniform_full_p17_a8,
)

function scan_architectures(config)
    [
        motion_architecture(
            :motion_compact_center_p12_native, collect(3:14), :compact_center,
        ),
        motion_architecture(
            :motion_uniform_full_p17_a8,
            collect(1:config.element_count),
            :uniform_full,
        ),
    ]
end

function pulse_whitening_transform(config)
    raw_weights = hann_burst_spectrum(config.frequency_ratios)
    weights = raw_weights / norm(raw_weights)
    lower = cholesky(Hermitian(
        electronics_frequency_covariance(config.frequency_ratios),
    )).L
    lower \ Diagonal(weights)
end

function monostatic_value(position, hypothesis, frequency, config)
    polarization = ComplexF64[1, 0]
    green = elastic_green_tensor(
        position,
        hypothesis.x_wavelengths,
        hypothesis.y_wavelengths,
        frequency,
        config,
    ) * polarization
    transpose(green) * hypothesis.scattering_tensor * green
end

"""Pulse-shaped, electronics-whitened monostatic blocks at arbitrary stops."""
function scan_features_at_positions(
    hypotheses,
    config,
    positions;
    log_amplitude_errors=zeros(length(positions)),
    delay_errors_cycles=zeros(length(positions)),
)
    stop_count = length(positions)
    length(log_amplitude_errors) == stop_count ||
        throw(DimensionMismatch("one log-amplitude error per stop is required"))
    length(delay_errors_cycles) == stop_count ||
        throw(DimensionMismatch("one delay error per stop is required"))
    frequency_count = length(config.frequency_ratios)
    raw = Matrix{ComplexF64}(undef,
        stop_count * frequency_count, length(hypotheses),
    )
    for (stop, position) in enumerate(positions)
        rows = ((stop - 1) * frequency_count + 1):(stop * frequency_count)
        for (frequency_index, frequency) in enumerate(config.frequency_ratios)
            contact = exp(log_amplitude_errors[stop]) * exp(
                -2pi * im * frequency * delay_errors_cycles[stop],
            )
            row = first(rows) + frequency_index - 1
            for (hypothesis_index, hypothesis) in enumerate(hypotheses)
                raw[row, hypothesis_index] = contact * monostatic_value(
                    position, hypothesis, frequency, config,
                )
            end
        end
    end
    transform = pulse_whitening_transform(config)
    features = similar(raw)
    for stop in 1:stop_count
        rows = ((stop - 1) * frequency_count + 1):(stop * frequency_count)
        features[rows, :] .= transform * view(raw, rows, :)
    end
    features
end

function linear_sample(grid, values, point)
    length(grid) == length(values) || throw(DimensionMismatch("grid mismatch"))
    length(grid) >= 2 || throw(ArgumentError("at least two grid points are required"))
    lower = if point <= first(grid)
        1
    elseif point >= last(grid)
        length(grid) - 1
    else
        min(searchsortedlast(grid, point), length(grid) - 1)
    end
    fraction = (point - grid[lower]) / (grid[lower + 1] - grid[lower])
    (1 - fraction) * values[lower] + fraction * values[lower + 1]
end

function ringing_features_at_positions(
    library,
    config,
    positions;
    log_amplitude_errors=zeros(length(positions)),
    delay_errors_cycles=zeros(length(positions)),
)
    native_positions = element_positions(config)
    frequency_count = length(config.frequency_ratios)
    atom_count = size(library, 4)
    raw = Matrix{ComplexF64}(undef,
        length(positions) * frequency_count, atom_count,
    )
    diagonals = Matrix{ComplexF64}(undef, length(native_positions), atom_count)
    for (frequency_index, frequency) in enumerate(config.frequency_ratios)
        for atom in 1:atom_count, index in eachindex(native_positions)
            diagonals[index, atom] = library[index, index, frequency_index, atom]
        end
        for (stop, position) in enumerate(positions)
            contact = exp(log_amplitude_errors[stop]) * exp(
                -2pi * im * frequency * delay_errors_cycles[stop],
            )
            row = (stop - 1) * frequency_count + frequency_index
            for atom in 1:atom_count
                raw[row, atom] = contact * linear_sample(
                    native_positions, view(diagonals, :, atom), position,
                )
            end
        end
    end
    transform = pulse_whitening_transform(config)
    features = similar(raw)
    for stop in eachindex(positions)
        rows = ((stop - 1) * frequency_count + 1):(stop * frequency_count)
        features[rows, :] .= transform * view(raw, rows, :)
    end
    features
end

function select_stops(full_features, indices, frequency_count)
    rows = reduce(vcat, [
        collect(((index - 1) * frequency_count + 1):(index * frequency_count))
        for index in indices
    ])
    Matrix(view(full_features, rows, :)) ./ sqrt(length(indices))
end

function build_frozen_models(
    architectures, config, hypotheses, split_libraries, scales, saved,
)
    nominal_library = response_library(hypotheses, config)
    models = Dict{Symbol,NamedTuple}()
    for architecture in architectures
        dictionary = architecture_dictionary(nominal_library, architecture, config)
        train_clutter = architecture_clutter(
            split_libraries[:train], architecture, config, scales,
        )
        test_clutter = architecture_clutter(
            split_libraries[:test], architecture, config, scales,
        )
        precision_applied = precision_dictionary(dictionary, train_clutter)
        normalization = max.(real.(diag(matched_mean(
            precision_applied, dictionary,
        ))), eps(Float64))
        threshold_index = findfirst(row ->
            row.architecture_id == architecture.id, saved["thresholds"])
        isnothing(threshold_index) && error("saved threshold is missing")
        threshold_row = saved["thresholds"][threshold_index]
        test_factor = covariance_factor(matched_output_covariance(
            precision_applied, test_clutter,
        ))
        models[architecture.id] = (;
            dictionary,
            precision_applied,
            normalization,
            threshold=threshold_row.threshold,
            saved_threshold_row=threshold_row,
            saved_architecture_index=threshold_index,
            nominal_test_factor=test_factor,
        )
    end
    models
end

function add_batch!(counts, architecture_id, scale, batch)
    for row in batch
        key = (architecture_id, scale, row.amplitude)
        previous = get(counts, key, (trials=0, detections=0, localized=0))
        counts[key] = (
            trials=previous.trials + row.trials,
            detections=previous.detections + row.detections,
            localized=previous.localized + row.localized,
        )
    end
end

function aggregate_curve_rows(counts, architectures, scales, amplitudes)
    rows = NamedTuple[]
    for architecture in architectures, scale in scales, amplitude in amplitudes
        count_row = counts[(architecture.id, scale, amplitude)]
        push!(rows, (
            architecture_id=architecture.id,
            scale,
            amplitude,
            trials=count_row.trials,
            detections=count_row.detections,
            localized=count_row.localized,
            pd=count_row.detections / count_row.trials,
            pdl=count_row.localized / count_row.trials,
        ))
    end
    rows
end

function curve_alpha90(rows, architecture_id, scale)
    selected = sort(filter(row ->
        row.architecture_id == architecture_id && row.scale == scale, rows,
    ); by=row -> row.amplitude)
    crossing_amplitude(
        getproperty.(selected, :amplitude), getproperty.(selected, :pd), 0.90,
    )
end

function curve_reference(rows, architecture_id, scale, reference_alpha)
    only(filter(row ->
        row.architecture_id == architecture_id && row.scale == scale &&
        isapprox(row.amplitude, reference_alpha; rtol=1e-12), rows,
    ))
end

function nominal_replay(
    architectures, models, hypotheses, config, amplitudes, reference_alpha, saved,
)
    counts = Dict{Tuple{Symbol,Float64,Float64},NamedTuple}()
    contexts = NamedTuple[]
    rng = Xoshiro(MOTION_SPLIT_SEEDS[:signal_test])
    for context_index in 1:MOTION_TEST_CONTEXTS
        context = fresh_signal_context(
            rng, hypotheses, config, 100_000 + context_index,
        )
        library = response_library(context.hypotheses, context.config)
        latent = complex_noise(rng, length(hypotheses), length(hypotheses))
        push!(contexts, context)
        for architecture in architectures
            model = models[architecture.id]
            dictionary = architecture_dictionary(library, architecture, context.config)
            batch = evaluate_with_common_noise(
                model.nominal_test_factor * latent,
                matched_mean(model.precision_applied, dictionary),
                model.normalization,
                model.threshold,
                amplitudes,
                hypotheses,
                context.hypotheses,
            )
            add_batch!(counts, architecture.id, 0.0, batch)
        end
    end
    curves = aggregate_curve_rows(counts, architectures, (0.0,), amplitudes)
    rows = NamedTuple[]
    for architecture in architectures
        id = architecture.id
        alpha90 = curve_alpha90(curves, id, 0.0)
        reference = curve_reference(curves, id, 0.0, reference_alpha)
        saved_summary = only(filter(row -> row.architecture_id == id,
            saved["summary_rows"]))
        model = models[id]
        h0_rng = Xoshiro(
            MOTION_SPLIT_SEEDS[:h0_test] +
            UInt64(10_000 * model.saved_architecture_index),
        )
        h0_scores = no_defect_scores(
            h0_rng, model.nominal_test_factor, model.normalization,
            MOTION_TEST_H0_TRIALS,
        )
        false_alarms = count(>(model.threshold), h0_scores)
        push!(rows, (
            architecture_id=id,
            saved_alpha90=saved_summary.alpha90,
            replay_alpha90=alpha90,
            alpha90_relative_error=abs(alpha90 - saved_summary.alpha90) /
                saved_summary.alpha90,
            saved_pd_reference=saved_summary.pd_reference,
            replay_pd_reference=reference.pd,
            saved_pdl_reference=saved_summary.pdl_reference,
            replay_pdl_reference=reference.pdl,
            saved_false_alarms=model.saved_threshold_row.false_alarms,
            replay_false_alarms=false_alarms,
            exact_count_reproduction=
                reference.detections == round(Int,
                    saved_summary.pd_reference * reference.trials) &&
                reference.localized == round(Int,
                    saved_summary.pdl_reference * reference.trials) &&
                false_alarms == model.saved_threshold_row.false_alarms,
        ))
    end
    all(row -> row.exact_count_reproduction && row.alpha90_relative_error <= 1e-12,
        rows) || error("frozen P12/P17 detector replay does not reproduce baseline")
    contexts, rows
end

function physical_realization_features(
    context, config, test_hypotheses, ringing,
    position_standard, amplitude_standard, delay_standard, scale,
)
    nominal_positions = element_positions(config)
    positions = nominal_positions .+
        scale * POSITION_RMS_LAMBDA .* position_standard
    log_amplitude = scale * CONTACT_LOG_AMPLITUDE_RMS .* amplitude_standard
    delay = scale * CONTACT_DELAY_RMS_CYCLES .* delay_standard
    signal = scan_features_at_positions(
        context.hypotheses, context.config, positions;
        log_amplitude_errors=log_amplitude,
        delay_errors_cycles=delay,
    )
    material = scan_features_at_positions(
        test_hypotheses.material, config, positions;
        log_amplitude_errors=log_amplitude,
        delay_errors_cycles=delay,
    )
    backwall = scan_features_at_positions(
        test_hypotheses.backwall, config, positions;
        log_amplitude_errors=log_amplitude,
        delay_errors_cycles=delay,
    )
    ringing_features = ringing_features_at_positions(
        ringing, config, positions;
        log_amplitude_errors=log_amplitude,
        delay_errors_cycles=delay,
    )
    (; signal, material, backwall, ringing=ringing_features)
end

function run_physical_campaign(
    architectures, models, contexts, hypotheses, config,
    split_libraries, scales, amplitudes, reference_alpha,
)
    physical_rng = Xoshiro(O6F_PHYSICAL_SEED)
    position_standard = randn(
        physical_rng, config.element_count, MOTION_TEST_CONTEXTS,
    )
    amplitude_standard = randn(
        physical_rng, config.element_count, MOTION_TEST_CONTEXTS,
    )
    delay_standard = randn(
        physical_rng, config.element_count, MOTION_TEST_CONTEXTS,
    )
    signal_noise_rng = Xoshiro(O6F_SIGNAL_NOISE_SEED)
    h0_noise_rng = Xoshiro(O6F_H0_NOISE_SEED)
    counts = Dict{Tuple{Symbol,Float64,Float64},NamedTuple}()
    false_alarms = Dict((architecture.id, scale) => 0
        for architecture in architectures for scale in O6F_SCALES)
    h0_trials = Dict((architecture.id, scale) => 0
        for architecture in architectures for scale in O6F_SCALES)
    test_hypotheses = (
        material=material_clutter_hypotheses(MOTION_SPLIT_SEEDS[:nuisance_test]),
        backwall=backwall_hypotheses(MOTION_SPLIT_SEEDS[:nuisance_test]),
    )
    ringing = split_libraries[:test][:ringing]
    frequency_count = length(config.frequency_ratios)
    for (context_index, context) in enumerate(contexts)
        signal_latent = complex_noise(
            signal_noise_rng, length(hypotheses), length(hypotheses),
        )
        h0_latent = complex_noise(
            h0_noise_rng, length(hypotheses), H0_TRIALS_PER_CONTEXT,
        )
        for scale in O6F_SCALES
            full = physical_realization_features(
                context,
                config,
                test_hypotheses,
                ringing,
                view(position_standard, :, context_index),
                view(amplitude_standard, :, context_index),
                view(delay_standard, :, context_index),
                scale,
            )
            for architecture in architectures
                id = architecture.id
                model = models[id]
                indices = architecture.mechanical_indices
                signal_dictionary = select_stops(
                    full.signal, indices, frequency_count,
                )
                clutter = hcat(
                    scales[:backwall] .* select_stops(
                        full.backwall, indices, frequency_count,
                    ),
                    scales[:ringing] .* select_stops(
                        full.ringing, indices, frequency_count,
                    ),
                    scales[:material] .* select_stops(
                        full.material, indices, frequency_count,
                    ),
                )
                factor = covariance_factor(matched_output_covariance(
                    model.precision_applied, clutter,
                ))
                batch = evaluate_with_common_noise(
                    factor * signal_latent,
                    matched_mean(model.precision_applied, signal_dictionary),
                    model.normalization,
                    model.threshold,
                    amplitudes,
                    hypotheses,
                    context.hypotheses,
                )
                add_batch!(counts, id, scale, batch)
                scores, _ = normalized_scores(
                    factor * h0_latent, model.normalization,
                )
                false_alarms[(id, scale)] += count(>(model.threshold), scores)
                h0_trials[(id, scale)] += length(scores)
            end
        end
        context_index % 5 == 0 && println(
            "[+] O6f-B physical contexts: $context_index / $(length(contexts))",
        )
    end
    curves = aggregate_curve_rows(counts, architectures, O6F_SCALES, amplitudes)
    baseline_alpha = Dict(architecture.id =>
        curve_alpha90(curves, architecture.id, 0.0)
        for architecture in architectures)
    summary = NamedTuple[]
    for architecture in architectures, scale in O6F_SCALES
        id = architecture.id
        alpha90 = curve_alpha90(curves, id, scale)
        reference = curve_reference(curves, id, scale, reference_alpha)
        failures = false_alarms[(id, scale)]
        trials = h0_trials[(id, scale)]
        pfa_upper = wilson_upper_bound(failures, trials)
        ratio = alpha90 / baseline_alpha[id]
        gates = (
            pfa=pfa_upper <= PFA_GATE,
            alpha90=isfinite(ratio) && ratio <= ALPHA90_RATIO_GATE,
            pdl=reference.pdl >= PDL_REFERENCE_GATE,
        )
        push!(summary, (
            architecture_id=id,
            scale,
            position_rms_lambda=scale * POSITION_RMS_LAMBDA,
            contact_log_amplitude_rms=scale * CONTACT_LOG_AMPLITUDE_RMS,
            contact_delay_rms_cycles=scale * CONTACT_DELAY_RMS_CYCLES,
            h0_trials=trials,
            false_alarms=failures,
            observed_pfa=failures / trials,
            pfa_upper_95=pfa_upper,
            alpha90,
            nominal_replay_alpha90=baseline_alpha[id],
            alpha90_ratio_to_nominal=ratio,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            pfa_gate=gates.pfa,
            alpha90_gate=gates.alpha90,
            pdl_gate=gates.pdl,
            joint_pass=all(values(gates)),
        ))
    end
    curves, summary, position_standard, amplitude_standard, delay_standard
end

function write_rows(path, rows)
    write_csv(path, propertynames(first(rows)), rows)
end

function self_check_arbitrary_scan(architectures, config, hypotheses)
    library = response_library(hypotheses, config)
    full = scan_features_at_positions(hypotheses, config, element_positions(config))
    for architecture in architectures
        expected = architecture_dictionary(library, architecture, config)
        actual = select_stops(
            full, architecture.mechanical_indices, length(config.frequency_ratios),
        )
        relative_error = norm(actual - expected) / norm(expected)
        relative_error <= 2e-12 || error(
            "arbitrary-position scan disagrees with frozen compressor: $relative_error",
        )
    end
    true
end

function plot_physical_frontier(summary)
    plot_object = plot(
        xlabel="physical-error scale",
        ylabel="alpha90 / nominal",
        title="O6f-B frozen motion physical-risk frontier",
        legend=:topleft,
    )
    for (id, marker, color) in zip(
        ARCHITECTURE_IDS, (:circle, :diamond), (:steelblue, :darkorange),
    )
        selected = sort(filter(row -> row.architecture_id == id, summary);
            by=row -> row.scale)
        plot!(
            plot_object,
            getproperty.(selected, :scale),
            getproperty.(selected, :alpha90_ratio_to_nominal);
            marker,
            color,
            linewidth=2,
            label=String(id),
        )
    end
    hline!(plot_object, [ALPHA90_RATIO_GATE];
        color=:black, linestyle=:dash, label="1.10 gate")
    savefig(plot_object, joinpath(
        O6F_MOTION_OUTPUT_ROOT, "motion_physical_risk_frontier.png",
    ))
end

function run_o6f_motion_risk()
    isfile(O6F_MOTION_SOURCE) || error("motion source result is missing")
    mkpath(O6F_MOTION_OUTPUT_ROOT)
    saved = load(O6F_MOTION_SOURCE)
    config = SurrogateConfig()
    hypotheses = localization_hypotheses()
    architectures = scan_architectures(config)
    self_check_arbitrary_scan(architectures, config, hypotheses)
    split_libraries = Dict(
        split => fresh_nuisance_libraries(config, split)
        for split in (:train, :test)
    )
    oracle = only(filter(architecture -> architecture.kind == :fmc,
        fixed_architectures(config)))
    scales, scale_rows = fresh_nuisance_scales(
        config, split_libraries, oracle,
    )
    for row in scale_rows
        saved_scale = only(filter(other -> other.category == row.category,
            saved["scale_rows"]))
        isapprox(row.physical_scale, saved_scale.physical_scale; rtol=1e-13) ||
            error("nuisance scale reconstruction failed")
    end
    models = build_frozen_models(
        architectures, config, hypotheses, split_libraries, scales, saved,
    )
    amplitudes = Float64.(saved["amplitudes"])
    reference_alpha = Float64(saved["reference_alpha"])
    println("[>] O6f-B exact frozen-detector replay")
    contexts, replay_rows = nominal_replay(
        architectures, models, hypotheses, config, amplitudes,
        reference_alpha, saved,
    )
    println("[+] exact P12/P17 signal and H0 counts reproduced")
    println("[>] O6f-B paired encoder/contact campaign")
    curves, summary, position_standard, amplitude_standard, delay_standard =
        run_physical_campaign(
            architectures,
            models,
            contexts,
            hypotheses,
            config,
            split_libraries,
            scales,
            amplitudes,
            reference_alpha,
        )
    threshold_rows = NamedTuple[]
    for architecture in architectures
        selected = filter(row ->
            row.architecture_id == architecture.id && row.joint_pass, summary)
        maximum_scale = isempty(selected) ? NaN : maximum(getproperty.(selected, :scale))
        push!(threshold_rows, (
            architecture_id=architecture.id,
            maximum_tested_passing_scale=maximum_scale,
            position_rms_lambda=maximum_scale * POSITION_RMS_LAMBDA,
            contact_log_amplitude_rms=maximum_scale * CONTACT_LOG_AMPLITUDE_RMS,
            contact_delay_rms_cycles=maximum_scale * CONTACT_DELAY_RMS_CYCLES,
            all_tested_scales_pass=isapprox(maximum_scale, maximum(O6F_SCALES)),
        ))
    end
    write_rows(joinpath(O6F_MOTION_OUTPUT_ROOT, "baseline_replay.csv"), replay_rows)
    write_rows(joinpath(O6F_MOTION_OUTPUT_ROOT, "physical_detection_curves.csv"), curves)
    write_rows(joinpath(O6F_MOTION_OUTPUT_ROOT, "physical_risk_summary.csv"), summary)
    write_rows(joinpath(O6F_MOTION_OUTPUT_ROOT, "physical_risk_thresholds.csv"), threshold_rows)
    seed_rows = [
        (role="physical errors", seed=string("0x", string(O6F_PHYSICAL_SEED; base=16))),
        (role="signal noise", seed=string("0x", string(O6F_SIGNAL_NOISE_SEED; base=16))),
        (role="H0 noise", seed=string("0x", string(O6F_H0_NOISE_SEED; base=16))),
    ]
    write_rows(joinpath(O6F_MOTION_OUTPUT_ROOT, "seed_contract.csv"), seed_rows)
    jldsave(
        joinpath(O6F_MOTION_OUTPUT_ROOT, "o6f_motion_risk.jld2");
        format_version=1,
        source_motion_path=abspath(O6F_MOTION_SOURCE),
        config,
        hypotheses,
        architecture_ids=ARCHITECTURE_IDS,
        amplitudes,
        reference_alpha,
        scale_rows,
        replay_rows,
        curves,
        summary,
        threshold_rows,
        physical_seed=O6F_PHYSICAL_SEED,
        signal_noise_seed=O6F_SIGNAL_NOISE_SEED,
        h0_noise_seed=O6F_H0_NOISE_SEED,
        h0_trials_per_context=H0_TRIALS_PER_CONTEXT,
        position_standard,
        amplitude_standard,
        delay_standard,
    )
    open(joinpath(O6F_MOTION_OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6f-B frozen motion physical-risk replay")
        println(io)
        println(io, "The saved P12/P17 GLRT thresholds, nominal templates and training nuisance model were reconstructed without retraining. Signal, localization and H0 false-alarm counts reproduce the saved baseline exactly.")
        println(io)
        println(io, "Each row uses 60 frozen signal contexts and 12,000 fresh aggregate H0 trials. Encoder, contact-amplitude and contact-delay errors are paired across architectures and nested across scales.")
        println(io)
        for threshold in threshold_rows
            @printf(io,
                "- `%s`: maximum tested joint-pass scale `%.2f` = position RMS `%.4f lambda`, contact log-amplitude RMS `%.2f%%`, contact delay RMS `%.4f cycles`; all tested scales pass: `%s`.\n",
                String(threshold.architecture_id),
                threshold.maximum_tested_passing_scale,
                threshold.position_rms_lambda,
                100threshold.contact_log_amplitude_rms,
                threshold.contact_delay_rms_cycles,
                threshold.all_tested_scales_pass,
            )
        end
        println(io)
        println(io, "This remains a reduced Born replay: it includes encoder and multiplicative/delay contact errors, but not wear, lift-off-dependent mode conversion, surface waves or scan-time economics.")
    end
    plot_physical_frontier(summary)
    for row in threshold_rows
        println("[+] $(row.architecture_id) maximum tested pass scale=$(row.maximum_tested_passing_scale)")
    end
    println("[+] wrote O6f-B to $O6F_MOTION_OUTPUT_ROOT")
    (; replay_rows, summary, threshold_rows)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6f_motion_risk()

end
