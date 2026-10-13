module RunOperatorO3F1CStochastic

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o3_fullwave_carrier.jl"))
include(joinpath(@__DIR__, "run_operator_o3_f1c_scenes.jl"))
include(joinpath(@__DIR__, "operator_o3_statistics.jl"))

using .RunOperatorO3FullWaveCarrier
using .RunOperatorO3F1CScenes
using .OperatorO3Statistics: wilson_upper_bound

export run_f1c_campaign

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_O3_F1C_OUTPUT",
    joinpath(PROJECT_ROOT, "results", "operator_space", "fullwave_f1c_stochastic_v0"),
)
const BASE_SEED = UInt64(0x051c2026)
const SEEDS = (
    train=BASE_SEED + 0x1001,
    validation_h0=BASE_SEED + 0x2002,
    test_h0=BASE_SEED + 0x3003,
    validation_signal=BASE_SEED + 0x4004,
    test_signal=BASE_SEED + 0x5005,
    validation_noise=BASE_SEED + 0x6006,
    test_noise=BASE_SEED + 0x7007,
)
const TARGET_Y_MM = RunOperatorO3F1CScenes.TARGET_Y_MM
const TARGET_IDS = [
    Symbol("target_y$(replace(string(y), '.' => 'p', '-' => 'm'))")
    for y in TARGET_Y_MM
]
const AMPLITUDES = [
    0.20, 0.30, 0.45, 0.65, 0.90, 1.20, 1.60, 2.20, 3.00, 4.20,
    6.00, 8.50, 12.0, 18.0, 27.0, 40.0, 60.0, 90.0, 135.0,
]
const TRAIN_H0_REPLICATES = 32
const TRAIN_SIGNAL_REPLICATES = 16
const VALIDATION_H0_REPLICATES = 256
const TEST_H0_REPLICATES = 512
const VALIDATION_SIGNAL_REPLICATES = 48
const TEST_SIGNAL_REPLICATES = 96
const NOISE_DRAWS_PER_SIGNAL_CONTEXT = 4
const H0_VALIDATION_TRIALS = 10_000
const H0_TEST_TRIALS = 10_000
const MAX_NUISANCE_RANK = 24
const NUISANCE_ENERGY_TARGET = 0.99
const PFA_TARGET = 0.01

architecture_offset(id) = UInt64(sum(Int, codeunits(string(id)))) * UInt64(0x10001)
positive_floor(value::Real) = max(Float64(value), floatmin(Float64))
train_context_ids() = [:train_1, :train_2]

function project_time_response(response, definition)
    time_count = size(response, 3)
    first_projection = RunOperatorO3FullWaveCarrier.project_operator(
        view(response, :, :, 1), definition,
    )
    projected = Matrix{Float64}(undef, length(first_projection), time_count)
    projected[:, 1] .= real.(first_projection)
    for time_index in 2:time_count
        projected[:, time_index] .= real.(
            RunOperatorO3FullWaveCarrier.project_operator(
                view(response, :, :, time_index), definition,
            ),
        )
    end
    projected
end

function load_architecture_scenes(definition, context_ids)
    scenes = Dict{Symbol,Any}()
    common_time = nothing
    for context_id in context_ids
        h0_data = JLD2.load(RunOperatorO3F1CScenes.response_path(context_id, :h0))
        time_s = Float64.(h0_data["time_s"])
        isnothing(common_time) ? (common_time = time_s) :
            common_time == time_s || error("F1c contexts use different time grids")
        h0 = project_time_response(Float64.(h0_data["normal_response_m2"]), definition)
        deltas = Matrix{Float64}[]
        for target_id in TARGET_IDS
            h1_data = JLD2.load(RunOperatorO3F1CScenes.response_path(context_id, target_id))
            Float64.(h1_data["time_s"]) == time_s || error("paired F1c time grids differ")
            h1 = project_time_response(Float64.(h1_data["normal_response_m2"]), definition)
            push!(deltas, h1 - h0)
        end
        scenes[context_id] = (; h0, deltas)
    end
    selected = [
        index for (index, time) in enumerate(common_time)
        if time >= 8.0e-6 && (index - 1) % 3 == 0
    ]
    (; time_s=common_time, selected, scenes)
end

function nuisance_specs(seed, context_ids, replicates; with_targets=false)
    specs = NamedTuple[]
    targets = with_targets ? eachindex(TARGET_IDS) : (0:0)
    for (context_index, context_id) in enumerate(context_ids), target_index in targets,
        replicate in 1:replicates
        stream = UInt64(100_000context_index + 1_000target_index + replicate)
        rng = Xoshiro(seed + stream)
        push!(specs, (
            context_id,
            target_index,
            global_gain=exp(0.05randn(rng) - 0.5 * 0.05^2),
            delay_s=0.20e-6 * randn(rng),
            time_stretch=0.005 * randn(rng),
            attenuation_slope=0.10 * randn(rng),
            ringing_scale=0.02 * randn(rng),
            gain_seed=rand(rng, UInt64),
        ))
    end
    specs
end

function sample_signal(signal, time_s, source_time)
    source_time <= first(time_s) && return copy(view(signal, :, 1))
    source_time >= last(time_s) && return zeros(Float64, size(signal, 1))
    upper = searchsortedfirst(time_s, source_time)
    lower = upper - 1
    weight = (source_time - time_s[lower]) / (time_s[upper] - time_s[lower])
    (1 - weight) .* view(signal, :, lower) .+ weight .* view(signal, :, upper)
end

function transform_context(scene, delta, time_s, selected, spec, architecture_id)
    feature_count = size(scene.h0, 1)
    selected_time = time_s[selected]
    core_h0 = Matrix{Float64}(undef, feature_count, length(selected))
    core_delta = similar(core_h0)
    ringing = similar(core_h0)
    for (column, time) in enumerate(selected_time)
        source_time = (time - spec.delay_s) * (1 + spec.time_stretch)
        ring_time = (time - spec.delay_s - 1.0e-6) * (1 + spec.time_stretch)
        core_h0[:, column] .= sample_signal(scene.h0, time_s, source_time)
        core_delta[:, column] .= sample_signal(delta, time_s, source_time)
        ringing[:, column] .= sample_signal(scene.h0, time_s, ring_time)
    end
    gain_rng = Xoshiro(spec.gain_seed + architecture_offset(architecture_id))
    feature_gain = max.(0.5, 1 .+ 0.02 .* randn(gain_rng, feature_count))
    attenuation = exp.(-spec.attenuation_slope .* selected_time ./ last(time_s))
    ring_decay = exp.(-selected_time ./ 40.0e-6)
    common_scale = feature_gain .* spec.global_gain
    h0 = common_scale .* core_h0 .* transpose(attenuation) .+
        feature_gain .* spec.ringing_scale .* ringing .* transpose(ring_decay)
    transformed_delta = common_scale .* core_delta .* transpose(attenuation)
    vec(h0), vec(transformed_delta)
end

project_nuisance(vector, basis) = vector - basis * (transpose(basis) * vector)

function build_architecture_model(definition; context_ids=train_context_ids())
    base = load_architecture_scenes(definition, context_ids)
    h0_specs = nuisance_specs(
        SEEDS.train, context_ids, TRAIN_H0_REPLICATES; with_targets=false,
    )
    signal_specs = nuisance_specs(
        SEEDS.train + 0x800, context_ids, TRAIN_SIGNAL_REPLICATES;
        with_targets=true,
    )
    dimension = size(first(values(base.scenes)).h0, 1) * length(base.selected)
    h0 = Matrix{Float64}(undef, dimension, length(h0_specs))
    zero_delta = zeros(Float64, size(first(values(base.scenes)).h0))
    for (index, spec) in enumerate(h0_specs)
        scene = base.scenes[spec.context_id]
        h0[:, index], _ = transform_context(
            scene, zero_delta, base.time_s, base.selected, spec, definition.id,
        )
    end
    calibration_mean = vec(mean(h0; dims=2))
    centered = h0 .- calibration_mean
    decomposition = svd(centered; full=false)
    cumulative = cumsum(abs2.(decomposition.S)) ./ positive_floor(sum(abs2, decomposition.S))
    energy_rank = something(findfirst(>=(NUISANCE_ENERGY_TARGET), cumulative), length(decomposition.S))
    nuisance_rank = min(MAX_NUISANCE_RANK, energy_rank, size(decomposition.U, 2))
    nuisance_basis = decomposition.U[:, 1:nuisance_rank]

    delta_sums = [zeros(Float64, dimension) for _ in TARGET_IDS]
    delta_counts = zeros(Int, length(TARGET_IDS))
    for spec in signal_specs
        scene = base.scenes[spec.context_id]
        _, delta = transform_context(
            scene, scene.deltas[spec.target_index], base.time_s, base.selected,
            spec, definition.id,
        )
        delta_sums[spec.target_index] .+= delta
        delta_counts[spec.target_index] += 1
    end
    raw_templates = hcat((delta_sums[index] ./ delta_counts[index] for index in eachindex(TARGET_IDS))...)
    projected_templates = similar(raw_templates)
    template_norms = zeros(Float64, length(TARGET_IDS))
    for index in eachindex(TARGET_IDS)
        projected_templates[:, index] .= project_nuisance(
            view(raw_templates, :, index), nuisance_basis,
        )
        template_norms[index] = norm(view(projected_templates, :, index))
        template_norms[index] > 0 || error("nuisance projection erased F1c target $index")
        projected_templates[:, index] ./= template_norms[index]
    end
    singular_values = svdvals(projected_templates)
    gram = abs.(transpose(projected_templates) * projected_templates)
    offdiagonal = [gram[i, j] for i in axes(gram, 1), j in axes(gram, 2) if i != j]
    (
        definition,
        base,
        calibration_mean,
        nuisance_basis,
        nuisance_rank,
        nuisance_energy=cumulative[nuisance_rank],
        templates=projected_templates,
        template_norms,
        maximum_template_coherence=maximum(offdiagonal),
        stable_rank=sum(abs2, singular_values) / first(singular_values)^2,
    )
end

function score_components(model, base, specs, electronics_sigma)
    h0_scores = Matrix{Float64}(undef, length(TARGET_IDS), length(specs))
    delta_scores = similar(h0_scores)
    zero_delta = zeros(Float64, size(first(values(base.scenes)).h0))
    for (index, spec) in enumerate(specs)
        scene = base.scenes[spec.context_id]
        physical_delta = spec.target_index == 0 ? zero_delta : scene.deltas[spec.target_index]
        h0, delta = transform_context(
            scene, physical_delta, base.time_s, base.selected,
            spec, model.definition.id,
        )
        projected_h0 = project_nuisance(h0 - model.calibration_mean, model.nuisance_basis)
        projected_delta = project_nuisance(delta, model.nuisance_basis)
        h0_scores[:, index] .= transpose(model.templates) * projected_h0 ./ electronics_sigma
        delta_scores[:, index] .= transpose(model.templates) * projected_delta ./ electronics_sigma
    end
    (; h0_scores, delta_scores, true_cases=getproperty.(specs, :target_index))
end

function noise_factor(model)
    covariance = Symmetric(transpose(model.templates) * model.templates)
    cholesky(covariance + 1.0e-10I).L
end

function h0_statistics(rng, components, factor, trial_count)
    statistics = Vector{Float64}(undef, trial_count)
    for trial in eachindex(statistics)
        context = rand(rng, axes(components.h0_scores, 2))
        scores = view(components.h0_scores, :, context) + factor * randn(rng, size(factor, 1))
        statistics[trial] = maximum(abs, scores)
    end
    statistics
end

function signal_metrics(rng, components, factor, threshold, amplitude)
    detections = 0
    localized = 0
    class_detections = zeros(Int, length(TARGET_IDS))
    class_localized = zeros(Int, length(TARGET_IDS))
    class_totals = zeros(Int, length(TARGET_IDS))
    for context in axes(components.h0_scores, 2), _ in 1:NOISE_DRAWS_PER_SIGNAL_CONTEXT
        truth = components.true_cases[context]
        scores = view(components.h0_scores, :, context) .+
            amplitude .* view(components.delta_scores, :, context) .+
            factor * randn(rng, size(factor, 1))
        statistic, predicted = findmax(abs.(scores))
        detected = statistic > threshold
        correct = detected && predicted == truth
        detections += detected
        localized += correct
        class_totals[truth] += 1
        class_detections[truth] += detected
        class_localized[truth] += correct
    end
    total = sum(class_totals)
    (
        pd=detections / total,
        pdl=localized / total,
        minimum_case_pd=minimum(class_detections ./ class_totals),
        minimum_case_pdl=minimum(class_localized ./ class_totals),
        case_pd=class_detections ./ class_totals,
        case_pdl=class_localized ./ class_totals,
    )
end

function crossing_amplitude(amplitudes, probabilities, target=0.90)
    index = findfirst(>=(target), probabilities)
    isnothing(index) && return Inf
    index == 1 && return amplitudes[1]
    low, high = index - 1, index
    p0, p1 = probabilities[low], probabilities[high]
    p1 == p0 && return amplitudes[high]
    fraction = clamp((target - p0) / (p1 - p0), 0.0, 1.0)
    exp(log(amplitudes[low]) + fraction * (log(amplitudes[high]) - log(amplitudes[low])))
end

function write_namedtuple_csv(path, rows)
    isempty(rows) && return path
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
    path
end

function run_f1c_campaign()
    mkpath(OUTPUT_ROOT)
    definitions = RunOperatorO3FullWaveCarrier.frozen_definitions()
    models = Dict{Symbol,Any}()
    println("[+] building train-only F1c nuisance models")
    for definition in definitions
        println("[+] model $(definition.id)")
        models[definition.id] = build_architecture_model(definition)
    end
    fmc_model = models[:paut_fmc_oracle]
    electronics_sigma = median(fmc_model.template_norms) / 5
    electronics_sigma > 0 || error("F1c electronics noise anchor vanished")
    println("[+] common electronics sigma = $electronics_sigma")

    validation_h0_specs = nuisance_specs(
        SEEDS.validation_h0, [:validation_1], VALIDATION_H0_REPLICATES;
        with_targets=false,
    )
    test_h0_specs = nuisance_specs(
        SEEDS.test_h0, [:test_1], TEST_H0_REPLICATES; with_targets=false,
    )
    validation_signal_specs = nuisance_specs(
        SEEDS.validation_signal, [:validation_1], VALIDATION_SIGNAL_REPLICATES;
        with_targets=true,
    )
    test_signal_specs = nuisance_specs(
        SEEDS.test_signal, [:test_1], TEST_SIGNAL_REPLICATES; with_targets=true,
    )
    evaluations = Dict{Symbol,Any}()
    for definition in definitions
        id = definition.id
        println("[+] stochastic evaluation $id")
        model = models[id]
        validation_base = load_architecture_scenes(definition, [:validation_1])
        test_base = load_architecture_scenes(definition, [:test_1])
        factor = noise_factor(model)
        validation_h0 = score_components(
            model, validation_base, validation_h0_specs, electronics_sigma,
        )
        test_h0 = score_components(model, test_base, test_h0_specs, electronics_sigma)
        validation_signal = score_components(
            model, validation_base, validation_signal_specs, electronics_sigma,
        )
        test_signal = score_components(model, test_base, test_signal_specs, electronics_sigma)
        validation_statistics = h0_statistics(
            Xoshiro(SEEDS.validation_noise + architecture_offset(id)),
            validation_h0, factor, H0_VALIDATION_TRIALS,
        )
        threshold = quantile(validation_statistics, 1 - PFA_TARGET)
        test_statistics = h0_statistics(
            Xoshiro(SEEDS.test_noise + architecture_offset(id)),
            test_h0, factor, H0_TEST_TRIALS,
        )
        false_alarms = count(>(threshold), test_statistics)
        validation_metrics = [
            signal_metrics(
                Xoshiro(SEEDS.validation_noise + architecture_offset(id) + UInt64(index)),
                validation_signal, factor, threshold, amplitude,
            )
            for (index, amplitude) in enumerate(AMPLITUDES)
        ]
        test_metrics = [
            signal_metrics(
                Xoshiro(SEEDS.test_noise + architecture_offset(id) + UInt64(100 + index)),
                test_signal, factor, threshold, amplitude,
            )
            for (index, amplitude) in enumerate(AMPLITUDES)
        ]
        evaluations[id] = (
            threshold,
            observed_pfa=false_alarms / H0_TEST_TRIALS,
            pfa_upper_95=wilson_upper_bound(false_alarms, H0_TEST_TRIALS),
            validation_metrics,
            test_metrics,
            alpha90=crossing_amplitude(AMPLITUDES, getproperty.(test_metrics, :pd)),
            alpha90_localized=crossing_amplitude(AMPLITUDES, getproperty.(test_metrics, :pdl)),
        )
    end

    fmc_validation_pd = getproperty.(evaluations[:paut_fmc_oracle].validation_metrics, :pd)
    reference_index = something(findfirst(>=(0.95), fmc_validation_pd), length(AMPLITUDES))
    reference_alpha = AMPLITUDES[reference_index]
    a_alpha90 = evaluations[:uniform_a_scan].alpha90
    fmc_alpha90 = evaluations[:paut_fmc_oracle].alpha90
    valid_gap = isfinite(a_alpha90) && isfinite(fmc_alpha90) && a_alpha90 > fmc_alpha90
    fmc_reference = evaluations[:paut_fmc_oracle].test_metrics[reference_index]
    positive_control_passed =
        fmc_model.stable_rank >= 3.0 &&
        fmc_model.maximum_template_coherence <= 0.80 &&
        isfinite(a_alpha90) && fmc_alpha90 <= 0.80a_alpha90 &&
        fmc_reference.pdl >= 0.90 &&
        evaluations[:paut_fmc_oracle].pfa_upper_95 <= 0.0125

    summary_rows = NamedTuple[]
    class_rows = NamedTuple[]
    curve_rows = NamedTuple[]
    for definition in definitions
        id = definition.id
        evaluation = evaluations[id]
        reference = evaluation.test_metrics[reference_index]
        gap_closure = valid_gap ?
            (a_alpha90 - evaluation.alpha90) / (a_alpha90 - fmc_alpha90) : NaN
        candidate_gate = definition.kind in (:static, :motion)
        gate_passed = positive_control_passed && candidate_gate &&
            evaluation.pfa_upper_95 <= 0.0125 &&
            reference.pd >= 0.90 &&
            reference.pdl >= 0.80 &&
            reference.minimum_case_pd >= 0.75 &&
            isfinite(gap_closure) && gap_closure >= 0.50
        push!(summary_rows, (
            architecture_id=id,
            kind=definition.kind,
            electronic_channels=definition.electronic_channels,
            physical_states=definition.physical_states,
            transmit_events=definition.transmit_events,
            nuisance_rank=models[id].nuisance_rank,
            nuisance_energy=models[id].nuisance_energy,
            maximum_template_coherence=models[id].maximum_template_coherence,
            template_stable_rank=models[id].stable_rank,
            observed_pfa=evaluation.observed_pfa,
            pfa_upper_95=evaluation.pfa_upper_95,
            reference_alpha,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            minimum_case_pd_reference=reference.minimum_case_pd,
            minimum_case_pdl_reference=reference.minimum_case_pdl,
            alpha90=evaluation.alpha90,
            alpha90_localized=evaluation.alpha90_localized,
            a_scan_to_fmc_gap_closure=gap_closure,
            positive_control_passed,
            f1c_gate_passed=gate_passed,
        ))
        for target_index in eachindex(TARGET_IDS)
            push!(class_rows, (
                architecture_id=id,
                target_id=TARGET_IDS[target_index],
                target_y_mm=TARGET_Y_MM[target_index],
                pd_reference=reference.case_pd[target_index],
                pdl_reference=reference.case_pdl[target_index],
            ))
        end
        for (index, amplitude) in enumerate(AMPLITUDES)
            validation = evaluation.validation_metrics[index]
            test = evaluation.test_metrics[index]
            push!(curve_rows, (
                architecture_id=id,
                amplitude,
                validation_pd=validation.pd,
                validation_pdl=validation.pdl,
                test_pd=test.pd,
                test_pdl=test.pdl,
                test_minimum_case_pd=test.minimum_case_pd,
                test_minimum_case_pdl=test.minimum_case_pdl,
            ))
        end
    end
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "summary.csv"), summary_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "pd_curves.csv"), curve_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "target_metrics.csv"), class_rows)
    jldsave(
        joinpath(OUTPUT_ROOT, "f1c_stochastic.jld2");
        format_version=1,
        contract="docs/operator_o3_f1c_structural_contract.md",
        base_seed=BASE_SEED,
        seeds=SEEDS,
        electronics_sigma,
        amplitudes=AMPLITUDES,
        reference_alpha,
        target_ids=TARGET_IDS,
        target_y_mm=TARGET_Y_MM,
        positive_control_passed,
        summary_rows,
        curve_rows,
        class_rows,
    )

    pd_plot = plot(
        xscale=:log10, xlabel="target amplitude multiplier", ylabel="probability",
        ylims=(0, 1.02), title="F1c structural-clutter test", legend=:bottomright,
    )
    for definition in definitions
        metrics = evaluations[definition.id].test_metrics
        plot!(pd_plot, AMPLITUDES, getproperty.(metrics, :pd);
            marker=:circle, lw=2, label="$(definition.id) Pd")
        plot!(pd_plot, AMPLITUDES, getproperty.(metrics, :pdl);
            ls=:dash, lw=1.5, label="$(definition.id) Pdl")
    end
    hline!(pd_plot, [0.90]; color=:black, ls=:dot, label="0.90 gate")
    vline!(pd_plot, [reference_alpha]; color=:gray40, ls=:dash, label="FMC reference")
    savefig(pd_plot, joinpath(OUTPUT_ROOT, "pd_pdl_curves.png"))

    readme_path = joinpath(OUTPUT_ROOT, "README.md")
    open(readme_path, "w") do io
        println(io, "# O3 F1c structural-clutter stochastic pilot")
        println(io)
        println(io, "> Four actual full-wave physical contexts; nuisance/noise replication is stochastic and is not an industrial POD sample size.")
        println(io)
        println(io, "Positive control: **$(positive_control_passed ? "PASS" : "FAIL")**. Common electronics sigma: `$(electronics_sigma)`; FMC validation reference amplitude: `$(reference_alpha)`.")
        println(io)
        println(io, @sprintf(
            "FMC stable rank `%.3f`, maximum template coherence `%.3f`, alpha90 ratio FMC/A `%.3f`, reference exact Pdl `%.3f`, test Pfa upper `%.3f%%`.",
            fmc_model.stable_rank, fmc_model.maximum_template_coherence,
            fmc_alpha90 / a_alpha90, fmc_reference.pdl,
            100evaluations[:paut_fmc_oracle].pfa_upper_95,
        ))
        println(io)
        println(io, "| Architecture | Pfa upper 95% | Pd ref | exact Pdl ref | min Pd | alpha90 | localized alpha90 | A->FMC closure | Gate |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|---|")
        for row in summary_rows
            closure = isfinite(row.a_scan_to_fmc_gap_closure) ? @sprintf("%.3f", row.a_scan_to_fmc_gap_closure) : "n/a"
            alpha = isfinite(row.alpha90) ? @sprintf("%.3f", row.alpha90) : "inf"
            alpha_local = isfinite(row.alpha90_localized) ? @sprintf("%.3f", row.alpha90_localized) : "inf"
            println(io, @sprintf(
                "| `%s` | %.3f%% | %.3f | %.3f | %.3f | %s | %s | %s | %s |",
                row.architecture_id, 100row.pfa_upper_95, row.pd_reference,
                row.pdl_reference, row.minimum_case_pd_reference, alpha, alpha_local,
                closure, row.f1c_gate_passed ? "pass" : "fail",
            ))
        end
        println(io)
        if !positive_control_passed
            println(io, "Frozen static/motion numbers are diagnostic only because the PAUT positive-control gate failed; they do not open O4.")
        end
        println(io, "All projectors/templates/nuisance bases use only `train_1` and `train_2`; threshold and reference amplitude use only `validation_1`; final metrics use only `test_1`.")
    end
    println("[+] F1c report: $readme_path")
    summary_rows
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_f1c_campaign()

end
