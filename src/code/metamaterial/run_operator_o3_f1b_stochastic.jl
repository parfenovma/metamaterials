module RunOperatorO3F1BStochastic

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o3_fullwave_carrier.jl"))
include(joinpath(@__DIR__, "run_operator_o3_f1_defect_ensemble.jl"))
include(joinpath(@__DIR__, "operator_o3_statistics.jl"))

using .RunOperatorO3FullWaveCarrier
using .RunOperatorO3F1DefectEnsemble
using .OperatorO3Statistics: wilson_upper_bound

export run_f1b_campaign

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_O3_F1B_OUTPUT",
    joinpath(PROJECT_ROOT, "results", "operator_space", "fullwave_f1b_pilot_v0"),
)
const BASE_SEED = UInt64(0x051b2026)
const SEEDS = (
    train=BASE_SEED + 0x101,
    validation_h0=BASE_SEED + 0x202,
    test_h0=BASE_SEED + 0x303,
    validation_signal=BASE_SEED + 0x404,
    test_signal=BASE_SEED + 0x505,
    validation_noise=BASE_SEED + 0x606,
    test_noise=BASE_SEED + 0x707,
)
const AMPLITUDES = [0.25, 0.4, 0.6, 0.8, 1.0, 1.3, 1.7, 2.2, 3.0, 4.0, 6.0, 8.0]
const TRAIN_CONTEXTS_PER_CASE = 16
const VALIDATION_H0_CONTEXTS_PER_CASE = 32
const TEST_H0_CONTEXTS_PER_CASE = 64
const VALIDATION_SIGNAL_CONTEXTS_PER_CASE = 30
const TEST_SIGNAL_CONTEXTS_PER_CASE = 60
const NOISE_DRAWS_PER_SIGNAL_CONTEXT = 4
const H0_VALIDATION_TRIALS = 10_000
const H0_TEST_TRIALS = 10_000
const MAX_NUISANCE_RANK = 24
const NUISANCE_ENERGY_TARGET = 0.99
const PFA_TARGET = 0.01

const CASE_ROWS = [
    (id=:void_r1_focus, x_mm=60.0, y_mm=0.0, size_class=:small_void),
    (id=:void_r2_near, x_mm=35.0, y_mm=0.0, size_class=:nominal_void),
    (id=:void_r2_focus, x_mm=60.0, y_mm=0.0, size_class=:nominal_void),
    (id=:void_r2_deep, x_mm=85.0, y_mm=0.0, size_class=:nominal_void),
    (id=:void_r2_offset, x_mm=60.0, y_mm=10.0, size_class=:nominal_void),
    (id=:slit_30_deg, x_mm=60.0, y_mm=0.0, size_class=:slit),
]

architecture_offset(id) = UInt64(sum(Int, codeunits(string(id)))) * UInt64(0x10001)
positive_floor(value::Real) = max(Float64(value), floatmin(Float64))

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

function load_architecture_signals(definition)
    signals = Dict{Symbol,Any}()
    common_time = nothing
    for row in CASE_ROWS
        free_data = JLD2.load(RunOperatorO3F1DefectEnsemble.response_path(row.id, :free))
        defect_data = JLD2.load(RunOperatorO3F1DefectEnsemble.response_path(row.id, :defect))
        time_s = Float64.(free_data["time_s"])
        time_s == Float64.(defect_data["time_s"]) || error("paired time grids differ for $(row.id)")
        isnothing(common_time) ? (common_time = time_s) :
            common_time == time_s || error("defect cases use different time grids")
        free = project_time_response(Float64.(free_data["normal_response_m2"]), definition)
        defect = project_time_response(Float64.(defect_data["normal_response_m2"]), definition)
        signals[row.id] = (; free, delta=defect - free)
    end
    selected = [
        index for (index, time) in enumerate(common_time)
        if time >= 8.0e-6 && (index - 1) % 3 == 0
    ]
    (; time_s=common_time, selected, signals)
end

function context_specs(seed, contexts_per_case)
    specs = NamedTuple[]
    for (case_index, row) in enumerate(CASE_ROWS), replicate in 1:contexts_per_case
        rng = Xoshiro(seed + UInt64(10_000case_index + replicate))
        push!(specs, (
            case_index,
            case_id=row.id,
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

function sample_signal(signal, time_s, selected, source_time)
    source_time <= first(time_s) && return copy(view(signal, :, 1))
    source_time >= last(time_s) && return zeros(Float64, size(signal, 1))
    upper = searchsortedfirst(time_s, source_time)
    lower = upper - 1
    weight = (source_time - time_s[lower]) / (time_s[upper] - time_s[lower])
    (1 - weight) .* view(signal, :, lower) .+ weight .* view(signal, :, upper)
end

function transform_context(base, time_s, selected, spec, architecture_id)
    feature_count = size(base.free, 1)
    selected_time = time_s[selected]
    core_free = Matrix{Float64}(undef, feature_count, length(selected))
    core_delta = similar(core_free)
    ringing = similar(core_free)
    for (column, time) in enumerate(selected_time)
        source_time = (time - spec.delay_s) * (1 + spec.time_stretch)
        ring_time = (time - spec.delay_s - 1.0e-6) * (1 + spec.time_stretch)
        core_free[:, column] .= sample_signal(base.free, time_s, selected, source_time)
        core_delta[:, column] .= sample_signal(base.delta, time_s, selected, source_time)
        ringing[:, column] .= sample_signal(base.free, time_s, selected, ring_time)
    end
    gain_rng = Xoshiro(spec.gain_seed + architecture_offset(architecture_id))
    feature_gain = max.(0.5, 1 .+ 0.02 .* randn(gain_rng, feature_count))
    attenuation = exp.(-spec.attenuation_slope .* selected_time ./ last(time_s))
    ring_decay = exp.(-selected_time ./ 40.0e-6)
    common_scale = feature_gain .* spec.global_gain
    free = common_scale .* core_free .* transpose(attenuation) .+
        feature_gain .* spec.ringing_scale .* ringing .* transpose(ring_decay)
    delta = common_scale .* core_delta .* transpose(attenuation)
    vec(free), vec(delta)
end

project_nuisance(vector, basis) = vector - basis * (transpose(basis) * vector)

function build_architecture_model(definition, specs)
    base = load_architecture_signals(definition)
    dimension = size(first(values(base.signals)).free, 1) * length(base.selected)
    h0 = Matrix{Float64}(undef, dimension, length(specs))
    delta_sums = [zeros(Float64, dimension) for _ in CASE_ROWS]
    delta_counts = zeros(Int, length(CASE_ROWS))
    for (index, spec) in enumerate(specs)
        free_vector, delta_vector = transform_context(
            base.signals[spec.case_id], base.time_s, base.selected, spec, definition.id,
        )
        h0[:, index] .= free_vector
        delta_sums[spec.case_index] .+= delta_vector
        delta_counts[spec.case_index] += 1
    end
    calibration_mean = vec(mean(h0; dims=2))
    centered = h0 .- calibration_mean
    decomposition = svd(centered; full=false)
    cumulative = cumsum(abs2.(decomposition.S)) ./ positive_floor(sum(abs2, decomposition.S))
    energy_rank = something(findfirst(>=(NUISANCE_ENERGY_TARGET), cumulative), length(decomposition.S))
    nuisance_rank = min(MAX_NUISANCE_RANK, energy_rank, size(decomposition.U, 2))
    nuisance_basis = decomposition.U[:, 1:nuisance_rank]
    raw_templates = hcat((
        delta_sums[index] ./ delta_counts[index]
        for index in eachindex(CASE_ROWS)
    )...)
    projected_templates = similar(raw_templates)
    template_norms = zeros(Float64, length(CASE_ROWS))
    for index in eachindex(CASE_ROWS)
        projected_templates[:, index] .= project_nuisance(
            view(raw_templates, :, index), nuisance_basis,
        )
        template_norms[index] = norm(view(projected_templates, :, index))
        template_norms[index] > 0 || error("nuisance projection erased template $(CASE_ROWS[index].id)")
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

function context_score_components(model, specs, electronics_sigma)
    h0_scores = Matrix{Float64}(undef, length(CASE_ROWS), length(specs))
    delta_scores = similar(h0_scores)
    for (index, spec) in enumerate(specs)
        free_vector, delta_vector = transform_context(
            model.base.signals[spec.case_id], model.base.time_s,
            model.base.selected, spec, model.definition.id,
        )
        h0_projected = project_nuisance(
            free_vector - model.calibration_mean, model.nuisance_basis,
        )
        delta_projected = project_nuisance(delta_vector, model.nuisance_basis)
        h0_scores[:, index] .= transpose(model.templates) * h0_projected ./ electronics_sigma
        delta_scores[:, index] .= transpose(model.templates) * delta_projected ./ electronics_sigma
    end
    (; h0_scores, delta_scores, true_cases=getproperty.(specs, :case_index))
end

function noise_factor(model)
    covariance = Symmetric(transpose(model.templates) * model.templates)
    cholesky(covariance + 1.0e-10I).L
end

function h0_statistics(rng, components, factor, trial_count)
    statistics = Vector{Float64}(undef, trial_count)
    context_count = size(components.h0_scores, 2)
    for trial in 1:trial_count
        context = rand(rng, 1:context_count)
        scores = view(components.h0_scores, :, context) + factor * randn(rng, size(factor, 1))
        statistics[trial] = maximum(abs, scores)
    end
    statistics
end

function localization_ok(predicted, truth)
    wavelength_mm = 6122.102437409232 / 242.0e3 * 1.0e3
    hypot(
        CASE_ROWS[predicted].x_mm - CASE_ROWS[truth].x_mm,
        CASE_ROWS[predicted].y_mm - CASE_ROWS[truth].y_mm,
    ) <= 0.5wavelength_mm
end

function signal_metrics(rng, components, factor, threshold, amplitude)
    detections = 0
    localized = 0
    class_detections = zeros(Int, length(CASE_ROWS))
    class_totals = zeros(Int, length(CASE_ROWS))
    context_count = size(components.h0_scores, 2)
    for context in 1:context_count, _ in 1:NOISE_DRAWS_PER_SIGNAL_CONTEXT
        truth = components.true_cases[context]
        scores = view(components.h0_scores, :, context) .+
            amplitude .* view(components.delta_scores, :, context) .+
            factor * randn(rng, size(factor, 1))
        statistic, predicted = findmax(abs.(scores))
        detected = statistic > threshold
        detections += detected
        class_totals[truth] += 1
        class_detections[truth] += detected
        localized += detected && localization_ok(predicted, truth)
    end
    total = sum(class_totals)
    (
        pd=detections / total,
        pdl=localized / total,
        minimum_case_pd=minimum(class_detections ./ class_totals),
        case_pd=class_detections ./ class_totals,
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

function run_f1b_campaign()
    mkpath(OUTPUT_ROOT)
    definitions = RunOperatorO3FullWaveCarrier.frozen_definitions()
    train_specs = context_specs(SEEDS.train, TRAIN_CONTEXTS_PER_CASE)
    models = Dict{Symbol,Any}()
    println("[+] building train-only F1b nuisance models")
    for definition in definitions
        println("[+] model $(definition.id)")
        models[definition.id] = build_architecture_model(definition, train_specs)
    end
    fmc_model = models[:paut_fmc_oracle]
    electronics_sigma = median(fmc_model.template_norms) / 5
    electronics_sigma > 0 || error("electronics noise anchor vanished")
    println("[+] common electronics sigma = $electronics_sigma")

    validation_h0_specs = context_specs(SEEDS.validation_h0, VALIDATION_H0_CONTEXTS_PER_CASE)
    test_h0_specs = context_specs(SEEDS.test_h0, TEST_H0_CONTEXTS_PER_CASE)
    validation_signal_specs = context_specs(SEEDS.validation_signal, VALIDATION_SIGNAL_CONTEXTS_PER_CASE)
    test_signal_specs = context_specs(SEEDS.test_signal, TEST_SIGNAL_CONTEXTS_PER_CASE)
    evaluations = Dict{Symbol,Any}()
    validation_curves = Dict{Symbol,Vector{Float64}}()
    test_curves = Dict{Symbol,Vector{Float64}}()

    for definition in definitions
        id = definition.id
        println("[+] stochastic evaluation $id")
        model = models[id]
        factor = noise_factor(model)
        validation_h0 = context_score_components(model, validation_h0_specs, electronics_sigma)
        test_h0 = context_score_components(model, test_h0_specs, electronics_sigma)
        validation_signal = context_score_components(model, validation_signal_specs, electronics_sigma)
        test_signal = context_score_components(model, test_signal_specs, electronics_sigma)
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
        pfa = false_alarms / H0_TEST_TRIALS
        pfa_upper_95 = wilson_upper_bound(false_alarms, H0_TEST_TRIALS)
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
        validation_curves[id] = getproperty.(validation_metrics, :pd)
        test_curves[id] = getproperty.(test_metrics, :pd)
        evaluations[id] = (
            threshold,
            observed_pfa=pfa,
            pfa_upper_95,
            validation_metrics,
            test_metrics,
            alpha90=crossing_amplitude(AMPLITUDES, test_curves[id]),
        )
    end

    fmc_validation = validation_curves[:paut_fmc_oracle]
    reference_index = something(findfirst(>=(0.95), fmc_validation), length(AMPLITUDES))
    reference_alpha = AMPLITUDES[reference_index]
    a_alpha90 = evaluations[:uniform_a_scan].alpha90
    fmc_alpha90 = evaluations[:paut_fmc_oracle].alpha90
    valid_gap = isfinite(a_alpha90) && isfinite(fmc_alpha90) && a_alpha90 > fmc_alpha90

    summary_rows = NamedTuple[]
    case_pd_rows = NamedTuple[]
    for definition in definitions
        id = definition.id
        evaluation = evaluations[id]
        reference = evaluation.test_metrics[reference_index]
        worst_case_index = argmin(reference.case_pd)
        gap_closure = valid_gap ?
            (a_alpha90 - evaluation.alpha90) / (a_alpha90 - fmc_alpha90) : NaN
        candidate_gate = definition.kind in (:static, :motion)
        gate_passed = candidate_gate &&
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
            worst_case_id=CASE_ROWS[worst_case_index].id,
            alpha90=evaluation.alpha90,
            a_scan_to_fmc_gap_closure=gap_closure,
            f1b_pilot_gate_passed=gate_passed,
        ))
        for (case_index, case) in enumerate(CASE_ROWS)
            push!(case_pd_rows, (
                architecture_id=id,
                case_id=case.id,
                x_mm=case.x_mm,
                y_mm=case.y_mm,
                pd_reference=reference.case_pd[case_index],
            ))
        end
    end

    curve_rows = NamedTuple[]
    for definition in definitions, (index, amplitude) in enumerate(AMPLITUDES)
        validation = evaluations[definition.id].validation_metrics[index]
        test = evaluations[definition.id].test_metrics[index]
        push!(curve_rows, (
            architecture_id=definition.id,
            amplitude,
            validation_pd=validation.pd,
            validation_pdl=validation.pdl,
            test_pd=test.pd,
            test_pdl=test.pdl,
            test_minimum_case_pd=test.minimum_case_pd,
        ))
    end
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "summary.csv"), summary_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "pd_curves.csv"), curve_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "case_pd_reference.csv"), case_pd_rows)
    jldsave(
        joinpath(OUTPUT_ROOT, "f1b_pilot.jld2");
        format_version=1,
        contract="docs/operator_o3_f1_contract.md",
        base_seed=BASE_SEED,
        seeds=SEEDS,
        electronics_sigma,
        amplitudes=AMPLITUDES,
        reference_alpha,
        case_rows=CASE_ROWS,
        summary_rows,
        curve_rows,
        case_pd_rows,
    )

    pd_plot = plot(
        xscale=:log10, xlabel="defect amplitude multiplier", ylabel="Pd",
        ylims=(0, 1.02), title="F1b full-wave finite-ROI pilot", legend=:bottomright,
    )
    for definition in definitions
        plot!(pd_plot, AMPLITUDES, test_curves[definition.id];
            marker=:circle, lw=2, label=string(definition.id))
    end
    hline!(pd_plot, [0.90]; color=:black, ls=:dash, label="Pd gate")
    vline!(pd_plot, [reference_alpha]; color=:gray40, ls=:dot, label="FMC reference alpha")
    savefig(pd_plot, joinpath(OUTPUT_ROOT, "pd_curves.png"))

    readme_path = joinpath(OUTPUT_ROOT, "README.md")
    open(readme_path, "w") do io
        println(io, "# O3 F1b finite-defect stochastic pilot")
        println(io)
        println(io, "> This is a six-class full-wave pilot, not industrial POD or continuous-ROI validation.")
        println(io)
        println(io, "Common electronics sigma: `$(electronics_sigma)`; FMC validation reference amplitude: `$(reference_alpha)`.")
        println(io)
        println(io, "| Architecture | Pfa upper 95% | Pd ref | Pdl ref | Min class Pd | Worst class | alpha90 | A->FMC closure | Gate |")
        println(io, "|---|---:|---:|---:|---:|---|---:|---:|---|")
        for row in summary_rows
            closure = isfinite(row.a_scan_to_fmc_gap_closure) ? @sprintf("%.3f", row.a_scan_to_fmc_gap_closure) : "n/a"
            alpha = isfinite(row.alpha90) ? @sprintf("%.3f", row.alpha90) : "inf"
            println(io, @sprintf(
                "| `%s` | %.3f%% | %.3f | %.3f | %.3f | `%s` | %s | %s | %s |",
                row.architecture_id, 100row.pfa_upper_95, row.pd_reference,
                row.pdl_reference, row.minimum_case_pd_reference, row.worst_case_id, alpha, closure,
                row.f1b_pilot_gate_passed ? "pass" : "fail",
            ))
        end
        println(io)
        if !valid_gap
            println(io, "The A-scan -> FMC closure is undefined because the uniform A-scan alpha90 is lower than the FMC alpha90 on this six-class task. This is a failed PAUT positive-control condition, not evidence that A-scan dominates general NDT.")
            println(io)
        end
        println(io, "Thresholds and nuisance PCs were frozen on train/validation splits. The detector never receives a per-contact defect-free reference; `free` signals enter only the fixed calibration/nuisance model and raw H0 trials.")
    end
    println("[+] F1b report: $readme_path")
    summary_rows
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_f1b_campaign()

end
