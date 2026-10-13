module RunOperatorO6HBVDLoading

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o6g_network_budget.jl"))
using .RunOperatorO6GNetworkBudget

export scattering_completion,
       bvd_impedance,
       loaded_transfer,
       run_o6h_bvd_loading

const N = RunOperatorO6GNetworkBudget
const F = N.F
const D = N.D
const O6H_BVD_OUTPUT_ROOT = joinpath(
    F.PROJECT_ROOT, "results", "operator_space", "o6h_bvd_loading_v0",
)
const O6H_BVD_SEED = UInt64(0x06f0d117)
const O6H_BVD_DRAWS = 20_000
const Q_GRID = (0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 35.0)
const RESONANCE_RMS = 0.01
const RESISTANCE_LOG_RMS = 0.02
const Q_LOG_RMS = 0.05
const RETURN_MAGNITUDE_GATE = 10.0^(-10 / 20)
const PASS_PROBABILITY_GATE = 0.99

function positive_sqrt(matrix)
    decomposition = eigen(Hermitian((matrix + matrix') / 2))
    minimum(decomposition.values) >= -1e-10 || error("matrix is not PSD")
    decomposition.vectors * Diagonal(sqrt.(max.(decomposition.values, 0.0))) *
        decomposition.vectors'
end

"""Exact matched-dump isometry and symmetric reciprocal completion."""
function scattering_completion(transfer)
    size(transfer, 2) == 3 || throw(DimensionMismatch("q8 must have 3 inputs"))
    deficit = Matrix{Float64}(I, 3, 3) - transpose(transfer) * transfer
    dump = positive_sqrt(deficit)
    forward = vcat(transfer, dump)
    isapprox(forward' * forward, Matrix{Float64}(I, 3, 3); atol=1e-11) ||
        error("matched-dump embedding is not isometric")
    complement = Matrix{Float64}(I, size(forward, 1), size(forward, 1)) -
        forward * forward'
    full_scattering = [zeros(3, 3) forward'; forward complement]
    isapprox(full_scattering, full_scattering'; atol=1e-12) ||
        error("completion lost reciprocity")
    isapprox(full_scattering' * full_scattering,
        Matrix{Float64}(I, size(full_scattering, 1), size(full_scattering, 1));
        atol=1e-10) || error("completion is not lossless")
    (; forward, dump, complement, full_scattering)
end

bvd_impedance(frequency_ratio, quality_factor, resonance_ratio, resistance) =
    resistance * (1 + im * quality_factor *
        (frequency_ratio / resonance_ratio - resonance_ratio / frequency_ratio))

function loaded_transfer(
    completion,
    frequency_ratio,
    quality_factors,
    resonance_ratios,
    resistances,
)
    physical_count = length(quality_factors)
    all(length(values) == physical_count for values in
        (resonance_ratios, resistances)) || throw(DimensionMismatch("load arrays"))
    dump_count = size(completion.forward, 1) - physical_count
    impedances = [bvd_impedance(
        frequency_ratio,
        quality_factors[index],
        resonance_ratios[index],
        resistances[index],
    ) for index in 1:physical_count]
    reflections = (impedances .- 1) ./ (impedances .+ 1)
    gamma = Diagonal(vcat(reflections, zeros(ComplexF64, dump_count)))
    incident = (Matrix{ComplexF64}(I, size(gamma, 1), size(gamma, 1)) -
        completion.complement * gamma) \ completion.forward
    public_reflection = completion.forward' * gamma * incident
    # Power-wave transmission into a positive-resistance BVD load. This is
    # bounded by one and retains the physically relevant load-current phase.
    delivered = 2 .* sqrt.(resistances) ./ (1 .+ impedances)
    mechanical = Diagonal(delivered) * incident[1:physical_count, :]
    (; mechanical, public_reflection, impedances, reflections, delivered)
end

function loaded_band_metrics(
    completion,
    frequencies_hz,
    target,
    quality_factors,
    resonance_ratios,
    resistances,
)
    reference_frequency = 242.0e3
    loads = [loaded_transfer(
        completion,
        frequency / reference_frequency,
        quality_factors,
        resonance_ratios,
        resistances,
    ) for frequency in frequencies_hz]
    transfer = cat(getproperty.(loads, :mechanical)...; dims=3)
    rows = [D.transfer_metrics(transfer[:, :, index], target)
        for index in eachindex(frequencies_hz)]
    delay = D.differential_delay(frequencies_hz, transfer, target)
    maximum_return = maximum(opnorm(load.public_reflection) for load in loads)
    summary = (
        minimum_fidelity=minimum(getproperty.(rows, :fidelity)),
        maximum_coherence=maximum(getproperty.(rows, :coherence)),
        maximum_imbalance_db=maximum(getproperty.(rows, :imbalance_db)),
        maximum_amplitude_rms=maximum(getproperty.(rows, :amplitude_rms)),
        maximum_delay_rms_ns=1e9delay.maximum_rms,
        minimum_useful_amplitude=minimum(getproperty.(rows, :useful_amplitude)),
        maximum_sigma=maximum(getproperty.(rows, :sigma_max)),
        maximum_public_reflection=maximum_return,
        minimum_public_return_loss_db=-20log10(max(maximum_return, eps(Float64))),
    )
    gates = (
        transfer_fidelity=summary.minimum_fidelity >= D.FIDELITY_GATE,
        port_coherence=summary.maximum_coherence <= D.COHERENCE_GATE,
        column_imbalance=summary.maximum_imbalance_db <= D.IMBALANCE_GATE_DB,
        amplitude_rms=summary.maximum_amplitude_rms <= D.AMPLITUDE_RMS_GATE,
        differential_delay=delay.maximum_rms <= D.DELAY_RMS_GATE_S,
        useful_amplitude=summary.minimum_useful_amplitude >= D.USEFUL_AMPLITUDE_GATE,
        passivity=summary.maximum_sigma <= 1 + 1e-10,
        reciprocity=true,
        public_return_loss=maximum_return <= RETURN_MAGNITUDE_GATE,
    )
    (; summary, summary..., gates, passed=all(values(gates)), transfer, loads, delay)
end

function q_result(completion, frequencies_hz, target, quality_factor)
    count = size(target, 1)
    loaded_band_metrics(
        completion,
        frequencies_hz,
        target,
        fill(quality_factor, count),
        ones(count),
        ones(count),
    )
end

function q_ceiling(completion, frequencies_hz, target)
    lower = 0.0
    upper = 1.0
    q_result(completion, frequencies_hz, target, lower).passed ||
        error("matched loads do not pass")
    while q_result(completion, frequencies_hz, target, upper).passed
        lower = upper
        upper *= 2
    end
    for _ in 1:50
        middle = (lower + upper) / 2
        if q_result(completion, frequencies_hz, target, middle).passed
            lower = middle
        else
            upper = middle
        end
    end
    lower, q_result(completion, frequencies_hz, target, lower)
end

function limiting_gate(result)
    margins = (
        transfer_fidelity=result.minimum_fidelity / D.FIDELITY_GATE,
        port_coherence=D.COHERENCE_GATE / max(result.maximum_coherence, 1e-15),
        column_imbalance=D.IMBALANCE_GATE_DB /
            max(result.maximum_imbalance_db, 1e-15),
        amplitude_rms=D.AMPLITUDE_RMS_GATE /
            max(result.maximum_amplitude_rms, 1e-15),
        differential_delay=1e9D.DELAY_RMS_GATE_S /
            max(result.maximum_delay_rms_ns, 1e-15),
        useful_amplitude=result.minimum_useful_amplitude / D.USEFUL_AMPLITUDE_GATE,
        passivity=result.maximum_sigma <= 1 + 1e-9 ? Inf :
            1 / max(result.maximum_sigma, 1e-15),
        public_return_loss=RETURN_MAGNITUDE_GATE /
            max(result.maximum_public_reflection, 1e-15),
    )
    names = propertynames(margins)
    String(names[argmin(collect(values(margins)))])
end

function monte_carlo_summary(
    scale,
    quality_factor,
    completion,
    frequencies_hz,
    target,
    standards,
)
    draws = size(standards.resonance, 2)
    passed = falses(draws)
    useful = zeros(draws)
    reflection = zeros(draws)
    amplitude_rms = zeros(draws)
    for draw in 1:draws
        resonance = exp.(scale * RESONANCE_RMS .* view(
            standards.resonance, :, draw,
        ))
        resistance = exp.(scale * RESISTANCE_LOG_RMS .* view(
            standards.resistance, :, draw,
        ))
        qualities = quality_factor .* exp.(scale * Q_LOG_RMS .* view(
            standards.quality, :, draw,
        ))
        result = loaded_band_metrics(
            completion,
            frequencies_hz,
            target,
            qualities,
            resonance,
            resistance,
        )
        passed[draw] = result.passed
        useful[draw] = result.minimum_useful_amplitude
        reflection[draw] = result.maximum_public_reflection
        amplitude_rms[draw] = result.maximum_amplitude_rms
    end
    successes = count(passed)
    (
        scale,
        effective_quality_factor=quality_factor,
        resonance_frequency_rms=scale * RESONANCE_RMS,
        resistance_log_rms=scale * RESISTANCE_LOG_RMS,
        quality_factor_log_rms=scale * Q_LOG_RMS,
        draws,
        successes,
        observed_joint_pass_probability=successes / draws,
        wilson_lower_95=F.wilson_lower_bound(successes, draws),
        q01_minimum_useful_amplitude=quantile(useful, 0.01),
        q99_maximum_public_reflection=quantile(reflection, 0.99),
        q99_maximum_amplitude_rms=quantile(amplitude_rms, 0.99),
        probability_gate_passed=
            F.wilson_lower_bound(successes, draws) >= PASS_PROBABILITY_GATE,
    )
end

function tolerance_ceiling(args...)
    zero = monte_carlo_summary(0.0, args...)
    zero.probability_gate_passed || return zero
    lower = 0.0
    upper = 1.0
    while monte_carlo_summary(upper, args...).probability_gate_passed
        lower = upper
        upper *= 2
    end
    for _ in 1:12
        middle = (lower + upper) / 2
        if monte_carlo_summary(middle, args...).probability_gate_passed
            lower = middle
        else
            upper = middle
        end
    end
    monte_carlo_summary(lower, args...)
end

function write_rows(path, rows)
    F.write_namedtuple_csv(path, rows)
end

function run_o6h_bvd_loading()
    mkpath(O6H_BVD_OUTPUT_ROOT)
    o6d = load(F.O6D_PATH)
    target = Float64.(o6d["target"])
    nominal = Float64.(o6d["quantized_transfers"][8])
    frequencies_hz = Float64.(o6d["frequencies_hz"])
    completion = scattering_completion(nominal)
    println("[>] O6h-A deterministic loaded-BVD Q frontier")
    q_rows = NamedTuple[]
    for quality_factor in Q_GRID
        result = q_result(completion, frequencies_hz, target, quality_factor)
        push!(q_rows, merge((; quality_factor), result.summary,
            (passed=result.passed, active_gate=limiting_gate(result))))
    end
    quality_ceiling, quality_boundary = q_ceiling(
        completion, frequencies_hz, target,
    )
    audit_quality = 0.5quality_ceiling
    rng = Xoshiro(O6H_BVD_SEED)
    standards = (
        resonance=randn(rng, 17, O6H_BVD_DRAWS),
        resistance=randn(rng, 17, O6H_BVD_DRAWS),
        quality=randn(rng, 17, O6H_BVD_DRAWS),
    )
    println("[>] O6h-A loaded-BVD mismatch tolerance")
    tolerance = tolerance_ceiling(
        audit_quality, completion, frequencies_hz, target, standards,
    )
    q35_result = q_result(completion, frequencies_hz, target, 35.0)
    branch_passed = quality_ceiling > 0 && tolerance.probability_gate_passed
    write_rows(joinpath(O6H_BVD_OUTPUT_ROOT, "loaded_q_frontier.csv"), q_rows)
    write_rows(joinpath(O6H_BVD_OUTPUT_ROOT, "loaded_q_ceiling.csv"), [(
        maximum_passing_effective_q=quality_ceiling,
        audit_effective_q=audit_quality,
        active_gate=limiting_gate(quality_boundary),
        q35_passed=q35_result.passed,
        q35_return_loss_db=q35_result.minimum_public_return_loss_db,
        q35_useful_amplitude=q35_result.minimum_useful_amplitude,
    )])
    write_rows(joinpath(O6H_BVD_OUTPUT_ROOT, "loaded_mismatch_threshold.csv"), [tolerance])
    jldsave(
        joinpath(O6H_BVD_OUTPUT_ROOT, "o6h_bvd_loading.jld2");
        format_version=1,
        source_o6d_path=abspath(F.O6D_PATH),
        target,
        nominal,
        frequencies_hz,
        completion,
        q_rows,
        quality_ceiling,
        quality_boundary_summary=quality_boundary.summary,
        audit_quality,
        tolerance,
        q35_summary=q35_result.summary,
        q35_gates=q35_result.gates,
        random_seed=O6H_BVD_SEED,
        draws=O6H_BVD_DRAWS,
        branch_passed,
    )
    open(joinpath(O6H_BVD_OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6h-A reciprocal q8 network with loaded BVD terminals")
        println(io)
        println(io, branch_passed ?
            "**PASS for a finite broadband-matched BVD region.**" :
            "**FAIL: no finite loaded-BVD region survives the frozen gates.**")
        println(io)
        @printf(io,
            "The q8 contraction is embedded in an exact `%d x %d` reciprocal lossless scattering matrix with three matched dump ports; isometry error is `%.3e`.\n\n",
            size(completion.full_scattering)...,
            norm(completion.forward' * completion.forward - Matrix{Float64}(I, 3, 3)),
        )
        @printf(io,
            "With all 17 series-BVD loads matched at f0, the maximum effective Q across the 0.8--1.2 band is `%.5f`; active gate `%s`. The mismatch audit is run at half this value, Q=`%.5f`.\n\n",
            quality_ceiling, limiting_gate(quality_boundary), audit_quality,
        )
        @printf(io,
            "At audit Q, the 99%%-Wilson load-error threshold is scale `%.4f`: resonance RMS `%.3f%%`, resistance log-RMS `%.3f%%`, Q log-RMS `%.3f%%`; Wilson lower-95 `%.3f%%`.\n\n",
            tolerance.scale, 100tolerance.resonance_frequency_rms,
            100tolerance.resistance_log_rms,
            100tolerance.quality_factor_log_rms,
            100tolerance.wilson_lower_95,
        )
        @printf(io,
            "The frozen intrinsic-Q=35 load is `%s`: worst-band public return loss `%.3f dB`, useful amplitude `%.4f`. It therefore requires broadband/multiresonant matching to reach the effective-Q ceiling; this is not a global passive-network no-go.\n\n",
            q35_result.passed ? "PASS" : "FAIL",
            q35_result.minimum_public_return_loss_db,
            q35_result.minimum_useful_amplitude,
        )
        println(io, "The model includes reciprocal multiple reflections and matched dump ports, but remains normalized: a concrete matching circuit, transformer leakage and electrical-noise model are still required before fabrication.")
    end
    q_plot = plot(
        getproperty.(q_rows, :quality_factor),
        getproperty.(q_rows, :minimum_useful_amplitude);
        marker=:circle,
        linewidth=2,
        xlabel="effective series-BVD Q",
        ylabel="minimum useful amplitude",
        title="O6h-A loaded q8 bandwidth",
        label="loaded q8",
    )
    hline!(q_plot, [D.USEFUL_AMPLITUDE_GATE];
        linestyle=:dash, color=:black, label="O4 gate")
    vline!(q_plot, [quality_ceiling];
        linestyle=:dot, color=:red, label="Q ceiling")
    savefig(q_plot, joinpath(O6H_BVD_OUTPUT_ROOT, "loaded_q_frontier.png"))
    println("[+] O6h-A effective Q ceiling=$(round(quality_ceiling; digits=5))")
    println("[+] wrote O6h-A to $O6H_BVD_OUTPUT_ROOT")
    (; branch_passed, quality_ceiling, tolerance)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6h_bvd_loading()

end
