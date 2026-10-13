module RunOperatorO6GNetworkBudget

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o6f_beamformer_network.jl"))
using .RunOperatorO6FBeamformerNetwork

export factor_orthogonal_transform,
       complete_network_transfer,
       run_o6g_network_budget

const F = RunOperatorO6FBeamformerNetwork
const D = F.RunOperatorO6DIsolatedBusBound
const O6G_NETWORK_OUTPUT_ROOT = joinpath(
    F.PROJECT_ROOT, "results", "operator_space", "o6g_network_budget_v0",
)
const O6G_NETWORK_SEED = UInt64(0x06f0c117)
const O6G_MONTE_CARLO_DRAWS = 20_000
const O6G_ERROR_SCALES = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
const COUPLER_ANGLE_RMS_RAD = deg2rad(0.1)
const COUPLER_DELAY_RMS_S = 2.0e-9
const ATTENUATOR_LOG_AMPLITUDE_RMS = 0.01
const ATTENUATOR_DELAY_RMS_S = 2.0e-9
const PASS_PROBABILITY_GATE = 0.99

"""Three-Givens signed factorization of a real 3 x 3 orthogonal map."""
function factor_orthogonal_transform(transform)
    size(transform) == (3, 3) || throw(DimensionMismatch("expected 3 x 3 map"))
    factor = F.synthesize_isometry_lattice(transform)
    length(factor.stages) == 3 || error("3 x 3 map did not use three rotations")
    factor
end

function apply_rotation_factor(
    field,
    factor,
    frequency_hz;
    insertion_loss_db=0.0,
    angle_errors=zeros(length(factor.stages)),
    delay_errors_s=zeros(length(factor.stages)),
)
    result = ComplexF64.(factor.core_rotation) * field
    amplitude = 10.0^(-insertion_loss_db / 20)
    for (index, stage) in enumerate(factor.stages)
        phase = exp(-2pi * im * frequency_hz * delay_errors_s[index])
        block = amplitude * phase .* F.rotation_block(
            stage.angle + angle_errors[index],
        )
        result[[stage.upper, stage.lower], :] .=
            block * result[[stage.upper, stage.lower], :]
    end
    result
end

function complete_network_transfer(
    output_factor,
    first_input_factor,
    second_input_factor,
    singular_values,
    frequency_hz;
    coupler_insertion_loss_db=0.0,
    attenuator_excess_loss_db=0.0,
    coupler_angle_errors=zeros(51),
    coupler_delay_errors_s=zeros(51),
    attenuator_log_amplitude_errors=zeros(3),
    attenuator_delay_errors_s=zeros(3),
)
    length(coupler_angle_errors) == 51 ||
        throw(DimensionMismatch("complete network has 51 couplers"))
    length(coupler_delay_errors_s) == 51 ||
        throw(DimensionMismatch("complete network has 51 couplers"))
    field = Matrix{ComplexF64}(I, 3, 3)
    field = apply_rotation_factor(
        field,
        first_input_factor,
        frequency_hz;
        insertion_loss_db=coupler_insertion_loss_db,
        angle_errors=view(coupler_angle_errors, 1:3),
        delay_errors_s=view(coupler_delay_errors_s, 1:3),
    )
    # A passive attenuator cannot realize gain above unity. Clipping each
    # perturbed local transmission at one is a component constraint, not a
    # global post-error contraction or rebalance.
    excess_amplitude = 10.0^(-attenuator_excess_loss_db / 20)
    gains = min.(1.0,
        singular_values .* excess_amplitude .* exp.(attenuator_log_amplitude_errors),
    )
    phases = exp.(-2pi * im * frequency_hz .* attenuator_delay_errors_s)
    field = Diagonal(gains .* phases) * field
    field = apply_rotation_factor(
        field,
        second_input_factor,
        frequency_hz;
        insertion_loss_db=coupler_insertion_loss_db,
        angle_errors=view(coupler_angle_errors, 4:6),
        delay_errors_s=view(coupler_delay_errors_s, 4:6),
    )
    embedded = zeros(ComplexF64, 17, 3)
    embedded[1:3, :] .= field
    amplitude = 10.0^(-coupler_insertion_loss_db / 20)
    for (local_index, stage) in enumerate(output_factor.stages)
        global_index = local_index + 6
        phase = exp(-2pi * im * frequency_hz *
            coupler_delay_errors_s[global_index])
        block = amplitude * phase .* F.rotation_block(
            stage.angle + coupler_angle_errors[global_index],
        )
        embedded[[stage.upper, stage.lower], :] .=
            block * embedded[[stage.upper, stage.lower], :]
    end
    embedded
end

function complete_band_metrics(
    output_factor,
    first_input_factor,
    second_input_factor,
    singular_values,
    frequencies_hz,
    target;
    coupler_insertion_loss_db=0.0,
    attenuator_excess_loss_db=0.0,
    coupler_angle_errors=zeros(51),
    coupler_delay_errors_s=zeros(51),
    attenuator_log_amplitude_errors=zeros(3),
    attenuator_delay_errors_s=zeros(3),
)
    transfer = cat((complete_network_transfer(
        output_factor,
        first_input_factor,
        second_input_factor,
        singular_values,
        frequency;
        coupler_insertion_loss_db,
        attenuator_excess_loss_db,
        coupler_angle_errors,
        coupler_delay_errors_s,
        attenuator_log_amplitude_errors,
        attenuator_delay_errors_s,
    ) for frequency in frequencies_hz)...; dims=3)
    rows = [D.transfer_metrics(transfer[:, :, index], target)
        for index in eachindex(frequencies_hz)]
    delay = D.differential_delay(frequencies_hz, transfer, target)
    summary = (
        minimum_fidelity=minimum(getproperty.(rows, :fidelity)),
        maximum_coherence=maximum(getproperty.(rows, :coherence)),
        maximum_imbalance_db=maximum(getproperty.(rows, :imbalance_db)),
        maximum_amplitude_rms=maximum(getproperty.(rows, :amplitude_rms)),
        maximum_delay_rms_ns=1e9delay.maximum_rms,
        minimum_useful_amplitude=minimum(getproperty.(rows, :useful_amplitude)),
        maximum_sigma=maximum(getproperty.(rows, :sigma_max)),
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
    )
    (; summary, summary..., gates, passed=all(values(gates)), transfer, delay)
end

function scalar_ceiling(evaluator; initial_upper)
    lower = 0.0
    evaluator(lower).passed || error("zero-error complete network does not pass")
    upper = initial_upper
    while evaluator(upper).passed
        lower = upper
        upper *= 2
    end
    for _ in 1:50
        middle = (lower + upper) / 2
        if evaluator(middle).passed
            lower = middle
        else
            upper = middle
        end
    end
    lower, evaluator(lower)
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
    )
    names = propertynames(margins)
    String(names[argmin(collect(values(margins)))])
end

function monte_carlo_summary(
    scale,
    output_factor,
    first_input_factor,
    second_input_factor,
    singular_values,
    frequencies_hz,
    target,
    standards;
    coupler_insertion_loss_db=0.0,
    attenuator_excess_loss_db=0.0,
)
    draws = size(standards.coupler_angle, 2)
    passes = falses(draws)
    amplitude_rms = zeros(draws)
    delay_rms_ns = zeros(draws)
    useful_amplitude = zeros(draws)
    for draw in 1:draws
        result = complete_band_metrics(
            output_factor,
            first_input_factor,
            second_input_factor,
            singular_values,
            frequencies_hz,
            target;
            coupler_insertion_loss_db,
            attenuator_excess_loss_db,
            coupler_angle_errors=scale * COUPLER_ANGLE_RMS_RAD .* view(
                standards.coupler_angle, :, draw,
            ),
            coupler_delay_errors_s=scale * COUPLER_DELAY_RMS_S .* view(
                standards.coupler_delay, :, draw,
            ),
            attenuator_log_amplitude_errors=
                scale * ATTENUATOR_LOG_AMPLITUDE_RMS .* view(
                    standards.attenuator_amplitude, :, draw,
                ),
            attenuator_delay_errors_s=scale * ATTENUATOR_DELAY_RMS_S .* view(
                standards.attenuator_delay, :, draw,
            ),
        )
        passes[draw] = result.passed
        amplitude_rms[draw] = result.maximum_amplitude_rms
        delay_rms_ns[draw] = result.maximum_delay_rms_ns
        useful_amplitude[draw] = result.minimum_useful_amplitude
    end
    successes = count(passes)
    (
        scale,
        coupler_angle_rms_deg=rad2deg(scale * COUPLER_ANGLE_RMS_RAD),
        coupler_delay_rms_ns=1e9scale * COUPLER_DELAY_RMS_S,
        attenuator_log_amplitude_rms=scale * ATTENUATOR_LOG_AMPLITUDE_RMS,
        attenuator_delay_rms_ns=1e9scale * ATTENUATOR_DELAY_RMS_S,
        coupler_insertion_loss_db,
        attenuator_excess_loss_db,
        draws,
        successes,
        observed_joint_pass_probability=successes / draws,
        wilson_lower_95=F.wilson_lower_bound(successes, draws),
        q99_amplitude_rms=quantile(amplitude_rms, 0.99),
        q99_delay_rms_ns=quantile(delay_rms_ns, 0.99),
        q01_useful_amplitude=quantile(useful_amplitude, 0.01),
        probability_gate_passed=
            F.wilson_lower_bound(successes, draws) >= PASS_PROBABILITY_GATE,
    )
end

function monte_carlo_threshold(args...; kwargs...)
    zero = monte_carlo_summary(0.0, args...; kwargs...)
    zero.probability_gate_passed || return zero
    lower = 0.0
    upper = 1.0
    while monte_carlo_summary(upper, args...; kwargs...).probability_gate_passed
        lower = upper
        upper *= 2
    end
    for _ in 1:12
        middle = (lower + upper) / 2
        result = monte_carlo_summary(middle, args...; kwargs...)
        if result.probability_gate_passed
            lower = middle
        else
            upper = middle
        end
    end
    monte_carlo_summary(lower, args...; kwargs...)
end

function write_rows(path, rows)
    F.write_namedtuple_csv(path, rows)
end

function run_o6g_network_budget()
    isfile(F.O6D_PATH) || error("O6d input is missing")
    mkpath(O6G_NETWORK_OUTPUT_ROOT)
    o6d = load(F.O6D_PATH)
    target = Float64.(o6d["target"])
    nominal = Float64.(o6d["quantized_transfers"][8])
    frequencies_hz = Float64.(o6d["frequencies_hz"])
    decomposition = svd(nominal; full=false)
    output_factor = F.synthesize_isometry_lattice(decomposition.U)
    first_input_factor = factor_orthogonal_transform(Matrix(decomposition.Vt))
    second_input_factor = factor_orthogonal_transform(output_factor.core_rotation)
    nominal_transfer = complete_network_transfer(
        output_factor,
        first_input_factor,
        second_input_factor,
        decomposition.S,
        242.0e3,
    )
    reconstruction_error = norm(nominal_transfer - nominal) / norm(nominal)
    reconstruction_error <= 1e-11 || error("complete-network reconstruction failed")
    total_couplers = length(output_factor.stages) +
        length(first_input_factor.stages) + length(second_input_factor.stages)
    total_couplers == 51 || error("complete q8 network should have 51 couplers")
    nominal_result = complete_band_metrics(
        output_factor,
        first_input_factor,
        second_input_factor,
        decomposition.S,
        frequencies_hz,
        target,
    )
    nominal_result.passed || error("nominal complete q8 network fails O4")

    println("[>] O6g-A complete-network scalar loss ceilings")
    coupler_ceiling, coupler_boundary = scalar_ceiling(
        loss -> complete_band_metrics(
            output_factor,
            first_input_factor,
            second_input_factor,
            decomposition.S,
            frequencies_hz,
            target;
            coupler_insertion_loss_db=loss,
        ); initial_upper=0.05,
    )
    attenuator_ceiling, attenuator_boundary = scalar_ceiling(
        loss -> complete_band_metrics(
            output_factor,
            first_input_factor,
            second_input_factor,
            decomposition.S,
            frequencies_hz,
            target;
            attenuator_excess_loss_db=loss,
        ); initial_upper=0.5,
    )
    half_coupler_loss = 0.5coupler_ceiling
    half_attenuator_loss = 0.5attenuator_ceiling
    half_result = complete_band_metrics(
        output_factor,
        first_input_factor,
        second_input_factor,
        decomposition.S,
        frequencies_hz,
        target;
        coupler_insertion_loss_db=half_coupler_loss,
        attenuator_excess_loss_db=half_attenuator_loss,
    )
    ray_ceiling, ray_boundary = scalar_ceiling(
        scale -> complete_band_metrics(
            output_factor,
            first_input_factor,
            second_input_factor,
            decomposition.S,
            frequencies_hz,
            target;
            coupler_insertion_loss_db=scale * half_coupler_loss,
            attenuator_excess_loss_db=scale * half_attenuator_loss,
        ); initial_upper=1.0,
    )
    loss_rows = NamedTuple[]
    for coupler_fraction in range(0.0, 1.1; length=12),
        attenuator_fraction in range(0.0, 1.1; length=12)
        result = complete_band_metrics(
            output_factor,
            first_input_factor,
            second_input_factor,
            decomposition.S,
            frequencies_hz,
            target;
            coupler_insertion_loss_db=coupler_fraction * coupler_ceiling,
            attenuator_excess_loss_db=attenuator_fraction * attenuator_ceiling,
        )
        push!(loss_rows, (
            coupler_fraction,
            attenuator_fraction,
            coupler_insertion_loss_db=coupler_fraction * coupler_ceiling,
            attenuator_excess_loss_db=attenuator_fraction * attenuator_ceiling,
            minimum_fidelity=result.minimum_fidelity,
            maximum_amplitude_rms=result.maximum_amplitude_rms,
            minimum_useful_amplitude=result.minimum_useful_amplitude,
            passed=result.passed,
        ))
    end

    rng = Xoshiro(O6G_NETWORK_SEED)
    standards = (
        coupler_angle=randn(rng, 51, O6G_MONTE_CARLO_DRAWS),
        coupler_delay=randn(rng, 51, O6G_MONTE_CARLO_DRAWS),
        attenuator_amplitude=randn(rng, 3, O6G_MONTE_CARLO_DRAWS),
        attenuator_delay=randn(rng, 3, O6G_MONTE_CARLO_DRAWS),
    )
    println("[>] O6g-A zero-loss joint manufacturing threshold")
    zero_loss_threshold = monte_carlo_threshold(
        output_factor,
        first_input_factor,
        second_input_factor,
        decomposition.S,
        frequencies_hz,
        target,
        standards,
    )
    println("[>] O6g-A half-scalar-ceiling loss allocation")
    allocated_threshold = monte_carlo_threshold(
        output_factor,
        first_input_factor,
        second_input_factor,
        decomposition.S,
        frequencies_hz,
        target,
        standards;
        coupler_insertion_loss_db=half_coupler_loss,
        attenuator_excess_loss_db=half_attenuator_loss,
    )
    ceiling_rows = [
        (
            parameter="uniform insertion loss per all 51 rotation couplers",
            ceiling=coupler_ceiling,
            unit="dB",
            active_gate=limiting_gate(coupler_boundary),
        ),
        (
            parameter="uniform excess loss per three attenuator paths",
            ceiling=attenuator_ceiling,
            unit="dB",
            active_gate=limiting_gate(attenuator_boundary),
        ),
        (
            parameter="half-ceiling joint loss ray multiplier",
            ceiling=ray_ceiling,
            unit="scale",
            active_gate=limiting_gate(ray_boundary),
        ),
    ]
    threshold_rows = [merge((allocation="zero_loss",), zero_loss_threshold),
        merge((allocation="half_each_scalar_ceiling",), allocated_threshold)]
    branch_passed = half_result.passed && allocated_threshold.probability_gate_passed
    write_rows(joinpath(O6G_NETWORK_OUTPUT_ROOT, "component_loss_ceilings.csv"), ceiling_rows)
    write_rows(joinpath(O6G_NETWORK_OUTPUT_ROOT, "joint_loss_frontier.csv"), loss_rows)
    write_rows(joinpath(O6G_NETWORK_OUTPUT_ROOT, "joint_tolerance_thresholds.csv"), threshold_rows)
    jldsave(
        joinpath(O6G_NETWORK_OUTPUT_ROOT, "o6g_network_budget.jld2");
        format_version=1,
        source_o6d_path=abspath(F.O6D_PATH),
        target,
        nominal,
        frequencies_hz,
        singular_values=decomposition.S,
        first_input_factor,
        second_input_factor,
        output_factor,
        total_couplers,
        reconstruction_error,
        nominal_summary=nominal_result.summary,
        coupler_ceiling,
        attenuator_ceiling,
        half_coupler_loss,
        half_attenuator_loss,
        half_summary=half_result.summary,
        ray_ceiling,
        loss_rows,
        zero_loss_threshold,
        allocated_threshold,
        random_seed=O6G_NETWORK_SEED,
        monte_carlo_draws=O6G_MONTE_CARLO_DRAWS,
        branch_passed,
    )
    open(joinpath(O6G_NETWORK_OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6g-A complete q8 passive-network budget")
        println(io)
        println(io, branch_passed ?
            "**PASS at the preregistered half-each scalar-loss allocation.**" :
            "**FAIL at the preregistered half-each scalar-loss allocation.**")
        println(io)
        @printf(io,
            "The complete signed factorization uses `%d` rotation couplers plus three singular attenuator/dump paths and reconstructs q8 with relative error `%.3e`. The singular transmissions are `%.9f / %.9f / %.9f`.\n\n",
            total_couplers, reconstruction_error, decomposition.S...,
        )
        @printf(io,
            "Zero-error scalar ceilings are `%.6f dB` per rotation coupler (%s) and `%.6f dB` excess per attenuator path (%s). The ray through half of both scalar ceilings reaches multiplier `%.4f`.\n\n",
            coupler_ceiling, limiting_gate(coupler_boundary),
            attenuator_ceiling, limiting_gate(attenuator_boundary), ray_ceiling,
        )
        @printf(io,
            "At zero loss, the 99%%-Wilson joint manufacturing threshold is scale `%.4f` (`%.4f deg`, `%.4f ns`, attenuator log-amplitude `%.3f%%`, attenuator delay `%.4f ns`).\n\n",
            zero_loss_threshold.scale,
            zero_loss_threshold.coupler_angle_rms_deg,
            zero_loss_threshold.coupler_delay_rms_ns,
            100zero_loss_threshold.attenuator_log_amplitude_rms,
            zero_loss_threshold.attenuator_delay_rms_ns,
        )
        @printf(io,
            "At the half-each loss allocation (`%.6f dB` per coupler, `%.6f dB` per attenuator path), the joint manufacturing threshold is scale `%.4f`, with Wilson lower-95 `%.3f%%`.\n\n",
            half_coupler_loss, half_attenuator_loss,
            allocated_threshold.scale, 100allocated_threshold.wilson_lower_95,
        )
        println(io, "This remains an ideal matched network. Loaded piezo impedance, transformer leakage, electrical noise and component-specific dispersion are promoted to O6h and are not claimed here.")
    end
    frontier_plot = scatter(
        getproperty.(loss_rows, :coupler_insertion_loss_db),
        getproperty.(loss_rows, :attenuator_excess_loss_db);
        marker_z=Float64.(getproperty.(loss_rows, :passed)),
        color=:viridis,
        xlabel="loss per rotation coupler, dB",
        ylabel="excess loss per attenuator path, dB",
        title="O6g-A complete q8 joint loss frontier",
        label="grid",
    )
    savefig(frontier_plot, joinpath(O6G_NETWORK_OUTPUT_ROOT, "joint_loss_frontier.png"))
    println(branch_passed ? "[+] O6g-A joint-budget PASS" : "[-] O6g-A joint-budget FAIL")
    println("[+] wrote O6g-A to $O6G_NETWORK_OUTPUT_ROOT")
    (; branch_passed, ceiling_rows, threshold_rows)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6g_network_budget()

end
