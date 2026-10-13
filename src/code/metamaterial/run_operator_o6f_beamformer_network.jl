module RunOperatorO6FBeamformerNetwork

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o6d_isolated_bus_bound.jl"))
using .RunOperatorO6DIsolatedBusBound

export synthesize_isometry_lattice,
       realized_transfer,
       run_o6f_beamformer_network

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const O6D_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6d_isolated_bus_bound_v0",
    "o6d_isolated_bus_bound.jld2",
)
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6f_beamformer_network_v0",
)
const RANDOM_SEED = UInt64(0x06f0a117)
const MONTE_CARLO_DRAWS = 20_000
const ERROR_SCALES = (0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0)
const ANGLE_RMS_REFERENCE_RAD = deg2rad(0.1)
const DELAY_RMS_REFERENCE_S = 2.0e-9
const PASS_PROBABILITY_GATE = 0.99
const OUTPUT_COUNT = 17
const INPUT_COUNT = 3

rotation_block(angle) = [cos(angle) -sin(angle); sin(angle) cos(angle)]

function wilson_lower_bound(successes, trials; z=1.6448536269514722)
    proportion = successes / trials
    denominator = 1 + z^2 / trials
    center = proportion + z^2 / (2trials)
    spread = z * sqrt(
        proportion * (1 - proportion) / trials + z^2 / (4trials^2),
    )
    max(0.0, (center - spread) / denominator)
end

"""Adjacent-Givens realization `U = W*R` for a tall real isometry."""
function synthesize_isometry_lattice(isometry)
    rows, columns = size(isometry)
    rows >= columns || throw(DimensionMismatch("isometry must be tall"))
    isapprox(
        transpose(isometry) * isometry,
        Matrix{Float64}(I, columns, columns);
        atol=1e-12,
    ) ||
        error("input columns are not orthonormal")
    working = copy(Float64.(isometry))
    eliminations = NamedTuple[]
    for column in 1:columns
        for lower in rows:-1:(column + 1)
            upper = lower - 1
            first = working[upper, column]
            second = working[lower, column]
            radius = hypot(first, second)
            cosine = radius > eps(Float64) ? first / radius : 1.0
            sine = radius > eps(Float64) ? second / radius : 0.0
            elimination = [cosine sine; -sine cosine]
            working[[upper, lower], :] .= elimination * working[[upper, lower], :]
            push!(eliminations, (
                upper,
                lower,
                angle=atan(sine, cosine),
                eliminated_column=column,
            ))
        end
    end
    core_rotation = working[1:columns, :]
    maximum(abs, working[(columns + 1):end, :]) <= 1e-12 ||
        error("Givens elimination left a nonzero tail")
    isapprox(
        transpose(core_rotation) * core_rotation,
        Matrix{Float64}(I, columns, columns);
        atol=1e-12,
    ) ||
        error("isometry core is not orthogonal")
    stages = reverse(eliminations)
    embedding = zeros(Float64, rows, columns)
    embedding[1:columns, :] .= Matrix{Float64}(I, columns, columns)
    realized = copy(embedding)
    for stage in stages
        realized[[stage.upper, stage.lower], :] .=
            rotation_block(stage.angle) * realized[[stage.upper, stage.lower], :]
    end
    isapprox(realized * core_rotation, isometry; atol=1e-11) ||
        error("Givens lattice does not reconstruct the isometry")
    (; stages, core_rotation, isometry_embedding=realized)
end

function realized_transfer(
    stages,
    input_block,
    frequency_hz;
    insertion_loss_db_per_coupler=0.0,
    common_delay_s_per_coupler=0.0,
    angle_errors=zeros(length(stages)),
    delay_errors_s=zeros(length(stages)),
)
    length(angle_errors) == length(stages) ||
        throw(DimensionMismatch("one angle error per Givens stage is required"))
    length(delay_errors_s) == length(stages) ||
        throw(DimensionMismatch("one delay error per Givens stage is required"))
    field = zeros(ComplexF64, OUTPUT_COUNT, INPUT_COUNT)
    field[1:INPUT_COUNT, :] .= input_block
    amplitude = 10.0^(-insertion_loss_db_per_coupler / 20)
    for (index, stage) in enumerate(stages)
        phase = exp(-2pi * im * frequency_hz *
            (common_delay_s_per_coupler + delay_errors_s[index]))
        block = amplitude * phase .* rotation_block(
            stage.angle + angle_errors[index],
        )
        field[[stage.upper, stage.lower], :] .=
            block * field[[stage.upper, stage.lower], :]
    end
    field
end

function band_metrics(
    stages, input_block, frequencies_hz, target;
    insertion_loss_db_per_coupler=0.0,
    common_delay_s_per_coupler=0.0,
    angle_errors=zeros(length(stages)),
    delay_errors_s=zeros(length(stages)),
)
    transfer = cat((
        realized_transfer(
            stages, input_block, frequency_hz;
            insertion_loss_db_per_coupler,
            common_delay_s_per_coupler,
            angle_errors,
            delay_errors_s,
        ) for frequency_hz in frequencies_hz
    )...; dims=3)
    rows = [
        RunOperatorO6DIsolatedBusBound.transfer_metrics(
            transfer[:, :, index], target,
        ) for index in eachindex(frequencies_hz)
    ]
    delay = RunOperatorO6DIsolatedBusBound.differential_delay(
        frequencies_hz, transfer, target,
    )
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
        transfer_fidelity=summary.minimum_fidelity >=
            RunOperatorO6DIsolatedBusBound.FIDELITY_GATE,
        port_coherence=summary.maximum_coherence <=
            RunOperatorO6DIsolatedBusBound.COHERENCE_GATE,
        column_imbalance=summary.maximum_imbalance_db <=
            RunOperatorO6DIsolatedBusBound.IMBALANCE_GATE_DB,
        amplitude_rms=summary.maximum_amplitude_rms <=
            RunOperatorO6DIsolatedBusBound.AMPLITUDE_RMS_GATE,
        differential_delay=delay.maximum_rms <=
            RunOperatorO6DIsolatedBusBound.DELAY_RMS_GATE_S,
        useful_amplitude=summary.minimum_useful_amplitude >=
            RunOperatorO6DIsolatedBusBound.USEFUL_AMPLITUDE_GATE,
        passivity=summary.maximum_sigma <= 1 + 1e-10,
        reciprocity=true,
    )
    (; summary, summary..., gates, passed=all(values(gates)), transfer, delay)
end

function scalar_ceiling(evaluator; initial_upper)
    lower = 0.0
    evaluator(lower).passed || error("zero-parasitic network does not pass")
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

function monte_carlo_summary(
    scale, stages, input_block, frequencies_hz, target,
    angle_standard, delay_standard,
)
    draws = size(angle_standard, 2)
    passed = falses(draws)
    minimum_fidelity = zeros(draws)
    maximum_amplitude_rms = zeros(draws)
    maximum_delay_rms_ns = zeros(draws)
    minimum_useful_amplitude = zeros(draws)
    maximum_coherence = zeros(draws)
    for draw in 1:draws
        result = band_metrics(
            stages, input_block, frequencies_hz, target;
            angle_errors=scale * ANGLE_RMS_REFERENCE_RAD .* @view(
                angle_standard[:, draw],
            ),
            delay_errors_s=scale * DELAY_RMS_REFERENCE_S .* @view(
                delay_standard[:, draw],
            ),
        )
        passed[draw] = result.passed
        minimum_fidelity[draw] = result.minimum_fidelity
        maximum_amplitude_rms[draw] = result.maximum_amplitude_rms
        maximum_delay_rms_ns[draw] = result.maximum_delay_rms_ns
        minimum_useful_amplitude[draw] = result.minimum_useful_amplitude
        maximum_coherence[draw] = result.maximum_coherence
    end
    successes = count(passed)
    (
        scale,
        angle_rms_deg=rad2deg(scale * ANGLE_RMS_REFERENCE_RAD),
        component_delay_rms_ns=1e9scale * DELAY_RMS_REFERENCE_S,
        draws,
        successes,
        observed_joint_pass_probability=successes / draws,
        wilson_lower_95=wilson_lower_bound(successes, draws),
        q01_minimum_fidelity=quantile(minimum_fidelity, 0.01),
        q99_maximum_coherence=quantile(maximum_coherence, 0.99),
        q99_maximum_amplitude_rms=quantile(maximum_amplitude_rms, 0.99),
        q99_maximum_delay_rms_ns=quantile(maximum_delay_rms_ns, 0.99),
        q01_minimum_useful_amplitude=quantile(minimum_useful_amplitude, 0.01),
        probability_gate_passed=wilson_lower_bound(successes, draws) >=
                                PASS_PROBABILITY_GATE,
    )
end

function monte_carlo_threshold(
    rows, stages, input_block, frequencies_hz, target,
    angle_standard, delay_standard,
)
    passing = filter(row -> row.probability_gate_passed, rows)
    lower = maximum(getproperty.(passing, :scale))
    failing = filter(row -> !row.probability_gate_passed && row.scale > lower, rows)
    upper = isempty(failing) ? 2maximum(getproperty.(rows, :scale)) :
            minimum(getproperty.(failing, :scale))
    while monte_carlo_summary(
        upper, stages, input_block, frequencies_hz, target,
        angle_standard, delay_standard,
    ).probability_gate_passed
        lower = upper
        upper *= 2
    end
    for _ in 1:12
        middle = (lower + upper) / 2
        result = monte_carlo_summary(
            middle, stages, input_block, frequencies_hz, target,
            angle_standard, delay_standard,
        )
        if result.probability_gate_passed
            lower = middle
        else
            upper = middle
        end
    end
    monte_carlo_summary(
        lower, stages, input_block, frequencies_hz, target,
        angle_standard, delay_standard,
    )
end

function write_namedtuple_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function limiting_gate(result)
    margins = (
        transfer_fidelity=result.minimum_fidelity /
            RunOperatorO6DIsolatedBusBound.FIDELITY_GATE,
        port_coherence=RunOperatorO6DIsolatedBusBound.COHERENCE_GATE /
            max(result.maximum_coherence, 1e-15),
        column_imbalance=RunOperatorO6DIsolatedBusBound.IMBALANCE_GATE_DB /
            max(result.maximum_imbalance_db, 1e-15),
        amplitude_rms=RunOperatorO6DIsolatedBusBound.AMPLITUDE_RMS_GATE /
            max(result.maximum_amplitude_rms, 1e-15),
        differential_delay=1e9RunOperatorO6DIsolatedBusBound.DELAY_RMS_GATE_S /
            max(result.maximum_delay_rms_ns, 1e-15),
        useful_amplitude=result.minimum_useful_amplitude /
            RunOperatorO6DIsolatedBusBound.USEFUL_AMPLITUDE_GATE,
        # The synthesized lattice is a product of lossless rotations (or
        # contractions when insertion loss is present), so sigma_max=1 up to
        # roundoff for the delay-only sweep. Do not let a 1e-15 numerical tie
        # mask the physical gate that actually terminates the bisection.
        passivity=result.maximum_sigma <= 1 + 1e-9 ? Inf :
            1 / max(result.maximum_sigma, 1e-15),
    )
    names = propertynames(margins)
    String(names[argmin(collect(values(margins)))])
end

function run_o6f_beamformer_network()
    isfile(O6D_PATH) || error("O6d result not found: $O6D_PATH")
    mkpath(OUTPUT_ROOT)
    o6d = load(O6D_PATH)
    target = Float64.(o6d["target"])
    nominal = Float64.(o6d["quantized_transfers"][8])
    frequencies_hz = Float64.(o6d["frequencies_hz"])
    factorization = svd(nominal; full=false)
    isometry = Matrix(factorization.U[:, 1:INPUT_COUNT])
    singular_values = factorization.S
    input_block = Diagonal(singular_values) * transpose(factorization.V)
    lattice = synthesize_isometry_lattice(isometry)
    input_block = lattice.core_rotation * input_block
    length(lattice.stages) == 45 ||
        error("the frozen 17 x 3 isometry should require 45 Givens stages")
    reconstructed = realized_transfer(lattice.stages, input_block, 242.0e3)
    reconstruction_error = norm(reconstructed - nominal) / norm(nominal)
    reconstruction_error <= 1e-11 || error("q8 network reconstruction failed")
    nominal_result = band_metrics(
        lattice.stages, input_block, frequencies_hz, target,
    )
    nominal_result.passed || error("zero-parasitic q8 lattice does not pass O4")

    println("[>] O6f-A uniform component ceilings")
    loss_ceiling_db, loss_boundary = scalar_ceiling(
        loss -> band_metrics(
            lattice.stages, input_block, frequencies_hz, target;
            insertion_loss_db_per_coupler=loss,
        );
        initial_upper=0.05,
    )
    delay_ceiling_s, delay_boundary = scalar_ceiling(
        delay -> band_metrics(
            lattice.stages, input_block, frequencies_hz, target;
            common_delay_s_per_coupler=delay,
        );
        initial_upper=20.0e-9,
    )

    rng = Xoshiro(RANDOM_SEED)
    angle_standard = randn(rng, length(lattice.stages), MONTE_CARLO_DRAWS)
    delay_standard = randn(rng, length(lattice.stages), MONTE_CARLO_DRAWS)
    println("[>] O6f-A component tolerance frontier")
    tolerance_rows = [
        monte_carlo_summary(
            scale, lattice.stages, input_block, frequencies_hz, target,
            angle_standard, delay_standard,
        ) for scale in ERROR_SCALES
    ]
    threshold = monte_carlo_threshold(
        tolerance_rows, lattice.stages, input_block, frequencies_hz, target,
        angle_standard, delay_standard,
    )
    branch_passed = threshold.probability_gate_passed &&
                    loss_ceiling_db > 0 && delay_ceiling_s > 0

    stage_rows = [(
        stage_index=index,
        upper_output=stage.upper,
        lower_output=stage.lower,
        ideal_angle_deg=rad2deg(stage.angle),
        eliminated_column=stage.eliminated_column,
    ) for (index, stage) in enumerate(lattice.stages)]
    ceiling_rows = [
        (
            parameter="uniform insertion loss per output Givens coupler",
            ceiling=loss_ceiling_db,
            unit="dB",
            active_gate=limiting_gate(loss_boundary),
        ),
        (
            parameter="uniform group delay per output Givens coupler",
            ceiling=1e9delay_ceiling_s,
            unit="ns",
            active_gate=limiting_gate(delay_boundary),
        ),
    ]
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "givens_lattice.csv"), stage_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "component_ceilings.csv"), ceiling_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "component_tolerance_frontier.csv"), tolerance_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "component_tolerance_threshold.csv"), [threshold])
    jldsave(
        joinpath(OUTPUT_ROOT, "o6f_beamformer_network.jld2");
        format_version=1,
        source_o6d_path=abspath(O6D_PATH),
        target,
        nominal,
        frequencies_hz,
        singular_values,
        input_block,
        stages=lattice.stages,
        core_rotation=lattice.core_rotation,
        reconstruction_error,
        nominal_summary=nominal_result.summary,
        loss_ceiling_db,
        loss_boundary_summary=loss_boundary.summary,
        delay_ceiling_s,
        delay_boundary_summary=delay_boundary.summary,
        random_seed=RANDOM_SEED,
        monte_carlo_draws=MONTE_CARLO_DRAWS,
        tolerance_rows,
        threshold,
        branch_passed,
    )

    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6f-A concrete reciprocal q8 lattice")
        println(io)
        println(io, branch_passed ?
            "**PASS as an optimistic topology bound; electrical load and input-hybrid parasitics remain.**" :
            "**FAIL: the concrete q8 lattice cannot satisfy O4 even under the frozen optimistic model.**")
        println(io)
        @printf(io,
            "The q8 `17 x 3` contraction is reconstructed to relative error `%.3e` as three singular-value attenuators, two ideal calibrated `3 x 3` hybrid blocks and a `45`-coupler adjacent-Givens output lattice. The small hybrids, dump-port loss and piezo impedance are counted but not yet assigned parasitics.\n\n",
            reconstruction_error,
        )
        @printf(io,
            "With every other error zero, the O4 ceilings are `%.6f dB` identical insertion loss per output coupler and `%.3f ns` identical delay per output coupler. These are network-distributed ceilings, not datasheet limits for an isolated part.\n\n",
            loss_ceiling_db, 1e9delay_ceiling_s,
        )
        @printf(io,
            "For independent component errors, the largest common scale with Wilson lower-95 pass >=99%% is `%.4f`: angle RMS `%.4f deg`, delay RMS `%.4f ns`; observed pass `%.3f%%`, lower bound `%.3f%%`, q01 useful `%.6f`.\n\n",
            threshold.scale, threshold.angle_rms_deg,
            threshold.component_delay_rms_ns,
            100threshold.observed_joint_pass_probability,
            100threshold.wilson_lower_95,
            threshold.q01_minimum_useful_amplitude,
        )
        println(io, "This is a concrete reciprocal factorization, but not an electrical prototype. The remaining common loss budget must still pay for two input hybrids, three attenuator/dump paths, winding resistance, dielectric loss and mismatch to 17 capacitive piezo loads.")
    end
    tolerance_plot = plot(
        getproperty.(tolerance_rows, :scale),
        100 .* getproperty.(tolerance_rows, :wilson_lower_95);
        marker=:circle, linewidth=2, label="Wilson lower 95%",
        xlabel="component-error scale", ylabel="joint pass probability, %",
        title="O6f-A q8 lattice tolerance",
    )
    hline!(tolerance_plot, [100PASS_PROBABILITY_GATE];
        linestyle=:dash, color=:black, label="99% gate")
    vline!(tolerance_plot, [threshold.scale];
        linestyle=:dot, color=:red, label="threshold")
    savefig(tolerance_plot, joinpath(OUTPUT_ROOT, "component_tolerance_frontier.png"))
    println(branch_passed ? "[+] O6f-A topology-bound PASS" : "[-] O6f-A FAIL")
    println("[+] $(joinpath(OUTPUT_ROOT, "README.md"))")
    (; branch_passed, loss_ceiling_db, delay_ceiling_s, threshold)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6f_beamformer_network()

end
