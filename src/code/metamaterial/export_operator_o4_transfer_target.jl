module ExportOperatorO4TransferTarget

using JLD2
using LinearAlgebra
using Printf
using Random
using Statistics

export export_o4_target

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const STATIC_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "static_frontier_v0",
    "static_frontier.jld2",
)
const F1C_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "fullwave_f1c_stochastic_v0",
    "f1c_stochastic.jld2",
)
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o4_k4q2_transfer_target_v0",
)
const ARCHITECTURE_ID = :trained_k4_q2_unconstrained
const FREQUENCY_HZ = 242.0e3 .* collect(range(0.8, 1.2; length=7))
const ELEMENT_COUNT = 17
const PORT_COUNT = 4
const STATE_COUNT = 2
const APERTURE_HEIGHT_MM = 71.6
const TOLERANCE_SEED = UInt64(0x04a42026)
const TOLERANCE_DRAWS = 10_000
const COEFFICIENT_AMPLITUDE_RMS = 0.03
const DIFFERENTIAL_DELAY_RMS_CYCLES = 0.02

function write_namedtuple_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function tolerance_envelope(states)
    rng = Xoshiro(TOLERANCE_SEED)
    rows = NamedTuple[]
    for _ in 1:TOLERANCE_DRAWS, state_index in eachindex(states), ratio in (0.8, 1.0, 1.2)
        amplitude = exp.(COEFFICIENT_AMPLITUDE_RMS .* randn(rng, size(states[state_index])))
        delay = DIFFERENTIAL_DELAY_RMS_CYCLES .* randn(rng, size(states[state_index]))
        raw = complex.(states[state_index]) .* amplitude .* exp.(-2pi * im * ratio .* delay)
        transfer = raw ./ max(1.0, opnorm(raw))
        gram = adjoint(transfer) * transfer
        gains = sqrt.(real.(diag(gram)))
        coherence = maximum(
            abs(gram[first, second]) / sqrt(real(gram[first, first] * gram[second, second]))
            for first in 1:PORT_COUNT, second in 1:PORT_COUNT if first != second
        )
        fidelity = abs(dot(vec(complex.(states[state_index])), vec(transfer))) /
            (norm(states[state_index]) * norm(transfer))
        push!(rows, (
            state_index,
            frequency_ratio=ratio,
            maximum_port_coherence=coherence,
            normalized_transfer_fidelity=fidelity,
            minimum_column_gain=minimum(gains),
            maximum_column_gain=maximum(gains),
            column_imbalance_db=20log10(maximum(gains) / minimum(gains)),
            sigma_before_contraction=opnorm(raw),
            sigma_after_contraction=opnorm(transfer),
        ))
    end
    summary = NamedTuple[]
    for metric in (
        :maximum_port_coherence,
        :normalized_transfer_fidelity,
        :minimum_column_gain,
        :maximum_column_gain,
        :column_imbalance_db,
        :sigma_before_contraction,
    )
        values = getproperty.(rows, metric)
        push!(summary, (
            metric,
            minimum=minimum(values),
            q01=quantile(values, 0.01),
            median=median(values),
            q99=quantile(values, 0.99),
            maximum=maximum(values),
        ))
    end
    rows, summary
end

function manufacturability_diagnostics(states)
    aperture_count = size(first(states), 1)
    basis = zeros(Float64, aperture_count, aperture_count)
    basis[:, 1] .= inv(sqrt(aperture_count))
    for mode in 2:aperture_count
        basis[:, mode] .= sqrt(2 / aperture_count) .* cos.(
            (mode - 1) .* ((1:aperture_count) .- 0.5) .* pi / aperture_count,
        )
    end
    rows = NamedTuple[]
    for state_index in eachindex(states), port_index in 1:size(states[state_index], 2)
        port = states[state_index][:, port_index]
        modal_energy = abs2.(transpose(basis) * port)
        push!(rows, (
            state_index,
            port_index,
            modes_for_95_percent=something(
                findfirst(>=(0.95), cumsum(modal_energy)), aperture_count,
            ),
            energy_first_4_modes=sum(modal_energy[1:4]),
            energy_first_8_modes=sum(modal_energy[1:8]),
            adjacent_roughness=norm(diff(port)) / norm(port),
            sign_flips=count(!=(0), diff(sign.(port))),
        ))
    end
    principal_cosines = svdvals(transpose(states[1]) * states[2])
    rows, principal_cosines, acosd.(clamp.(principal_cosines, -1, 1))
end

function export_o4_target()
    mkpath(OUTPUT_ROOT)
    states = JLD2.load(STATIC_PATH, "architecture_states")[ARCHITECTURE_ID]
    length(states) == STATE_COUNT || error("O4 state count changed")
    all(size(state) == (ELEMENT_COUNT, PORT_COUNT) for state in states) ||
        error("O4 state dimensions changed")
    all(isapprox(adjoint(state) * state, Matrix{Float64}(I, PORT_COUNT, PORT_COUNT); atol=1e-10)
        for state in states) || error("O4 state ports are no longer orthonormal")

    f1c = JLD2.load(F1C_PATH)
    summary_rows = f1c["summary_rows"]
    candidate = only(filter(row -> row.architecture_id == ARCHITECTURE_ID, summary_rows))
    candidate.f1c_gate_passed || error("O4 handoff requires a passed O3 candidate")
    reference_alpha = Float64(f1c["reference_alpha"])
    roundtrip_gain_min = candidate.alpha90 / reference_alpha
    one_way_gain_min = sqrt(roundtrip_gain_min)
    one_way_insertion_loss_db_max = -20log10(one_way_gain_min)
    roundtrip_insertion_loss_db_max = -20log10(roundtrip_gain_min)
    missing_useful_power_max = 1 - one_way_gain_min^2

    ideal_transfer = Array{ComplexF64}(undef, ELEMENT_COUNT, PORT_COUNT, length(FREQUENCY_HZ), STATE_COUNT)
    for state_index in 1:STATE_COUNT, frequency_index in eachindex(FREQUENCY_HZ)
        # The common scalar phase/delay is a calibrated gauge and is stored as zero.
        ideal_transfer[:, :, frequency_index, state_index] .= states[state_index]
    end
    pitch_mm = APERTURE_HEIGHT_MM / ELEMENT_COUNT
    y_mm = collect(range(
        -APERTURE_HEIGHT_MM / 2 + pitch_mm / 2,
        APERTURE_HEIGHT_MM / 2 - pitch_mm / 2;
        length=ELEMENT_COUNT,
    ))
    coefficient_rows = [
        (
            state_index,
            element_index,
            y_mm=y_mm[element_index],
            port_1=states[state_index][element_index, 1],
            port_2=states[state_index][element_index, 2],
            port_3=states[state_index][element_index, 3],
            port_4=states[state_index][element_index, 4],
        )
        for state_index in 1:STATE_COUNT, element_index in 1:ELEMENT_COUNT
    ]
    frequency_rows = [
        (
            anchor_index=index,
            frequency_hz=frequency,
            frequency_ratio=frequency / 242.0e3,
            ideal_common_gain=1.0,
            ideal_common_phase_deg=0.0,
            ideal_differential_group_delay_s=0.0,
        )
        for (index, frequency) in enumerate(FREQUENCY_HZ)
    ]
    delay_rms_s = DIFFERENTIAL_DELAY_RMS_CYCLES / 242.0e3
    tolerance_rows = [
        (quantity="coefficient log-amplitude error", gate="RMS <= 0.03", provenance="passed O2c/static C3 ensemble"),
        (quantity="differential coefficient group delay", gate=@sprintf("RMS <= %.6g s (0.02 cycles)", delay_rms_s), provenance="passed O2c/static C3 ensemble"),
        (quantity="carrier differential phase", gate="RMS <= 7.2 deg", provenance="0.02 cycles at 242 kHz"),
        (quantity="normalized transfer fidelity", gate=">= 0.98", provenance="Q01 of frozen 3%/0.02-cycle tolerance ensemble"),
        (quantity="maximum normalized port coherence", gate="<= 0.15 (isolation >= 16.5 dB)", provenance="above Q99=0.139 tolerance envelope"),
        (quantity="column gain imbalance", gate="<= 0.75 dB", provenance="above Q99 tolerance envelope"),
        (quantity="largest singular value", gate="<= 1 + numerical tolerance", provenance="passivity/contractivity invariant"),
        (quantity="one-way useful amplitude", gate=@sprintf(">= %.6f", one_way_gain_min), provenance="F1c alpha90/reference margin"),
        (quantity="one-way insertion loss", gate=@sprintf("<= %.3f dB", one_way_insertion_loss_db_max), provenance="F1c alpha90/reference margin"),
        (quantity="round-trip insertion loss", gate=@sprintf("<= %.3f dB", roundtrip_insertion_loss_db_max), provenance="F1c alpha90/reference margin"),
        (quantity="reflection + parasitic power", gate=@sprintf("<= %.2f%% one way", 100missing_useful_power_max), provenance="passivity plus useful-gain floor"),
        (quantity="reciprocity residual", gate="<= 3%", provenance="frozen reciprocal operator contract"),
        (quantity="defect-dictionary calibration fidelity", gate=">= 0.90", provenance="passed O2c gate"),
        (quantity="state/event energy", gate="1/sqrt(2) per state; 1/sqrt(4) per Tx basis", provenance="frozen common-energy rule"),
    ]
    envelope_rows, envelope_summary = tolerance_envelope(states)
    manufacturability_rows, principal_cosines, principal_angles_deg =
        manufacturability_diagnostics(states)

    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "state_coefficients.csv"), vec(coefficient_rows))
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "frequency_contract.csv"), frequency_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "tolerance_budget.csv"), tolerance_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "tolerance_envelope_summary.csv"), envelope_summary)
    write_namedtuple_csv(
        joinpath(OUTPUT_ROOT, "manufacturability_diagnostics.csv"),
        manufacturability_rows,
    )
    jldsave(
        joinpath(OUTPUT_ROOT, "o4_transfer_target.jld2");
        format_version=1,
        architecture_id=ARCHITECTURE_ID,
        source_static_path=abspath(STATIC_PATH),
        source_f1c_path=abspath(F1C_PATH),
        frequency_hz=FREQUENCY_HZ,
        state_matrices=states,
        ideal_transfer,
        aperture_y_mm=y_mm,
        coefficient_amplitude_rms=COEFFICIENT_AMPLITUDE_RMS,
        differential_delay_rms_cycles=DIFFERENTIAL_DELAY_RMS_CYCLES,
        differential_delay_rms_s=delay_rms_s,
        reference_alpha,
        candidate_alpha90=candidate.alpha90,
        one_way_gain_min,
        one_way_insertion_loss_db_max,
        roundtrip_insertion_loss_db_max,
        missing_useful_power_max,
        tolerance_seed=TOLERANCE_SEED,
        tolerance_draws=TOLERANCE_DRAWS,
        envelope_summary,
        manufacturability_rows,
        state_principal_cosines=principal_cosines,
        state_principal_angles_deg=principal_angles_deg,
    )

    readme_path = joinpath(OUTPUT_ROOT, "README.md")
    open(readme_path, "w") do io
        println(io, "# O4 transfer target: static K4,Q2")
        println(io)
        println(io, "> This is a realizable transfer-function contract, not a geometry and not a 17-channel electronics prescription.")
        println(io)
        println(io, "The selected O3 operator is `$(ARCHITECTURE_ID)`: four electronic ports, two sequential passive states and eight transmit events. Each state is a real orthonormal `17 x 4` map from four device ports to a discretized continuous contact aperture.")
        println(io)
        println(io, "The target is")
        println(io)
        println(io, "```text")
        println(io, "H_q(omega) = g_q(omega) exp(-i omega tau_q) [S_q + E_q(omega)], q=1,2")
        println(io, "reverse transfer = transpose(H_q)             (reciprocal)")
        println(io, "sigma_max(H_q) <= 1                           (passive)")
        println(io, "```")
        println(io)
        println(io, "`S_q` is stored in `state_coefficients.csv`; the common scalar gain and delay are calibrated gauges. Only differential errors change the operator shape. The qualified band is `193.6--290.4 kHz` (seven anchors, five-cycle 242 kHz pulse).")
        println(io)
        println(io, "## Hard acceptance numbers")
        println(io)
        println(io, @sprintf("- coefficient amplitude error: <= 3%% RMS; differential group-delay error: <= %.1f ns RMS (7.2 deg at carrier);", 1e9delay_rms_s))
        println(io, "- normalized transfer fidelity >= 0.98; maximum port coherence <= 0.15; column imbalance <= 0.75 dB;")
        println(io, @sprintf("- useful one-way amplitude >= %.4f, hence one-way insertion loss <= %.3f dB and round-trip loss <= %.3f dB;", one_way_gain_min, one_way_insertion_loss_db_max, roundtrip_insertion_loss_db_max))
        println(io, @sprintf("- reflected plus parasitic one-way power <= %.2f%%; reciprocity residual <= 3%%;", 100missing_useful_power_max))
        println(io, "- calibrated defect-dictionary fidelity >= 0.90 and the frozen F1c acceptance replay must still pass.")
        println(io)
        println(io, "The loss bound is derived from the observed O3 margin `alpha90=$(round(candidate.alpha90; digits=4))` at `reference alpha=$(reference_alpha)`. Because pulse echo traverses the network twice, the one-way amplitude floor is the square root of the round-trip margin.")
        println(io)
        println(io, "The 10,000-draw tolerance envelope is deterministic (`seed=0x04a42026`). It is a consistency envelope for the already-tested 3%/0.02-cycle perturbation law, not additional evidence or a replacement for the F1c replay.")
        println(io)
        println(io, "## Manufacturability diagnostic")
        println(io)
        println(io, "The exact unconstrained target is not a low-order lens field. Its four state-to-state principal angles are `$(join(round.(principal_angles_deg; digits=1), ", ")) deg`; individual ports need `$(minimum(getproperty.(manufacturability_rows, :modes_for_95_percent)))--$(maximum(getproperty.(manufacturability_rows, :modes_for_95_percent)))` DCT aperture modes for 95% energy and contain `$(minimum(getproperty.(manufacturability_rows, :sign_flips)))--$(maximum(getproperty.(manufacturability_rows, :sign_flips)))` sign changes. O5 must therefore test a geometry-aware smooth projection before any dense geometry optimization.")
    end
    println("[+] O4 transfer target: $readme_path")
    readme_path
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_o4_target()

end
