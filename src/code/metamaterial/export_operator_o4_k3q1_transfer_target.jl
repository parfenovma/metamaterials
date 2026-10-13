module ExportOperatorO4K3Q1TransferTarget

using JLD2
using LinearAlgebra
using Printf
using Random
using Statistics

export export_o4_target

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const F1D_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "fullwave_f1d_expanded_v0",
    "f1d_expanded.jld2",
)
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o4_k3q1_transfer_target_v0",
)
const ARCHITECTURE_ID = :o5a_dct_k3_q1
const CARRIER_HZ = 242.0e3
const FREQUENCY_HZ = CARRIER_HZ .* collect(range(0.8, 1.2; length=7))
const ELEMENT_COUNT = 17
const PORT_COUNT = 3
const APERTURE_HEIGHT_MM = 71.6
const COEFFICIENT_AMPLITUDE_RMS = 0.03
const DIFFERENTIAL_DELAY_RMS_CYCLES = 0.02
const TOLERANCE_SEED = UInt64(0x04a43101)
const TOLERANCE_DRAWS = 10_000

function dct_basis(count)
    basis = zeros(Float64, count, count)
    basis[:, 1] .= inv(sqrt(count))
    for mode in 2:count
        basis[:, mode] .= sqrt(2 / count) .* cos.(
            (mode - 1) .* ((1:count) .- 0.5) .* pi / count,
        )
    end
    basis
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

function tolerance_envelope(state)
    rng = Xoshiro(TOLERANCE_SEED)
    rows = NamedTuple[]
    for _ in 1:TOLERANCE_DRAWS, ratio in (0.8, 1.0, 1.2)
        amplitude = exp.(COEFFICIENT_AMPLITUDE_RMS .* randn(rng, size(state)))
        delay = DIFFERENTIAL_DELAY_RMS_CYCLES .* randn(rng, size(state))
        raw = complex.(state) .* amplitude .* exp.(-2pi * im * ratio .* delay)
        transfer = raw ./ max(1.0, opnorm(raw))
        gram = adjoint(transfer) * transfer
        gains = sqrt.(real.(diag(gram)))
        coherence = maximum(
            abs(gram[first, second]) / sqrt(real(gram[first, first] * gram[second, second]))
            for first in 1:PORT_COUNT, second in 1:PORT_COUNT if first != second
        )
        fidelity = abs(dot(vec(complex.(state)), vec(transfer))) /
            (norm(state) * norm(transfer))
        push!(rows, (
            frequency_ratio=ratio,
            maximum_port_coherence=coherence,
            normalized_transfer_fidelity=fidelity,
            minimum_column_gain=minimum(gains),
            maximum_column_gain=maximum(gains),
            column_imbalance_db=20log10(maximum(gains) / minimum(gains)),
            sigma_before_contraction=opnorm(raw),
        ))
    end
    summary = NamedTuple[]
    for metric in (
        :maximum_port_coherence, :normalized_transfer_fidelity,
        :minimum_column_gain, :maximum_column_gain, :column_imbalance_db,
        :sigma_before_contraction,
    )
        values = getproperty.(rows, metric)
        push!(summary, (
            metric,
            minimum=minimum(values), q01=quantile(values, 0.01),
            median=median(values), q99=quantile(values, 0.99),
            maximum=maximum(values),
        ))
    end
    summary
end

function export_o4_target()
    mkpath(OUTPUT_ROOT)
    state = dct_basis(ELEMENT_COUNT)[:, 1:PORT_COUNT]
    isapprox(transpose(state) * state, Matrix{Float64}(I, PORT_COUNT, PORT_COUNT); atol=1e-12) ||
        error("K3,Q1 target lost orthonormality")
    f1d = JLD2.load(F1D_PATH)
    selected_id = f1d["selected_id"]
    selected_id == ARCHITECTURE_ID || error("F1d did not select K3,Q1")
    candidate = only(filter(
        row -> row.architecture_id == ARCHITECTURE_ID, f1d["summary_rows"],
    ))
    candidate.f1d_gate_passed || error("O4 requires passed F1d candidate")
    reference_alpha = Float64(f1d["reference_alpha"])
    roundtrip_gain_min = candidate.alpha90 / reference_alpha
    one_way_gain_min = sqrt(roundtrip_gain_min)
    one_way_insertion_loss_db_max = -20log10(one_way_gain_min)
    roundtrip_insertion_loss_db_max = -20log10(roundtrip_gain_min)
    missing_useful_power_max = 1 - one_way_gain_min^2
    delay_rms_s = DIFFERENTIAL_DELAY_RMS_CYCLES / CARRIER_HZ
    envelope_summary = tolerance_envelope(state)

    pitch_mm = APERTURE_HEIGHT_MM / ELEMENT_COUNT
    y_mm = collect(range(
        -APERTURE_HEIGHT_MM / 2 + pitch_mm / 2,
        APERTURE_HEIGHT_MM / 2 - pitch_mm / 2;
        length=ELEMENT_COUNT,
    ))
    coefficient_rows = [
        (
            element_index,
            y_mm=y_mm[element_index],
            uniform_even=state[element_index, 1],
            first_odd=state[element_index, 2],
            second_even=state[element_index, 3],
        )
        for element_index in 1:ELEMENT_COUNT
    ]
    frequency_rows = [
        (
            anchor_index=index,
            frequency_hz=frequency,
            frequency_ratio=frequency / CARRIER_HZ,
            ideal_common_gain=1.0,
            ideal_common_phase_deg=0.0,
            ideal_differential_group_delay_s=0.0,
        )
        for (index, frequency) in enumerate(FREQUENCY_HZ)
    ]
    tolerance_rows = [
        (quantity="coefficient log-amplitude error", gate="RMS <= 0.03", provenance="frozen O4 acceptance law"),
        (quantity="differential coefficient group delay", gate=@sprintf("RMS <= %.6g s (0.02 cycles)", delay_rms_s), provenance="causal robustness law"),
        (quantity="carrier differential phase", gate="RMS <= 7.2 deg", provenance="0.02 cycles at 242 kHz"),
        (quantity="normalized transfer fidelity", gate=">= 0.98", provenance="below Q01 tolerance envelope"),
        (quantity="maximum normalized port coherence", gate="<= 0.15 (isolation >= 16.5 dB)", provenance="above Q99 tolerance envelope"),
        (quantity="column gain imbalance", gate="<= 0.75 dB", provenance="above Q99 tolerance envelope"),
        (quantity="largest singular value", gate="<= 1 + numerical tolerance", provenance="passivity/contractivity invariant"),
        (quantity="one-way useful amplitude", gate=@sprintf(">= %.6f", one_way_gain_min), provenance="F1d alpha90/reference margin"),
        (quantity="one-way insertion loss without drive compensation", gate=@sprintf("<= %.3f dB", one_way_insertion_loss_db_max), provenance="F1d alpha90/reference margin"),
        (quantity="round-trip insertion loss without drive compensation", gate=@sprintf("<= %.3f dB", roundtrip_insertion_loss_db_max), provenance="F1d alpha90/reference margin"),
        (quantity="reflection + parasitic power without compensation", gate=@sprintf("<= %.2f%% one way", 100missing_useful_power_max), provenance="passivity plus useful-gain floor"),
        (quantity="reciprocity residual", gate="<= 3%", provenance="frozen reciprocal operator contract"),
        (quantity="defect-dictionary calibration fidelity", gate=">= 0.90", provenance="causal O2c gate"),
        (quantity="state/event energy", gate="one state; 1/sqrt(3) per Tx basis", provenance="frozen common-energy rule"),
    ]
    ideal_transfer = repeat(
        reshape(complex.(state), ELEMENT_COUNT, PORT_COUNT, 1),
        1, 1, length(FREQUENCY_HZ),
    )
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "state_coefficients.csv"), coefficient_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "frequency_contract.csv"), frequency_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "tolerance_budget.csv"), tolerance_rows)
    write_namedtuple_csv(
        joinpath(OUTPUT_ROOT, "tolerance_envelope_summary.csv"), envelope_summary,
    )
    jldsave(
        joinpath(OUTPUT_ROOT, "o4_transfer_target.jld2");
        format_version=1,
        architecture_id=ARCHITECTURE_ID,
        source_f1d_path=abspath(F1D_PATH),
        frequency_hz=FREQUENCY_HZ,
        state_matrix=state,
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
    )
    readme_path = joinpath(OUTPUT_ROOT, "README.md")
    open(readme_path, "w") do io
        println(io, "# O4 transfer target: smooth K3,Q1")
        println(io)
        println(io, "> Final operator handoff after expanded F1d physical-context calibration.")
        println(io)
        println(io, "The target has three fixed reciprocal modal ports and one passive state. The aperture maps are the first three orthonormal DCT modes: uniform even, first odd and second even. There is no element-wise phase control and no state switch.")
        println(io)
        println(io, "```text")
        println(io, "H(omega) = g(omega) exp(-i omega tau0) [S_DCT3 + E(omega)]")
        println(io, "reverse transfer = transpose(H); sigma_max(H) <= 1")
        println(io, "```")
        println(io)
        println(io, "Qualified band: `193.6--290.4 kHz`. Common gain/delay are calibrated gauges; differential shape errors are acceptance variables.")
        println(io)
        println(io, "## Hard acceptance numbers")
        println(io)
        println(io, @sprintf("- coefficient amplitude RMS <= 3%%; differential delay RMS <= %.1f ns (7.2 deg at carrier);", 1e9delay_rms_s))
        println(io, "- transfer fidelity >= 0.98; port coherence <= 0.15; column imbalance <= 0.75 dB; reciprocity residual <= 3%;")
        println(io, @sprintf("- with no drive/noise compensation: one-way useful amplitude >= %.4f, insertion loss <= %.3f dB one way and <= %.3f dB round trip;", one_way_gain_min, one_way_insertion_loss_db_max, roundtrip_insertion_loss_db_max))
        println(io, "- as-built calibrated H(omega) must replay the expanded F1d gates; coefficient tolerances alone are not a substitute.")
        println(io)
        println(io, "The strict loss floor assumes equal device-side energy. Extra pulser headroom is allowed only when its cost and output-referred noise penalty are included explicitly in the resource vector.")
    end
    println("[+] final O4 transfer target: $readme_path")
    readme_path
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && export_o4_target()

end
