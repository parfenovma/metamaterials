module RunOperatorO6BStackTwin

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Statistics

include(joinpath(@__DIR__, "operator_o6b_stack_mesher.jl"))
include(joinpath(@__DIR__, "operator_o6b_stack_harmonic_solver.jl"))
include(joinpath(@__DIR__, "run_operator_o6b_direct_coupler_bound.jl"))

using .OperatorO6BStackMesher
using .OperatorO6BStackHarmonicSolver
using .RunOperatorO6BDirectCouplerBound

export run_o6b_stack_twin

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6b_stack_twin_v0",
)
const FREQUENCIES_HZ = 242.0e3 .* collect(range(0.8, 1.2; length=13))
const CARRIER_INDEX = 7
const APERTURE_COUNT = 17
const PORT_COUNT = 3

mesh_path() = joinpath(OUTPUT_ROOT, "meshes", "pzt6_bond0p1_al35.msh")

function dct_basis(count=APERTURE_COUNT)
    basis = zeros(Float64, count, PORT_COUNT)
    basis[:, 1] .= inv(sqrt(count))
    for mode in 2:PORT_COUNT
        basis[:, mode] .= sqrt(2 / count) .* cos.(
            (mode - 1) .* ((1:count) .- 0.5) .* pi / count,
        )
    end
    basis
end

function unwrap_phase(values)
    result = angle.(values)
    for index in 2:length(result)
        result[index:end] .-= 2pi * round((result[index] - result[index - 1]) / (2pi))
    end
    result
end

function derivative(values, coordinates)
    result = similar(values, Float64)
    result[1] = (values[2] - values[1]) / (coordinates[2] - coordinates[1])
    result[end] = (values[end] - values[end - 1]) /
                  (coordinates[end] - coordinates[end - 1])
    for index in 2:(length(values) - 1)
        result[index] = (values[index + 1] - values[index - 1]) /
                        (coordinates[index + 1] - coordinates[index - 1])
    end
    result
end

function transfer_metrics(transfer, target)
    gram = adjoint(transfer) * transfer
    gains = sqrt.(real.(diag(gram)))
    common = dot(vec(target), vec(transfer)) / sum(abs2, target)
    coherence = maximum(
        abs(gram[first, second]) / sqrt(real(gram[first, first] * gram[second, second]))
        for first in 1:PORT_COUNT, second in 1:PORT_COUNT if first != second
    )
    (
        common,
        useful_amplitude=abs(common),
        fidelity=abs(dot(vec(target), vec(transfer))) / (norm(target) * norm(transfer)),
        coherence,
        imbalance_db=20log10(maximum(gains) / minimum(gains)),
        amplitude_rms=norm(abs.(transfer) .- abs.(common .* target)) /
                      norm(abs.(common .* target)),
        sigma_max=opnorm(transfer),
    )
end

function differential_delay(frequencies_hz, transfer, target)
    omega = 2pi .* frequencies_hz
    active = abs.(target) .>= 0.10maximum(abs, target)
    delays = fill(NaN, size(transfer))
    for port in 1:PORT_COUNT, bin in 1:APERTURE_COUNT
        active[bin, port] || continue
        relative = transfer[bin, port, :] ./ target[bin, port]
        delays[bin, port, :] .= -derivative(unwrap_phase(relative), omega)
    end
    weights = abs2.(target)
    rms = zeros(Float64, length(frequencies_hz))
    common = similar(rms)
    for frequency_index in eachindex(frequencies_hz)
        valid = active .& isfinite.(delays[:, :, frequency_index])
        local_weights = weights[valid]
        local_delays = delays[:, :, frequency_index][valid]
        common[frequency_index] = sum(local_weights .* local_delays) / sum(local_weights)
        rms[frequency_index] = sqrt(
            sum(local_weights .* (local_delays .- common[frequency_index]) .^ 2) /
            sum(local_weights),
        )
    end
    (; delays, common, rms, maximum_rms=maximum(rms))
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

function run_o6b_stack_twin()
    mkpath(OUTPUT_ROOT)
    geometry = O6BStackGeometry()
    mesh_metadata = if isfile(mesh_path())
        (
            output_path=abspath(mesh_path()), node_count=-1, element_count=-1,
            pzt_back_x_m=-(geometry.pzt_thickness_mm + geometry.bond_thickness_mm) * 1e-3,
            contact_x_m=0.0,
            aluminium_end_x_m=geometry.aluminium_length_mm * 1e-3,
        )
    else
        build_o6b_stack_mesh(mesh_path(); config=geometry)
    end

    raw_transfer = zeros(ComplexF64, APERTURE_COUNT, PORT_COUNT, length(FREQUENCIES_HZ))
    port_work = zeros(ComplexF64, PORT_COUNT, PORT_COUNT, length(FREQUENCIES_HZ))
    reciprocity = zeros(Float64, length(FREQUENCIES_HZ))
    matrix_size = (0, 0)
    matrix_nonzeros = 0
    for (frequency_index, frequency_hz) in enumerate(FREQUENCIES_HZ)
        println("[>] O6b finite twin $(frequency_index)/$(length(FREQUENCIES_HZ)): $(frequency_hz / 1e3) kHz")
        result = solve_o6b_stack_harmonic(
            mesh_path();
            config=O6BStackHarmonicConfig(frequency_hz=frequency_hz),
        )
        raw_transfer[:, :, frequency_index] .= result.output_displacement_m
        port_work[:, :, frequency_index] .= result.port_work_compliance
        reciprocity[frequency_index] = result.reciprocity_residual
        matrix_size = result.matrix_size
        matrix_nonzeros = result.matrix_nonzeros
    end

    target = dct_basis()
    carrier_gains = [
        abs(dot(target[:, port], raw_transfer[:, port, CARRIER_INDEX]))
        for port in 1:PORT_COUNT
    ]
    reference_gain = exp(mean(log.(carrier_gains)))
    drive_scale = reference_gain ./ carrier_gains
    calibrated = raw_transfer .* reshape(drive_scale, 1, PORT_COUNT, 1)

    bound = evaluate_o6b_bound()
    bound_indices = collect(1:4:length(bound.frequencies_hz))
    maximum(abs.(bound.frequencies_hz[bound_indices] .- FREQUENCIES_HZ)) <= 1e-6 ||
        error("O6b bound and finite-twin frequency grids do not align")
    transfer = similar(calibrated)
    pzt_transit_loss = zeros(Float64, length(FREQUENCIES_HZ))
    common_raw_response = zeros(ComplexF64, length(FREQUENCIES_HZ))
    for (frequency_index, bound_index) in enumerate(bound_indices)
        common_raw = dot(
            vec(target), vec(calibrated[:, :, frequency_index]),
        ) / sum(abs2, target)
        common_raw_response[frequency_index] = common_raw
        omega = 2pi * FREQUENCIES_HZ[frequency_index]
        pzt_loss = exp(-omega * geometry.pzt_thickness_mm * 1e-3 /
            (2OperatorO6BStackHarmonicSolver.PZT.quality_factor *
             OperatorO6BStackHarmonicSolver.PZT.pressure_speed_m_s))
        bond_loss = exp(-omega * geometry.bond_thickness_mm * 1e-3 /
            (2OperatorO6BStackHarmonicSolver.BOND.quality_factor *
             OperatorO6BStackHarmonicSolver.BOND.pressure_speed_m_s))
        pzt_transit_loss[frequency_index] = pzt_loss * bond_loss
        common_power_wave = mean(bound.modal_h[:, bound_index]) *
                            pzt_transit_loss[frequency_index]
        transfer[:, :, frequency_index] .=
            calibrated[:, :, frequency_index] ./ common_raw .* common_power_wave
    end

    delay = differential_delay(FREQUENCIES_HZ, transfer, target)
    rows = NamedTuple[]
    for (index, frequency_hz) in enumerate(FREQUENCIES_HZ)
        metrics = transfer_metrics(transfer[:, :, index], target)
        push!(rows, (
            frequency_hz,
            frequency_ratio=frequency_hz / 242.0e3,
            normalized_transfer_fidelity=metrics.fidelity,
            maximum_port_coherence=metrics.coherence,
            column_imbalance_db=metrics.imbalance_db,
            differential_amplitude_rms=metrics.amplitude_rms,
            differential_delay_rms_ns=1e9delay.rms[index],
            useful_amplitude=metrics.useful_amplitude,
            sigma_max=metrics.sigma_max,
            mechanical_reciprocity_residual=reciprocity[index],
            pzt_bond_transit_loss=pzt_transit_loss[index],
        ))
    end
    gates = (
        transfer_fidelity=minimum(getproperty.(rows, :normalized_transfer_fidelity)) >= 0.98,
        port_coherence=maximum(getproperty.(rows, :maximum_port_coherence)) <= 0.15,
        column_imbalance=maximum(getproperty.(rows, :column_imbalance_db)) <= 0.75,
        amplitude_rms=maximum(getproperty.(rows, :differential_amplitude_rms)) <= 0.03,
        differential_delay=delay.maximum_rms <= 82.6e-9,
        useful_amplitude=minimum(getproperty.(rows, :useful_amplitude)) >= 0.8994,
        passivity=maximum(getproperty.(rows, :sigma_max)) <= 1 + 1e-10,
        reciprocity=maximum(reciprocity) <= 0.03,
    )
    passed = all(values(gates))
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "dense_band_metrics.csv"), rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "gates.csv"), [
        (gate=String(name), passed=value) for (name, value) in pairs(gates)
    ])
    jldsave(
        joinpath(OUTPUT_ROOT, "o6b_stack_twin.jld2");
        format_version=1,
        frequencies_hz=FREQUENCIES_HZ,
        geometry,
        mesh_metadata,
        target,
        raw_transfer,
        calibrated_transfer=calibrated,
        as_built_transfer=transfer,
        port_work_compliance=port_work,
        reciprocity,
        drive_scale,
        common_raw_response,
        pzt_transit_loss,
        differential_delay_s=delay.delays,
        rows,
        gates,
        passed,
        matrix_size,
        matrix_nonzeros,
    )
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6b finite PZT/bond/Al mechanical twin")
        println(io)
        println(io, passed ?
            "**PASS: the as-built H(omega) is authorized for frozen F1d replay.**" :
            "**FAIL: the finite direct-coupler twin does not satisfy O4.**")
        println(io)
        println(io, "One frozen `6.0 mm PZT-5A mechanical surrogate -> 0.10 mm epoxy -> 35 mm matched Al load` geometry was solved at 13 frequencies. The three right-hand sides are continuous DCT0--2 tractions on the common PZT back face; the 17-bin output is sampled at the bonded Al interface. No length, thickness or profile sweep was performed.")
        println(io)
        println(io, "## Worst dense-band metrics")
        println(io)
        @printf(io, "- transfer fidelity: `%.6f`;\n", minimum(getproperty.(rows, :normalized_transfer_fidelity)))
        @printf(io, "- maximum port coherence: `%.6f`;\n", maximum(getproperty.(rows, :maximum_port_coherence)))
        @printf(io, "- column imbalance: `%.4f dB`;\n", maximum(getproperty.(rows, :column_imbalance_db)))
        @printf(io, "- differential amplitude RMS: `%.3f%%`;\n", 100maximum(getproperty.(rows, :differential_amplitude_rms)))
        @printf(io, "- differential-delay RMS: `%.3f ns`;\n", 1e9delay.maximum_rms)
        @printf(io, "- minimum useful amplitude: `%.6f`;\n", minimum(getproperty.(rows, :useful_amplitude)))
        @printf(io, "- maximum sigma: `%.6f`;\n", maximum(getproperty.(rows, :sigma_max)))
        @printf(io, "- mechanical reciprocity residual: `%.3e`.\n\n", maximum(reciprocity))
        println(io, "The common power-wave scale is inherited from the preregistered O6b PZT/bond/Al bound and includes nominal PZT/bond propagation loss. The FEM contributes only the finite-aperture differential shape, phase and cross-port terms; common displacement compliance is a calibrated gauge.")
        println(io)
        println(io, "This remains a mechanical distributed-port twin, not an electrical terminal model. A pass authorizes F1d replay first; electrode topology, capacitance and receive electronics remain blocked until that replay.")
    end
    metric_plot = plot(
        FREQUENCIES_HZ ./ 1e3,
        getproperty.(rows, :normalized_transfer_fidelity);
        linewidth=2, marker=:circle, label="fidelity",
        xlabel="frequency, kHz", ylabel="metric", title="O6b finite twin",
    )
    plot!(metric_plot, FREQUENCIES_HZ ./ 1e3,
        getproperty.(rows, :useful_amplitude); linewidth=2, label="useful amplitude")
    hline!(metric_plot, [0.98, 0.8994]; linestyle=:dash, color=:black, label="gates")
    savefig(metric_plot, joinpath(OUTPUT_ROOT, "finite_twin_metrics.png"))
    println(passed ? "[+] O6b finite twin PASS" : "[-] O6b finite twin FAIL")
    println("[+] $(joinpath(OUTPUT_ROOT, "README.md"))")
    (; passed, gates, rows, transfer)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6b_stack_twin()

end
