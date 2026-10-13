module RunOperatorO6CCompositeTwin

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Statistics

include(joinpath(@__DIR__, "operator_o6b_stack_harmonic_solver.jl"))
include(joinpath(@__DIR__, "run_operator_o6b_direct_coupler_bound.jl"))

using .OperatorO6BStackHarmonicSolver
using .RunOperatorO6BDirectCouplerBound

export composite_material, composite_metadata, run_o6c_composite_twin

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6c_composite_twin_v0",
)
const O6B_MESH_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6b_stack_twin_v0",
    "meshes", "pzt6_bond0p1_al35.msh",
)
const FREQUENCIES_HZ = 242.0e3 .* collect(range(0.8, 1.2; length=13))
const CARRIER_INDEX = 7
const APERTURE_COUNT = 17
const PORT_COUNT = 3
const CERAMIC_VOLUME_FRACTION = 0.65
const PZT_THICKNESS_M = 6.0e-3
const BOND_THICKNESS_M = 100.0e-6
const AXIAL_IMPEDANCE_RAYL = 17.0e6
const COMPOSITE_QUALITY_FACTOR = 35.0

const CERAMIC = (
    density_kg_m3=7750.0,
    pressure_speed_m_s=4350.0,
    shear_speed_m_s=2260.0,
)
const POLYMER = (
    density_kg_m3=1210.0,
    pressure_speed_m_s=2340.0,
    shear_speed_m_s=1170.0,
)

function isotropic_stiffness(material)
    shear = material.density_kg_m3 * material.shear_speed_m_s^2
    lambda = material.density_kg_m3 * material.pressure_speed_m_s^2 - 2shear
    (; lambda, shear, longitudinal=lambda + 2shear)
end

harmonic_mix(first, second, fraction) = inv(fraction / first + (1 - fraction) / second)

function composite_metadata()
    fraction = CERAMIC_VOLUME_FRACTION
    ceramic_stiffness = isotropic_stiffness(CERAMIC)
    polymer_stiffness = isotropic_stiffness(POLYMER)
    density = fraction * CERAMIC.density_kg_m3 +
              (1 - fraction) * POLYMER.density_kg_m3
    axial = AXIAL_IMPEDANCE_RAYL^2 / density
    lateral = harmonic_mix(
        ceramic_stiffness.longitudinal,
        polymer_stiffness.longitudinal,
        fraction,
    )
    cross = harmonic_mix(
        ceramic_stiffness.lambda,
        polymer_stiffness.lambda,
        fraction,
    )
    shear = harmonic_mix(
        ceramic_stiffness.shear,
        polymer_stiffness.shear,
        fraction,
    )
    axial_speed = sqrt(axial / density)
    lateral_speed = sqrt(lateral / density)
    shear_speed = sqrt(shear / density)
    (
        ceramic_volume_fraction=fraction,
        polymer_volume_fraction=1 - fraction,
        density_kg_m3=density,
        axial_stiffness_pa=axial,
        lateral_stiffness_pa=lateral,
        cross_stiffness_pa=cross,
        shear_stiffness_pa=shear,
        axial_impedance_rayl=AXIAL_IMPEDANCE_RAYL,
        axial_speed_m_s=axial_speed,
        lateral_speed_m_s=lateral_speed,
        shear_speed_m_s=shear_speed,
        quality_factor=COMPOSITE_QUALITY_FACTOR,
        positive_definite=axial > 0 && lateral > 0 && shear > 0 &&
                          axial * lateral - cross^2 > 0,
        plane_strain_determinant_pa2=axial * lateral - cross^2,
        nominal_half_wave_frequency_hz=axial_speed / (2PZT_THICKNESS_M),
    )
end

function composite_material()
    metadata = composite_metadata()
    metadata.positive_definite || error("frozen O6c tensor is not positive definite")
    O6COrthotropicMaterial(
        density_kg_m3=metadata.density_kg_m3,
        axial_stiffness_pa=metadata.axial_stiffness_pa,
        lateral_stiffness_pa=metadata.lateral_stiffness_pa,
        cross_stiffness_pa=metadata.cross_stiffness_pa,
        shear_stiffness_pa=metadata.shear_stiffness_pa,
        quality_factor=metadata.quality_factor,
    )
end

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
        result[index:end] .-= 2pi * round(
            (result[index] - result[index - 1]) / (2pi),
        )
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
        abs(gram[first, second]) /
        sqrt(real(gram[first, first] * gram[second, second]))
        for first in 1:PORT_COUNT, second in 1:PORT_COUNT if first != second
    )
    (
        common,
        useful_amplitude=abs(common),
        fidelity=abs(dot(vec(target), vec(transfer))) /
                 (norm(target) * norm(transfer)),
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
        common[frequency_index] = sum(local_weights .* local_delays) /
                                  sum(local_weights)
        rms[frequency_index] = sqrt(
            sum(local_weights .* (local_delays .- common[frequency_index]) .^ 2) /
            sum(local_weights),
        )
    end
    (; delays, common, rms, maximum_rms=maximum(rms))
end

function per_port_diagnostics(transfer, target, differential_delays)
    rows = NamedTuple[]
    active = abs.(target) .>= 0.10maximum(abs, target)
    for port in 1:PORT_COUNT
        fidelities = Float64[]
        amplitude_rms = Float64[]
        delay_rms = Float64[]
        for frequency_index in eachindex(FREQUENCIES_HZ)
            column = transfer[:, port, frequency_index]
            target_column = target[:, port]
            common = dot(target_column, column) / sum(abs2, target_column)
            push!(fidelities,
                abs(dot(target_column, column)) / (norm(target_column) * norm(column)))
            push!(amplitude_rms,
                norm(abs.(column) .- abs.(common .* target_column)) /
                norm(abs.(common .* target_column)))
            valid = active[:, port]
            weights = abs2.(target[valid, port])
            local_delays = differential_delays[valid, port, frequency_index]
            local_common = sum(weights .* local_delays) / sum(weights)
            push!(delay_rms, sqrt(
                sum(weights .* (local_delays .- local_common) .^ 2) / sum(weights),
            ))
        end
        push!(rows, (
            port="DCT$(port - 1)",
            worst_fidelity=minimum(fidelities),
            worst_amplitude_rms=maximum(amplitude_rms),
            worst_delay_rms_ns=1e9maximum(delay_rms),
        ))
    end
    rows
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

function write_material_csv(path, metadata)
    open(path, "w") do io
        println(io, "quantity,value")
        for name in propertynames(metadata)
            println(io, name, ',', getproperty(metadata, name))
        end
    end
end

function evaluate_as_built(raw_transfer, reciprocity, metadata)
    target = dct_basis()
    carrier_gains = [
        abs(dot(target[:, port], raw_transfer[:, port, CARRIER_INDEX]))
        for port in 1:PORT_COUNT
    ]
    reference_gain = exp(mean(log.(carrier_gains)))
    drive_scale = reference_gain ./ carrier_gains
    calibrated = raw_transfer .* reshape(drive_scale, 1, PORT_COUNT, 1)

    bound_config = O6BConfig(
        pzt_density_kg_m3=metadata.density_kg_m3,
        pzt_pressure_speed_m_s=metadata.axial_speed_m_s,
    )
    bound = evaluate_o6b_bound(bound_config)
    bound_indices = collect(1:4:length(bound.frequencies_hz))
    maximum(abs.(bound.frequencies_hz[bound_indices] .- FREQUENCIES_HZ)) <= 1e-6 ||
        error("O6c bound and finite-twin frequency grids do not align")
    transfer = similar(calibrated)
    transit_loss = zeros(Float64, length(FREQUENCIES_HZ))
    common_raw_response = zeros(ComplexF64, length(FREQUENCIES_HZ))
    for (frequency_index, bound_index) in enumerate(bound_indices)
        common_raw = dot(vec(target), vec(calibrated[:, :, frequency_index])) /
                     sum(abs2, target)
        common_raw_response[frequency_index] = common_raw
        omega = 2pi * FREQUENCIES_HZ[frequency_index]
        composite_loss = exp(-omega * PZT_THICKNESS_M /
            (2metadata.quality_factor * metadata.axial_speed_m_s))
        bond_loss = exp(-omega * BOND_THICKNESS_M /
            (2OperatorO6BStackHarmonicSolver.BOND.quality_factor *
             OperatorO6BStackHarmonicSolver.BOND.pressure_speed_m_s))
        transit_loss[frequency_index] = composite_loss * bond_loss
        common_power_wave = mean(bound.modal_h[:, bound_index]) *
                            transit_loss[frequency_index]
        transfer[:, :, frequency_index] .= calibrated[:, :, frequency_index] ./
                                            common_raw .* common_power_wave
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
            composite_bond_transit_loss=transit_loss[index],
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
    (
        target, calibrated, transfer, drive_scale, common_raw_response,
        transit_loss, delay, rows, gates, passed=all(values(gates)), bound,
    )
end

function run_o6c_composite_twin()
    isfile(O6B_MESH_PATH) || error("frozen O6b mesh not found: $O6B_MESH_PATH")
    mkpath(OUTPUT_ROOT)
    metadata = composite_metadata()
    material = composite_material()
    raw_transfer = zeros(
        ComplexF64, APERTURE_COUNT, PORT_COUNT, length(FREQUENCIES_HZ),
    )
    port_work = zeros(ComplexF64, PORT_COUNT, PORT_COUNT, length(FREQUENCIES_HZ))
    reciprocity = zeros(Float64, length(FREQUENCIES_HZ))
    matrix_size = (0, 0)
    matrix_nonzeros = 0
    for (frequency_index, frequency_hz) in enumerate(FREQUENCIES_HZ)
        println(
            "[>] O6c composite twin $(frequency_index)/$(length(FREQUENCIES_HZ)): " *
            "$(frequency_hz / 1e3) kHz",
        )
        result = solve_o6b_stack_harmonic(
            O6B_MESH_PATH;
            config=O6BStackHarmonicConfig(frequency_hz=frequency_hz),
            pzt_material=material,
        )
        raw_transfer[:, :, frequency_index] .= result.output_displacement_m
        port_work[:, :, frequency_index] .= result.port_work_compliance
        reciprocity[frequency_index] = result.reciprocity_residual
        matrix_size = result.matrix_size
        matrix_nonzeros = result.matrix_nonzeros
    end

    evaluated = evaluate_as_built(raw_transfer, reciprocity, metadata)
    port_diagnostics = per_port_diagnostics(
        evaluated.transfer, evaluated.target, evaluated.delay.delays,
    )
    write_material_csv(joinpath(OUTPUT_ROOT, "frozen_composite_tensor.csv"), metadata)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "dense_band_metrics.csv"), evaluated.rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "per_port_diagnostics.csv"), port_diagnostics)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "gates.csv"), [
        (gate=String(name), passed=value) for (name, value) in pairs(evaluated.gates)
    ])
    jldsave(
        joinpath(OUTPUT_ROOT, "o6c_composite_twin.jld2");
        format_version=1,
        frequencies_hz=FREQUENCIES_HZ,
        mesh_path=abspath(O6B_MESH_PATH),
        pzt_thickness_m=PZT_THICKNESS_M,
        bond_thickness_m=BOND_THICKNESS_M,
        composite_metadata=metadata,
        target=evaluated.target,
        raw_transfer,
        calibrated_transfer=evaluated.calibrated,
        as_built_transfer=evaluated.transfer,
        port_work_compliance=port_work,
        reciprocity,
        drive_scale=evaluated.drive_scale,
        common_raw_response=evaluated.common_raw_response,
        composite_bond_transit_loss=evaluated.transit_loss,
        differential_delay_s=evaluated.delay.delays,
        rows=evaluated.rows,
        port_diagnostics,
        gates=evaluated.gates,
        passed=evaluated.passed,
        matrix_size,
        matrix_nonzeros,
    )

    rows = evaluated.rows
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6c homogenized 1-3 composite finite twin")
        println(io)
        println(io, evaluated.passed ?
            "**PASS: mesh anchors and one explicit-cell plausibility check are authorized.**" :
            "**FAIL: the fixed-thickness homogenized 1-3 mechanism is stopped.**")
        println(io)
        @printf(io,
            "One frozen `65%% ceramic` orthotropic effective medium was substituted into the unchanged `6.0 mm composite -> 0.10 mm epoxy -> 35 mm Al` O6b mesh. Its density is `%.1f kg/m^3`, axial impedance `%.2f MRayl`, and tensor `(C11,C22,C12,C66)=(%.3f, %.3f, %.3f, %.3f) GPa`. No composition, tensor, thickness or geometry sweep was performed.\n\n",
            metadata.density_kg_m3, metadata.axial_impedance_rayl / 1e6,
            metadata.axial_stiffness_pa / 1e9, metadata.lateral_stiffness_pa / 1e9,
            metadata.cross_stiffness_pa / 1e9, metadata.shear_stiffness_pa / 1e9,
        )
        println(io, "## Worst dense-band metrics")
        println(io)
        @printf(io, "- transfer fidelity: `%.6f`;\n", minimum(getproperty.(rows, :normalized_transfer_fidelity)))
        @printf(io, "- maximum port coherence: `%.6f`;\n", maximum(getproperty.(rows, :maximum_port_coherence)))
        @printf(io, "- column imbalance: `%.4f dB`;\n", maximum(getproperty.(rows, :column_imbalance_db)))
        @printf(io, "- differential amplitude RMS: `%.3f%%`;\n", 100maximum(getproperty.(rows, :differential_amplitude_rms)))
        @printf(io, "- differential-delay RMS: `%.3f ns`;\n", 1e9evaluated.delay.maximum_rms)
        @printf(io, "- minimum useful amplitude: `%.6f`;\n", minimum(getproperty.(rows, :useful_amplitude)))
        @printf(io, "- maximum sigma: `%.6f`;\n", maximum(getproperty.(rows, :sigma_max)))
        @printf(io, "- mechanical reciprocity residual: `%.3e`.\n\n", maximum(reciprocity))
        println(io, "Compared with homogeneous O6b, the composite raises worst fidelity from `0.97331` to `0.99383`, reduces differential amplitude RMS from `19.59%` to `9.02%`, and reduces differential-delay RMS from `846 ns` to `641 ns`. The direction is correct, but the remaining misses are still factors of `3.01` and `7.76` beyond their gates. DCT2 is the worst individual port (`10.65%`, `784 ns`).")
        println(io)
        @printf(io,
            "The unchanged `6 mm` layer has a nominal axial half-wave frequency of `%.2f kHz`, inside the qualified band; that resonance was not retuned away. The common compliance and common delay remain calibrated gauges exactly as in O6b; FEM supplies the differential finite-aperture transfer and reciprocity.\n\n",
            metadata.nominal_half_wave_frequency_hz / 1e3,
        )
        println(io, evaluated.passed ?
            "This is not yet a realizable device claim. The next required checks are three frequency anchors on the preregistered refined mesh, then one pillar/pitch homogenization plausibility point, and only then frozen F1d replay." :
            "Per the preregistered stopping rule, no fill-fraction/thickness sweep, explicit-cell geometry, electrical terminal model or F1d replay is authorized from this result.")
    end
    metric_plot = plot(
        FREQUENCIES_HZ ./ 1e3,
        getproperty.(rows, :normalized_transfer_fidelity);
        linewidth=2, marker=:circle, label="fidelity",
        xlabel="frequency, kHz", ylabel="metric", title="O6c composite twin",
    )
    plot!(metric_plot, FREQUENCIES_HZ ./ 1e3,
        getproperty.(rows, :useful_amplitude); linewidth=2, label="useful amplitude")
    hline!(metric_plot, [0.98, 0.8994]; linestyle=:dash, color=:black, label="gates")
    savefig(metric_plot, joinpath(OUTPUT_ROOT, "finite_twin_metrics.png"))
    println(evaluated.passed ? "[+] O6c composite twin PASS" : "[-] O6c composite twin FAIL")
    println("[+] $(joinpath(OUTPUT_ROOT, "README.md"))")
    (; evaluated..., metadata, raw_transfer, port_work, reciprocity)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6c_composite_twin()

end
