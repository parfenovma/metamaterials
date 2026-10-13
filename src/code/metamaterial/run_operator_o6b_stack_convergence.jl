module RunOperatorO6BStackConvergence

using JLD2
using LinearAlgebra
using Printf
using Statistics

include(joinpath(@__DIR__, "run_operator_o6b_stack_twin.jl"))
using .RunOperatorO6BStackTwin

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6b_stack_twin_v0",
)
const ANCHORS_HZ = [193.6e3, 242.0e3, 290.4e3]

function write_namedtuple_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function normalized_transfer(raw, geometry)
    target = RunOperatorO6BStackTwin.dct_basis()
    carrier_gains = [
        abs(dot(target[:, port], raw[:, port, 2])) for port in 1:3
    ]
    reference_gain = exp(mean(log.(carrier_gains)))
    drive_scale = reference_gain ./ carrier_gains
    calibrated = raw .* reshape(drive_scale, 1, 3, 1)
    bound = RunOperatorO6BStackTwin.evaluate_o6b_bound()
    bound_indices = [1, 25, 49]
    transfer = similar(raw)
    for (frequency_index, bound_index) in enumerate(bound_indices)
        common_raw = dot(vec(target), vec(calibrated[:, :, frequency_index])) /
                     sum(abs2, target)
        omega = 2pi * ANCHORS_HZ[frequency_index]
        pzt = RunOperatorO6BStackTwin.OperatorO6BStackHarmonicSolver.PZT
        bond = RunOperatorO6BStackTwin.OperatorO6BStackHarmonicSolver.BOND
        transit_loss = exp(
            -omega * geometry.pzt_thickness_mm * 1e-3 /
                (2pzt.quality_factor * pzt.pressure_speed_m_s) -
            omega * geometry.bond_thickness_mm * 1e-3 /
                (2bond.quality_factor * bond.pressure_speed_m_s),
        )
        common_power_wave = mean(bound.modal_h[:, bound_index]) * transit_loss
        transfer[:, :, frequency_index] .= calibrated[:, :, frequency_index] ./
            common_raw .* common_power_wave
    end
    transfer, target
end

function resolution_rows(label, transfer, target, reciprocity)
    delay = RunOperatorO6BStackTwin.differential_delay(ANCHORS_HZ, transfer, target)
    [
        begin
            metrics = RunOperatorO6BStackTwin.transfer_metrics(
                transfer[:, :, index], target,
            )
            (
                resolution=label,
                frequency_hz=ANCHORS_HZ[index],
                normalized_transfer_fidelity=metrics.fidelity,
                maximum_port_coherence=metrics.coherence,
                column_imbalance_db=metrics.imbalance_db,
                differential_amplitude_rms=metrics.amplitude_rms,
                differential_delay_rms_ns=1e9delay.rms[index],
                useful_amplitude=metrics.useful_amplitude,
                sigma_max=metrics.sigma_max,
                reciprocity_residual=reciprocity[index],
            )
        end for index in eachindex(ANCHORS_HZ)
    ]
end

function run_convergence()
    coarse = JLD2.load(joinpath(OUTPUT_ROOT, "o6b_stack_twin.jld2"))
    coarse_indices = [1, 7, 13]
    coarse_transfer = coarse["as_built_transfer"][:, :, coarse_indices]
    coarse_reciprocity = coarse["reciprocity"][coarse_indices]
    target = coarse["target"]
    rows = resolution_rows("base_pzt0p60_al1p00", coarse_transfer, target, coarse_reciprocity)

    refined_geometry = RunOperatorO6BStackTwin.O6BStackGeometry(
        pzt_mesh_mm=0.45,
        bond_mesh_mm=0.10,
        aluminium_mesh_mm=0.75,
    )
    refined_path = joinpath(OUTPUT_ROOT, "meshes", "pzt6_bond0p1_al35_refined.msh")
    mesh_metadata = RunOperatorO6BStackTwin.build_o6b_stack_mesh(
        refined_path; config=refined_geometry,
    )
    raw = zeros(ComplexF64, 17, 3, 3)
    reciprocity = zeros(Float64, 3)
    for (index, frequency_hz) in enumerate(ANCHORS_HZ)
        println("[>] O6b refined convergence $(index)/3: $(frequency_hz / 1e3) kHz")
        result = RunOperatorO6BStackTwin.solve_o6b_stack_harmonic(
            refined_path;
            config=RunOperatorO6BStackTwin.O6BStackHarmonicConfig(
                frequency_hz=frequency_hz,
            ),
        )
        raw[:, :, index] .= result.output_displacement_m
        reciprocity[index] = result.reciprocity_residual
    end
    refined_transfer, refined_target = normalized_transfer(raw, refined_geometry)
    append!(rows, resolution_rows(
        "refined_pzt0p45_al0p75", refined_transfer, refined_target, reciprocity,
    ))
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "convergence_anchors.csv"), rows)
    jldsave(
        joinpath(OUTPUT_ROOT, "o6b_stack_convergence.jld2");
        format_version=1,
        anchors_hz=ANCHORS_HZ,
        mesh_metadata,
        raw_transfer=raw,
        as_built_transfer=refined_transfer,
        reciprocity,
        rows,
    )
    base_worst_fidelity = minimum(
        row.normalized_transfer_fidelity for row in rows if startswith(row.resolution, "base"),
    )
    refined_worst_fidelity = minimum(
        row.normalized_transfer_fidelity for row in rows if startswith(row.resolution, "refined"),
    )
    base_worst_amplitude = maximum(
        row.differential_amplitude_rms for row in rows if startswith(row.resolution, "base"),
    )
    refined_worst_amplitude = maximum(
        row.differential_amplitude_rms for row in rows if startswith(row.resolution, "refined"),
    )
    open(joinpath(OUTPUT_ROOT, "convergence_verdict.md"), "w") do io
        println(io, "# O6b finite-twin convergence")
        println(io)
        @printf(io, "Base/refined worst anchor fidelity: `%.6f / %.6f` (gate `0.98`).\n\n",
            base_worst_fidelity, refined_worst_fidelity)
        @printf(io, "Base/refined worst anchor amplitude RMS: `%.3f%% / %.3f%%` (gate `3%%`).\n\n",
            100base_worst_amplitude, 100refined_worst_amplitude)
        println(io, "The refinement preserves the sign and scale of the failure. The homogeneous PZT slab is rejected without a geometry sweep; the dense base-grid delay result remains diagnostic, while fidelity and amplitude provide the mesh-checked no-go.")
    end
    println("[+] O6b convergence complete")
    rows
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_convergence()

end
