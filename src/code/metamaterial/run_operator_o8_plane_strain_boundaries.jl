module RunOperatorO8PlaneStrainBoundaries

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Statistics

include(joinpath(@__DIR__, "operator_o8_boundary_mesher.jl"))
include(joinpath(@__DIR__, "operator_o8_plane_strain_solver.jl"))

using .OperatorO8BoundaryMesher
using .OperatorO8PlaneStrainSolver

export run_o8_plane_strain_boundaries

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o8_plane_strain_v0",
)
const O7_BUNDLE = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o7_global_state_v0",
    "o7_global_states.jld2",
)
const CARRIER_HZ = 242.0e3
const DENSITY_KG_M3 = 2700.0
const PRESSURE_SPEED_M_S = 6122.102437409232
const SHEAR_SPEED_M_S = 3083.810277185563
const PRESSURE_WAVELENGTH_M = PRESSURE_SPEED_M_S / CARRIER_HZ
const SHEAR_WAVELENGTH_M = SHEAR_SPEED_M_S / CARRIER_HZ
const APERTURE_HEIGHT_M = 8PRESSURE_WAVELENGTH_M
const BLOCK_LENGTH_M = 2PRESSURE_WAVELENGTH_M
const FREQUENCY_RATIOS = collect(range(0.8, 1.2; length=7))
const FREQUENCIES_HZ = CARRIER_HZ .* FREQUENCY_RATIOS
const CARRIER_INDEX = 4
const REFINED_INDICES = (1, 4, 7)
const BASE_MESH_SIZE_M = SHEAR_WAVELENGTH_M / 6
const REFINED_MESH_SIZE_M = SHEAR_WAVELENGTH_M / 8
const BITS = ((0, 0), (1, 0), (0, 1), (1, 1))

function write_csv(path, rows)
    isempty(rows) && return
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function geometry(parameters, mesh_size)
    O8BoundaryGeometry(
        aperture_height_m=APERTURE_HEIGHT_M,
        block_length_m=BLOCK_LENGTH_M,
        position_a=parameters.position_a,
        width_a=parameters.width_a,
        position_b=parameters.position_b,
        width_b=parameters.width_b,
        mesh_size_m=mesh_size,
    )
end

mesh_path(mechanism, level) = joinpath(
    OUTPUT_ROOT, "meshes", "$(mechanism)_$(level).msh",
)

function load_physical_definitions()
    isfile(O7_BUNDLE) || error("missing O7 bundle: $O7_BUNDLE")
    bundle = JLD2.load(O7_BUNDLE)
    pzt_parameters = bundle["pzt_parameters"]
    mechanical_parameters = bundle["mechanical_parameters"]
    pzt = O8PZTBoundary()
    load_a = O8Load(
        inductance_h=94.40989561054282e-6,
        resistance_ohm=24.07157769563397,
    )
    load_b = O8Load(
        inductance_h=98.6120227995537e-6,
        resistance_ohm=37.099652264801485,
    )
    natural_surface_stiffness = DENSITY_KG_M3 * PRESSURE_SPEED_M_S^2 /
        PRESSURE_WAVELENGTH_M
    mechanical = O8MechanicalBoundary(
        stiffness_a_n_m3=mechanical_parameters.strength_a * natural_surface_stiffness,
        stiffness_b_n_m3=mechanical_parameters.strength_b * natural_surface_stiffness,
        loss_a=mechanical_parameters.loss_a,
        loss_b=mechanical_parameters.loss_b,
    )
    (;
        pzt_parameters, mechanical_parameters, pzt, load_a, load_b,
        mechanical, natural_surface_stiffness,
    )
end

function solve_branch(
    mechanism, mesh, frequency_indices, definition;
    zero_coupling=false,
)
    results = NamedTuple[]
    for (run_index, frequency_index) in enumerate(frequency_indices)
        frequency_hz = FREQUENCIES_HZ[frequency_index]
        println(
            "[>] O8 $(mechanism) $(run_index)/$(length(frequency_indices)): " *
            "$(frequency_hz / 1e3) kHz",
        )
        config = O8PlaneStrainConfig(
            frequency_hz=frequency_hz,
            aperture_height_m=APERTURE_HEIGHT_M,
            block_length_m=BLOCK_LENGTH_M,
        )
        result = if mechanism == :pzt
            solve_o8_frequency(
                mesh; mechanism, config, pzt=definition.pzt,
                load_a=definition.load_a, load_b=definition.load_b,
                zero_coupling,
            )
        else
            solve_o8_frequency(
                mesh; mechanism, config, mechanical=definition.mechanical,
                zero_coupling,
            )
        end
        push!(results, merge((; frequency_index), result))
    end
    results
end

function state_arrays(results)
    outputs = zeros(ComplexF64, 17, 3, 4, length(results))
    raw_outputs = similar(outputs)
    references = zeros(ComplexF64, 17, 3, length(results))
    for (index, result) in enumerate(results)
        outputs[:, :, :, index] .= result.state_outputs
        raw_outputs[:, :, :, index] .= result.raw_state_outputs
        references[:, :, index] .= result.reference_output
    end
    outputs, raw_outputs, references
end

function state_rows(mechanism, level, results)
    rows = NamedTuple[]
    for result in results, state in 1:4
        push!(rows, (
            mechanism, mesh_level=level, frequency_hz=result.frequency_hz,
            frequency_ratio=result.frequency_hz / CARRIER_HZ,
            state="q$(BITS[state][1])$(BITS[state][2])",
            raw_reflection_error=result.raw_mirror_errors[state],
            projected_reflection_error=result.projected_mirror_errors[state],
            reciprocity_error=result.reciprocity_errors[state],
            amplitude_reference_ratio=result.amplitude_ratios[state],
            minimum_accepted_work=minimum(result.accepted_work[:, state]),
        ))
    end
    rows
end

function distance_rows(mechanism, level, results)
    pairs = [(first, second) for first in 1:4 for second in (first + 1):4]
    rows = NamedTuple[]
    for result in results, (pair_index, pair) in enumerate(pairs)
        push!(rows, (
            mechanism, mesh_level=level, frequency_hz=result.frequency_hz,
            frequency_ratio=result.frequency_hz / CARRIER_HZ,
            first_state="q$(BITS[pair[1]][1])$(BITS[pair[1]][2])",
            second_state="q$(BITS[pair[2]][1])$(BITS[pair[2]][2])",
            projector_distance=result.state_distances[pair_index],
        ))
    end
    rows
end

function zero_collapse(result)
    reference = result.state_outputs[:, :, 1]
    maximum(
        norm(result.state_outputs[:, :, state] - reference) /
            max(norm(reference), eps(Float64))
        for state in 2:4
    )
end

function convergence_rows(mechanism, base_results, refined_results)
    base = Dict(result.frequency_index => result for result in base_results)
    rows = NamedTuple[]
    for refined in refined_results, state in 1:4
        coarse = base[refined.frequency_index]
        projector_error = state_projector_distance(
            coarse.state_outputs[:, :, state], refined.state_outputs[:, :, state],
        )
        amplitude_relative_error = abs(
            refined.amplitude_ratios[state] - coarse.amplitude_ratios[state],
        ) / max(coarse.amplitude_ratios[state], eps(Float64))
        push!(rows, (
            mechanism, frequency_hz=refined.frequency_hz,
            state="q$(BITS[state][1])$(BITS[state][2])",
            projector_error, amplitude_relative_error,
            base_raw_reflection=coarse.raw_mirror_errors[state],
            refined_raw_reflection=refined.raw_mirror_errors[state],
        ))
    end
    rows
end

function physical_gates(
    mechanism, base_results, refined_results, convergence, collapse,
    definition,
)
    carrier = only(filter(result -> result.frequency_index == CARRIER_INDEX, base_results))
    all_state_rows = state_rows(mechanism, :base, base_results)
    component_positive = if mechanism == :pzt
        definition.load_a.inductance_h > 0 && definition.load_a.resistance_ohm > 0 &&
        definition.load_b.inductance_h > 0 && definition.load_b.resistance_ohm > 0 &&
        definition.pzt.capacitance_per_patch_f > 0 &&
        definition.pzt.c33e_pa > 0 && definition.pzt.mechanical_quality_factor > 0
    else
        definition.mechanical.stiffness_a_n_m3 > 0 &&
        definition.mechanical.stiffness_b_n_m3 > 0 &&
        definition.mechanical.loss_a >= 0 && definition.mechanical.loss_b >= 0
    end
    gates = (
        reciprocity=maximum(getproperty.(all_state_rows, :reciprocity_error)) <= 1.0e-6,
        raw_mesh_reflection=maximum(getproperty.(all_state_rows, :raw_reflection_error)) <= 0.02,
        projected_reflection=maximum(getproperty.(all_state_rows, :projected_reflection_error)) <= 1.0e-8,
        component_positive,
        carrier_state_distance=minimum(carrier.state_distances) >= 0.05,
        zero_coupling_collapse=collapse <= 1.0e-10,
        amplitude_floor=minimum(getproperty.(all_state_rows, :amplitude_reference_ratio)) >= 0.70,
        amplitude_ceiling=maximum(getproperty.(all_state_rows, :amplitude_reference_ratio)) <= 1.50,
        convergence_projector=maximum(getproperty.(convergence, :projector_error)) <= 0.10,
        convergence_amplitude=maximum(getproperty.(convergence, :amplitude_relative_error)) <= 0.10,
        refinement_reflection=maximum(getproperty.(convergence, :refined_raw_reflection)) <= 0.02,
    )
    gates, all(values(gates))
end

function normalized_frequency_states(outputs)
    raw_maximum_sigma = maximum(
        opnorm(outputs[:, :, state, frequency])
        for state in 1:4 for frequency in axes(outputs, 4)
    )
    common_scale = 0.98 / raw_maximum_sigma
    states = [
        [common_scale .* outputs[:, :, state, frequency] for frequency in axes(outputs, 4)]
        for state in 1:4
    ]
    states, common_scale, raw_maximum_sigma
end

function plot_branch(mechanism, state_rows_value, distance_rows_value)
    amplitude = plot(
        xlabel="frequency / f0", ylabel="output/reference amplitude",
        title="O8 $(mechanism): physical amplitude", legend=:best,
        ylim=(0, 1.55), size=(900, 550),
    )
    for state in ("q00", "q10", "q01", "q11")
        selected = filter(row -> row.state == state, state_rows_value)
        plot!(amplitude, getproperty.(selected, :frequency_ratio),
            getproperty.(selected, :amplitude_reference_ratio); marker=:circle,
            label=state)
    end
    hline!(amplitude, [0.70, 1.50]; color=:black, linestyle=:dash, label="gate")
    savefig(amplitude, joinpath(OUTPUT_ROOT, "$(mechanism)_amplitude_band.png"))

    distance = plot(
        xlabel="frequency / f0", ylabel="projector distance",
        title="O8 $(mechanism): four-state separation", legend=false,
        size=(900, 550),
    )
    for first in ("q00", "q10", "q01")
        for second in ("q10", "q01", "q11")
            selected = filter(row -> row.first_state == first && row.second_state == second,
                distance_rows_value)
            isempty(selected) && continue
            plot!(distance, getproperty.(selected, :frequency_ratio),
                getproperty.(selected, :projector_distance); marker=:circle)
        end
    end
    hline!(distance, [0.05]; color=:black, linestyle=:dash)
    savefig(distance, joinpath(OUTPUT_ROOT, "$(mechanism)_state_distance.png"))
end

function plot_pzt_schematic(parameters)
    config = geometry(parameters, BASE_MESH_SIZE_M)
    intervals = control_intervals(config)
    scale = 1.0e3
    length_mm = BLOCK_LENGTH_M * scale
    half_mm = APERTURE_HEIGHT_M * scale / 2
    figure = plot(
        xlim=(-30, length_mm + 8), ylim=(-half_mm - 8, half_mm + 8),
        aspect_ratio=:equal, xlabel="propagation direction x / mm",
        ylabel="aperture coordinate y / mm",
        title="O8 frozen stationary PZT head (plane-strain section)",
        legend=:topright, size=(1050, 760),
    )
    body = Shape(
        [0.0, length_mm, length_mm, 0.0],
        [-half_mm, -half_mm, half_mm, half_mm],
    )
    plot!(figure, body; color=:lightgray, linecolor=:black, label="Al block")
    for (zone, interval) in pairs(intervals)
        lower, upper = scale .* collect(interval)
        color = startswith(String(zone), "a_") ? :dodgerblue : :darkorange
        patch = Shape([-3.0, 0.0, 0.0, -3.0], [lower, lower, upper, upper])
        plot!(figure, patch; color, linecolor=:black,
            label=zone in (:a_minus, :b_minus) ?
                (startswith(String(zone), "a_") ? "pair A: even RL" : "pair B: odd RL") : "")
    end
    plot!(figure, [0.0, length_mm], [0.0, 0.0];
        color=:black, linestyle=:dash, linewidth=1, label="mirror plane")
    annotate!(figure, length_mm / 2, 0, text("3 fixed ports\nDCT0 / DCT2 / DCT1", 10))
    annotate!(figure, length_mm + 2, 0, text("17-bin\ncontact face", 9, :left))
    annotate!(figure, -12, 52, text("bit A", 9, :blue))
    annotate!(figure, -12, 92, text("bit B", 9, :darkorange))
    savefig(figure, joinpath(OUTPUT_ROOT, "pzt_head_schematic.png"))
end

function write_readme(branch_rows, gate_rows)
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O8 frozen plane-strain boundary realization")
        println(io)
        println(io, "Evidence: plane-strain elastic FEM with finite control zones. PZT electrodes/load are represented by an exact one-port constitutive Schur complement; the patch electric field is not domain-meshed.")
        println(io)
        println(io, "| Branch | min amplitude | max amplitude | carrier min state distance | max raw mirror error | max reciprocity error | Physical gate |")
        println(io, "|---|---:|---:|---:|---:|---:|---|")
        for branch in branch_rows
            println(io, "| $(branch.mechanism) | $(round(branch.minimum_amplitude; digits=4)) | " *
                "$(round(branch.maximum_amplitude; digits=4)) | $(round(branch.carrier_minimum_state_distance; digits=4)) | " *
                "$(round(branch.maximum_raw_reflection_error; digits=5)) | $(branch.maximum_reciprocity_error) | $(branch.physical_pass) |")
        end
        println(io)
        for branch in branch_rows
            failed = [row.gate for row in gate_rows if row.mechanism == branch.mechanism && !row.passed]
            println(io, "- `$(branch.mechanism)`: failed gates = $(isempty(failed) ? "none" : join(failed, ", ")).")
        end
        println(io)
        println(io, "Only a physically passing branch may enter the fresh O8 detector replay. The O7 final-test split was not used.")
    end
end

function run_o8_plane_strain_boundaries()
    mkpath(OUTPUT_ROOT)
    definition = load_physical_definitions()
    plot_pzt_schematic(definition.pzt_parameters)
    for mechanism in (:pzt, :mechanical)
        parameters = mechanism == :pzt ? definition.pzt_parameters :
            definition.mechanical_parameters
        for (level, mesh_size) in ((:base, BASE_MESH_SIZE_M), (:refined, REFINED_MESH_SIZE_M))
            path = mesh_path(mechanism, level)
            metadata = build_o8_boundary_mesh(
                path; config=geometry(parameters, mesh_size),
            )
            println("[+] O8 $(mechanism) $(level) mesh: $(metadata.node_count) nodes")
        end
    end

    all_state_rows = NamedTuple[]
    all_distance_rows = NamedTuple[]
    all_convergence_rows = NamedTuple[]
    gate_rows = NamedTuple[]
    branch_rows = NamedTuple[]
    bundles = Dict{Symbol,Any}()

    for mechanism in (:pzt, :mechanical)
        base_results = solve_branch(
            mechanism, mesh_path(mechanism, :base), collect(1:7), definition,
        )
        refined_results = solve_branch(
            mechanism, mesh_path(mechanism, :refined), collect(REFINED_INDICES), definition,
        )
        zero_result = only(solve_branch(
            mechanism, mesh_path(mechanism, :base), [CARRIER_INDEX], definition;
            zero_coupling=true,
        ))
        collapse = zero_collapse(zero_result)
        convergence = convergence_rows(mechanism, base_results, refined_results)
        states = state_rows(mechanism, :base, base_results)
        distances = distance_rows(mechanism, :base, base_results)
        refined_states = state_rows(mechanism, :refined, refined_results)
        refined_distances = distance_rows(mechanism, :refined, refined_results)
        append!(all_state_rows, states)
        append!(all_state_rows, refined_states)
        append!(all_distance_rows, distances)
        append!(all_distance_rows, refined_distances)
        append!(all_convergence_rows, convergence)
        gates, passed = physical_gates(
            mechanism, base_results, refined_results, convergence, collapse, definition,
        )
        append!(gate_rows, [
            (mechanism, gate=String(name), passed=value) for (name, value) in pairs(gates)
        ])
        carrier = base_results[CARRIER_INDEX]
        outputs, raw_outputs, references = state_arrays(base_results)
        frequency_states, common_scale, raw_maximum_sigma =
            normalized_frequency_states(outputs)
        push!(branch_rows, (
            mechanism,
            minimum_amplitude=minimum(getproperty.(states, :amplitude_reference_ratio)),
            maximum_amplitude=maximum(getproperty.(states, :amplitude_reference_ratio)),
            carrier_minimum_state_distance=minimum(carrier.state_distances),
            maximum_raw_reflection_error=maximum(getproperty.(states, :raw_reflection_error)),
            maximum_reciprocity_error=maximum(getproperty.(states, :reciprocity_error)),
            maximum_convergence_projector=maximum(getproperty.(convergence, :projector_error)),
            maximum_convergence_amplitude=maximum(getproperty.(convergence, :amplitude_relative_error)),
            zero_collapse=collapse, common_scale, raw_maximum_sigma,
            physical_pass=passed,
        ))
        bundles[mechanism] = (
            base_results, refined_results, zero_result, outputs, raw_outputs,
            references, frequency_states, gates, passed, common_scale,
            raw_maximum_sigma,
        )
        plot_branch(mechanism, states, distances)
        jldsave(
            joinpath(OUTPUT_ROOT, "$(mechanism)_physical_bundle.jld2");
            mechanism, base_results, refined_results, zero_result, outputs,
            raw_outputs, references, frequency_states, gates, passed,
            common_scale, raw_maximum_sigma,
        )
        println("[+] O8 $(mechanism) physical pass=$(passed)")
    end

    write_csv(joinpath(OUTPUT_ROOT, "state_band_metrics.csv"), all_state_rows)
    write_csv(joinpath(OUTPUT_ROOT, "state_distances.csv"), all_distance_rows)
    write_csv(joinpath(OUTPUT_ROOT, "mesh_convergence.csv"), all_convergence_rows)
    write_csv(joinpath(OUTPUT_ROOT, "physical_gates.csv"), gate_rows)
    write_csv(joinpath(OUTPUT_ROOT, "branch_summary.csv"), branch_rows)
    write_readme(branch_rows, gate_rows)
    jldsave(
        joinpath(OUTPUT_ROOT, "o8_plane_strain_boundaries.jld2");
        frequencies_hz=FREQUENCIES_HZ, frequency_ratios=FREQUENCY_RATIOS,
        base_mesh_size_m=BASE_MESH_SIZE_M,
        refined_mesh_size_m=REFINED_MESH_SIZE_M,
        definition, branch_rows, gate_rows,
    )
    branch_rows
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_o8_plane_strain_boundaries()
end

end
