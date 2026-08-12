module MultimodePortHarmonicSolver

using Gridap
using GridapGmsh
using LinearAlgebra
using SparseArrays

if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
end
if !isdefined(parentmodule(@__MODULE__), :ElasticPortModes)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "port_mode_solver.jl"))
end
if !isdefined(parentmodule(@__MODULE__), :ElasticModalProjection)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "modal_projection.jl"))
end
using ..SinusoidalMaterialLens
using ..ElasticPortModes
using ..ElasticModalProjection

export MultimodePortConfig, solve_multimode_port_harmonic

Base.@kwdef struct MultimodePortConfig
    frequency_hz::Float64 = 242.0e3
    element_order::Int = 2
    quadrature_degree::Int = 4
    port_element_count::Int = 120
    evanescent_modes_per_port::Int = 8
end

function elastic_stress(material::ElasticMaterial, strain)
    rho = material.density_kg_m3
    mu = rho * material.shear_wave_speed_m_s^2
    lambda = rho * material.pressure_wave_speed_m_s^2 - 2mu
    identity = one(TensorValue{2, 2, ComplexF64})
    lambda * tr(strain) * identity + 2mu * strain
end

function interpolate_mode(mode, physical_y, center_y)
    y = physical_y - center_y
    displacement = VectorValue(
        ElasticModalProjection.linear_interpolate(mode.y_m, mode.displacement_x, y),
        ElasticModalProjection.linear_interpolate(mode.y_m, mode.displacement_y, y),
    )
    traction_x = VectorValue(
        ElasticModalProjection.linear_interpolate(mode.y_m, mode.traction_xx_pa, y),
        ElasticModalProjection.linear_interpolate(mode.y_m, mode.traction_xy_pa, y),
    )
    displacement, traction_x
end

function mode_cell_fields(boundary, mode, center_y, outward_sign)
    coordinate = get_physical_coordinate(boundary)
    displacement = (point -> first(interpolate_mode(mode, point[2], center_y))) ∘ coordinate
    outward_traction = (
        point -> outward_sign * last(interpolate_mode(mode, point[2], center_y))
    ) ∘ coordinate
    displacement, outward_traction
end

function port_modes(material, frequency_hz, width_m, element_count)
    solve_port_modes(PortModeConfig(
        height_m=width_m,
        density=material.density_kg_m3,
        pressure_wave_speed=material.pressure_wave_speed_m_s,
        shear_wave_speed=material.shear_wave_speed_m_s,
        frequency_hz=frequency_hz,
        element_count=element_count,
    ))
end

function outgoing_basis(modes, direction, evanescent_count)
    propagating = filter(
        mode -> mode.kind == :propagating && mode.direction == direction,
        modes,
    )
    evanescent = filter(
        mode -> mode.kind == :evanescent && mode.direction == direction,
        modes,
    )
    sort!(evanescent; by=mode -> abs(real(mode.gamma_per_m)))
    vcat(propagating, first(evanescent, min(evanescent_count, length(evanescent))))
end

function boundary_width_and_center(boundary)
    data = only(Gridap.Visualization.visualization_data(
        boundary, "multimode_port_coordinates"; order=1,
    ))
    coordinates = collect(Gridap.Geometry.get_node_coordinates(data.grid))
    y = getindex.(coordinates, 2)
    maximum(y) - minimum(y), (minimum(y) + maximum(y)) / 2
end

function port_operator_data(
    test_space,
    model,
    tag,
    material,
    config,
    direction,
    outward_sign,
)
    boundary = BoundaryTriangulation(model; tags=[String(tag)])
    measure = Measure(boundary, config.quadrature_degree)
    width_m, center_y = boundary_width_and_center(boundary)
    modes = port_modes(
        material, config.frequency_hz, width_m, config.port_element_count,
    )
    basis_modes = outgoing_basis(
        modes, direction, config.evanescent_modes_per_port,
    )
    displacement_fields = Any[]
    traction_fields = Any[]
    for mode in basis_modes
        displacement, traction = mode_cell_fields(
            boundary, mode, center_y, outward_sign,
        )
        push!(displacement_fields, displacement)
        push!(traction_fields, traction)
    end
    count = length(basis_modes)
    gram = Matrix{ComplexF64}(undef, count, count)
    conjugated(field) = conj(field)
    for row in 1:count, column in 1:count
        conjugate_displacement = conjugated ∘ displacement_fields[row]
        gram[row, column] = sum(
            ∫(conjugate_displacement ⋅ displacement_fields[column])measure
        )
    end
    gram_inverse = inv(gram)
    displacement_vectors = [
        assemble_vector(
            v -> begin
                conjugate_displacement = conjugated ∘ displacement_fields[index]
                ∫(conjugate_displacement ⋅ v)measure
            end,
            test_space,
        ) for index in 1:count
    ]
    traction_vectors = [
        assemble_vector(v -> ∫(traction_fields[index] ⋅ v)measure, test_space)
        for index in 1:count
    ]
    (
        tag=String(tag), boundary, measure, width_m, center_y, modes, basis_modes,
        displacement_fields, traction_fields, gram, gram_inverse,
        displacement_vectors, traction_vectors,
        gram_condition_number=cond(gram),
    )
end

function add_port_operator!(matrix, port)
    count = length(port.basis_modes)
    for output_index in 1:count, input_index in 1:count
        coefficient = port.gram_inverse[output_index, input_index]
        abs(coefficient) == 0 && continue
        output_vector = sparse(port.traction_vectors[output_index])
        input_vector = sparse(port.displacement_vectors[input_index])
        matrix .-= coefficient .* (output_vector * transpose(input_vector))
    end
    matrix
end

function incident_load_vector(test_space, source_port, incident_mode)
    incident_displacement, incident_traction = mode_cell_fields(
        source_port.boundary,
        incident_mode,
        source_port.center_y,
        -1.0,
    )
    conjugated(field) = conj(field)
    overlap = ComplexF64[
        sum(∫(
            (conjugated ∘ source_port.displacement_fields[index]) ⋅
            incident_displacement
        )source_port.measure)
        for index in eachindex(source_port.basis_modes)
    ]
    outgoing_coefficients = source_port.gram_inverse * overlap
    load = assemble_vector(
        v -> ∫(incident_traction ⋅ v)source_port.measure,
        test_space,
    )
    for index in eachindex(outgoing_coefficients)
        load .-= outgoing_coefficients[index] .* source_port.traction_vectors[index]
    end
    load
end

function propagating_amplitudes(state, tag, modes)
    right_modes = filter(
        mode -> mode.kind == :propagating && mode.direction == :right,
        modes,
    )
    [
        begin
            left_mode = counterpropagating_mode(right_mode, modes)
            projection = project_boundary_mode_pair(
                state, tag, right_mode, left_mode;
                quadrature_degree=state.quadrature_degree,
            )
            (
                parity=right_mode.parity,
                wavenumber_per_m=right_mode.wavenumber_per_m,
                right_amplitude=projection.right_amplitude,
                left_amplitude=projection.left_amplitude,
            )
        end
        for right_mode in right_modes
    ]
end

function solve_multimode_port_harmonic(
    mesh_path::AbstractString,
    output_tags::AbstractVector{<:AbstractString},
    material::ElasticMaterial=aluminium_6061();
    config::MultimodePortConfig=MultimodePortConfig(),
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    isempty(output_tags) && throw(ArgumentError("at least one output port is required"))
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    rho = material.density_kg_m3
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model, reference; conformity=:H1, vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)
    domain = Triangulation(model; tags=["Domain"])
    measure = Measure(domain, config.quadrature_degree)
    stiffness(u, v) = elastic_stress(material, ε(u)) ⊙ ε(v)
    bilinear(u, v) = ∫(
        stiffness(u, v) - rho * omega^2 * (u ⋅ v) +
        im * omega * (
            material.rayleigh_alpha_s_inv * rho * (u ⋅ v) +
            material.rayleigh_beta_s * stiffness(u, v)
        )
    )measure
    matrix = assemble_matrix(bilinear, trial_space, test_space)
    source_port = port_operator_data(
        test_space, model, "Source", material, config, :left, -1.0,
    )
    output_ports = [
        port_operator_data(test_space, model, tag, material, config, :right, 1.0)
        for tag in output_tags
    ]
    add_port_operator!(matrix, source_port)
    for port in output_ports
        add_port_operator!(matrix, port)
    end
    incident_mode = fundamental_quasi_longitudinal(source_port.modes)
    load = incident_load_vector(test_space, source_port, incident_mode)
    coefficients = matrix \ load
    displacement = FEFunction(trial_space, coefficients)
    stress = elastic_stress(material, ε(displacement))
    state = (
        model, displacement, stress, quadrature_degree=config.quadrature_degree,
    )
    source_amplitudes = propagating_amplitudes(state, "Source", source_port.modes)
    output_amplitudes = [
        propagating_amplitudes(state, tag, port.modes)
        for (tag, port) in zip(output_tags, output_ports)
    ]
    incident_power = sum(abs2(item.right_amplitude) for item in source_amplitudes)
    reflected_power = sum(abs2(item.left_amplitude) for item in source_amplitudes)
    transmitted_power = [
        sum(abs2(item.right_amplitude) for item in amplitudes)
        for amplitudes in output_amplitudes
    ]
    incoming_output_power = [
        sum(abs2(item.left_amplitude) for item in amplitudes)
        for amplitudes in output_amplitudes
    ]
    (
        frequency_hz=config.frequency_hz,
        incident_mode,
        incident_power_w_per_m=incident_power,
        reflected_power_w_per_m=reflected_power,
        transmitted_power_w_per_m=transmitted_power,
        incoming_output_power_w_per_m=incoming_output_power,
        source_amplitudes,
        output_amplitudes,
        port_basis_sizes=[length(source_port.basis_modes); length.(getproperty.(output_ports, :basis_modes))],
        port_gram_condition_numbers=[source_port.gram_condition_number; getproperty.(output_ports, :gram_condition_number)],
        displacement,
        model,
    )
end

end
