module ModalHarmonicElasticity

using Gridap

export ModalBoundaryMetrics, solve_symmetry_reduced_modal_state

if !isdefined(@__MODULE__, :HarmonicElasticity)
    include(joinpath(@__DIR__, "harmonic_solver.jl"))
end
using .HarmonicElasticity: HarmonicConfig

if !isdefined(@__MODULE__, :ElasticPortModes)
    include(joinpath(@__DIR__, "port_mode_solver.jl"))
end
using .ElasticPortModes

struct ModalBoundaryMetrics
    frequency_hz::Float64
    source_work_w_per_m::Float64
    left_boundary_power_w_per_m::Float64
    right_boundary_power_w_per_m::Float64
    internal_dissipated_power_w_per_m::Float64
    balance_residual_w_per_m::Float64
end

function linear_interpolate(x, values, query)
    query <= first(x) && return values[1]
    query >= last(x) && return values[end]
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function boundary_center_y(boundary)
    data = only(Gridap.Visualization.visualization_data(boundary, "port_coordinates"; order=1))
    coordinates = collect(Gridap.Geometry.get_node_coordinates(data.grid))
    y = getindex.(coordinates, 2)
    (minimum(y) + maximum(y)) / 2
end

function mode_values(mode, physical_y, center_y)
    y = physical_y - center_y
    displacement = VectorValue(
        linear_interpolate(mode.y_m, mode.displacement_x, y),
        linear_interpolate(mode.y_m, mode.displacement_y, y),
    )
    traction_x = VectorValue(
        linear_interpolate(mode.y_m, mode.traction_xx_pa, y),
        linear_interpolate(mode.y_m, mode.traction_xy_pa, y),
    )
    displacement, traction_x
end

"""Return the rank-one tensor Z satisfying Z*u_mode = t_mode pointwise."""
function rank_one_traction_map(displacement, traction)
    denominator = abs2(displacement[1]) + abs2(displacement[2])
    denominator > floatmin(Float64) || error("modal displacement vanishes on the port")
    z11 = traction[1] * conj(displacement[1]) / denominator
    z12 = traction[1] * conj(displacement[2]) / denominator
    z21 = traction[2] * conj(displacement[1]) / denominator
    z22 = traction[2] * conj(displacement[2]) / denominator
    # TensorValue arguments are stored column-major.
    TensorValue(z11, z21, z12, z22)
end

function modal_boundary_fields(boundary, outgoing_mode, outward_sign)
    center_y = boundary_center_y(boundary)
    coordinate = get_physical_coordinate(boundary)
    operator_function = point -> begin
        displacement, traction_x = mode_values(outgoing_mode, point[2], center_y)
        rank_one_traction_map(displacement, outward_sign * traction_x)
    end
    operator_function ∘ coordinate
end

function incident_source_field(boundary, right_mode, left_operator)
    center_y = boundary_center_y(boundary)
    coordinate = get_physical_coordinate(boundary)
    incident_displacement_function = point -> first(mode_values(right_mode, point[2], center_y))
    incident_outward_traction_function = point -> begin
        _, traction_x = mode_values(right_mode, point[2], center_y)
        -traction_x
    end
    incident_displacement = incident_displacement_function ∘ coordinate
    incident_outward_traction = incident_outward_traction_function ∘ coordinate
    incident_outward_traction - left_operator ⋅ incident_displacement
end

function matched_port_modes(config, frequency_hz, port_height_m, port_element_count)
    port_config = PortModeConfig(
        height_m=Float64(port_height_m),
        density=config.density,
        pressure_wave_speed=config.pressure_wave_speed,
        shear_wave_speed=config.shear_wave_speed,
        frequency_hz=Float64(frequency_hz),
        element_count=Int(port_element_count),
    )
    modes = solve_port_modes(port_config)
    symmetric_right = filter(
        mode -> mode.kind == :propagating &&
                mode.direction == :right &&
                mode.parity == :symmetric,
        modes,
    )
    length(symmetric_right) == 1 || error(
        "symmetry-reduced boundary requires exactly one right-going symmetric mode; " *
        "found $(length(symmetric_right))",
    )
    right_mode = only(symmetric_right)
    symmetric_left = filter(
        mode -> mode.kind == :propagating &&
                mode.direction == :left &&
                mode.parity == :symmetric,
        modes,
    )
    isempty(symmetric_left) && error("left-going symmetric partner was not found")
    left_mode = first(sort(symmetric_left; by=mode -> abs(
        mode.wavenumber_per_m + right_mode.wavenumber_per_m,
    )))
    modes, right_mode, left_mode
end

"""
Solve a y-symmetric cell with modal injection and a rank-one matched boundary.

This boundary is exact for the sole propagating symmetric mode. It is not a
full multimode DtN: odd modes must remain symmetry-forbidden and straight lead
sections must be long enough for symmetric evanescent modes to decay.
"""
function solve_symmetry_reduced_modal_state(
    model_path::AbstractString,
    frequency_hz::Real;
    config::HarmonicConfig=HarmonicConfig(rayleigh_alpha=0.0, rayleigh_beta=0.0),
    port_height_m::Real=4.2e-3,
    port_element_count::Integer=80,
)
    isfile(model_path) || error("Gridap model not found: $model_path")
    model = DiscreteModelFromFile(model_path)
    omega = 2pi * Float64(frequency_hz)
    rho = config.density
    mu = rho * config.shear_wave_speed^2
    lambda = rho * config.pressure_wave_speed^2 - 2mu
    identity_tensor = one(TensorValue{2, 2, ComplexF64})
    sigma(strain) = lambda * tr(strain) * identity_tensor + 2mu * strain

    modes, right_mode, left_mode = matched_port_modes(
        config,
        frequency_hz,
        port_height_m,
        port_element_count,
    )

    reference_element = ReferenceFE(
        lagrangian,
        VectorValue{2, Float64},
        config.element_order,
    )
    test_space = TestFESpace(
        model,
        reference_element;
        conformity=:H1,
        vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)
    domain = Triangulation(model)
    domain_measure = Measure(domain, config.quadrature_degree)
    left_port = BoundaryTriangulation(model; tags=["Source"])
    right_port = BoundaryTriangulation(model; tags=["Microphone"])
    left_measure = Measure(left_port, config.quadrature_degree)
    right_measure = Measure(right_port, config.quadrature_degree)

    # At the right port outward traction is +sigma*e_x for the right-going
    # mode. At the left it is -sigma*e_x for the left-going mode.
    left_operator = modal_boundary_fields(left_port, left_mode, -1.0)
    right_operator = modal_boundary_fields(right_port, right_mode, 1.0)
    source_traction = incident_source_field(left_port, right_mode, left_operator)

    stiffness(u, v) = sigma(ε(u)) ⊙ ε(v)
    bilinear(u, v) =
        ∫(
            stiffness(u, v) - rho * omega^2 * (u ⋅ v) +
            im * omega * (
                config.rayleigh_alpha * rho * (u ⋅ v) +
                config.rayleigh_beta * stiffness(u, v)
            ),
        )domain_measure -
        ∫((left_operator ⋅ u) ⋅ v)left_measure -
        ∫((right_operator ⋅ u) ⋅ v)right_measure
    linear(v) = ∫(source_traction ⋅ v)left_measure

    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))
    stress = sigma ∘ ε(displacement)
    velocity = im * omega * displacement
    conjugated(value) = conj(value)
    conjugate_velocity = conjugated ∘ velocity
    conjugate_displacement = conjugated ∘ displacement
    conjugate_strain = conjugated ∘ ε(displacement)
    left_outgoing_traction = left_operator ⋅ displacement
    right_outgoing_traction = right_operator ⋅ displacement
    source_work = 0.5 * real(sum(∫(source_traction ⋅ conjugate_velocity)left_measure))
    left_power = -0.5 * real(sum(∫(left_outgoing_traction ⋅ conjugate_velocity)left_measure))
    right_power = -0.5 * real(sum(∫(right_outgoing_traction ⋅ conjugate_velocity)right_measure))
    internal_power = 0.5 * omega^2 * real(sum(
        ∫(
            config.rayleigh_alpha * rho * (displacement ⋅ conjugate_displacement) +
            config.rayleigh_beta * (sigma(ε(displacement)) ⊙ conjugate_strain)
        )domain_measure,
    ))
    balance_residual = source_work - left_power - right_power - internal_power
    metrics = ModalBoundaryMetrics(
        Float64(frequency_hz),
        Float64(source_work),
        Float64(left_power),
        Float64(right_power),
        Float64(internal_power),
        Float64(balance_residual),
    )
    (; metrics, model, domain, displacement, stress, omega, modes, right_mode, left_mode)
end

end
