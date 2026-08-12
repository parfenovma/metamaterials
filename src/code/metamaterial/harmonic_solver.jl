module HarmonicElasticity

using Gridap
using JLD2

export HarmonicConfig,
       HarmonicPoint,
       HarmonicField,
       solve_harmonic_point,
       solve_harmonic_field,
       solve_harmonic_sweep

if !isdefined(@__MODULE__, :MetamaterialProfiles)
    include(joinpath(@__DIR__, "profiles.jl"))
end
using .MetamaterialProfiles

Base.@kwdef struct HarmonicConfig
    pressure_amplitude_pa::Float64 = 1.0e6
    density::Float64 = 1210.0
    pressure_wave_speed::Float64 = 2340.0
    shear_wave_speed::Float64 = 1170.0
    rayleigh_alpha::Float64 = 79560.0
    rayleigh_beta::Float64 = 2.5e-9
    element_order::Int = 1
    quadrature_degree::Int = 2
end

struct HarmonicPoint
    frequency_hz::Float64
    left_displacement::ComplexF64
    right_displacement::ComplexF64
    left_traction_pa::ComplexF64
    right_traction_pa::ComplexF64
    source_power_w_per_m::Float64
    left_absorbed_power_w_per_m::Float64
    right_absorbed_power_w_per_m::Float64
    internal_dissipated_power_w_per_m::Float64
end

struct HarmonicField
    point::HarmonicPoint
    node_x_m::Vector{Float64}
    node_y_m::Vector{Float64}
    cell_node_ids::Vector{Vector{Int}}
    displacement_x_m::Vector{ComplexF64}
    displacement_y_m::Vector{ComplexF64}
end

function solve_harmonic_state(
    model_path::AbstractString,
    frequency_hz::Real;
    config::HarmonicConfig=HarmonicConfig(),
)
    isfile(model_path) || error("Gridap model not found: $model_path")
    model = DiscreteModelFromFile(model_path)
    omega = 2pi * Float64(frequency_hz)
    rho = config.density
    mu = rho * config.shear_wave_speed^2
    lambda = rho * config.pressure_wave_speed^2 - 2mu
    identity_tensor = one(TensorValue{2, 2, ComplexF64})
    sigma(strain) = lambda * tr(strain) * identity_tensor + 2mu * strain

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
    normal_left = VectorValue(-1.0, 0.0)
    normal_right = VectorValue(1.0, 0.0)

    function impedance_displacement(displacement, normal)
        normal_part = (displacement ⋅ normal) * normal
        tangent_part = displacement - normal_part
        rho * config.pressure_wave_speed * normal_part +
            rho * config.shear_wave_speed * tangent_part
    end

    stiffness(u, v) = sigma(ε(u)) ⊙ ε(v)
    bilinear(u, v) =
        ∫(
            stiffness(u, v) - rho * omega^2 * (u ⋅ v) +
            im * omega * (
                config.rayleigh_alpha * rho * (u ⋅ v) +
                config.rayleigh_beta * stiffness(u, v)
            ),
        )domain_measure +
        ∫(im * omega * impedance_displacement(u, normal_left) ⋅ v)left_measure +
        ∫(im * omega * impedance_displacement(u, normal_right) ⋅ v)right_measure
    linear(v) =
        ∫(-config.pressure_amplitude_pa * (normal_left ⋅ v))left_measure

    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))
    stress = sigma ∘ ε(displacement)
    traction_left = stress ⋅ normal_left
    traction_right = stress ⋅ normal_right
    velocity = im * omega * displacement
    conjugated(value) = conj(value)
    conjugate_velocity = conjugated ∘ velocity
    conjugate_displacement = conjugated ∘ displacement
    conjugate_strain = conjugated ∘ ε(displacement)
    left_length = sum(∫(1.0)left_measure)
    right_length = sum(∫(1.0)right_measure)
    source_traction = -config.pressure_amplitude_pa * normal_left
    source_power = 0.5 * real(sum(∫(source_traction ⋅ conjugate_velocity)left_measure))
    left_absorbed_power = 0.5 * omega^2 * real(sum(
        ∫(impedance_displacement(displacement, normal_left) ⋅ conjugate_displacement)left_measure,
    ))
    right_absorbed_power = 0.5 * omega^2 * real(sum(
        ∫(impedance_displacement(displacement, normal_right) ⋅ conjugate_displacement)right_measure,
    ))
    internal_dissipated_power = 0.5 * omega^2 * real(sum(
        ∫(
            config.rayleigh_alpha * rho * (displacement ⋅ conjugate_displacement) +
            config.rayleigh_beta * (sigma(ε(displacement)) ⊙ conjugate_strain)
        )domain_measure,
    ))
    point = HarmonicPoint(
        Float64(frequency_hz),
        ComplexF64(sum(∫(displacement ⋅ normal_left)left_measure) / left_length),
        ComplexF64(sum(∫(displacement ⋅ normal_right)right_measure) / right_length),
        ComplexF64(sum(∫(normal_left ⋅ traction_left)left_measure) / left_length),
        ComplexF64(sum(∫(normal_right ⋅ traction_right)right_measure) / right_length),
        Float64(source_power),
        Float64(left_absorbed_power),
        Float64(right_absorbed_power),
        Float64(internal_dissipated_power),
    )
    # The additional fields are intentionally kept in this internal state so
    # modal post-processing can use distributed boundary displacement and
    # traction instead of the legacy cross-section averages.
    (; point, model, domain, displacement, stress, omega)
end

function solve_harmonic_point(
    model_path::AbstractString,
    frequency_hz::Real;
    config::HarmonicConfig=HarmonicConfig(),
)
    solve_harmonic_state(model_path, frequency_hz; config).point
end

function solve_harmonic_field(
    model_path::AbstractString,
    frequency_hz::Real;
    config::HarmonicConfig=HarmonicConfig(),
    visualization_order::Integer=1,
)
    state = solve_harmonic_state(model_path, frequency_hz; config)
    data = only(Gridap.Visualization.visualization_data(
        state.domain,
        "harmonic_field";
        order=visualization_order,
        cellfields=Dict("displacement" => state.displacement),
    ))
    coordinates = collect(Gridap.Geometry.get_node_coordinates(data.grid))
    displacement = collect(data.nodaldata["displacement"])
    HarmonicField(
        state.point,
        Float64[coordinate[1] for coordinate in coordinates],
        Float64[coordinate[2] for coordinate in coordinates],
        [Int.(collect(ids)) for ids in Gridap.Geometry.get_cell_node_ids(data.grid)],
        ComplexF64[value[1] for value in displacement],
        ComplexF64[value[2] for value in displacement],
    )
end

function solve_harmonic_sweep(
    model_path::AbstractString,
    frequencies_hz::AbstractVector{<:Real};
    config::HarmonicConfig=HarmonicConfig(),
)
    [solve_harmonic_point(model_path, frequency_hz; config) for frequency_hz in frequencies_hz]
end

end
