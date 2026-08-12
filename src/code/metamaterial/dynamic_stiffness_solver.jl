module DynamicStiffnessSolver

using Gridap

export DynamicStiffnessConfig,
       DynamicStiffnessPoint,
       solve_dynamic_stiffness,
       solve_dynamic_stiffness_sweep

Base.@kwdef struct DynamicStiffnessConfig
    density::Float64 = 1210.0
    pressure_wave_speed::Float64 = 2340.0
    shear_wave_speed::Float64 = 1170.0
    rayleigh_alpha::Float64 = 0.0
    rayleigh_beta::Float64 = 0.0
    imposed_displacement_m::Float64 = 1.0e-6
    element_order::Int = 2
    quadrature_degree::Int = 4
end

struct DynamicStiffnessPoint
    frequency_hz::Float64
    reaction_force_n_per_m::ComplexF64
    stiffness_n_per_m2::ComplexF64
    compliance_m2_per_n::ComplexF64
end

"""
Return the generalized axial reaction to a uniform imposed right-interface
displacement. Both displacement components are fixed at the two interfaces;
all other boundaries retain their natural free condition.

The reaction is evaluated with a lifting field. This avoids differentiating
boundary tractions at rounded slot ends and is exactly the work-conjugate
reaction of the prescribed displacement.
"""
function solve_dynamic_stiffness(
    model_path::AbstractString,
    frequency_hz::Real;
    component_length_m::Real,
    config::DynamicStiffnessConfig=DynamicStiffnessConfig(),
)
    isfile(model_path) || error("Gridap model not found: $model_path")
    component_length_m > 0 || throw(ArgumentError("component length must be positive"))
    config.imposed_displacement_m != 0 ||
        throw(ArgumentError("imposed displacement must be nonzero"))
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
        dirichlet_tags=["LeftInterface", "RightInterface"],
    )
    zero_displacement(_) = VectorValue(0.0 + 0im, 0.0 + 0im)
    right_displacement(_) = VectorValue(
        ComplexF64(config.imposed_displacement_m),
        0.0 + 0im,
    )
    trial_space = TrialFESpace(test_space, [zero_displacement, right_displacement])
    domain = Triangulation(model)
    measure = Measure(domain, config.quadrature_degree)

    stiffness(u, v) = sigma(ε(u)) ⊙ ε(v)
    bilinear(u, v) = ∫(
        stiffness(u, v) - rho * omega^2 * (u ⋅ v) +
        im * omega * (
            config.rayleigh_alpha * rho * (u ⋅ v) +
            config.rayleigh_beta * stiffness(u, v)
        ),
    )measure
    zero_vector = VectorValue(0.0 + 0im, 0.0 + 0im)
    linear(v) = ∫(zero_vector ⋅ v)measure
    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))

    coordinate = get_physical_coordinate(domain)
    lifting = (point -> VectorValue(
        ComplexF64(point[1] / component_length_m),
        0.0 + 0im,
    )) ∘ coordinate
    lifting_strain = TensorValue(
        ComplexF64(1 / component_length_m), 0.0 + 0im,
        0.0 + 0im, 0.0 + 0im,
    )
    reaction = sum(∫(
        sigma(ε(displacement)) ⊙ lifting_strain -
        rho * omega^2 * (displacement ⋅ lifting) +
        im * omega * (
            config.rayleigh_alpha * rho * (displacement ⋅ lifting) +
            config.rayleigh_beta * (sigma(ε(displacement)) ⊙ lifting_strain)
        ),
    )measure)
    stiffness_value = ComplexF64(reaction / config.imposed_displacement_m)
    compliance = inv(stiffness_value)
    DynamicStiffnessPoint(
        Float64(frequency_hz),
        ComplexF64(reaction),
        stiffness_value,
        compliance,
    )
end

function solve_dynamic_stiffness_sweep(
    model_path::AbstractString,
    frequencies_hz::AbstractVector{<:Real};
    component_length_m::Real,
    config::DynamicStiffnessConfig=DynamicStiffnessConfig(),
)
    [
        solve_dynamic_stiffness(
            model_path,
            frequency_hz;
            component_length_m,
            config,
        ) for frequency_hz in frequencies_hz
    ]
end

end
