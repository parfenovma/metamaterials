module HornPointRadiatorHarmonicSolver

using Gridap
using GridapGmsh

if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "sinusoidal_material_lens.jl"),
    )
end
using ..SinusoidalMaterialLens

export HornPointRadiatorHarmonicConfig, solve_horn_point_radiator_harmonic

Base.@kwdef struct HornPointRadiatorHarmonicConfig
    frequency_hz::Float64 = 242.0e3
    pressure_amplitude_pa::Float64 = 1.0e6
    element_order::Int = 1
    quadrature_degree::Int = 2
end

function stress(material::ElasticMaterial, strain)
    rho = material.density_kg_m3
    mu = rho * material.shear_wave_speed_m_s^2
    lambda = rho * material.pressure_wave_speed_m_s^2 - 2mu
    identity = one(TensorValue{2, 2, ComplexF64})
    lambda * tr(strain) * identity + 2mu * strain
end

function damping_form(material, u, v)
    rho = material.density_kg_m3
    sigma_u = stress(material, ε(u))
    material.rayleigh_alpha_s_inv * rho * (u ⋅ v) +
    material.rayleigh_beta_s * (sigma_u ⊙ ε(v))
end

function absorbing_displacement(displacement, normal, material)
    normal_part = (displacement ⋅ normal) * normal
    tangent_part = displacement - normal_part
    material.density_kg_m3 * material.pressure_wave_speed_m_s * normal_part +
    material.density_kg_m3 * material.shear_wave_speed_m_s * tangent_part
end

function solve_horn_point_radiator_harmonic(
    mesh_path::AbstractString,
    target_x_mm::Real,
    target_y_mm::Real,
    material::ElasticMaterial=aluminium_6061();
    config::HornPointRadiatorHarmonicConfig=HornPointRadiatorHarmonicConfig(),
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model, reference; conformity=:H1, vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)
    domain_measure = Measure(Triangulation(model; tags=["Domain"]), config.quadrature_degree)
    source = BoundaryTriangulation(model; tags=["Source"])
    radiation = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    source_measure = Measure(source, config.quadrature_degree)
    radiation_measure = Measure(radiation, config.quadrature_degree)
    normal_source = get_normal_vector(source)
    normal_radiation = get_normal_vector(radiation)
    rho = material.density_kg_m3
    material_integrand(u, v) =
        stress(material, ε(u)) ⊙ ε(v) - rho * omega^2 * (u ⋅ v) +
        im * omega * damping_form(material, u, v)
    bilinear(u, v) =
        ∫(material_integrand(u, v))domain_measure +
        ∫(im * omega * absorbing_displacement(u, normal_source, material) ⋅ v)source_measure +
        ∫(im * omega * absorbing_displacement(u, normal_radiation, material) ⋅ v)radiation_measure
    linear(v) = ∫(-config.pressure_amplitude_pa * (normal_source ⋅ v))source_measure
    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))
    target = displacement(Point(Float64(target_x_mm) * 1e-3, Float64(target_y_mm) * 1e-3))
    source_integral_m2 = sum(∫(displacement ⋅ normal_source)source_measure)
    source_velocity_integral_m2_s = im * omega * source_integral_m2
    complex_input_power_w_per_m =
        0.5 * (-config.pressure_amplitude_pa) * conj(source_velocity_integral_m2_s)
    (
        target_displacement_m=ComplexF64[target[1], target[2]],
        source_normal_displacement_integral_m2=source_integral_m2,
        active_input_power_w_per_m=real(complex_input_power_w_per_m),
        reactive_input_power_var_per_m=imag(complex_input_power_w_per_m),
        displacement,
        model,
    )
end

end
