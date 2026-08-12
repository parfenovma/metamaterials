module PassiveFeedSplitterHarmonicSolver

using Gridap
using GridapGmsh

if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "sinusoidal_material_lens.jl"),
    )
end
using ..SinusoidalMaterialLens

export PassiveFeedSplitterHarmonicConfig, solve_passive_feed_splitter_harmonic

Base.@kwdef struct PassiveFeedSplitterHarmonicConfig
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

function absorbing_displacement(displacement, normal, material)
    normal_part = (displacement ⋅ normal) * normal
    tangent_part = displacement - normal_part
    material.density_kg_m3 * material.pressure_wave_speed_m_s * normal_part +
    material.density_kg_m3 * material.shear_wave_speed_m_s * tangent_part
end

function solve_passive_feed_splitter_harmonic(
    mesh_path::AbstractString,
    material::ElasticMaterial=aluminium_6061();
    config::PassiveFeedSplitterHarmonicConfig=PassiveFeedSplitterHarmonicConfig(),
)
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model, reference; conformity=:H1, vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)
    domain_measure = Measure(Triangulation(model; tags=["Domain"]), config.quadrature_degree)
    source = BoundaryTriangulation(model; tags=["Source"])
    small = BoundaryTriangulation(model; tags=["OutputSmall"])
    large = BoundaryTriangulation(model; tags=["OutputLarge"])
    source_measure = Measure(source, config.quadrature_degree)
    small_measure = Measure(small, config.quadrature_degree)
    large_measure = Measure(large, config.quadrature_degree)
    normal_source = get_normal_vector(source)
    normal_small = get_normal_vector(small)
    normal_large = get_normal_vector(large)
    rho = material.density_kg_m3
    function damping_form(u, v)
        sigma_u = stress(material, ε(u))
        material.rayleigh_alpha_s_inv * rho * (u ⋅ v) +
        material.rayleigh_beta_s * (sigma_u ⊙ ε(v))
    end
    integrand(u, v) =
        stress(material, ε(u)) ⊙ ε(v) - rho * omega^2 * (u ⋅ v) +
        im * omega * damping_form(u, v)
    bilinear(u, v) =
        ∫(integrand(u, v))domain_measure +
        ∫(im * omega * absorbing_displacement(u, normal_source, material) ⋅ v)source_measure +
        ∫(im * omega * absorbing_displacement(u, normal_small, material) ⋅ v)small_measure +
        ∫(im * omega * absorbing_displacement(u, normal_large, material) ⋅ v)large_measure
    linear(v) = ∫(-config.pressure_amplitude_pa * (normal_source ⋅ v))source_measure
    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))
    source_length_m = sum(∫(1.0)source_measure)
    small_length_m = sum(∫(1.0)small_measure)
    large_length_m = sum(∫(1.0)large_measure)
    source_integral = sum(∫(displacement ⋅ normal_source)source_measure)
    axis_x = VectorValue(1.0, 0.0)
    axis_y = VectorValue(0.0, 1.0)
    small_ux = sum(∫(displacement ⋅ axis_x)small_measure) / small_length_m
    small_uy = sum(∫(displacement ⋅ axis_y)small_measure) / small_length_m
    large_ux = sum(∫(displacement ⋅ axis_x)large_measure) / large_length_m
    large_uy = sum(∫(displacement ⋅ axis_y)large_measure) / large_length_m
    input_power = abs(real(
        0.5 * (-config.pressure_amplitude_pa) *
        conj(im * omega * source_integral),
    ))
    function absorbed_boundary_power(boundary_normal, boundary_measure)
        impedance_displacement = absorbing_displacement(
            displacement, boundary_normal, material,
        )
        0.5 * omega^2 * sum(
            ∫(real(conj(displacement) ⋅ impedance_displacement))boundary_measure
        )
    end
    small_absorbed_power = absorbed_boundary_power(normal_small, small_measure)
    large_absorbed_power = absorbed_boundary_power(normal_large, large_measure)
    function modal_power_proxy(ux, uy, width_m)
        0.5 * rho * omega^2 * width_m * (
            material.pressure_wave_speed_m_s * abs2(ux) +
            material.shear_wave_speed_m_s * abs2(uy)
        )
    end
    small_power = modal_power_proxy(small_ux, small_uy, small_length_m)
    large_power = modal_power_proxy(large_ux, large_uy, large_length_m)
    (
        source_length_m,
        small_length_m,
        large_length_m,
        small_displacement_m=ComplexF64[small_ux, small_uy],
        large_displacement_m=ComplexF64[large_ux, large_uy],
        active_input_power_w_per_m=input_power,
        small_absorbed_power_w_per_m=small_absorbed_power,
        large_absorbed_power_w_per_m=large_absorbed_power,
        small_modal_power_proxy_w_per_m=small_power,
        large_modal_power_proxy_w_per_m=large_power,
        displacement,
        model,
    )
end

end
