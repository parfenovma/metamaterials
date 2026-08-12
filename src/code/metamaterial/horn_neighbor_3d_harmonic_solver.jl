module HornNeighbor3DHarmonicSolver

using Gridap
using GridapGmsh

if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
end
using ..SinusoidalMaterialLens

export HornNeighbor3DHarmonicConfig, solve_horn_neighbor_3d_harmonic

Base.@kwdef struct HornNeighbor3DHarmonicConfig
    frequency_hz::Float64 = 242.0e3
    pressure_amplitude_pa::Float64 = 1.0e6
    element_order::Int = 1
    quadrature_degree::Int = 2
end

function stress(material::ElasticMaterial, strain)
    rho = material.density_kg_m3
    mu = rho * material.shear_wave_speed_m_s^2
    lambda = rho * material.pressure_wave_speed_m_s^2 - 2mu
    identity = one(TensorValue{3, 3, ComplexF64})
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

function solve_horn_neighbor_3d_harmonic(
    mesh_path::AbstractString,
    source_count::Integer,
    probe_points_mm,
    material::ElasticMaterial=photopolymer();
    config::HornNeighbor3DHarmonicConfig=HornNeighbor3DHarmonicConfig(),
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    source_count > 0 || throw(ArgumentError("source count must be positive"))
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{3, Float64}, config.element_order)
    test_space = TestFESpace(model, reference; conformity=:H1, vector_type=Vector{ComplexF64})
    trial_space = TrialFESpace(test_space)
    domain = Triangulation(model; tags=["Domain"])
    domain_measure = Measure(domain, config.quadrature_degree)
    radiation = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    radiation_measure = Measure(radiation, config.quadrature_degree)
    normal_radiation = get_normal_vector(radiation)
    source_boundaries = [
        BoundaryTriangulation(model; tags=["Source_$index"])
        for index in 1:source_count
    ]
    source_measures = [Measure(boundary, config.quadrature_degree) for boundary in source_boundaries]
    source_normals = get_normal_vector.(source_boundaries)

    rho = material.density_kg_m3
    material_integrand(u, v) =
        stress(material, ε(u)) ⊙ ε(v) - rho * omega^2 * (u ⋅ v) +
        im * omega * damping_form(material, u, v)
    bilinear(u, v) =
        ∫(material_integrand(u, v))domain_measure +
        ∫(im * omega * absorbing_displacement(u, normal_radiation, material) ⋅ v)radiation_measure
    function linear(v)
        sum(1:source_count) do index
            ∫(-config.pressure_amplitude_pa *
              (source_normals[index] ⋅ v))source_measures[index]
        end
    end
    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))
    values = Matrix{ComplexF64}(undef, length(probe_points_mm.names), 3)
    for index in eachindex(probe_points_mm.names)
        value = displacement(Point(
            probe_points_mm.x_mm[index] * 1e-3,
            probe_points_mm.y_mm[index] * 1e-3,
            probe_points_mm.z_mm[index] * 1e-3,
        ))
        values[index, :] .= (value[1], value[2], value[3])
    end
    (
        probe_names=probe_points_mm.names,
        probe_x_mm=probe_points_mm.x_mm,
        probe_y_mm=probe_points_mm.y_mm,
        probe_z_mm=probe_points_mm.z_mm,
        probe_displacement_m=values,
        displacement,
        model,
    )
end

end
