module MeasuredApertureHarmonicSolver

using Gridap
using GridapGmsh

if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
end
using ..SinusoidalMaterialLens

export MeasuredApertureHarmonicConfig, solve_measured_aperture_harmonic

Base.@kwdef struct MeasuredApertureHarmonicConfig
    frequency_hz::Float64 = 242.0e3
    pressure_amplitude_pa::Float64 = 1.0e6
    focal_distance_mm::Float64 = 35.0
    element_order::Int = 1
    quadrature_degree::Int = 2
    scan_x_min_mm::Float64 = 15.0
    scan_x_max_mm::Float64 = 55.0
    scan_y_half_width_mm::Float64 = 20.0
    scan_step_mm::Float64 = 0.5
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

function solve_measured_aperture_harmonic(
    mesh_path::AbstractString,
    weights::AbstractVector{<:Complex},
    material::ElasticMaterial=photopolymer();
    config::MeasuredApertureHarmonicConfig=MeasuredApertureHarmonicConfig(),
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(model, reference; conformity=:H1, vector_type=Vector{ComplexF64})
    trial_space = TrialFESpace(test_space)
    domain = Triangulation(model; tags=["Domain"])
    domain_measure = Measure(domain, config.quadrature_degree)
    radiation = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    radiation_measure = Measure(radiation, config.quadrature_degree)
    normal_radiation = get_normal_vector(radiation)

    radiator_boundaries = [
        BoundaryTriangulation(model; tags=["Radiator_$index"])
        for index in eachindex(weights)
    ]
    radiator_measures = [Measure(boundary, config.quadrature_degree) for boundary in radiator_boundaries]
    radiator_normals = get_normal_vector.(radiator_boundaries)

    rho = material.density_kg_m3
    material_integrand(u, v) =
        stress(material, ε(u)) ⊙ ε(v) - rho * omega^2 * (u ⋅ v) +
        im * omega * damping_form(material, u, v)
    bilinear(u, v) =
        ∫(material_integrand(u, v))domain_measure +
        ∫(im * omega * absorbing_displacement(u, normal_radiation, material) ⋅ v)radiation_measure
    function linear(v)
        sum(eachindex(weights)) do index
            weight = ComplexF64(weights[index])
            ∫(-config.pressure_amplitude_pa * weight *
              (radiator_normals[index] ⋅ v))radiator_measures[index]
        end
    end
    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))

    focus_value = displacement(Point(config.focal_distance_mm * 1e-3, 0.0))
    x_mm = collect(config.scan_x_min_mm:config.scan_step_mm:config.scan_x_max_mm)
    y_mm = collect(-config.scan_y_half_width_mm:config.scan_step_mm:config.scan_y_half_width_mm)
    ux_m = fill(ComplexF64(NaN, NaN), length(y_mm), length(x_mm))
    uy_m = similar(ux_m)
    for (x_index, x) in enumerate(x_mm), (y_index, y) in enumerate(y_mm)
        value = displacement(Point(x * 1e-3, y * 1e-3))
        ux_m[y_index, x_index] = value[1]
        uy_m[y_index, x_index] = value[2]
    end
    (
        focus_displacement=ComplexF64[focus_value[1], focus_value[2]],
        focus_longitudinal_amplitude_m=abs(focus_value[1]),
        focus_total_amplitude_m=hypot(abs(focus_value[1]), abs(focus_value[2])),
        scan_x_mm=x_mm,
        scan_y_mm=y_mm,
        scan_ux_m=ux_m,
        scan_uy_m=uy_m,
        scan_total_amplitude_m=sqrt.(abs2.(ux_m) .+ abs2.(uy_m)),
        displacement,
        model,
    )
end

end
