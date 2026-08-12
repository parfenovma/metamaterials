module MaterialLensHarmonicSolver

using Gridap
using GridapGmsh

if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
end
using ..SinusoidalMaterialLens

export MaterialLensHarmonicConfig, solve_material_lens_harmonic

Base.@kwdef struct MaterialLensHarmonicConfig
    frequency_hz::Float64 = 242.0e3
    pressure_amplitude_pa::Float64 = 1.0e6
    element_order::Int = 2
    quadrature_degree::Int = 4
    focus_x_mm::Float64
    focus_y_mm::Float64 = 0.0
    scan_before_focus_mm::Float64 = 40.0
    scan_after_focus_mm::Float64 = 30.0
    scan_transverse_half_width_mm::Float64 = 30.0
    scan_step_mm::Float64 = 2.0
    source_excitation_tags::Vector{String} = String["Source"]
end

function stress(material::ElasticMaterial, strain)
    rho = material.density_kg_m3
    mu = rho * material.shear_wave_speed_m_s^2
    lambda = rho * material.pressure_wave_speed_m_s^2 - 2mu
    identity = one(TensorValue{2, 2, ComplexF64})
    lambda * tr(strain) * identity + 2mu * strain
end

function damping_form(material, omega, u, v)
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

function solve_material_lens_harmonic(
    mesh_path::AbstractString,
    lens_material::ElasticMaterial,
    output_material::ElasticMaterial;
    matching_material::Union{Nothing, ElasticMaterial}=nothing,
    config::MaterialLensHarmonicConfig,
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model,
        reference;
        conformity=:H1,
        vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)

    lens_domain = Triangulation(model; tags=["Lens"])
    aluminium_domain = Triangulation(model; tags=["Aluminium"])
    lens_measure = Measure(lens_domain, config.quadrature_degree)
    aluminium_measure = Measure(aluminium_domain, config.quadrature_degree)
    matching_measure = isnothing(matching_material) ? nothing :
                       Measure(Triangulation(model; tags=["MatchingLayer"]), config.quadrature_degree)
    source_boundary = BoundaryTriangulation(model; tags=["Source"])
    excitation_boundary = BoundaryTriangulation(model; tags=config.source_excitation_tags)
    radiation_boundary = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    source_measure = Measure(source_boundary, config.quadrature_degree)
    excitation_measure = Measure(excitation_boundary, config.quadrature_degree)
    radiation_measure = Measure(radiation_boundary, config.quadrature_degree)
    normal_source = get_normal_vector(source_boundary)
    normal_excitation = get_normal_vector(excitation_boundary)
    normal_radiation = get_normal_vector(radiation_boundary)

    function material_integrand(material, u, v)
        rho = material.density_kg_m3
        stress(material, ε(u)) ⊙ ε(v) - rho * omega^2 * (u ⋅ v) +
        im * omega * damping_form(material, omega, u, v)
    end
    function bilinear(u, v)
        result = ∫(material_integrand(lens_material, u, v))lens_measure +
                 ∫(material_integrand(output_material, u, v))aluminium_measure +
                 ∫(im * omega * absorbing_displacement(u, normal_source, lens_material) ⋅ v)source_measure +
                 ∫(im * omega * absorbing_displacement(u, normal_radiation, output_material) ⋅ v)radiation_measure
        if !isnothing(matching_material)
            result += ∫(material_integrand(matching_material, u, v))matching_measure
        end
        result
    end
    linear(v) = ∫(-config.pressure_amplitude_pa * (normal_excitation ⋅ v))excitation_measure
    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))

    invalid_displacement = ComplexF64[ComplexF64(NaN, NaN), ComplexF64(NaN, NaN)]
    function point_displacement(x_mm, y_mm)
        try
            value = displacement(Point(x_mm * 1.0e-3, y_mm * 1.0e-3))
            ComplexF64[value[1], value[2]]
        catch
            copy(invalid_displacement)
        end
    end
    focus_displacement = point_displacement(config.focus_x_mm, config.focus_y_mm)
    x_values_mm = collect(
        (config.focus_x_mm - config.scan_before_focus_mm):config.scan_step_mm:
        (config.focus_x_mm + config.scan_after_focus_mm),
    )
    y_values_mm = collect(
        (config.focus_y_mm - config.scan_transverse_half_width_mm):config.scan_step_mm:
        (config.focus_y_mm + config.scan_transverse_half_width_mm),
    )
    best_amplitude = -Inf
    best_x_mm, best_y_mm = config.focus_x_mm, config.focus_y_mm
    best_displacement = focus_displacement
    scan_amplitude_m = fill(NaN, length(y_values_mm), length(x_values_mm))
    scan_ux_m = fill(ComplexF64(NaN, NaN), length(y_values_mm), length(x_values_mm))
    scan_uy_m = fill(ComplexF64(NaN, NaN), length(y_values_mm), length(x_values_mm))
    for (x_index, x_mm) in enumerate(x_values_mm), (y_index, y_mm) in enumerate(y_values_mm)
        value = point_displacement(x_mm, y_mm)
        all(isfinite, real.(value)) && all(isfinite, imag.(value)) || continue
        amplitude = sqrt(abs2(value[1]) + abs2(value[2]))
        scan_amplitude_m[y_index, x_index] = amplitude
        scan_ux_m[y_index, x_index] = value[1]
        scan_uy_m[y_index, x_index] = value[2]
        if amplitude > best_amplitude
            best_amplitude = amplitude
            best_x_mm, best_y_mm = x_mm, y_mm
            best_displacement = value
        end
    end
    (
        focus_displacement,
        focus_longitudinal_amplitude_m=abs(focus_displacement[1]),
        focus_total_amplitude_m=sqrt(abs2(focus_displacement[1]) + abs2(focus_displacement[2])),
        local_peak_amplitude_m=best_amplitude,
        local_peak_x_mm=best_x_mm,
        local_peak_y_mm=best_y_mm,
        local_peak_displacement=ComplexF64[best_displacement[1], best_displacement[2]],
        scan_x_mm=x_values_mm,
        scan_y_mm=y_values_mm,
        scan_amplitude_m,
        scan_ux_m,
        scan_uy_m,
        displacement,
        model,
        lens_domain,
        aluminium_domain,
    )
end

end
