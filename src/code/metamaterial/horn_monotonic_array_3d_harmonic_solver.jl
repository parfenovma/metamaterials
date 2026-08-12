module HornMonotonicArray3DHarmonicSolver

using Gridap
using GridapGmsh

if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
end
using ..SinusoidalMaterialLens

export HornMonotonicArray3DHarmonicConfig,
       solve_horn_monotonic_array_3d_harmonic,
       solve_horn_monotonic_array_3d_response_matrix

Base.@kwdef struct HornMonotonicArray3DHarmonicConfig
    frequency_hz::Float64 = 242.0e3
    pressure_amplitude_pa::Float64 = 1.0e6
    element_order::Int = 1
    quadrature_degree::Int = 2
    focal_distance_mm::Float64 = 35.0
    scan_x_min_mm::Float64 = 18.0
    scan_x_max_mm::Float64 = 52.0
    scan_y_half_width_mm::Float64 = 20.0
    scan_y_min_mm::Union{Nothing, Float64} = nothing
    scan_step_mm::Float64 = 1.0
    symmetry_y::Bool = false
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

function evaluate_inside(displacement, x_m::Real, y_m::Real, z_m::Real)
    offsets_m = (
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0e-8),
        (0.0, 1.0e-8, 0.0),
        (1.0e-8, 0.0, 0.0),
        (0.0, 0.0, -1.0e-8),
        (0.0, -1.0e-8, 0.0),
        (1.0e-6, 0.0, 0.0),
        (-1.0e-6, 0.0, 0.0),
        (0.0, 1.0e-6, 0.0),
        (0.0, -1.0e-6, 0.0),
        (0.0, 0.0, 1.0e-6),
        (0.0, 0.0, -1.0e-6),
        (1.0e-6, 1.0e-6, 1.0e-6),
        (1.0e-5, 0.0, 0.0),
        (0.0, 1.0e-5, 0.0),
        (0.0, -1.0e-5, 0.0),
        (0.0, 0.0, 1.0e-5),
        (0.0, 0.0, -1.0e-5),
        (1.0e-4, 0.0, 0.0),
        (-1.0e-4, 0.0, 0.0),
        (0.0, 1.0e-4, 0.0),
        (0.0, -1.0e-4, 0.0),
        (0.0, 0.0, 1.0e-4),
        (0.0, 0.0, -1.0e-4),
        (3.0e-4, 0.0, 0.0),
    )
    last_error = nothing
    for (dx_m, dy_m, dz_m) in offsets_m
        try
            return displacement(Point(x_m + dx_m, y_m + dy_m, z_m + dz_m))
        catch error
            error isa AssertionError || rethrow()
            last_error = error
        end
    end
    throw(last_error)
end

function solve_horn_monotonic_array_3d_harmonic(
    mesh_path::AbstractString,
    source_count::Integer,
    outlet_x_mm::Real,
    probe_points_mm,
    material::ElasticMaterial=photopolymer();
    config::HornMonotonicArray3DHarmonicConfig=HornMonotonicArray3DHarmonicConfig(),
    source_weights=nothing,
    compute_source_power::Bool=false,
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    source_count > 0 || throw(ArgumentError("source count must be positive"))
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{3, Float64}, config.element_order)
    test_space = if config.symmetry_y
        TestFESpace(
            model,
            reference;
            conformity=:H1,
            vector_type=Vector{ComplexF64},
            dirichlet_tags=["SymmetryBoundary"],
            dirichlet_masks=[(false, true, false)],
        )
    else
        TestFESpace(model, reference; conformity=:H1, vector_type=Vector{ComplexF64})
    end
    trial_space = TrialFESpace(test_space)
    domain_measure = Measure(Triangulation(model; tags=["Domain"]), config.quadrature_degree)
    radiation = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    radiation_measure = Measure(radiation, config.quadrature_degree)
    normal_radiation = get_normal_vector(radiation)
    source_boundaries = [
        BoundaryTriangulation(model; tags=["Source_$index"])
        for index in 1:source_count
    ]
    source_measures = [Measure(boundary, config.quadrature_degree) for boundary in source_boundaries]
    source_normals = get_normal_vector.(source_boundaries)
    weights = isnothing(source_weights) ? ones(ComplexF64, source_count) :
              ComplexF64.(source_weights)
    length(weights) == source_count || throw(DimensionMismatch("one weight per source"))

    rho = material.density_kg_m3
    material_integrand(u, v) =
        stress(material, ε(u)) ⊙ ε(v) - rho * omega^2 * (u ⋅ v) +
        im * omega * damping_form(material, u, v)
    bilinear(u, v) =
        ∫(material_integrand(u, v))domain_measure +
        ∫(im * omega * absorbing_displacement(u, normal_radiation, material) ⋅ v)radiation_measure
    function linear(v)
        sum(1:source_count) do index
            ∫(-config.pressure_amplitude_pa * weights[index] *
              (source_normals[index] ⋅ v))source_measures[index]
        end
    end
    println("[+] assembling $source_count-source 3D carrier operator")
    operator = AffineFEOperator(bilinear, linear, trial_space, test_space)
    println("[+] solving $source_count-source 3D carrier operator")
    displacement = solve(operator)

    probe_values = Matrix{ComplexF64}(undef, length(probe_points_mm.names), 3)
    for index in eachindex(probe_points_mm.names)
        value = evaluate_inside(
            displacement,
            probe_points_mm.x_mm[index] * 1e-3,
            probe_points_mm.y_mm[index] * 1e-3,
            probe_points_mm.z_mm[index] * 1e-3,
        )
        probe_values[index, :] .= (value[1], value[2], value[3])
    end

    scan_x_mm = collect(config.scan_x_min_mm:config.scan_step_mm:config.scan_x_max_mm)
    scan_y_min_mm = isnothing(config.scan_y_min_mm) ?
                    -config.scan_y_half_width_mm : config.scan_y_min_mm
    scan_y_mm = collect(scan_y_min_mm:config.scan_step_mm:config.scan_y_half_width_mm)
    invalid = ComplexF64(NaN, NaN)
    scan_ux_m = fill(invalid, length(scan_y_mm), length(scan_x_mm))
    scan_uy_m = fill(invalid, length(scan_y_mm), length(scan_x_mm))
    scan_uz_m = fill(invalid, length(scan_y_mm), length(scan_x_mm))
    for (x_index, relative_x_mm) in enumerate(scan_x_mm),
        (y_index, y_mm) in enumerate(scan_y_mm)
        value = try
            evaluate_inside(
                displacement,
                (Float64(outlet_x_mm) + relative_x_mm) * 1e-3,
                y_mm * 1e-3,
                0.0,
            )
        catch error
            error isa AssertionError || rethrow()
            continue
        end
        scan_ux_m[y_index, x_index] = value[1]
        scan_uy_m[y_index, x_index] = value[2]
        scan_uz_m[y_index, x_index] = value[3]
    end
    focus_value = evaluate_inside(
        displacement,
        (Float64(outlet_x_mm) + config.focal_distance_mm) * 1e-3,
        0.0,
        0.0,
    )
    source_normal_displacement_integral_m3 = if compute_source_power
        ComplexF64[
            sum(∫(displacement ⋅ source_normals[index])source_measures[index])
            for index in 1:source_count
        ]
    else
        fill(ComplexF64(NaN, NaN), source_count)
    end
    complex_input_power_va = if compute_source_power
        source_velocity_integral_m3_s =
            im * omega .* source_normal_displacement_integral_m3
        0.5sum(
            (-config.pressure_amplitude_pa .* weights) .*
            conj.(source_velocity_integral_m3_s),
        )
    else
        ComplexF64(NaN, NaN)
    end
    (
        probe_names=probe_points_mm.names,
        probe_x_mm=probe_points_mm.x_mm,
        probe_y_mm=probe_points_mm.y_mm,
        probe_z_mm=probe_points_mm.z_mm,
        probe_displacement_m=probe_values,
        focus_displacement_m=ComplexF64[focus_value[1], focus_value[2], focus_value[3]],
        scan_x_mm,
        scan_y_mm,
        scan_ux_m,
        scan_uy_m,
        scan_uz_m,
        scan_total_amplitude_m=sqrt.(abs2.(scan_ux_m) .+ abs2.(scan_uy_m) .+ abs2.(scan_uz_m)),
        source_weights=weights,
        source_normal_displacement_integral_m3,
        complex_input_power_va,
        active_input_power_w=real(complex_input_power_va),
        reactive_input_power_var=imag(complex_input_power_va),
        displacement,
        model,
    )
end

"""
Factor the frozen carrier operator once and solve one right-hand side per
source group. In a symmetry-constrained half-domain, source 1 is the centre
channel and sources 2:end represent mirrored channel pairs.
"""
function solve_horn_monotonic_array_3d_response_matrix(
    mesh_path::AbstractString,
    source_count::Integer,
    outlet_x_mm::Real,
    probe_points_mm,
    material::ElasticMaterial=photopolymer();
    config::HornMonotonicArray3DHarmonicConfig=HornMonotonicArray3DHarmonicConfig(
        scan_y_min_mm=0.0,
        symmetry_y=true,
    ),
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    source_count > 0 || throw(ArgumentError("source count must be positive"))
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{3, Float64}, config.element_order)
    test_space = if config.symmetry_y
        TestFESpace(
            model,
            reference;
            conformity=:H1,
            vector_type=Vector{ComplexF64},
            dirichlet_tags=["SymmetryBoundary"],
            dirichlet_masks=[(false, true, false)],
        )
    else
        TestFESpace(model, reference; conformity=:H1, vector_type=Vector{ComplexF64})
    end
    trial_space = TrialFESpace(test_space)
    domain_measure = Measure(Triangulation(model; tags=["Domain"]), config.quadrature_degree)
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
    source_linear(index) = v -> ∫(-config.pressure_amplitude_pa *
        (source_normals[index] ⋅ v))source_measures[index]

    println("[+] assembling shared $source_count-source 3D carrier matrix")
    matrix = assemble_matrix(bilinear, trial_space, test_space)
    rhs = [assemble_vector(source_linear(index), test_space) for index in 1:source_count]
    free_values = similar(first(rhs))
    probe_response_m = Array{ComplexF64}(undef, length(probe_points_mm.names), 3, source_count)
    focus_response_m = Matrix{ComplexF64}(undef, 3, source_count)
    scan_y_min_mm = isnothing(config.scan_y_min_mm) ? 0.0 : config.scan_y_min_mm
    scan_y_mm = collect(scan_y_min_mm:config.scan_step_mm:config.scan_y_half_width_mm)
    scan_ux_m = Matrix{ComplexF64}(undef, length(scan_y_mm), source_count)
    scan_uy_m = similar(scan_ux_m)
    scan_uz_m = similar(scan_ux_m)
    source_normal_displacement_integral_m3 = Matrix{ComplexF64}(
        undef,
        source_count,
        source_count,
    )
    println("[+] factoring shared $source_count-source 3D carrier matrix")
    linear_solver = LUSolver()
    numerical = numerical_setup(symbolic_setup(linear_solver, matrix), matrix)

    for source_index in 1:source_count
        println("[+] solving source group $source_index/$source_count")
        fill!(free_values, 0)
        solve!(free_values, numerical, rhs[source_index])
        displacement = FEFunction(trial_space, free_values)
        for probe_index in eachindex(probe_points_mm.names)
            value = evaluate_inside(
                displacement,
                probe_points_mm.x_mm[probe_index] * 1e-3,
                probe_points_mm.y_mm[probe_index] * 1e-3,
                probe_points_mm.z_mm[probe_index] * 1e-3,
            )
            probe_response_m[probe_index, :, source_index] .= (value[1], value[2], value[3])
        end
        focus_value = evaluate_inside(
            displacement,
            (Float64(outlet_x_mm) + config.focal_distance_mm) * 1e-3,
            0.0,
            0.0,
        )
        focus_response_m[:, source_index] .= (focus_value[1], focus_value[2], focus_value[3])
        for (y_index, y_mm) in enumerate(scan_y_mm)
            value = evaluate_inside(
                displacement,
                (Float64(outlet_x_mm) + config.focal_distance_mm) * 1e-3,
                y_mm * 1e-3,
                0.0,
            )
            scan_ux_m[y_index, source_index] = value[1]
            scan_uy_m[y_index, source_index] = value[2]
            scan_uz_m[y_index, source_index] = value[3]
        end
        for receiver_source in 1:source_count
            source_normal_displacement_integral_m3[receiver_source, source_index] =
                sum(∫(displacement ⋅ source_normals[receiver_source])source_measures[receiver_source])
        end
    end

    (
        probe_names=probe_points_mm.names,
        probe_x_mm=probe_points_mm.x_mm,
        probe_y_mm=probe_points_mm.y_mm,
        probe_z_mm=probe_points_mm.z_mm,
        probe_response_m,
        focus_response_m,
        scan_y_mm,
        scan_ux_m,
        scan_uy_m,
        scan_uz_m,
        source_normal_displacement_integral_m3,
    )
end

end
