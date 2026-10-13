module NDTFullWaveHarmonicSolver

using Gridap
using GridapGmsh
using LinearAlgebra
using SparseArrays: nnz

export NDTFullWaveHarmonicConfig, solve_ndt_element_response_matrix

Base.@kwdef struct NDTFullWaveHarmonicConfig
    frequency_hz::Float64 = 242.0e3
    pressure_amplitude_pa::Float64 = 1.0
    density_kg_m3::Float64 = 2700.0
    pressure_wave_speed_m_s::Float64 = 6122.102437409232
    shear_wave_speed_m_s::Float64 = 3083.810277185563
    rayleigh_alpha_s_inv::Float64 = 0.0
    rayleigh_beta_s::Float64 = 0.0
    element_count::Int = 17
    element_order::Int = 1
    quadrature_degree::Int = 2
end

function validate(config::NDTFullWaveHarmonicConfig)
    config.frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    config.pressure_amplitude_pa > 0 || throw(ArgumentError("pressure must be positive"))
    config.density_kg_m3 > 0 || throw(ArgumentError("density must be positive"))
    config.pressure_wave_speed_m_s > config.shear_wave_speed_m_s > 0 ||
        throw(ArgumentError("wave speeds must satisfy cp > cs > 0"))
    config.rayleigh_alpha_s_inv >= 0 || throw(ArgumentError("alpha must be nonnegative"))
    config.rayleigh_beta_s >= 0 || throw(ArgumentError("beta must be nonnegative"))
    config.element_count >= 3 || throw(ArgumentError("at least three elements are required"))
    config
end

function stress(config, strain)
    rho = config.density_kg_m3
    mu = rho * config.shear_wave_speed_m_s^2
    lambda = rho * config.pressure_wave_speed_m_s^2 - 2mu
    identity = one(TensorValue{2, 2, ComplexF64})
    lambda * tr(strain) * identity + 2mu * strain
end

function absorbing_displacement(displacement, normal, config)
    normal_part = (displacement ⋅ normal) * normal
    tangent_part = displacement - normal_part
    config.density_kg_m3 * config.pressure_wave_speed_m_s * normal_part +
        config.density_kg_m3 * config.shear_wave_speed_m_s * tangent_part
end

"""
Solve all elemental transmit right-hand sides after factoring one shared
complex elastic FEM matrix. Integrated receive displacement is used so the
normal-normal matrix obeys the same work-conjugate reciprocity convention as
the applied unit traction.
"""
function solve_ndt_element_response_matrix(
    mesh_path::AbstractString;
    config::NDTFullWaveHarmonicConfig=NDTFullWaveHarmonicConfig(),
)
    validate(config)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model, reference; conformity=:H1, vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)

    domain = Triangulation(model; tags=["Specimen"])
    domain_measure = Measure(domain, config.quadrature_degree)
    source = BoundaryTriangulation(model; tags=["Source"])
    radiation = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    source_measure = Measure(source, config.quadrature_degree)
    radiation_measure = Measure(radiation, config.quadrature_degree)
    normal_source = get_normal_vector(source)
    normal_radiation = get_normal_vector(radiation)
    elements = [
        BoundaryTriangulation(model; tags=["SourceElement$(index)"])
        for index in 1:config.element_count
    ]
    element_measures = [Measure(boundary, config.quadrature_degree) for boundary in elements]
    element_normals = get_normal_vector.(elements)
    tangent = VectorValue(0.0, 1.0)

    function material_integrand(u, v)
        sigma_u = stress(config, ε(u))
        sigma_u ⊙ ε(v) - config.density_kg_m3 * omega^2 * (u ⋅ v) +
            im * omega * (
                config.rayleigh_alpha_s_inv * config.density_kg_m3 * (u ⋅ v) +
                config.rayleigh_beta_s * (sigma_u ⊙ ε(v))
            )
    end
    bilinear(u, v) =
        ∫(material_integrand(u, v))domain_measure +
        ∫(im * omega * absorbing_displacement(u, normal_source, config) ⋅ v)source_measure +
        ∫(im * omega * absorbing_displacement(u, normal_radiation, config) ⋅ v)radiation_measure
    source_linear(index) = v -> ∫(
        -config.pressure_amplitude_pa * (element_normals[index] ⋅ v)
    )element_measures[index]

    matrix = assemble_matrix(bilinear, trial_space, test_space)
    right_hand_sides = [
        assemble_vector(source_linear(index), test_space)
        for index in 1:config.element_count
    ]
    linear_solver = LUSolver()
    numerical = numerical_setup(symbolic_setup(linear_solver, matrix), matrix)
    free_values = similar(first(right_hand_sides))
    normal_response_m2 = Matrix{ComplexF64}(
        undef, config.element_count, config.element_count,
    )
    tangent_response_m2 = similar(normal_response_m2)
    element_lengths_m = Float64[
        sum(∫(1.0)measure) for measure in element_measures
    ]

    for transmit_index in 1:config.element_count
        fill!(free_values, 0)
        solve!(free_values, numerical, right_hand_sides[transmit_index])
        displacement = FEFunction(trial_space, free_values)
        for receive_index in 1:config.element_count
            normal_response_m2[receive_index, transmit_index] = ComplexF64(sum(
                ∫(displacement ⋅ element_normals[receive_index])element_measures[receive_index],
            ))
            tangent_response_m2[receive_index, transmit_index] = ComplexF64(sum(
                ∫(displacement ⋅ tangent)element_measures[receive_index],
            ))
        end
    end

    (
        frequency_hz=config.frequency_hz,
        normal_response_m2,
        tangent_response_m2,
        element_lengths_m,
        matrix_size=size(matrix),
        matrix_nonzeros=nnz(matrix),
    )
end

end

