module OperatorO6BStackHarmonicSolver

using Gridap
using GridapGmsh
using LinearAlgebra
using SparseArrays: nnz
using Statistics

export O6BElasticMaterial,
       O6COrthotropicMaterial,
       O6BStackHarmonicConfig,
       solve_o6b_stack_harmonic

Base.@kwdef struct O6BElasticMaterial
    density_kg_m3::Float64
    pressure_speed_m_s::Float64
    shear_speed_m_s::Float64
    quality_factor::Float64
end

Base.@kwdef struct O6COrthotropicMaterial
    density_kg_m3::Float64
    axial_stiffness_pa::Float64
    lateral_stiffness_pa::Float64
    cross_stiffness_pa::Float64
    shear_stiffness_pa::Float64
    quality_factor::Float64
end

Base.@kwdef struct O6BStackHarmonicConfig
    frequency_hz::Float64 = 242.0e3
    aperture_height_m::Float64 = 71.6e-3
    pressure_amplitude_pa::Float64 = 1.0
    aperture_sample_count::Int = 17
    sample_points_per_bin::Int = 7
    element_order::Int = 1
    quadrature_degree::Int = 2
end

const ALUMINIUM = O6BElasticMaterial(
    density_kg_m3=2700.0,
    pressure_speed_m_s=6122.102437409232,
    shear_speed_m_s=3083.810277185563,
    quality_factor=200.0,
)

const PZT = O6BElasticMaterial(
    density_kg_m3=7750.0,
    pressure_speed_m_s=4350.0,
    shear_speed_m_s=2260.0,
    quality_factor=55.0,
)

const BOND = O6BElasticMaterial(
    density_kg_m3=1150.0,
    pressure_speed_m_s=2500.0,
    shear_speed_m_s=1100.0,
    quality_factor=20.0,
)

function stress(material, strain)
    rho = material.density_kg_m3
    mu = rho * material.shear_speed_m_s^2
    lambda = rho * material.pressure_speed_m_s^2 - 2mu
    identity = one(TensorValue{2, 2, ComplexF64})
    lambda * tr(strain) * identity + 2mu * strain
end

function stress(material::O6COrthotropicMaterial, strain)
    function pointwise_orthotropic_stress(local_strain)
        exx = local_strain[1, 1]
        eyy = local_strain[2, 2]
        exy = local_strain[1, 2]
        sxx = material.axial_stiffness_pa * exx +
              material.cross_stiffness_pa * eyy
        syy = material.cross_stiffness_pa * exx +
              material.lateral_stiffness_pa * eyy
        sxy = 2material.shear_stiffness_pa * exy
        TensorValue(sxx, sxy, sxy, syy)
    end
    pointwise_orthotropic_stress ∘ strain
end

pressure_speed(material::O6BElasticMaterial) = material.pressure_speed_m_s
shear_speed(material::O6BElasticMaterial) = material.shear_speed_m_s
pressure_speed(material::O6COrthotropicMaterial) = sqrt(
    material.axial_stiffness_pa / material.density_kg_m3,
)
shear_speed(material::O6COrthotropicMaterial) = sqrt(
    material.shear_stiffness_pa / material.density_kg_m3,
)

function absorbing_displacement(displacement, normal, material)
    normal_part = (displacement ⋅ normal) * normal
    tangent_part = displacement - normal_part
    material.density_kg_m3 * pressure_speed(material) * normal_part +
        material.density_kg_m3 * shear_speed(material) * tangent_part
end

function profile(mode_index, y, height)
    mode_index == 0 && return 1.0
    mode_index == 1 && return -sqrt(2.0) * sin(pi * y / height)
    mode_index == 2 && return -sqrt(2.0) * cos(2pi * y / height)
    throw(ArgumentError("only DCT modes 0, 1 and 2 are allowed"))
end

function solve_o6b_stack_harmonic(
    mesh_path::AbstractString;
    config::O6BStackHarmonicConfig=O6BStackHarmonicConfig(),
    pzt_material::Union{O6BElasticMaterial,O6COrthotropicMaterial}=PZT,
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    model = GmshDiscreteModel(mesh_path)
    println("[+] O6b model loaded")
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model, reference; conformity=:H1, vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)
    println("[+] O6b FE spaces ready")
    pzt = Triangulation(model; tags=["PZT"])
    bond = Triangulation(model; tags=["Bond"])
    aluminium = Triangulation(model; tags=["Aluminium"])
    d_pzt = Measure(pzt, config.quadrature_degree)
    d_bond = Measure(bond, config.quadrature_degree)
    d_aluminium = Measure(aluminium, config.quadrature_degree)
    source = BoundaryTriangulation(model; tags=["Source"])
    backing = BoundaryTriangulation(model; tags=["BackingBoundary"])
    radiation = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    d_source = Measure(source, config.quadrature_degree)
    d_backing = Measure(backing, config.quadrature_degree)
    d_radiation = Measure(radiation, config.quadrature_degree)
    n_source = get_normal_vector(source)
    n_backing = get_normal_vector(backing)
    n_radiation = get_normal_vector(radiation)
    println("[+] O6b domains and boundaries ready")

    function material_form(material, measure, u, v)
        sigma = stress(material, ε(u))
        ∫(
            (1 + im / material.quality_factor) * (sigma ⊙ ε(v)) -
            material.density_kg_m3 * omega^2 * (u ⋅ v)
        )measure
    end
    bilinear(u, v) =
        material_form(pzt_material, d_pzt, u, v) +
        material_form(BOND, d_bond, u, v) +
        material_form(ALUMINIUM, d_aluminium, u, v) +
        ∫(im * omega * absorbing_displacement(u, n_backing, pzt_material) ⋅ v)d_backing +
        ∫(im * omega * absorbing_displacement(u, n_radiation, ALUMINIUM) ⋅ v)d_radiation

    coordinate = get_physical_coordinate(source)
    source_profile(mode_index) = (point -> profile(
        mode_index, point[2], config.aperture_height_m,
    )) ∘ coordinate
    source_linear(mode_index) = v -> ∫(
        -config.pressure_amplitude_pa * source_profile(mode_index) * (n_source ⋅ v)
    )d_source

    println("[+] O6b assembling matrix")
    matrix = assemble_matrix(bilinear, trial_space, test_space)
    println("[+] O6b assembling right-hand sides")
    right_hand_sides = [
        assemble_vector(source_linear(mode_index), test_space)
        for mode_index in 0:2
    ]
    println("[+] O6b factor $(size(matrix, 1)) x $(size(matrix, 2)) at $(config.frequency_hz / 1e3) kHz")
    solver = LUSolver()
    numerical = numerical_setup(symbolic_setup(solver, matrix), matrix)
    coefficients = similar(first(right_hand_sides))
    output = zeros(ComplexF64, config.aperture_sample_count, 3)
    port_work = zeros(ComplexF64, 3, 3)
    sample_height = config.aperture_height_m / config.aperture_sample_count
    for transmit in 1:3
        println("[+] O6b solve port $transmit/3")
        fill!(coefficients, 0)
        solve!(coefficients, numerical, right_hand_sides[transmit])
        displacement = FEFunction(trial_space, coefficients)
        for bin in 1:config.aperture_sample_count
            lower = -config.aperture_height_m / 2 + (bin - 1) * sample_height
            samples = range(
                lower + sample_height / (2config.sample_points_per_bin),
                lower + sample_height - sample_height / (2config.sample_points_per_bin);
                length=config.sample_points_per_bin,
            )
            values = ComplexF64[
                displacement(Point(1.0e-9, y))[1] for y in samples
            ]
            output[bin, transmit] = mean(values)
        end
        for receive in 1:3
            port_work[receive, transmit] = dot(
                right_hand_sides[receive], coefficients,
            )
        end
    end
    (
        frequency_hz=config.frequency_hz,
        output_displacement_m=output,
        port_work_compliance=port_work,
        reciprocity_residual=norm(port_work - transpose(port_work)) / norm(port_work),
        matrix_size=size(matrix),
        matrix_nonzeros=nnz(matrix),
    )
end

end
