module NDTFullWaveTransientSolver

using Gridap
using GridapGmsh
using LinearAlgebra
using SparseArrays: nnz

export NDTFullWaveTransientConfig,
       validate,
       hann_burst,
       solve_ndt_element_transient_matrix

Base.@kwdef struct NDTFullWaveTransientConfig
    frequency_hz::Float64 = 242.0e3
    pulse_cycles::Float64 = 5.0
    pressure_amplitude_pa::Float64 = 1.0e6
    final_time_s::Float64 = 70.0e-6
    samples_per_cycle::Int = 30
    density_kg_m3::Float64 = 2700.0
    pressure_wave_speed_m_s::Float64 = 6122.102437409232
    shear_wave_speed_m_s::Float64 = 3083.810277185563
    rayleigh_alpha_s_inv::Float64 = 0.0
    rayleigh_beta_s::Float64 = 0.0
    element_count::Int = 17
    element_order::Int = 1
    quadrature_degree::Int = 2
    newmark_beta::Float64 = 0.25
    newmark_gamma::Float64 = 0.5
end

function validate(config::NDTFullWaveTransientConfig)
    config.frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    config.pulse_cycles > 0 || throw(ArgumentError("pulse cycle count must be positive"))
    config.pressure_amplitude_pa > 0 || throw(ArgumentError("pressure must be positive"))
    config.final_time_s > config.pulse_cycles / config.frequency_hz ||
        throw(ArgumentError("record must be longer than the source burst"))
    config.samples_per_cycle >= 10 || throw(ArgumentError("at least ten samples per cycle are required"))
    config.density_kg_m3 > 0 || throw(ArgumentError("density must be positive"))
    config.pressure_wave_speed_m_s > config.shear_wave_speed_m_s > 0 ||
        throw(ArgumentError("wave speeds must satisfy cp > cs > 0"))
    config.rayleigh_alpha_s_inv >= 0 || throw(ArgumentError("alpha must be nonnegative"))
    config.rayleigh_beta_s >= 0 || throw(ArgumentError("beta must be nonnegative"))
    config.element_count >= 3 || throw(ArgumentError("at least three elements are required"))
    config.newmark_beta == 0.25 && config.newmark_gamma == 0.5 ||
        throw(ArgumentError("F1 freezes average-acceleration Newmark parameters"))
    config
end

function hann_burst(time_s::Real, frequency_hz::Real, cycles::Real)
    duration_s = Float64(cycles) / Float64(frequency_hz)
    time = Float64(time_s)
    0 <= time <= duration_s || return 0.0
    0.5 * (1 - cos(2pi * time / duration_s)) * sin(2pi * frequency_hz * time)
end

function stress(config, strain)
    rho = config.density_kg_m3
    mu = rho * config.shear_wave_speed_m_s^2
    lambda = rho * config.pressure_wave_speed_m_s^2 - 2mu
    identity = one(TensorValue{2, 2, Float64})
    lambda * tr(strain) * identity + 2mu * strain
end

function impedance_velocity(velocity, normal, config)
    normal_part = (velocity ⋅ normal) * normal
    tangent_part = velocity - normal_part
    config.density_kg_m3 * config.pressure_wave_speed_m_s * normal_part +
        config.density_kg_m3 * config.shear_wave_speed_m_s * tangent_part
end

"""
Direct linear transient solve with one sparse factorization and all elemental
transmit experiments carried as simultaneous right-hand sides.
"""
function solve_ndt_element_transient_matrix(
    mesh_path::AbstractString;
    config::NDTFullWaveTransientConfig=NDTFullWaveTransientConfig(),
)
    validate(config)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    model = GmshDiscreteModel(mesh_path)
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(model, reference; conformity=:H1)
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

    mass_form(u, v) = ∫(config.density_kg_m3 * (u ⋅ v))domain_measure
    stiffness_form(u, v) = ∫(stress(config, ε(u)) ⊙ ε(v))domain_measure
    damping_form(u, v) =
        ∫(
            config.rayleigh_alpha_s_inv * config.density_kg_m3 * (u ⋅ v) +
            config.rayleigh_beta_s * (stress(config, ε(u)) ⊙ ε(v))
        )domain_measure +
        ∫(impedance_velocity(u, normal_source, config) ⋅ v)source_measure +
        ∫(impedance_velocity(u, normal_radiation, config) ⋅ v)radiation_measure
    source_linear(index) = v -> ∫(
        -config.pressure_amplitude_pa * (element_normals[index] ⋅ v)
    )element_measures[index]
    tangent_linear(index) = v -> ∫(tangent ⋅ v)element_measures[index]

    println("[+] assembling transient M/C/K matrices")
    flush(stdout)
    mass = assemble_matrix(mass_form, trial_space, test_space)
    damping = assemble_matrix(damping_form, trial_space, test_space)
    stiffness = assemble_matrix(stiffness_form, trial_space, test_space)
    load_matrix = hcat((
        assemble_vector(source_linear(index), test_space)
        for index in 1:config.element_count
    )...)
    tangent_observer = hcat((
        assemble_vector(tangent_linear(index), test_space)
        for index in 1:config.element_count
    )...)

    nominal_dt = inv(config.frequency_hz * config.samples_per_cycle)
    step_count = ceil(Int, config.final_time_s / nominal_dt)
    dt = config.final_time_s / step_count
    time_s = collect(0:step_count) .* dt
    drive = hann_burst.(time_s, config.frequency_hz, config.pulse_cycles)
    beta = config.newmark_beta
    gamma = config.newmark_gamma
    effective = mass + gamma * dt * damping + beta * dt^2 * stiffness
    println("[+] factoring Newmark matrix $(size(effective)) with $(nnz(effective)) nonzeros")
    flush(stdout)
    factor = cholesky(Symmetric(effective))

    dof_count = size(mass, 1)
    source_count = config.element_count
    displacement = zeros(Float64, dof_count, source_count)
    velocity = zeros(Float64, dof_count, source_count)
    acceleration = zeros(Float64, dof_count, source_count)
    predicted_displacement = zeros(Float64, dof_count, source_count)
    predicted_velocity = zeros(Float64, dof_count, source_count)
    right_hand_side = zeros(Float64, dof_count, source_count)
    normal_response_m2 = zeros(Float64, source_count, source_count, step_count + 1)
    tangent_response_m2 = zeros(Float64, source_count, source_count, step_count + 1)

    progress_stride = max(div(step_count, 5), 1)
    for step in 1:step_count
        @. predicted_displacement =
            displacement + dt * velocity + dt^2 * (0.5 - beta) * acceleration
        @. predicted_velocity = velocity + dt * (1 - gamma) * acceleration
        mul!(right_hand_side, damping, predicted_velocity, -1.0, 0.0)
        mul!(right_hand_side, stiffness, predicted_displacement, -1.0, 1.0)
        right_hand_side .+= drive[step + 1] .* load_matrix
        acceleration .= factor \ right_hand_side
        @. displacement = predicted_displacement + beta * dt^2 * acceleration
        @. velocity = predicted_velocity + gamma * dt * acceleration
        normal_response_m2[:, :, step + 1] .=
            -(transpose(load_matrix) * displacement) ./ config.pressure_amplitude_pa
        tangent_response_m2[:, :, step + 1] .= transpose(tangent_observer) * displacement
        if step % progress_stride == 0 || step == step_count
            println("[+] transient step $step/$step_count, t=$(round(1.0e6time_s[step + 1]; digits=2)) us")
            flush(stdout)
        end
    end

    (
        time_s,
        source_drive=drive,
        normal_response_m2,
        tangent_response_m2,
        dt_s=dt,
        matrix_size=size(effective),
        mass_nonzeros=nnz(mass),
        damping_nonzeros=nnz(damping),
        stiffness_nonzeros=nnz(stiffness),
        effective_nonzeros=nnz(effective),
        element_lengths_m=Float64[
            sum(∫(1.0)measure) for measure in element_measures
        ],
    )
end

end
