module TopologyHybridMonostaticTransient

using Gridap
using JLD2
using LinearAlgebra
using SparseArrays

include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
using .TopologyHybridTTDGRINDesign

export HybridMonostaticTransientConfig,
       hann_burst,
       solve_hybrid_monostatic_transient,
       nominal_incident_energy_j_per_m

Base.@kwdef struct HybridMonostaticTransientConfig
    frequency_hz::Float64 = 242.0e3
    pulse_cycles::Float64 = 5.0
    final_time_s::Float64 = 80.0e-6
    samples_per_period::Int = 18
    pressure_amplitude_pa::Float64 = 1.0e6
    reference_aperture_m::Float64 = 71.6e-3
    element_order::Int = 1
    quadrature_degree::Int = 2
    rayleigh_alpha_per_s::Float64 = 0.0
    rayleigh_beta_s::Float64 = 0.0
end

function validate(config::HybridMonostaticTransientConfig)
    config.frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    config.pulse_cycles > 0 || throw(ArgumentError("cycles must be positive"))
    config.final_time_s > config.pulse_cycles / config.frequency_hz ||
        throw(ArgumentError("final time does not contain the drive"))
    config.samples_per_period >= 12 || throw(ArgumentError("time grid is too coarse"))
    config.pressure_amplitude_pa > 0 || throw(ArgumentError("pressure must be positive"))
    config.reference_aperture_m > 0 || throw(ArgumentError("aperture must be positive"))
    config.element_order in (1, 2) || throw(ArgumentError("unsupported element order"))
    config
end

function hann_burst(t_s::Real, config::HybridMonostaticTransientConfig)
    duration = config.pulse_cycles / config.frequency_hz
    0 <= t_s < duration || return 0.0
    0.5 * (1 - cospi(2Float64(t_s) / duration)) *
    sinpi(2config.frequency_hz * Float64(t_s))
end

function time_derivative(time, values)
    derivative = zeros(Float64, length(values))
    length(values) >= 2 || return derivative
    derivative[1] = (values[2] - values[1]) / (time[2] - time[1])
    for index in 2:(length(values) - 1)
        derivative[index] = (values[index + 1] - values[index - 1]) /
                            (time[index + 1] - time[index - 1])
    end
    derivative[end] = (values[end] - values[end - 1]) / (time[end] - time[end - 1])
    derivative
end

function elastic_material(fraction, design)
    guided = guided_material(fraction, design)
    young = (1 - fraction) * aluminium_6061().young_modulus_pa +
            fraction * copper_c110().young_modulus_pa
    poisson = (1 - fraction) * aluminium_6061().poisson_ratio +
              fraction * copper_c110().poisson_ratio
    rho = guided.density_kg_m3
    mu = young / (2(1 + poisson))
    lambda = young * poisson / ((1 + poisson) * (1 - 2poisson))
    cp = sqrt((lambda + 2mu) / rho)
    cs = sqrt(mu / rho)
    identity = one(TensorValue{2, 2, Float64})
    sigma(strain) = lambda * tr(strain) * identity + 2mu * strain
    (; fraction=Float64(fraction), rho, young, poisson, mu, lambda, cp, cs, sigma)
end

function absorbing_traction(velocity, normal, material)
    normal_velocity = (velocity ⋅ normal) * normal
    tangent_velocity = velocity - normal_velocity
    material.rho * material.cp * normal_velocity +
    material.rho * material.cs * tangent_velocity
end

function nominal_incident_energy_j_per_m(
    source_length_m::Real,
    pressure_scale::Real,
    config::HybridMonostaticTransientConfig=HybridMonostaticTransientConfig(),
)
    validate(config)
    al = elastic_material(0.0, HybridTTDGRINConfig())
    duration = config.pulse_cycles / config.frequency_hz
    dt = duration / 10000
    time = 0:dt:duration
    # The characteristic boundary forcing equals twice the incident pressure.
    incident_pressure = 0.5config.pressure_amplitude_pa * Float64(pressure_scale)
    Float64(source_length_m) / (al.rho * al.cp) *
        sum(abs2, incident_pressure .* hann_burst.(time, Ref(config))) * dt
end

function solve_hybrid_monostatic_transient(
    model_path::AbstractString,
    output_path::AbstractString;
    topology::Symbol,
    defect_id::Symbol,
    material_fractions::AbstractVector{<:Real}=[0.0],
    design::HybridTTDGRINConfig=HybridTTDGRINConfig(),
    config::HybridMonostaticTransientConfig=HybridMonostaticTransientConfig(),
)
    isfile(model_path) || error("model not found: $model_path")
    validate(config)
    mkpath(dirname(output_path))
    model = DiscreteModelFromFile(model_path)
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(model, reference; conformity=:H1)
    trial_space = TrialFESpace(test_space)

    contributions = NamedTuple[]
    for fraction in sort(unique(Float64.(material_fractions)))
        key = round(Int, 1000fraction)
        tag = "MaterialF$(lpad(string(key), 4, '0'))"
        triangulation = Triangulation(model; tags=[tag])
        measure = Measure(triangulation, config.quadrature_degree)
        push!(contributions, (material=elastic_material(fraction, design), measure, tag))
    end
    isempty(contributions) && error("no material domains")

    source = BoundaryTriangulation(model; tags=["Source"])
    radiation = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    source_measure = Measure(source, config.quadrature_degree)
    radiation_measure = Measure(radiation, config.quadrature_degree)
    source_normal = VectorValue(-1.0, 0.0)
    source_tangent = VectorValue(0.0, 1.0)
    radiation_normal = get_normal_vector(radiation)
    al = elastic_material(0.0, design)
    source_length_m = sum(∫(1.0)source_measure)
    pressure_scale = sqrt(config.reference_aperture_m / source_length_m)
    drive(t) = hann_burst(t, config)

    function stiffness_form(u, v)
        first_part = first(contributions)
        material = first_part.material
        result = ∫(material.sigma(ε(u)) ⊙ ε(v))first_part.measure
        for part in Iterators.drop(contributions, 1)
            material = part.material
            result = result + ∫(material.sigma(ε(u)) ⊙ ε(v))part.measure
        end
        result
    end
    function mass_form(u, v)
        first_part = first(contributions)
        result = ∫(first_part.material.rho * u ⋅ v)first_part.measure
        for part in Iterators.drop(contributions, 1)
            result = result + ∫(part.material.rho * u ⋅ v)part.measure
        end
        result
    end
    function damping_form(u, v)
        first_part = first(contributions)
        material = first_part.material
        result = ∫(
            config.rayleigh_alpha_per_s * material.rho * u ⋅ v +
            config.rayleigh_beta_s * (material.sigma(ε(u)) ⊙ ε(v)),
        )first_part.measure
        for part in Iterators.drop(contributions, 1)
            material = part.material
            result = result + ∫(
                config.rayleigh_alpha_per_s * material.rho * u ⋅ v +
                config.rayleigh_beta_s * (material.sigma(ε(u)) ⊙ ε(v)),
            )part.measure
        end
        result +
        ∫(absorbing_traction(u, source_normal, al) ⋅ v)source_measure +
        ∫(absorbing_traction(u, radiation_normal, al) ⋅ v)radiation_measure
    end
    zero_field = VectorValue(0.0, 0.0)
    first_measure = first(contributions).measure
    zero_linear(v) = ∫(zero_field ⋅ v)first_measure
    stiffness_operator = AffineFEOperator(stiffness_form, zero_linear, trial_space, test_space)
    mass_operator = AffineFEOperator(mass_form, zero_linear, trial_space, test_space)
    damping_operator = AffineFEOperator(damping_form, zero_linear, trial_space, test_space)
    stiffness = sparse(get_matrix(stiffness_operator))
    mass = sparse(get_matrix(mass_operator))
    damping = sparse(get_matrix(damping_operator))

    # Assemble the spatial source and the two modal receiver row vectors once.
    source_linear(v) = ∫(source_normal ⋅ v)source_measure
    tangent_linear(v) = ∫(source_tangent ⋅ v)source_measure
    source_operator = AffineFEOperator(stiffness_form, source_linear, trial_space, test_space)
    tangent_operator = AffineFEOperator(stiffness_form, tangent_linear, trial_space, test_space)
    source_vector = get_vector(source_operator)
    tangent_vector = get_vector(tangent_operator)

    dt_s = 1 / (config.frequency_hz * config.samples_per_period)
    beta, gamma = 0.25, 0.5
    effective_acceleration = mass + gamma * dt_s * damping + beta * dt_s^2 * stiffness
    started_ns = time_ns()
    factorization = lu(effective_acceleration)
    dof_count = size(mass, 1)
    displacement = zeros(Float64, dof_count)
    velocity = zeros(Float64, dof_count)
    acceleration = zeros(Float64, dof_count)
    time_s = collect(0.0:dt_s:config.final_time_s)
    source_normal_displacement_m = Float64[]
    source_tangent_displacement_m = Float64[]
    source_normal_velocity_m_s = Float64[]
    source_tangent_velocity_m_s = Float64[]
    for (index, time) in enumerate(time_s)
        push!(source_normal_displacement_m, dot(source_vector, displacement) / source_length_m)
        push!(source_tangent_displacement_m, dot(tangent_vector, displacement) / source_length_m)
        push!(source_normal_velocity_m_s, dot(source_vector, velocity) / source_length_m)
        push!(source_tangent_velocity_m_s, dot(tangent_vector, velocity) / source_length_m)
        index == length(time_s) && break
        next_time = time_s[index + 1]
        predicted_displacement = displacement + dt_s * velocity +
                                 dt_s^2 * (0.5 - beta) * acceleration
        predicted_velocity = velocity + dt_s * (1 - gamma) * acceleration
        load = config.pressure_amplitude_pa * pressure_scale * drive(next_time) .* source_vector
        rhs = load - damping * predicted_velocity - stiffness * predicted_displacement
        next_acceleration = factorization \ rhs
        displacement = predicted_displacement + beta * dt_s^2 * next_acceleration
        velocity = predicted_velocity + gamma * dt_s * next_acceleration
        acceleration = next_acceleration
    end
    elapsed_s = (time_ns() - started_ns) * 1e-9
    source_drive_pa = config.pressure_amplitude_pa * pressure_scale .* drive.(time_s)
    incident_energy_j_per_m = nominal_incident_energy_j_per_m(
        source_length_m, pressure_scale, config,
    )
    jldsave(
        output_path;
        topology,
        defect_id,
        material_fractions=Float64.(material_fractions),
        time_s,
        source_drive_pa,
        source_normal_displacement_m,
        source_tangent_displacement_m,
        source_normal_velocity_m_s,
        source_tangent_velocity_m_s,
        source_length_m,
        pressure_scale,
        incident_energy_j_per_m,
        dt_s,
        dof_count,
        stiffness_nnz=nnz(stiffness),
        mass_nnz=nnz(mass),
        damping_nnz=nnz(damping),
        elapsed_s,
        config,
        design,
    )
    println("[+] $(topology)/$(defect_id): $(length(time_s)) steps, $(round(elapsed_s; digits=2)) s -> $output_path")
    output_path
end

end
