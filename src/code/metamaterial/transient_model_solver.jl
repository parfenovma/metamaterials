module TransientModelElasticity

using Gridap
using Gridap.ODEs
using JLD2

export TransientMaterialConfig,
       TransientModelConfig,
       pulse_signal,
       time_derivative,
       run_transient_model

Base.@kwdef struct TransientMaterialConfig
    density::Float64 = 1210.0
    pressure_wave_speed::Float64 = 2340.0
    shear_wave_speed::Float64 = 1170.0
    rayleigh_alpha::Float64 = 79560.0
    rayleigh_beta::Float64 = 2.5e-9
end

Base.@kwdef struct TransientModelConfig
    frequency_hz::Float64 = 242.0e3
    pulse_cycles::Float64 = 5.0
    final_time_s::Float64 = 80.0e-6
    samples_per_period::Int = 30
    pressure_amplitude_pa::Float64 = 1.0e6
    element_order::Int = 1
    quadrature_degree::Int = 2
end

function pulse_signal(t_s::Real, config::TransientModelConfig=TransientModelConfig())
    duration = config.pulse_cycles / config.frequency_hz
    0 <= t_s < duration || return 0.0
    0.5 * (1 - cospi(2t_s / duration)) * sinpi(2config.frequency_hz * t_s)
end

function time_derivative(time::AbstractVector, values::AbstractVector)
    length(time) == length(values) || throw(DimensionMismatch("time and values must have equal lengths"))
    length(time) >= 2 || return zeros(Float64, length(time))
    derivative = similar(values, Float64)
    derivative[1] = (values[2] - values[1]) / (time[2] - time[1])
    for index in 2:(length(time) - 1)
        derivative[index] = (values[index + 1] - values[index - 1]) /
                            (time[index + 1] - time[index - 1])
    end
    derivative[end] = (values[end] - values[end - 1]) / (time[end] - time[end - 1])
    derivative
end

function absorbing_traction(velocity, normal, material)
    normal_velocity = (velocity ⋅ normal) * normal
    tangent_velocity = velocity - normal_velocity
    material.density * material.pressure_wave_speed * normal_velocity +
    material.density * material.shear_wave_speed * tangent_velocity
end

"""Run a five-cycle, fully elastic transient on any tagged Gridap model."""
function run_transient_model(
    model_path::AbstractString,
    output_path::AbstractString;
    id::Symbol,
    material::TransientMaterialConfig=TransientMaterialConfig(),
    config::TransientModelConfig=TransientModelConfig(),
)
    isfile(model_path) || error("Gridap model not found: $model_path")
    config.frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    config.pulse_cycles > 0 || throw(ArgumentError("pulse cycles must be positive"))
    config.samples_per_period >= 12 || throw(ArgumentError("transient needs at least 12 samples per period"))
    mkpath(dirname(output_path))

    model = DiscreteModelFromFile(model_path)
    rho = material.density
    mu = rho * material.shear_wave_speed^2
    lambda = rho * material.pressure_wave_speed^2 - 2mu
    sigma(strain) = lambda * tr(strain) * one(strain) + 2mu * strain

    reference_element = ReferenceFE(
        lagrangian,
        VectorValue{2, Float64},
        config.element_order,
    )
    test_space = TestFESpace(model, reference_element; conformity=:H1)
    trial_space = TransientTrialFESpace(test_space)
    domain = Triangulation(model)
    domain_measure = Measure(domain, config.quadrature_degree)
    left_port = BoundaryTriangulation(model; tags=["Source"])
    right_port = BoundaryTriangulation(model; tags=["Microphone"])
    left_measure = Measure(left_port, config.quadrature_degree)
    right_measure = Measure(right_port, config.quadrature_degree)
    normal_left = VectorValue(-1.0, 0.0)
    normal_right = VectorValue(1.0, 0.0)
    tangent = VectorValue(0.0, 1.0)
    drive(t) = pulse_signal(t, config)

    residual(t, u, v) =
        ∫(
            rho * ∂tt(u) ⋅ v +
            material.rayleigh_alpha * rho * ∂t(u) ⋅ v +
            material.rayleigh_beta * (sigma ∘ ε(∂t(u)) ⊙ ε(v)) +
            sigma ∘ ε(u) ⊙ ε(v),
        )domain_measure +
        ∫(absorbing_traction(∂t(u), normal_right, material) ⋅ v)right_measure +
        ∫(absorbing_traction(∂t(u), normal_left, material) ⋅ v)left_measure -
        ∫((config.pressure_amplitude_pa * drive(t)) * (normal_left ⋅ v))left_measure

    jacobian(t, u, du, v) = ∫(sigma ∘ ε(du) ⊙ ε(v))domain_measure
    jacobian_t(t, u, dut, v) =
        ∫(
            material.rayleigh_alpha * rho * dut ⋅ v +
            material.rayleigh_beta * (sigma ∘ ε(dut) ⊙ ε(v)),
        )domain_measure +
        ∫(absorbing_traction(dut, normal_right, material) ⋅ v)right_measure +
        ∫(absorbing_traction(dut, normal_left, material) ⋅ v)left_measure
    jacobian_tt(t, u, dutt, v) = ∫(rho * dutt ⋅ v)domain_measure

    dt_s = 1 / (config.frequency_hz * config.samples_per_period)
    operator = TransientFEOperator(
        residual,
        (jacobian, jacobian_t, jacobian_tt),
        trial_space,
        test_space,
    )
    trial_at_zero = trial_space(0.0)
    displacement_zero = interpolate_everywhere(x -> VectorValue(0.0, 0.0), trial_at_zero)
    velocity_zero = interpolate_everywhere(x -> VectorValue(0.0, 0.0), trial_at_zero)
    time_solver = Newmark(NLSolver(show_trace=false, method=:newton), dt_s, 0.5, 0.25)
    solution = solve(
        time_solver,
        operator,
        0.0,
        config.final_time_s,
        (displacement_zero, velocity_zero),
    )

    time_s = Float64[]
    left_normal_displacement_m = Float64[]
    right_normal_displacement_m = Float64[]
    right_tangent_displacement_m = Float64[]
    left_normal_traction_mpa = Float64[]
    right_normal_traction_mpa = Float64[]
    right_tangent_traction_mpa = Float64[]
    left_length = sum(∫(1.0)left_measure)
    right_length = sum(∫(1.0)right_measure)

    for (current_time, displacement) in solution
        stress = sigma ∘ ε(displacement)
        traction_left = stress ⋅ normal_left
        traction_right = stress ⋅ normal_right
        push!(time_s, current_time)
        push!(left_normal_displacement_m, sum(∫(displacement ⋅ normal_left)left_measure) / left_length)
        push!(right_normal_displacement_m, sum(∫(displacement ⋅ normal_right)right_measure) / right_length)
        push!(right_tangent_displacement_m, sum(∫(displacement ⋅ tangent)right_measure) / right_length)
        push!(left_normal_traction_mpa, sum(∫(normal_left ⋅ traction_left)left_measure) / left_length / 1e6)
        push!(right_normal_traction_mpa, sum(∫(normal_right ⋅ traction_right)right_measure) / right_length / 1e6)
        push!(right_tangent_traction_mpa, sum(∫(tangent ⋅ traction_right)right_measure) / right_length / 1e6)
    end

    left_normal_velocity_m_per_s = time_derivative(time_s, left_normal_displacement_m)
    right_normal_velocity_m_per_s = time_derivative(time_s, right_normal_displacement_m)
    right_tangent_velocity_m_per_s = time_derivative(time_s, right_tangent_displacement_m)
    source_drive_mpa = (config.pressure_amplitude_pa / 1e6) .* drive.(time_s)
    jldsave(
        output_path;
        id,
        time_s,
        source_drive_mpa,
        left_normal_displacement_m,
        right_normal_displacement_m,
        right_tangent_displacement_m,
        left_normal_velocity_m_per_s,
        right_normal_velocity_m_per_s,
        right_tangent_velocity_m_per_s,
        left_normal_traction_mpa,
        right_normal_traction_mpa,
        right_tangent_traction_mpa,
        material,
        config,
    )
    println("[+] $id transient -> $output_path")
    output_path
end

end
