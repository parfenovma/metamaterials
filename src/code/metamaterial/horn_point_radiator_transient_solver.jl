module HornPointRadiatorTransientSolver

using Gridap
using Gridap.ODEs
using JLD2

if !isdefined(parentmodule(@__MODULE__), :TransientModelElasticity)
    include(joinpath(@__DIR__, "transient_model_solver.jl"))
end
using ..TransientModelElasticity

export run_horn_point_radiator_transient

function absorbing_traction(velocity, normal, material)
    normal_velocity = (velocity ⋅ normal) * normal
    tangent_velocity = velocity - normal_velocity
    material.density * material.pressure_wave_speed * normal_velocity +
    material.density * material.shear_wave_speed * tangent_velocity
end

function rowwise_derivative(time_s, values)
    rows = [
        permutedims(TransientModelElasticity.time_derivative(time_s, vec(values[index, :])))
        for index in axes(values, 1)
    ]
    reduce(vcat, rows)
end

"""Fully elastic transient with a wide driven collector and point probes in the receiver."""
function run_horn_point_radiator_transient(
    model_path::AbstractString,
    output_path::AbstractString;
    id::Symbol,
    probe_names::AbstractVector{<:AbstractString},
    probe_x_mm::AbstractVector{<:Real},
    probe_y_mm::AbstractVector{<:Real},
    outlet_x_mm::Real,
    material::TransientMaterialConfig=TransientMaterialConfig(),
    config::TransientModelConfig=TransientModelConfig(),
)
    isfile(model_path) || error("Gridap model not found: $model_path")
    length(probe_names) == length(probe_x_mm) == length(probe_y_mm) ||
        throw(DimensionMismatch("probe metadata must have equal lengths"))
    mkpath(dirname(output_path))

    model = DiscreteModelFromFile(model_path)
    rho = material.density
    mu = rho * material.shear_wave_speed^2
    lambda = rho * material.pressure_wave_speed^2 - 2mu
    sigma(strain) = lambda * tr(strain) * one(strain) + 2mu * strain

    reference_element = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(model, reference_element; conformity=:H1)
    trial_space = TransientTrialFESpace(test_space)
    domain_measure = Measure(Triangulation(model), config.quadrature_degree)
    source_boundary = BoundaryTriangulation(model; tags=["Source"])
    radiation_boundary = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    source_measure = Measure(source_boundary, config.quadrature_degree)
    radiation_measure = Measure(radiation_boundary, config.quadrature_degree)
    normal_source = get_normal_vector(source_boundary)
    normal_radiation = get_normal_vector(radiation_boundary)
    drive(t) = TransientModelElasticity.pulse_signal(t, config)

    residual(t, u, v) =
        ∫(
            rho * ∂tt(u) ⋅ v +
            material.rayleigh_alpha * rho * ∂t(u) ⋅ v +
            material.rayleigh_beta * (sigma ∘ ε(∂t(u)) ⊙ ε(v)) +
            sigma ∘ ε(u) ⊙ ε(v),
        )domain_measure +
        ∫(absorbing_traction(∂t(u), normal_source, material) ⋅ v)source_measure +
        ∫(absorbing_traction(∂t(u), normal_radiation, material) ⋅ v)radiation_measure -
        ∫((config.pressure_amplitude_pa * drive(t)) * (normal_source ⋅ v))source_measure

    jacobian(t, u, du, v) = ∫(sigma ∘ ε(du) ⊙ ε(v))domain_measure
    jacobian_t(t, u, dut, v) =
        ∫(
            material.rayleigh_alpha * rho * dut ⋅ v +
            material.rayleigh_beta * (sigma ∘ ε(dut) ⊙ ε(v)),
        )domain_measure +
        ∫(absorbing_traction(dut, normal_source, material) ⋅ v)source_measure +
        ∫(absorbing_traction(dut, normal_radiation, material) ⋅ v)radiation_measure
    jacobian_tt(t, u, dutt, v) = ∫(rho * dutt ⋅ v)domain_measure

    dt_s = 1 / (config.frequency_hz * config.samples_per_period)
    operator = TransientFEOperator(
        residual,
        (jacobian, jacobian_t, jacobian_tt),
        trial_space,
        test_space,
    )
    trial_at_zero = trial_space(0.0)
    zero = interpolate_everywhere(x -> VectorValue(0.0, 0.0), trial_at_zero)
    time_solver = Newmark(NLSolver(show_trace=false, method=:newton), dt_s, 0.5, 0.25)
    solution = solve(time_solver, operator, 0.0, config.final_time_s, (zero, zero))

    points = [Point(Float64(x) * 1e-3, Float64(y) * 1e-3) for (x, y) in zip(probe_x_mm, probe_y_mm)]
    time_s = Float64[]
    source_normal_displacement_m = Float64[]
    probe_x = [Float64[] for _ in points]
    probe_y = [Float64[] for _ in points]
    source_length = sum(∫(1.0)source_measure)

    for (current_time, displacement) in solution
        push!(time_s, current_time)
        push!(source_normal_displacement_m,
              sum(∫(displacement ⋅ normal_source)source_measure) / source_length)
        for (index, point) in enumerate(points)
            value = displacement(point)
            push!(probe_x[index], value[1])
            push!(probe_y[index], value[2])
        end
    end

    probe_displacement_x_m = reduce(vcat, permutedims.(probe_x))
    probe_displacement_y_m = reduce(vcat, permutedims.(probe_y))
    probe_velocity_x_m_per_s = rowwise_derivative(time_s, probe_displacement_x_m)
    probe_velocity_y_m_per_s = rowwise_derivative(time_s, probe_displacement_y_m)
    source_normal_velocity_m_per_s =
        TransientModelElasticity.time_derivative(time_s, source_normal_displacement_m)
    source_drive_mpa = (config.pressure_amplitude_pa / 1e6) .* drive.(time_s)
    names = String.(probe_names)
    x_mm = Float64.(probe_x_mm)
    y_mm = Float64.(probe_y_mm)

    jldsave(
        output_path;
        format_version=1,
        id,
        time_s,
        source_drive_mpa,
        source_normal_displacement_m,
        source_normal_velocity_m_per_s,
        probe_names=names,
        probe_x_mm=x_mm,
        probe_y_mm=y_mm,
        probe_displacement_x_m,
        probe_displacement_y_m,
        probe_velocity_x_m_per_s,
        probe_velocity_y_m_per_s,
        outlet_x_mm=Float64(outlet_x_mm),
        material,
        config,
    )
    println("[+] $id horn-point transient -> $output_path")
    output_path
end

end
