using Gridap
using Gridap.ODEs
using JLD2

if !isdefined(@__MODULE__, :MetamaterialProfiles)
    include(joinpath(@__DIR__, "profiles.jl"))
end
using .MetamaterialProfiles

const PRINT_LOCK = ReentrantLock()
const DATA_FORMAT_VERSION = 2

Base.@kwdef struct MaterialConfig
    density::Float64 = 1210.0
    pressure_wave_speed::Float64 = 2340.0
    shear_wave_speed::Float64 = 1170.0
    rayleigh_alpha::Float64 = 79560.0
    rayleigh_beta::Float64 = 2.5e-9
end

Base.@kwdef struct SimulationConfig
    frequencies_hz::Vector{Float64} = [100e3, 120e3, 150e3, 200e3, 220e3, 300e3, 500e3, 830e3]
    final_time_s::Float64 = 60.0e-6
    samples_per_period::Int = 30
    pulse_cycles::Float64 = 4.0
    pressure_amplitude_pa::Float64 = 1.0e6
    element_order::Int = 1
    quadrature_degree::Int = 2
    save_vtk::Bool = false
    skip_existing::Bool = false
    right_boundary_condition::Symbol = :absorbing
    model_dir::String = "1_models"
    vtk_dir::String = "2_vtks"
    signal_dir::String = "3_signals"
end

function safe_println(args...)
    lock(PRINT_LOCK) do
        println(args...)
        flush(stdout)
    end
end

function struct_description(value)
    fields = Dict{String, Any}()
    for name in fieldnames(typeof(value))
        fields[string(name)] = getfield(value, name)
    end
    string(nameof(typeof(value))), fields
end

profile_description(profile::WallProfile) = struct_description(profile)

function pulse_signal(t, frequency_hz, cycles)
    duration = cycles / frequency_hz
    t < duration ?
        0.5 * (1.0 - cos(2.0 * pi * t / duration)) * sin(2.0 * pi * frequency_hz * t) :
        0.0
end

function time_derivative(time::AbstractVector, values::AbstractVector)
    length(time) == length(values) || throw(DimensionMismatch("time and values must have equal lengths"))
    length(time) >= 2 || return zeros(Float64, length(time))

    derivative = similar(values, Float64)
    derivative[1] = (values[2] - values[1]) / (time[2] - time[1])
    for i in 2:(length(time) - 1)
        derivative[i] = (values[i + 1] - values[i - 1]) / (time[i + 1] - time[i - 1])
    end
    derivative[end] = (values[end] - values[end - 1]) / (time[end] - time[end - 1])
    derivative
end


function absorbing_traction(velocity, normal, material::MaterialConfig)
    velocity_normal = (velocity ⋅ normal) * normal
    velocity_tangent = velocity - velocity_normal
    material.density * material.pressure_wave_speed * velocity_normal +
        material.density * material.shear_wave_speed * velocity_tangent
end


function run_acoustic_simulation(
    profile::WallProfile,
    frequency_hz::Real;
    material::MaterialConfig=MaterialConfig(),
    simulation::SimulationConfig=SimulationConfig(),
)
    simulation.right_boundary_condition in (:absorbing, :free_reflecting) ||
        throw(ArgumentError(
            "right_boundary_condition must be :absorbing or :free_reflecting",
        ))
    slug = profile_slug(profile)
    model_file = joinpath(simulation.model_dir, "model_$(slug).json")
    if !isfile(model_file)
        safe_println("  [!] Skipped: Gridap model $model_file was not found; run step1b_convert_models.jl")
        return nothing
    end

    # `run_acoustic_simulation` is also a public single-job entry point used by
    # process-based sweeps, so it must not rely on `run_sweep` creating these.
    mkpath(simulation.signal_dir)
    simulation.save_vtk && mkpath(simulation.vtk_dir)

    frequency_hz = Float64(frequency_hz)
    frequency_khz = frequency_hz / 1000.0
    boundary_suffix = simulation.right_boundary_condition == :absorbing ?
                      "" : "_BC_$(simulation.right_boundary_condition)"
    save_path = joinpath(
        simulation.signal_dir,
        "data_$(slug)_F_$(frequency_khz)$(boundary_suffix).jld2",
    )
    if simulation.skip_existing && isfile(save_path)
        safe_println("  [=] Existing result: $save_path")
        return save_path
    end
    safe_println("-> Start: $slug, frequency=$(frequency_khz) kHz, thread=$(Threads.threadid())")

    model = DiscreteModelFromFile(model_file)

    rho = material.density
    cp = material.pressure_wave_speed
    cs = material.shear_wave_speed
    mu = rho * cs^2
    lambda = rho * cp^2 - 2.0 * mu

    sigma(strain) = lambda * tr(strain) * one(strain) + 2.0 * mu * strain
    drive(t) = pulse_signal(t, frequency_hz, simulation.pulse_cycles)

    reference_element = ReferenceFE(
        lagrangian,
        VectorValue{2, Float64},
        simulation.element_order,
    )
    test_space = TestFESpace(model, reference_element, conformity=:H1)
    trial_space = TransientTrialFESpace(test_space)

    domain = Triangulation(model)
    domain_measure = Measure(domain, simulation.quadrature_degree)
    left_port = BoundaryTriangulation(model, tags=["Source"])
    right_port = BoundaryTriangulation(model, tags=["Microphone"])
    left_measure = Measure(left_port, simulation.quadrature_degree)
    right_measure = Measure(right_port, simulation.quadrature_degree)

    # The current meshes have vertical ports. These explicit conventions are
    # recorded in the output and will later be replaced by modal port objects.
    normal_left = VectorValue(-1.0, 0.0)
    normal_right = VectorValue(1.0, 0.0)
    tangent = VectorValue(0.0, 1.0)
    right_absorption_weight = simulation.right_boundary_condition == :absorbing ? 1.0 : 0.0

    residual(t, u, v) =
        ∫(
            rho * ∂tt(u) ⋅ v +
            material.rayleigh_alpha * rho * ∂t(u) ⋅ v +
            material.rayleigh_beta * (sigma ∘ (ε(∂t(u))) ⊙ ε(v)) +
            sigma ∘ (ε(u)) ⊙ ε(v),
        )domain_measure +
        ∫(right_absorption_weight * absorbing_traction(∂t(u), normal_right, material) ⋅ v)right_measure +
        ∫(absorbing_traction(∂t(u), normal_left, material) ⋅ v)left_measure -
        ∫((simulation.pressure_amplitude_pa * drive(t)) * (normal_left ⋅ v))left_measure

    jacobian(t, u, du, v) = ∫(sigma ∘ (ε(du)) ⊙ ε(v))domain_measure
    jacobian_t(t, u, dut, v) =
        ∫(
            material.rayleigh_alpha * rho * dut ⋅ v +
            material.rayleigh_beta * (sigma ∘ (ε(dut)) ⊙ ε(v)),
        )domain_measure +
        ∫(right_absorption_weight * absorbing_traction(dut, normal_right, material) ⋅ v)right_measure +
        ∫(absorbing_traction(dut, normal_left, material) ⋅ v)left_measure
    jacobian_tt(t, u, dutt, v) = ∫(rho * dutt ⋅ v)domain_measure

    dt = (1.0 / frequency_hz) / simulation.samples_per_period
    operator = TransientFEOperator(
        residual,
        (jacobian, jacobian_t, jacobian_tt),
        trial_space,
        test_space,
    )

    trial_at_t0 = trial_space(0.0)
    displacement_0 = interpolate_everywhere(x -> VectorValue(0.0, 0.0), trial_at_t0)
    velocity_0 = interpolate_everywhere(x -> VectorValue(0.0, 0.0), trial_at_t0)

    nonlinear_solver = NLSolver(show_trace=false, method=:newton)
    time_solver = Newmark(nonlinear_solver, dt, 0.5, 0.25)
    solution = solve(
        time_solver,
        operator,
        0.0,
        simulation.final_time_s,
        (displacement_0, velocity_0),
    )

    time = Float64[]
    left_normal_displacement_m = Float64[]
    left_tangent_displacement_m = Float64[]
    right_normal_displacement_m = Float64[]
    right_tangent_displacement_m = Float64[]
    left_normal_traction_mpa = Float64[]
    left_tangent_traction_mpa = Float64[]
    right_normal_traction_mpa = Float64[]
    right_tangent_traction_mpa = Float64[]

    left_length = sum(∫(1.0)left_measure)
    right_length = sum(∫(1.0)right_measure)
    pvd = simulation.save_vtk ?
          createpvd(joinpath(simulation.vtk_dir, "anim_$(slug)_F_$(frequency_khz)")) :
          nothing

    for (step, (current_time, displacement)) in enumerate(solution)
        stress = sigma ∘ (ε(displacement))
        traction_left = stress ⋅ normal_left
        traction_right = stress ⋅ normal_right

        push!(time, current_time)
        push!(left_normal_displacement_m, sum(∫(displacement ⋅ normal_left)left_measure) / left_length)
        push!(left_tangent_displacement_m, sum(∫(displacement ⋅ tangent)left_measure) / left_length)
        push!(right_normal_displacement_m, sum(∫(displacement ⋅ normal_right)right_measure) / right_length)
        push!(right_tangent_displacement_m, sum(∫(displacement ⋅ tangent)right_measure) / right_length)
        push!(left_normal_traction_mpa, sum(∫(normal_left ⋅ traction_left)left_measure) / left_length / 1.0e6)
        push!(left_tangent_traction_mpa, sum(∫(tangent ⋅ traction_left)left_measure) / left_length / 1.0e6)
        push!(right_normal_traction_mpa, sum(∫(normal_right ⋅ traction_right)right_measure) / right_length / 1.0e6)
        push!(right_tangent_traction_mpa, sum(∫(tangent ⋅ traction_right)right_measure) / right_length / 1.0e6)

        if simulation.save_vtk && step % 3 == 0
            pvd[current_time] = createvtk(
                domain,
                joinpath(simulation.vtk_dir, "anim_$(slug)_F_$(frequency_khz)_$(step).vtu"),
                cellfields=["u" => displacement],
            )
        end
    end

    simulation.save_vtk && savepvd(pvd)

    source_drive_mpa = (simulation.pressure_amplitude_pa / 1.0e6) .* drive.(time)
    left_normal_velocity_m_per_s = time_derivative(time, left_normal_displacement_m)
    left_tangent_velocity_m_per_s = time_derivative(time, left_tangent_displacement_m)
    right_normal_velocity_m_per_s = time_derivative(time, right_normal_displacement_m)
    right_tangent_velocity_m_per_s = time_derivative(time, right_tangent_displacement_m)
    profile_type, profile_parameters = profile_description(profile)
    _, material_parameters = struct_description(material)
    _, simulation_parameters = struct_description(simulation)
    jldsave(
        save_path;
        data_format_version=DATA_FORMAT_VERSION,
        profile_type,
        profile_parameters,
        frequency_hz,
        time_s=time,
        source_drive_mpa,
        left_normal_displacement_m,
        left_tangent_displacement_m,
        right_normal_displacement_m,
        right_tangent_displacement_m,
        left_normal_velocity_m_per_s,
        left_tangent_velocity_m_per_s,
        right_normal_velocity_m_per_s,
        right_tangent_velocity_m_per_s,
        left_normal_traction_mpa,
        left_tangent_traction_mpa,
        right_normal_traction_mpa,
        right_tangent_traction_mpa,
        port_normal_left=Tuple(normal_left),
        port_normal_right=Tuple(normal_right),
        port_tangent=Tuple(tangent),
        material_parameters,
        simulation_parameters,
        # Compatibility aliases for the existing Pluto notebook.
        A=amplitude_mm(profile),
        freq=frequency_hz,
        time=time,
        signal_in=source_drive_mpa,
        signal_out=-right_normal_traction_mpa,
    )

    safe_println("  [+] Done: $slug, frequency=$(frequency_khz) kHz -> $save_path")
    save_path
end

function run_sweep(
    profiles::AbstractVector{<:WallProfile};
    material::MaterialConfig=MaterialConfig(),
    simulation::SimulationConfig=SimulationConfig(),
)
    mkpath(simulation.vtk_dir)
    mkpath(simulation.signal_dir)
    tasks = [(profile, frequency) for profile in profiles for frequency in simulation.frequencies_hz]

    Threads.@threads for task in tasks
        profile, frequency = task
        try
            run_acoustic_simulation(profile, frequency; material, simulation)
        catch error
            safe_println("  [ERROR] $(profile_slug(profile)), frequency=$frequency: ", sprint(showerror, error))
        end
    end
end

function main()
    amplitudes = [0.0, 0.6, 1.2, 1.6, 1.9, 2.2, 2.5, 2.6, 2.7, 2.8, 3.0, 3.2, 3.4, 3.6]
    profiles = WallProfile[LegacySinusoidalProfile(A) for A in amplitudes]
    simulation = SimulationConfig()
    println("=== Start evaluation: $(Threads.nthreads()) threads ===")
    run_sweep(profiles; simulation)
    println("=== Evaluation completed ===")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
