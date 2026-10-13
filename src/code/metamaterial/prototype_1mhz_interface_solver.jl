module Prototype1MHzInterfaceSolver

using Gridap
using JLD2
using LinearAlgebra

include(joinpath(@__DIR__, "elastic_materials.jl"))
using .ElasticMaterials

export PrototypeConfig, solve_sweep, analyze_sweep

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT = joinpath(ROOT, "tmp", "prototype_1mhz_interface")

Base.@kwdef struct PrototypeConfig
    frequencies_hz::Vector{Float64} = collect(0.8e6:0.1e6:1.2e6)
    pressure_amplitude_pa::Float64 = 1.0e6
    element_order::Int = 1
    quadrature_degree::Int = 2
    polymer_length_mm::Float64 = 17.0 * 0.14
    observation_offset_mm::Float64 = 1.0
end

function stress(material::ElasticMaterial, strain)
    rho = material.density_kg_m3
    mu = rho * material.shear_wave_speed_m_s^2
    lambda = rho * material.pressure_wave_speed_m_s^2 - 2mu
    identity = one(TensorValue{2, 2, ComplexF64})
    lambda * tr(strain) * identity + 2mu * strain
end

function damping(material, u, v)
    material.rayleigh_alpha_s_inv * material.density_kg_m3 * (u ⋅ v) +
    material.rayleigh_beta_s * (stress(material, ε(u)) ⊙ ε(v))
end

function absorbing_displacement(displacement, normal, material)
    normal_part = (displacement ⋅ normal) * normal
    tangent_part = displacement - normal_part
    material.density_kg_m3 * material.pressure_wave_speed_m_s * normal_part +
    material.density_kg_m3 * material.shear_wave_speed_m_s * tangent_part
end

function solve_frequency(model_path, frequency_hz; config=PrototypeConfig())
    model = DiscreteModelFromFile(model_path)
    polymer = photopolymer()
    aluminium = aluminium_6061()
    omega = 2pi * frequency_hz

    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model, reference; conformity=:H1, vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)
    polymer_domain = Triangulation(model; tags=["Polymer"])
    aluminium_domain = Triangulation(model; tags=["Aluminium"])
    polymer_measure = Measure(polymer_domain, config.quadrature_degree)
    aluminium_measure = Measure(aluminium_domain, config.quadrature_degree)
    source_boundary = BoundaryTriangulation(model; tags=["Source"])
    radiation_boundary = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    source_measure = Measure(source_boundary, config.quadrature_degree)
    radiation_measure = Measure(radiation_boundary, config.quadrature_degree)
    normal_source = get_normal_vector(source_boundary)
    normal_radiation = get_normal_vector(radiation_boundary)

    function material_form(material, u, v)
        stress(material, ε(u)) ⊙ ε(v) -
        material.density_kg_m3 * omega^2 * (u ⋅ v) +
        im * omega * damping(material, u, v)
    end
    bilinear(u, v) =
        ∫(material_form(polymer, u, v))polymer_measure +
        ∫(material_form(aluminium, u, v))aluminium_measure +
        ∫(im * omega * absorbing_displacement(u, normal_source, polymer) ⋅ v)source_measure +
        ∫(im * omega * absorbing_displacement(u, normal_radiation, aluminium) ⋅ v)radiation_measure
    linear(v) = ∫(-config.pressure_amplitude_pa * (normal_source ⋅ v))source_measure

    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))
    velocity = im * omega * displacement
    conjugated(value) = conj(value)
    conjugate_velocity = conjugated ∘ velocity
    conjugate_displacement = conjugated ∘ displacement
    conjugate_strain = conjugated ∘ ε(displacement)

    source_length_m = sum(∫(1.0)source_measure)
    x_direction = VectorValue(1.0, 0.0)
    source_displacement_x_m = ComplexF64(
        sum(∫(displacement ⋅ x_direction)source_measure) / source_length_m,
    )
    observation = displacement(Point(
        (config.polymer_length_mm + config.observation_offset_mm) * 1.0e-3,
        0.0,
    ))
    observation_ux_m = ComplexF64(observation[1])
    observation_uy_m = ComplexF64(observation[2])

    source_traction = -config.pressure_amplitude_pa * normal_source
    source_power = 0.5 * real(sum(
        ∫(source_traction ⋅ conjugate_velocity)source_measure,
    ))
    source_absorbed_power = 0.5 * omega^2 * real(sum(
        ∫(absorbing_displacement(
            displacement, normal_source, polymer,
        ) ⋅ conjugate_displacement)source_measure,
    ))
    radiation_absorbed_power = 0.5 * omega^2 * real(sum(
        ∫(absorbing_displacement(
            displacement, normal_radiation, aluminium,
        ) ⋅ conjugate_displacement)radiation_measure,
    ))
    function dissipated_power(material, measure)
        0.5 * omega^2 * real(sum(∫(
            material.rayleigh_alpha_s_inv * material.density_kg_m3 *
                (displacement ⋅ conjugate_displacement) +
            material.rayleigh_beta_s *
                (stress(material, ε(displacement)) ⊙ conjugate_strain)
        )measure))
    end
    polymer_dissipated_power = dissipated_power(polymer, polymer_measure)
    aluminium_dissipated_power = dissipated_power(aluminium, aluminium_measure)
    accounted_power = source_absorbed_power + radiation_absorbed_power +
                      polymer_dissipated_power + aluminium_dissipated_power
    balance_residual = (source_power - accounted_power) / max(abs(source_power), eps())

    (
        frequency_hz=Float64(frequency_hz),
        source_power_w_per_m=Float64(source_power),
        source_absorbed_power_w_per_m=Float64(source_absorbed_power),
        radiation_absorbed_power_w_per_m=Float64(radiation_absorbed_power),
        polymer_dissipated_power_w_per_m=Float64(polymer_dissipated_power),
        aluminium_dissipated_power_w_per_m=Float64(aluminium_dissipated_power),
        balance_residual=Float64(balance_residual),
        source_displacement_x_m,
        observation_ux_m,
        observation_uy_m,
    )
end

function model_path(role, resolution)
    joinpath(OUTPUT, "models", "model_$(role)_$(resolution).json")
end

function solve_sweep(; config=PrototypeConfig())
    rows = NamedTuple[]
    for role in (:straight, :stg), resolution in (:coarse, :fine)
        path = model_path(role, resolution)
        isfile(path) || error("missing model $path; generate and convert meshes first")
        for frequency_hz in config.frequencies_hz
            println("[>] $(role) $(resolution) $(frequency_hz / 1e6) MHz")
            row = solve_frequency(path, frequency_hz; config)
            push!(rows, merge((role=String(role), resolution=String(resolution)), row))
        end
    end
    mkpath(OUTPUT)
    jldsave(joinpath(OUTPUT, "sweep.jld2"); rows, config)
    rows
end

function unwrap_phase(values)
    phases = angle.(values)
    for index in 2:length(phases)
        delta = phases[index] - phases[index - 1]
        while delta > pi
            phases[index] -= 2pi
            delta -= 2pi
        end
        while delta < -pi
            phases[index] += 2pi
            delta += 2pi
        end
    end
    phases
end

function group_delay(rows)
    ordered = sort(rows; by=row -> row.frequency_hz)
    frequencies = getproperty.(ordered, :frequency_hz)
    phases = unwrap_phase(getproperty.(ordered, :observation_ux_m))
    center = argmin(abs.(frequencies .- 1.0e6))
    1 < center < length(ordered) || error("frequency grid must bracket 1 MHz")
    -(phases[center + 1] - phases[center - 1]) /
        (2pi * (frequencies[center + 1] - frequencies[center - 1]))
end

normalized_amplitude(row) =
    abs(row.observation_ux_m) / sqrt(max(row.source_power_w_per_m, eps()))

function relative_difference(a, b)
    abs(a - b) / max(abs(b), eps())
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join(
                (string(getproperty(row, column)) for column in columns),
                ',',
            ))
        end
    end
end

function analyze_sweep(rows=JLD2.load(joinpath(OUTPUT, "sweep.jld2"), "rows"))
    delays = Dict{Tuple{String, String}, Float64}()
    for role in ("straight", "stg"), resolution in ("coarse", "fine")
        subset = filter(row -> row.role == role && row.resolution == resolution, rows)
        delays[(role, resolution)] = group_delay(subset)
    end
    carrier = Dict(
        (row.role, row.resolution) => row for row in rows
        if row.frequency_hz == 1.0e6
    )
    convergence = Dict{String, NamedTuple}()
    for role in ("straight", "stg")
        coarse = carrier[(role, "coarse")]
        fine = carrier[(role, "fine")]
        convergence[role] = (
            normalized_amplitude_relative=relative_difference(
                normalized_amplitude(coarse), normalized_amplitude(fine),
            ),
            delay_difference_s=abs(delays[(role, "coarse")] - delays[(role, "fine")]),
        )
    end
    fine_straight = carrier[("straight", "fine")]
    fine_stg = carrier[("stg", "fine")]
    equal_energy_observation_ratio =
        normalized_amplitude(fine_stg) / normalized_amplitude(fine_straight)
    max_balance_residual = maximum(abs(row.balance_residual) for row in rows)
    passed = max_balance_residual < 0.01 &&
        all(row.source_power_w_per_m > 0 for row in rows) &&
        all(value > 0 for value in values(delays)) &&
        all(value.normalized_amplitude_relative < 0.05 for value in values(convergence)) &&
        all(value.delay_difference_s < 0.15e-6 for value in values(convergence))

    write_csv(joinpath(OUTPUT, "sweep.csv"), rows)
    report = joinpath(OUTPUT, "VERDICT.md")
    open(report, "w") do io
        println(io, "# Проверка интерфейса фотополимер--алюминий на 1 МГц\n")
        println(io, "Это проверка прямого прохождения без дефектов, а не результат по разрешению дефектоскопа.\n")
        println(io, "- результат численной проверки: **$(passed ? "PASS" : "FAIL")**")
        println(io, "- максимальная невязка баланса мощности: $(max_balance_residual)")
        println(io, "- отношение осевого перемещения STG / прямая полоса при одинаковой принятой энергии: $(equal_energy_observation_ratio)")
        println(io, "\nКритерии: невязка мощности < 1%, положительная задержка, изменение нормированной амплитуды < 5% и задержки < 0.15 мкс при измельчении сетки.\n")
        println(io, "## Результаты на несущей\n")
        println(io, "| форма | сетка | мощность источника, Вт/м | мощность на поглощающих границах алюминия, Вт/м | потери в полимере, Вт/м | невязка | нормированное осевое перемещение | групповая задержка, мкс |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---:|")
        for role in ("straight", "stg"), resolution in ("coarse", "fine")
            row = carrier[(role, resolution)]
            println(io, "| $role | $resolution | $(row.source_power_w_per_m) | $(row.radiation_absorbed_power_w_per_m) | $(row.polymer_dissipated_power_w_per_m) | $(row.balance_residual) | $(normalized_amplitude(row)) | $(1e6 * delays[(role, resolution)]) |")
        end
        println(io, "\n## Сходимость по сетке\n")
        for role in ("straight", "stg")
            value = convergence[role]
            println(io, "- $role: изменение нормированной амплитуды $(100 * value.normalized_amplitude_relative)%, задержки $(1e6 * value.delay_difference_s) мкс")
        end
        println(io, "\nОтношение STG / прямая полоса выше характеризует только поле при прямом прохождении. Оно не показывает разрешение дефектов и не должно выдаваться за такой результат.")
    end
    println("[+] $report")
    (; passed, max_balance_residual, equal_energy_observation_ratio, delays, convergence, report)
end

function main(args=ARGS)
    stage = isempty(args) ? "solve" : first(args)
    if stage == "solve"
        rows = solve_sweep()
        analyze_sweep(rows)
    elseif stage == "analyze"
        analyze_sweep()
    else
        error("usage: julia prototype_1mhz_interface_solver.jl [solve|analyze]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
