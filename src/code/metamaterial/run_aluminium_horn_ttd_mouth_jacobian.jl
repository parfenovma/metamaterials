module AluminiumHornTTDMouthJacobian

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_AL_HORN_TTD_MOUTH_JACOBIAN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_mouth_jacobian_242khz"),
)

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "analyze")
const CASE_NAME = argument("case", nothing)
const FREQUENCY_HZ = 242.0e3
const LOWER_FREQUENCY_HZ = 162.2e3
const ALUMINIUM_DENSITY_KG_M3 = 2700.0
const ALUMINIUM_CP_M_S = 6122.102437409232
const ALUMINIUM_CS_M_S = 3083.810277185563
const THROAT_HEIGHT_MM = 3.5
const BASE_INPUT_HEIGHT_MM = 7.0
const GUIDE_AXIAL_MM = 132.5
const MAXIMUM_PATH_MM = 155.53
const DIFFUSER_LENGTH_MM = 65.41
const DIFFUSER_OUTPUT_HEIGHT_MM = 7.0

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
using .HornPointRadiatorMesher

const WEIGHT_PATH = joinpath(
    PROJECT_ROOT,
    "tmp",
    "aluminium_horn_ttd_full_diffuser_aperture_242khz",
    "15_power_optimal_weights.csv",
)

function load_extreme_weights()
    rows = readlines(WEIGHT_PATH)[2:end]
    weights = [parse(Float64, split(row, ',')[3]) for row in rows]
    multiplicity = vcat(1.0, fill(2.0, length(weights) - 1))
    weighted_mean = sum(multiplicity .* weights) / sum(multiplicity)
    (; low=last(weights), high=first(weights), weighted_mean)
end

const EXTREMES = load_extreme_weights()
const LOW_INPUT_HEIGHT_MM = BASE_INPUT_HEIGHT_MM * EXTREMES.low / EXTREMES.weighted_mean
const HIGH_INPUT_HEIGHT_MM = BASE_INPUT_HEIGHT_MM * EXTREMES.high / EXTREMES.weighted_mean

function required_common_horn_length_mm()
    maximum_ratio = max(
        HIGH_INPUT_HEIGHT_MM / THROAT_HEIGHT_MM,
        THROAT_HEIGHT_MM / LOW_INPUT_HEIGHT_MM,
    )
    ALUMINIUM_CP_M_S * log(maximum_ratio) /
    (8 * LOWER_FREQUENCY_HZ * 0.05) * 1e3
end

const HORN_LENGTH_MM = required_common_horn_length_mm()

function physical_config(input_height_mm)
    HornPointRadiatorConfig(
        input_height_mm=Float64(input_height_mm),
        throat_height_mm=THROAT_HEIGHT_MM,
        horn_length_mm=HORN_LENGTH_MM,
        lower_band_frequency_hz=LOWER_FREQUENCY_HZ,
        pressure_wave_speed_m_s=ALUMINIUM_CP_M_S,
        straight_guide_length_mm=MAXIMUM_PATH_MM,
        rounded_axial_length_mm=GUIDE_AXIAL_MM,
        receiver_length_mm=85.0,
        receiver_half_height_mm=45.0,
        target_distance_mm=60.0,
        diffuser_output_height_mm=DIFFUSER_OUTPUT_HEIGHT_MM,
        diffuser_length_mm=DIFFUSER_LENGTH_MM,
        profile_samples=360,
        arc_integration_samples=8001,
        minimum_inner_radius_mm=1.0,
    )
end

const CASES = Dict(
    "baseline_edge" => (
        variant=:collector_straight_diffuser,
        config=physical_config(BASE_INPUT_HEIGHT_MM),
    ),
    "physical_edge" => (
        variant=:collector_straight_diffuser,
        config=physical_config(LOW_INPUT_HEIGHT_MM),
    ),
    "baseline_center" => (
        variant=:collector_smooth_diffuser,
        config=physical_config(BASE_INPUT_HEIGHT_MM),
    ),
    "physical_center" => (
        variant=:collector_smooth_diffuser,
        config=physical_config(HIGH_INPUT_HEIGHT_MM),
    ),
)

if STAGE == "mesh"
    using Gmsh: gmsh
elseif STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif STAGE == "solve"
    include(joinpath(@__DIR__, "transient_model_solver.jl"))
    include(joinpath(@__DIR__, "horn_point_radiator_transient_solver.jl"))
    using .TransientModelElasticity
    using .HornPointRadiatorTransientSolver
elseif STAGE in ("harmonic", "harmonic_analyze")
    using JLD2
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
    include(joinpath(@__DIR__, "horn_point_radiator_harmonic_solver.jl"))
    using .SinusoidalMaterialLens
    using .HornPointRadiatorHarmonicSolver
elseif STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
    include(joinpath(@__DIR__, "transient_model_solver.jl"))
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    include(joinpath(@__DIR__, "spectral_analysis.jl"))
    using .ImpulseRiskAnalysis
    using .SpectralAnalysis
end

mesh_path(case_name) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(case_name).msh")
model_path(case_name) = joinpath(OUTPUT_ROOT, "models", "model_$(case_name).json")
signal_path(case_name) = joinpath(OUTPUT_ROOT, "signals", "$(case_name).jld2")
harmonic_path(case_name) = joinpath(OUTPUT_ROOT, "harmonic", "$(case_name).jld2")

function selected_case()
    isnothing(CASE_NAME) && error("--case is required for stage $STAGE")
    haskey(CASES, CASE_NAME) || error("unknown case: $CASE_NAME")
    CASE_NAME, CASES[CASE_NAME]
end

function run_mesh_stage()
    case_name, case = selected_case()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        build_horn_point_radiator_mesh(
            mesh_path(case_name);
            config=case.config,
            variant=case.variant,
            size_path_mm=0.55,
            size_receiver_mm=1.20,
            size_radiator_mm=0.55,
        )
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    case_name, _ = selected_case()
    convert_mesh(mesh_path(case_name); output_dir=dirname(model_path(case_name)))
end

function run_solve_stage()
    case_name, case = selected_case()
    probes = probe_points_mm(case.config, case.variant)
    run_horn_point_radiator_transient(
        model_path(case_name), signal_path(case_name);
        id=Symbol(case_name),
        probe_names=probes.names,
        probe_x_mm=probes.x_mm,
        probe_y_mm=probes.y_mm,
        outlet_x_mm=outlet_x_mm(case.config, case.variant),
        material=TransientMaterialConfig(
            density=ALUMINIUM_DENSITY_KG_M3,
            pressure_wave_speed=ALUMINIUM_CP_M_S,
            shear_wave_speed=ALUMINIUM_CS_M_S,
            rayleigh_alpha=0.0,
            rayleigh_beta=0.0,
        ),
        config=TransientModelConfig(
            frequency_hz=FREQUENCY_HZ,
            pulse_cycles=5.0,
            final_time_s=145.0e-6,
            samples_per_period=40,
            pressure_amplitude_pa=1.0e6,
            element_order=1,
            quadrature_degree=2,
        ),
    )
end

function run_harmonic_stage()
    case_name, case = selected_case()
    probes = probe_points_mm(case.config, case.variant)
    target_index = findfirst(==("target_axis"), probes.names)
    isnothing(target_index) && error("target_axis probe is absent")
    result = solve_horn_point_radiator_harmonic(
        mesh_path(case_name), probes.x_mm[target_index], probes.y_mm[target_index],
        aluminium_6061();
        config=HornPointRadiatorHarmonicConfig(
            frequency_hz=FREQUENCY_HZ,
            pressure_amplitude_pa=1.0e6,
            element_order=1,
            quadrature_degree=2,
        ),
    )
    mkpath(dirname(harmonic_path(case_name)))
    JLD2.jldsave(
        harmonic_path(case_name);
        format_version=1,
        case_name,
        frequency_hz=FREQUENCY_HZ,
        input_height_mm=case.config.input_height_mm,
        target_displacement_m=result.target_displacement_m,
        source_normal_displacement_integral_m2=
            result.source_normal_displacement_integral_m2,
        active_input_power_w_per_m=result.active_input_power_w_per_m,
        reactive_input_power_var_per_m=result.reactive_input_power_var_per_m,
    )
    println("[+] $case_name target |ux|=" *
            "$(abs(result.target_displacement_m[1]) * 1e9) nm")
    println("[+] $case_name active input=$(result.active_input_power_w_per_m) W/m")
    println("[+] $(harmonic_path(case_name))")
end

wrap_phase_deg(value) = rad2deg(atan(sin(angle(value)), cos(angle(value))))

function run_harmonic_analyze_stage()
    data = Dict(name => JLD2.load(harmonic_path(name)) for name in keys(CASES))
    ux = Dict(name => ComplexF64(item["target_displacement_m"][1])
              for (name, item) in data)
    power = Dict(name => abs(Float64(item["active_input_power_w_per_m"]))
                 for (name, item) in data)
    edge_ratio = abs(ux["physical_edge"] / ux["baseline_edge"])
    center_ratio = abs(ux["physical_center"] / ux["baseline_center"])
    edge_phase_deg = wrap_phase_deg(ux["physical_edge"] / ux["baseline_edge"])
    center_phase_deg = wrap_phase_deg(ux["physical_center"] / ux["baseline_center"])
    target_edge_ratio = LOW_INPUT_HEIGHT_MM / BASE_INPUT_HEIGHT_MM
    target_center_ratio = HIGH_INPUT_HEIGHT_MM / BASE_INPUT_HEIGHT_MM
    extreme_ratio_error = (center_ratio / edge_ratio) /
                          (target_center_ratio / target_edge_ratio) - 1
    edge_efficiency_ratio = edge_ratio * sqrt(
        power["baseline_edge"] / power["physical_edge"],
    )
    center_efficiency_ratio = center_ratio * sqrt(
        power["baseline_center"] / power["physical_center"],
    )
    passed = abs(extreme_ratio_error) <= 0.20 &&
             abs(edge_phase_deg) <= 10.0 && abs(center_phase_deg) <= 10.0

    summary_path = joinpath(OUTPUT_ROOT, "mouth_jacobian_harmonic_summary.csv")
    open(summary_path, "w") do io
        println(io, "common_horn_length_mm,edge_input_height_mm,center_input_height_mm,target_edge_amplitude_ratio,target_center_amplitude_ratio,measured_edge_amplitude_ratio,measured_center_amplitude_ratio,extreme_ratio_relative_error,edge_phase_shift_deg,center_phase_shift_deg,edge_equal_power_efficiency_ratio,center_equal_power_efficiency_ratio,baseline_edge_power_w_per_m,physical_edge_power_w_per_m,baseline_center_power_w_per_m,physical_center_power_w_per_m,passed")
        println(io, join((
            HORN_LENGTH_MM, LOW_INPUT_HEIGHT_MM, HIGH_INPUT_HEIGHT_MM,
            target_edge_ratio, target_center_ratio, edge_ratio, center_ratio,
            extreme_ratio_error, edge_phase_deg, center_phase_deg,
            edge_efficiency_ratio, center_efficiency_ratio,
            power["baseline_edge"], power["physical_edge"],
            power["baseline_center"], power["physical_center"], passed,
        ), ','))
    end
    verdict_path = joinpath(OUTPUT_ROOT, "mouth_jacobian_harmonic_verdict.txt")
    open(verdict_path, "w") do io
        println(io, passed ?
            "PASS: physical input-mouth area is a usable amplitude regulator." :
            "STOP: physical input-mouth area does not reproduce the required weights.")
        println(io, "extreme_ratio_error=$extreme_ratio_error")
        println(io, "edge_phase_shift_deg=$edge_phase_deg")
        println(io, "center_phase_shift_deg=$center_phase_deg")
    end
    println("[+] harmonic physical edge ratio=$edge_ratio (target=$target_edge_ratio)")
    println("[+] harmonic physical centre ratio=$center_ratio (target=$target_center_ratio)")
    println("[+] extreme-ratio error=$(100extreme_ratio_error)%")
    println("[+] phase shifts edge/centre=$edge_phase_deg / $center_phase_deg deg")
    println("[+] passed=$passed")
    println("[+] $summary_path")
end

function probe_signal(data, name, component=:x)
    index = findfirst(==(name), data["probe_names"])
    isnothing(index) && error("probe $name is absent")
    key = component == :x ? "probe_velocity_x_m_per_s" : "probe_velocity_y_m_per_s"
    vec(data[key][index, :])
end

function source_work_j_per_m(data, source_height_mm)
    time_s = data["time_s"]
    pressure_pa = data["source_drive_mpa"] .* 1e6
    velocity = data["source_normal_velocity_m_per_s"]
    integrand = pressure_pa .* velocity .* (source_height_mm * 1e-3)
    abs(sum(
        (integrand[index] + integrand[index + 1]) / 2 *
        (time_s[index + 1] - time_s[index])
        for index in 1:(length(time_s) - 1)
    ))
end

envelope_peak(signal) = maximum(analytic_envelope(signal))

function run_analyze_stage()
    data = Dict(name => JLD2.load(signal_path(name)) for name in keys(CASES))
    time_s = data["baseline_edge"]["time_s"]
    all(item["time_s"] == time_s for item in values(data)) ||
        error("all mouth-Jacobian cases must use one time grid")
    dt_s = time_s[2] - time_s[1]
    signal = Dict(name => probe_signal(item, "target_axis") for (name, item) in data)
    transverse = Dict(name => probe_signal(item, "target_axis", :y)
                      for (name, item) in data)
    peak = Dict(name => envelope_peak(value) for (name, value) in signal)
    work = Dict(name => source_work_j_per_m(
        data[name], CASES[name].config.input_height_mm,
    ) for name in keys(CASES))

    edge_ratio = peak["physical_edge"] / peak["baseline_edge"]
    center_ratio = peak["physical_center"] / peak["baseline_center"]
    target_edge_ratio = LOW_INPUT_HEIGHT_MM / BASE_INPUT_HEIGHT_MM
    target_center_ratio = HIGH_INPUT_HEIGHT_MM / BASE_INPUT_HEIGHT_MM
    extreme_ratio_error = (center_ratio / edge_ratio) /
                          (target_center_ratio / target_edge_ratio) - 1
    edge_pulse = pulse_metrics(
        signal["physical_edge"], signal["baseline_edge"], signal["baseline_edge"], dt_s,
    )
    center_pulse = pulse_metrics(
        signal["physical_center"], signal["baseline_center"],
        signal["baseline_center"], dt_s,
    )
    edge_transverse = sum(abs2, transverse["physical_edge"]) /
                      sum(abs2, signal["physical_edge"])
    center_transverse = sum(abs2, transverse["physical_center"]) /
                        sum(abs2, signal["physical_center"])
    passed = abs(extreme_ratio_error) <= 0.20 &&
             edge_pulse.pulse_correlation >= 0.90 &&
             center_pulse.pulse_correlation >= 0.90 &&
             edge_pulse.broadening_ratio <= 1.25 &&
             center_pulse.broadening_ratio <= 1.25 &&
             edge_pulse.postcursor_ratio <= 0.12 &&
             center_pulse.postcursor_ratio <= 0.12 &&
             max(edge_transverse, center_transverse) <= 0.10

    mkpath(OUTPUT_ROOT)
    summary_path = joinpath(OUTPUT_ROOT, "mouth_jacobian_summary.csv")
    open(summary_path, "w") do io
        println(io, "common_horn_length_mm,edge_input_height_mm,center_input_height_mm,target_edge_amplitude_ratio,target_center_amplitude_ratio,measured_edge_amplitude_ratio,measured_center_amplitude_ratio,extreme_ratio_relative_error,edge_equal_work_efficiency_ratio,center_equal_work_efficiency_ratio,edge_broadening_ratio,center_broadening_ratio,edge_pulse_correlation,center_pulse_correlation,edge_postcursor_ratio,center_postcursor_ratio,edge_transverse_energy_ratio,center_transverse_energy_ratio,passed")
        println(io, join((
            HORN_LENGTH_MM, LOW_INPUT_HEIGHT_MM, HIGH_INPUT_HEIGHT_MM,
            target_edge_ratio, target_center_ratio, edge_ratio, center_ratio,
            extreme_ratio_error,
            edge_ratio * sqrt(work["baseline_edge"] / work["physical_edge"]),
            center_ratio * sqrt(work["baseline_center"] / work["physical_center"]),
            edge_pulse.broadening_ratio, center_pulse.broadening_ratio,
            edge_pulse.pulse_correlation, center_pulse.pulse_correlation,
            edge_pulse.postcursor_ratio, center_pulse.postcursor_ratio,
            edge_transverse, center_transverse, passed,
        ), ','))
    end
    cases_path = joinpath(OUTPUT_ROOT, "mouth_jacobian_cases.csv")
    open(cases_path, "w") do io
        println(io, "case,input_height_mm,guide_state,target_envelope_peak_m_per_s,source_work_j_per_m")
        for name in sort(collect(keys(CASES)))
            state = occursin("edge", name) ? "straight" : "maximum_delay"
            println(io, join((name, CASES[name].config.input_height_mm, state,
                              peak[name], work[name]), ','))
        end
    end
    panel = plot(
        time_s .* 1e6, analytic_envelope(signal["baseline_edge"]);
        label="edge: 7.00 mm", linewidth=2, color=:gray45,
        xlabel="time, μs", ylabel="|analytic vₓ|, m/s",
        title="Physical input-mouth Jacobian: extreme channels",
        gridalpha=0.25,
    )
    plot!(panel, time_s .* 1e6, analytic_envelope(signal["physical_edge"]);
          label="edge: $(round(LOW_INPUT_HEIGHT_MM; digits=2)) mm", linewidth=2.3)
    plot!(panel, time_s .* 1e6, analytic_envelope(signal["baseline_center"]);
          label="centre: 7.00 mm", linewidth=2, color=:gray65)
    plot!(panel, time_s .* 1e6, analytic_envelope(signal["physical_center"]);
          label="centre: $(round(HIGH_INPUT_HEIGHT_MM; digits=2)) mm", linewidth=2.3)
    figure_path = joinpath(OUTPUT_ROOT, "mouth_jacobian_transient.png")
    savefig(panel, figure_path)
    println("[+] physical edge ratio=$edge_ratio (target=$target_edge_ratio)")
    println("[+] physical centre ratio=$center_ratio (target=$target_center_ratio)")
    println("[+] extreme-ratio error=$(100extreme_ratio_error)%, passed=$passed")
    println("[+] common input-transition length=$HORN_LENGTH_MM mm")
    println("[+] $summary_path")
    println("[+] $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "convert" ? run_convert_stage() :
    STAGE == "solve" ? run_solve_stage() :
    STAGE == "harmonic" ? run_harmonic_stage() :
    STAGE == "harmonic_analyze" ? run_harmonic_analyze_stage() :
    STAGE == "analyze" ? run_analyze_stage() :
    error("unknown stage: $STAGE")
end

end
