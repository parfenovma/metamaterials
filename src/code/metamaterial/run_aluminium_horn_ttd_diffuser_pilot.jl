module AluminiumHornTTDDiffuserPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_AL_HORN_TTD_DIFFUSER_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_diffuser_pilot_242khz"),
)
const ABRUPT_ROOT = joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_pilot_242khz")

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "analyze")
const CASE_NAME = argument("case", nothing)
const FREQUENCY_HZ = 242.0e3
const ALUMINIUM_DENSITY_KG_M3 = 2700.0
const ALUMINIUM_CP_M_S = 6122.102437409232
const ALUMINIUM_CS_M_S = 3083.810277185563
const HORN_LENGTH_MM = 65.41
const GUIDE_AXIAL_MM = 132.5
const MAXIMUM_PATH_MM = 155.53
const DIFFUSER_LENGTH_MM = 65.41
const SOURCE_HEIGHT_M = 7.0e-3

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
using .HornPointRadiatorMesher

function physical_config()
    HornPointRadiatorConfig(
        input_height_mm=7.0,
        throat_height_mm=3.5,
        horn_length_mm=HORN_LENGTH_MM,
        lower_band_frequency_hz=162.2e3,
        pressure_wave_speed_m_s=ALUMINIUM_CP_M_S,
        straight_guide_length_mm=MAXIMUM_PATH_MM,
        rounded_axial_length_mm=GUIDE_AXIAL_MM,
        receiver_length_mm=85.0,
        receiver_half_height_mm=45.0,
        target_distance_mm=60.0,
        diffuser_output_height_mm=7.0,
        diffuser_length_mm=DIFFUSER_LENGTH_MM,
        profile_samples=300,
        arc_integration_samples=8001,
        minimum_inner_radius_mm=1.0,
    )
end

const CASES = Dict(
    "diffuser_equal_path" => (
        variant=:collector_straight_diffuser,
        config=physical_config(),
    ),
    "diffuser_maximum" => (
        variant=:collector_smooth_diffuser,
        config=physical_config(),
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
            final_time_s=130.0e-6,
            samples_per_period=40,
            pressure_amplitude_pa=1.0e6,
            element_order=1,
            quadrature_degree=2,
        ),
    )
end

function probe_signal(data, name, component=:x)
    index = findfirst(==(name), data["probe_names"])
    isnothing(index) && error("probe $name is absent")
    key = component == :x ? "probe_velocity_x_m_per_s" : "probe_velocity_y_m_per_s"
    vec(data[key][index, :])
end

function source_work_j_per_m(data)
    time_s = data["time_s"]
    pressure_pa = data["source_drive_mpa"] .* 1e6
    velocity = data["source_normal_velocity_m_per_s"]
    integrand = pressure_pa .* velocity .* SOURCE_HEIGHT_M
    work = sum(
        (integrand[index] + integrand[index + 1]) / 2 *
        (time_s[index + 1] - time_s[index])
        for index in 1:(length(time_s) - 1)
    )
    abs(work)
end

function run_analyze_stage()
    diffuser_equal = JLD2.load(signal_path("diffuser_equal_path"))
    diffuser_maximum = JLD2.load(signal_path("diffuser_maximum"))
    abrupt_equal = JLD2.load(joinpath(
        ABRUPT_ROOT, "signals", "equal_path_time_refined.jld2",
    ))
    abrupt_maximum = JLD2.load(joinpath(
        ABRUPT_ROOT, "signals", "maximum_time_refined.jld2",
    ))
    time_s = diffuser_maximum["time_s"]
    time_s == diffuser_equal["time_s"] == abrupt_equal["time_s"] ==
              abrupt_maximum["time_s"] || error("all pilot time grids must match")
    dt_s = time_s[2] - time_s[1]
    diffuser_equal_x = probe_signal(diffuser_equal, "target_axis")
    diffuser_maximum_x = probe_signal(diffuser_maximum, "target_axis")
    diffuser_maximum_y = probe_signal(diffuser_maximum, "target_axis", :y)
    abrupt_equal_x = probe_signal(abrupt_equal, "target_axis")
    abrupt_maximum_x = probe_signal(abrupt_maximum, "target_axis")
    pulse = pulse_metrics(
        diffuser_maximum_x, diffuser_equal_x, diffuser_equal_x, dt_s,
    )
    work_diffuser_equal = source_work_j_per_m(diffuser_equal)
    work_diffuser_maximum = source_work_j_per_m(diffuser_maximum)
    work_abrupt_equal = source_work_j_per_m(abrupt_equal)
    work_abrupt_maximum = source_work_j_per_m(abrupt_maximum)
    envelope_peak(signal) = maximum(analytic_envelope(signal))
    diffuser_equal_over_abrupt = envelope_peak(diffuser_equal_x) /
        envelope_peak(abrupt_equal_x) * sqrt(work_abrupt_equal / work_diffuser_equal)
    diffuser_maximum_over_abrupt = envelope_peak(diffuser_maximum_x) /
        envelope_peak(abrupt_maximum_x) * sqrt(work_abrupt_maximum / work_diffuser_maximum)
    energy_over_equal = sum(abs2, diffuser_maximum_x) /
                        sum(abs2, diffuser_equal_x)
    transverse_ratio = sum(abs2, diffuser_maximum_y) /
                       sum(abs2, diffuser_maximum_x)
    passed = pulse.gain_peak >= 0.8 && pulse.broadening_ratio <= 1.25 &&
             pulse.pulse_correlation >= 0.90 && pulse.postcursor_ratio <= 0.10 &&
             transverse_ratio <= 0.10 && diffuser_maximum_over_abrupt >= 1.10

    mkpath(OUTPUT_ROOT)
    summary_path = joinpath(OUTPUT_ROOT, "diffuser_pilot_summary.csv")
    open(summary_path, "w") do io
        println(io, "diffuser_length_mm,diffuser_output_height_mm,diffuser_epsilon_ad,peak_over_diffuser_equal,energy_over_diffuser_equal,broadening_ratio,pulse_correlation,postcursor_ratio,transverse_energy_ratio,diffuser_equal_over_abrupt_equal_power,diffuser_maximum_over_abrupt_equal_power,work_diffuser_equal_j_per_m,work_diffuser_maximum_j_per_m,work_abrupt_equal_j_per_m,work_abrupt_maximum_j_per_m,passed")
        println(io, join((
            DIFFUSER_LENGTH_MM, 7.0,
            diffuser_adiabatic_parameter(CASES["diffuser_maximum"].config),
            pulse.gain_peak, energy_over_equal, pulse.broadening_ratio,
            pulse.pulse_correlation, pulse.postcursor_ratio, transverse_ratio,
            diffuser_equal_over_abrupt, diffuser_maximum_over_abrupt,
            work_diffuser_equal, work_diffuser_maximum,
            work_abrupt_equal, work_abrupt_maximum, passed,
        ), ','))
    end
    verdict_path = joinpath(OUTPUT_ROOT, "verdict.txt")
    open(verdict_path, "w") do io
        println(io, passed ?
            "PASS: output diffuser improves equal-energy target amplitude." :
            "STOP: output diffuser does not provide the required +10% throughput.")
        println(io, "maximum_over_abrupt_equal_power=$diffuser_maximum_over_abrupt")
        println(io, "pulse_peak_over_equal=$(pulse.gain_peak)")
        println(io, "postcursor=$(pulse.postcursor_ratio)")
    end

    panel = plot(
        time_s .* 1e6, analytic_envelope(abrupt_maximum_x);
        label="abrupt maximum", linewidth=2, color=:gray40,
        xlabel="time, μs", ylabel="|analytic vₓ|, m/s",
        title="All-Al output diffuser: equal-source target envelopes",
        gridalpha=0.25,
    )
    plot!(panel, time_s .* 1e6, analytic_envelope(diffuser_maximum_x);
          label="diffuser maximum", linewidth=2.5, color=:royalblue)
    plot!(panel, time_s .* 1e6, analytic_envelope(diffuser_equal_x);
          label="diffuser equal path", linewidth=2, color=:seagreen)
    figure_path = joinpath(OUTPUT_ROOT, "diffuser_pilot_transient.png")
    savefig(panel, figure_path)
    println("[+] diffuser maximum/equal peak=$(pulse.gain_peak), " *
            "Bt=$(pulse.broadening_ratio), rho=$(pulse.pulse_correlation)")
    println("[+] diffuser maximum/abrupt equal-power=" *
            "$diffuser_maximum_over_abrupt, passed=$passed")
    println("[+] $summary_path")
    println("[+] $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "convert" ? run_convert_stage() :
    STAGE == "solve" ? run_solve_stage() :
    STAGE == "analyze" ? run_analyze_stage() :
    error("unknown stage: $STAGE")
end

end
