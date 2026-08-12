module HornLongSmoothPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_LONG_SMOOTH_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_long_smooth_pilot"),
)

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "analyze")
const CASE_NAME = argument("case", nothing)
const AXIAL_LENGTH_MM = 38.158058
const EXTRA_PATH_MM = 4.42
const PATH_LENGTH_MM = AXIAL_LENGTH_MM + EXTRA_PATH_MM

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
using .HornPointRadiatorMesher

const CASES = Dict(
    "long_device" => (
        variant=:collector_device,
        config=HornPointRadiatorConfig(straight_guide_length_mm=AXIAL_LENGTH_MM),
    ),
    "long_equal_path" => (
        variant=:collector_straight,
        config=HornPointRadiatorConfig(
            straight_guide_length_mm=PATH_LENGTH_MM,
            rounded_axial_length_mm=AXIAL_LENGTH_MM,
        ),
    ),
    "long_smooth" => (
        variant=:collector_smooth,
        config=HornPointRadiatorConfig(
            straight_guide_length_mm=PATH_LENGTH_MM,
            rounded_axial_length_mm=AXIAL_LENGTH_MM,
        ),
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
        model_path(case_name),
        signal_path(case_name);
        id=Symbol(case_name),
        probe_names=probes.names,
        probe_x_mm=probes.x_mm,
        probe_y_mm=probes.y_mm,
        outlet_x_mm=outlet_x_mm(case.config, case.variant),
        material=TransientMaterialConfig(),
        config=TransientModelConfig(
            frequency_hz=242.0e3,
            pulse_cycles=5.0,
            final_time_s=85.0e-6,
            samples_per_period=30,
            element_order=1,
            quadrature_degree=2,
        ),
    )
end

function probe_index(data, name)
    index = findfirst(==(name), data["probe_names"])
    isnothing(index) && error("probe $name is absent")
    index
end

function probe_signal(data, name, component=:x)
    key = component == :x ? "probe_velocity_x_m_per_s" : "probe_velocity_y_m_per_s"
    vec(data[key][probe_index(data, name), :])
end

function transfer_at_carrier(data, reference, probe_name)
    config = SpectrumConfig(window=:rectangular, zero_padding_factor=16, input_floor_relative=1e-3)
    transfer = analyze_transfer(
        data["time_s"], data["source_drive_mpa"], probe_signal(data, probe_name); config,
    )
    reference_transfer = analyze_transfer(
        reference["time_s"], reference["source_drive_mpa"], probe_signal(reference, probe_name); config,
    )
    value_at_frequency(relative_transfer(transfer, reference_transfer), 242.0e3)
end

function run_analyze_stage()
    data = Dict(name => JLD2.load(signal_path(name)) for name in keys(CASES))
    device = data["long_device"]
    equal_path = data["long_equal_path"]
    smooth = data["long_smooth"]
    time_s = smooth["time_s"]
    dt_s = time_s[2] - time_s[1]
    target_equal = probe_signal(equal_path, "target_axis")
    target_device = probe_signal(device, "target_axis")
    target_smooth = probe_signal(smooth, "target_axis")
    target_smooth_y = probe_signal(smooth, "target_axis", :y)
    pulse = pulse_metrics(target_smooth, target_equal, target_equal, dt_s)
    smooth_energy = sum(abs2, target_smooth) * dt_s
    equal_energy = sum(abs2, target_equal) * dt_s
    transverse_energy_ratio = sum(abs2, target_smooth_y) / sum(abs2, target_smooth)
    near_coefficient = transfer_at_carrier(smooth, device, "near_radiator")
    target_coefficient = transfer_at_carrier(smooth, device, "target_axis")
    device_peak = maximum(analytic_envelope(target_device))
    smooth_peak = maximum(analytic_envelope(target_smooth))
    passed = pulse.gain_peak >= 0.8 && pulse.broadening_ratio <= 1.25 &&
             pulse.pulse_correlation >= 0.90 && pulse.postcursor_ratio <= 0.10 &&
             transverse_energy_ratio <= 0.10

    summary_path = joinpath(OUTPUT_ROOT, "horn_long_smooth_summary.csv")
    mkpath(OUTPUT_ROOT)
    open(summary_path, "w") do io
        println(io, "axial_length_mm,path_length_mm,smooth_amplitude_mm,peak_over_equal_path,peak_over_device,energy_over_equal_path,broadening_ratio,pulse_correlation,postcursor_ratio,transverse_energy_ratio,near_amplitude_over_device,near_phase_deg_over_device,near_group_delay_us_over_device,target_amplitude_over_device,target_phase_deg_over_device,passed")
        println(io, join((
            AXIAL_LENGTH_MM,
            PATH_LENGTH_MM,
            smooth_guide_amplitude_mm(CASES["long_smooth"].config),
            pulse.gain_peak,
            smooth_peak / device_peak,
            smooth_energy / equal_energy,
            pulse.broadening_ratio,
            pulse.pulse_correlation,
            pulse.postcursor_ratio,
            transverse_energy_ratio,
            near_coefficient.amplitude,
            rad2deg(near_coefficient.phase_rad),
            near_coefficient.group_delay_s * 1e6,
            target_coefficient.amplitude,
            rad2deg(target_coefficient.phase_rad),
            passed,
        ), ','))
    end

    panel = plot(
        time_s .* 1e6, target_device;
        label="device, L=$(round(AXIAL_LENGTH_MM; digits=2)) mm",
        linewidth=2, xlabel="time, μs", ylabel="target vx, m/s",
        title="Curvature-limited sin⁴ delay: five-cycle target pulse",
        xlims=(25, 82), gridalpha=0.25,
    )
    plot!(panel, equal_path["time_s"] .* 1e6, target_equal;
          label="straight equal path", linewidth=2)
    plot!(panel, time_s .* 1e6, target_smooth;
          label="long sin⁴ delay", linewidth=2.5)
    figure_path = joinpath(OUTPUT_ROOT, "horn_long_smooth_transient.png")
    savefig(panel, figure_path)
    println("[+] long smooth peak/equal=$(pulse.gain_peak), peak/device=$(smooth_peak / device_peak)")
    println("[+] energy/equal=$(smooth_energy / equal_energy), Bt=$(pulse.broadening_ratio), rho=$(pulse.pulse_correlation)")
    println("[+] transverse energy=$transverse_energy_ratio, gate passed=$passed")
    println("[+] near smooth/device=$(near_coefficient.amplitude) ∠ $(rad2deg(near_coefficient.phase_rad)) deg")
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
