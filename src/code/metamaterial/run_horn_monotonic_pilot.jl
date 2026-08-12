module HornMonotonicPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_aperture")
const DESIGN_PATH = get(
    ENV,
    "METAMATERIALS_HORN_MONOTONIC_DESIGN",
    joinpath(DESIGN_ROOT, "monotonic_aperture_impulse.jld2"),
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_MONOTONIC_PILOT_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_pilot"),
)

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "analyze")
const CASE_NAME = argument("case", nothing)

using JLD2
const DESIGN = JLD2.load(DESIGN_PATH)
const AXIAL_LENGTH_MM = Float64(DESIGN["common_axial_length_mm"])
const EXTRA_PATH_MM = haskey(DESIGN, "maximum_extra_path_mm") ?
    Float64(DESIGN["maximum_extra_path_mm"]) :
    maximum(Float64.(DESIGN["extra_path_mm"]))
const PATH_LENGTH_MM = AXIAL_LENGTH_MM + EXTRA_PATH_MM

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
using .HornPointRadiatorMesher

const CASES = Dict(
    "monotonic_device" => (
        variant=:collector_device,
        config=HornPointRadiatorConfig(
            straight_guide_length_mm=AXIAL_LENGTH_MM,
            rounded_axial_length_mm=AXIAL_LENGTH_MM - 0.1,
        ),
    ),
    "monotonic_equal_path" => (
        variant=:collector_straight,
        config=HornPointRadiatorConfig(
            straight_guide_length_mm=PATH_LENGTH_MM,
            rounded_axial_length_mm=AXIAL_LENGTH_MM,
        ),
    ),
    "monotonic_max" => (
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
            final_time_s=100.0e-6,
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
        reference["time_s"], reference["source_drive_mpa"],
        probe_signal(reference, probe_name); config,
    )
    value_at_frequency(relative_transfer(transfer, reference_transfer), 242.0e3)
end

function run_analyze_stage()
    data = Dict(name => JLD2.load(signal_path(name)) for name in keys(CASES))
    device = data["monotonic_device"]
    equal_path = data["monotonic_equal_path"]
    maximum_state = data["monotonic_max"]
    time_s = maximum_state["time_s"]
    dt_s = time_s[2] - time_s[1]
    target_equal = probe_signal(equal_path, "target_axis")
    target_device = probe_signal(device, "target_axis")
    target_maximum = probe_signal(maximum_state, "target_axis")
    target_maximum_y = probe_signal(maximum_state, "target_axis", :y)
    pulse = pulse_metrics(target_maximum, target_equal, target_equal, dt_s)
    maximum_energy = sum(abs2, target_maximum) * dt_s
    equal_energy = sum(abs2, target_equal) * dt_s
    transverse_energy_ratio = sum(abs2, target_maximum_y) / sum(abs2, target_maximum)
    near_coefficient = transfer_at_carrier(maximum_state, device, "near_radiator")
    target_coefficient = transfer_at_carrier(maximum_state, device, "target_axis")
    device_peak = maximum(analytic_envelope(target_device))
    maximum_peak = maximum(analytic_envelope(target_maximum))
    passed = pulse.gain_peak >= 0.8 && pulse.broadening_ratio <= 1.25 &&
             pulse.pulse_correlation >= 0.90 && pulse.postcursor_ratio <= 0.10 &&
             transverse_energy_ratio <= 0.10

    summary_path = joinpath(OUTPUT_ROOT, "horn_monotonic_maximum_summary.csv")
    mkpath(OUTPUT_ROOT)
    open(summary_path, "w") do io
        println(io, "design_path,axial_length_mm,path_length_mm,smooth_amplitude_mm,peak_over_equal_path,peak_over_device,energy_over_equal_path,broadening_ratio,pulse_correlation,postcursor_ratio,transverse_energy_ratio,near_amplitude_over_device,near_phase_deg_over_device,near_group_delay_us_over_device,target_amplitude_over_device,target_phase_deg_over_device,passed")
        println(io, join((
            DESIGN_PATH,
            AXIAL_LENGTH_MM,
            PATH_LENGTH_MM,
            smooth_guide_amplitude_mm(CASES["monotonic_max"].config),
            pulse.gain_peak,
            maximum_peak / device_peak,
            maximum_energy / equal_energy,
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
        time_s .* 1e6,
        target_device;
        label="straight device",
        linewidth=2,
        xlabel="time, μs",
        ylabel="target vx, m/s",
        title="Monotonic-lens maximum delay: five-cycle target pulse",
        xlims=(30, 95),
        gridalpha=0.25,
    )
    plot!(panel, equal_path["time_s"] .* 1e6, target_equal;
          label="straight equal path", linewidth=2)
    plot!(panel, time_s .* 1e6, target_maximum;
          label="maximum sin⁴ delay", linewidth=2.5)
    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_maximum_transient.png")
    savefig(panel, figure_path)
    println("[+] monotonic maximum peak/equal=$(pulse.gain_peak), peak/device=$(maximum_peak / device_peak)")
    println("[+] energy/equal=$(maximum_energy / equal_energy), Bt=$(pulse.broadening_ratio), rho=$(pulse.pulse_correlation)")
    println("[+] postcursor=$(pulse.postcursor_ratio), transverse energy=$transverse_energy_ratio, gate passed=$passed")
    println("[+] near maximum/device=$(near_coefficient.amplitude) ∠ $(rad2deg(near_coefficient.phase_rad)) deg, delay=$(near_coefficient.group_delay_s * 1e6) μs")
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
