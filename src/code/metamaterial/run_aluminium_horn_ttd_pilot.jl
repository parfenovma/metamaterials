module AluminiumHornTTDPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_AL_HORN_TTD_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_pilot_242khz"),
)

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "analyze")
const CASE_NAME = argument("case", nothing)
const SIGNAL_TAG = argument(
    "signal-tag", get(ENV, "METAMATERIALS_AL_HORN_TTD_SIGNAL_TAG", ""),
)
const FINAL_TIME_US = parse(Float64, argument(
    "final-time-us", get(ENV, "METAMATERIALS_AL_HORN_TTD_FINAL_TIME_US", "100.0"),
))
const SAMPLES_PER_PERIOD = parse(Int, argument(
    "samples-per-period",
    get(ENV, "METAMATERIALS_AL_HORN_TTD_SAMPLES_PER_PERIOD", "30"),
))

const CENTER_FREQUENCY_HZ = 242.0e3
const LOWER_BAND_FREQUENCY_HZ = 162.2e3
const ALUMINIUM_DENSITY_KG_M3 = 2700.0
const ALUMINIUM_CP_M_S = 6122.102437409232
const ALUMINIUM_CS_M_S = 3083.810277185563
const HORN_LENGTH_MM = 65.41
const AXIAL_LENGTH_MM = 132.5
const MIDDLE_EXTRA_PATH_MM = 14.65
const MAXIMUM_EXTRA_PATH_MM = 23.03
const MAXIMUM_PATH_LENGTH_MM = AXIAL_LENGTH_MM + MAXIMUM_EXTRA_PATH_MM

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
using .HornPointRadiatorMesher

function physical_config(path_length_mm; axial_length_mm=AXIAL_LENGTH_MM)
    HornPointRadiatorConfig(
        input_height_mm=7.0,
        throat_height_mm=3.5,
        horn_length_mm=HORN_LENGTH_MM,
        lower_band_frequency_hz=LOWER_BAND_FREQUENCY_HZ,
        pressure_wave_speed_m_s=ALUMINIUM_CP_M_S,
        straight_guide_length_mm=Float64(path_length_mm),
        rounded_axial_length_mm=Float64(axial_length_mm),
        receiver_length_mm=85.0,
        receiver_half_height_mm=45.0,
        target_distance_mm=60.0,
        profile_samples=300,
        arc_integration_samples=8001,
        minimum_inner_radius_mm=1.0,
    )
end

const CASES = Dict(
    "zero" => (
        variant=:collector_device,
        config=physical_config(AXIAL_LENGTH_MM; axial_length_mm=AXIAL_LENGTH_MM - 0.1),
        extra_path_mm=0.0,
    ),
    "middle" => (
        variant=:collector_smooth,
        config=physical_config(AXIAL_LENGTH_MM + MIDDLE_EXTRA_PATH_MM),
        extra_path_mm=MIDDLE_EXTRA_PATH_MM,
    ),
    "maximum" => (
        variant=:collector_smooth,
        config=physical_config(MAXIMUM_PATH_LENGTH_MM),
        extra_path_mm=MAXIMUM_EXTRA_PATH_MM,
    ),
    "equal_path" => (
        variant=:collector_straight,
        config=physical_config(MAXIMUM_PATH_LENGTH_MM),
        extra_path_mm=MAXIMUM_EXTRA_PATH_MM,
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
signal_path(case_name) = joinpath(OUTPUT_ROOT, "signals", "$(case_name)$(SIGNAL_TAG).jld2")

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
    aluminium = TransientMaterialConfig(
        density=ALUMINIUM_DENSITY_KG_M3,
        pressure_wave_speed=ALUMINIUM_CP_M_S,
        shear_wave_speed=ALUMINIUM_CS_M_S,
        rayleigh_alpha=0.0,
        rayleigh_beta=0.0,
    )
    run_horn_point_radiator_transient(
        model_path(case_name),
        signal_path(case_name);
        id=Symbol("aluminium_horn_ttd_$(case_name)"),
        probe_names=probes.names,
        probe_x_mm=probes.x_mm,
        probe_y_mm=probes.y_mm,
        outlet_x_mm=outlet_x_mm(case.config, case.variant),
        material=aluminium,
        config=TransientModelConfig(
            frequency_hz=CENTER_FREQUENCY_HZ,
            pulse_cycles=5.0,
            final_time_s=FINAL_TIME_US * 1e-6,
            samples_per_period=SAMPLES_PER_PERIOD,
            pressure_amplitude_pa=1.0e6,
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

function carrier_relative(sample, reference, probe_name)
    spectrum_config = SpectrumConfig(
        window=:rectangular,
        zero_padding_factor=16,
        input_floor_relative=1e-3,
    )
    sample_transfer = analyze_transfer(
        sample["time_s"], sample["source_drive_mpa"], probe_signal(sample, probe_name);
        config=spectrum_config,
    )
    reference_transfer = analyze_transfer(
        reference["time_s"], reference["source_drive_mpa"],
        probe_signal(reference, probe_name); config=spectrum_config,
    )
    value_at_frequency(relative_transfer(sample_transfer, reference_transfer), CENTER_FREQUENCY_HZ)
end

function state_row(name, data, zero, equal_path)
    signal = probe_signal(data, "target_axis")
    transverse = probe_signal(data, "target_axis", :y)
    zero_signal = probe_signal(zero, "target_axis")
    equal_signal = probe_signal(equal_path, "target_axis")
    dt_s = data["time_s"][2] - data["time_s"][1]
    relative_to_zero = pulse_metrics(signal, zero_signal, zero_signal, dt_s)
    relative_to_equal = pulse_metrics(signal, equal_signal, equal_signal, dt_s)
    carrier = carrier_relative(data, zero, "target_axis")
    case = CASES[name]
    amplitude = case.variant == :collector_smooth ?
                smooth_guide_amplitude_mm(case.config) : 0.0
    inner_radius = case.variant == :collector_smooth ?
                   minimum_smooth_inner_radius_mm(case.config) : Inf
    (
        state=name,
        path_length_mm=case.config.straight_guide_length_mm,
        extra_path_mm=case.extra_path_mm,
        target_delay_us=case.extra_path_mm / ALUMINIUM_CP_M_S * 1e3,
        smooth_amplitude_mm=amplitude,
        minimum_inner_radius_mm=inner_radius,
        outlet_x_mm=outlet_x_mm(case.config, case.variant),
        peak_velocity_m_s=maximum(analytic_envelope(signal)),
        peak_over_zero=relative_to_zero.gain_peak,
        peak_over_equal_path=relative_to_equal.gain_peak,
        energy_over_zero=sum(abs2, signal) / sum(abs2, zero_signal),
        energy_over_equal_path=sum(abs2, signal) / sum(abs2, equal_signal),
        broadening_over_zero=relative_to_zero.broadening_ratio,
        broadening_over_equal_path=relative_to_equal.broadening_ratio,
        correlation_over_zero=relative_to_zero.pulse_correlation,
        correlation_over_equal_path=relative_to_equal.pulse_correlation,
        postcursor_ratio=relative_to_zero.postcursor_ratio,
        transverse_energy_ratio=sum(abs2, transverse) / sum(abs2, signal),
        carrier_amplitude_over_zero=carrier.amplitude,
        carrier_phase_deg_over_zero=rad2deg(carrier.phase_rad),
        carrier_group_delay_us_over_zero=carrier.group_delay_s * 1e6,
        envelope_delay_us_over_zero=envelope_delay(
            data["time_s"], zero_signal, signal,
        ) * 1e6,
        envelope_delay_us_over_equal_path=envelope_delay(
            data["time_s"], equal_signal, signal,
        ) * 1e6,
    )
end

function write_rows(path, rows)
    keys_order = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(keys_order), ','))
        for row in rows
            println(io, join((getproperty(row, key) for key in keys_order), ','))
        end
    end
end

function geometry_plot()
    maximum_case = CASES["maximum"]
    maximum_config = maximum_case.config
    maximum_outlet = outlet_x_mm(maximum_config, maximum_case.variant)
    local_x = collect(range(0.0, maximum_config.rounded_axial_length_mm; length=700))
    amplitude = smooth_guide_amplitude_mm(maximum_config)
    centre_y = HornPointRadiatorMesher.smooth_centerline_y_mm.(
        Ref(maximum_config), local_x, amplitude,
    )
    slope = HornPointRadiatorMesher.smooth_centerline_slope.(
        Ref(maximum_config), local_x, amplitude,
    )
    normalization = hypot.(1.0, slope)
    normal_x = .-slope ./ normalization
    normal_y = 1.0 ./ normalization
    half_width = maximum_config.throat_height_mm / 2
    guide_x = HORN_LENGTH_MM .+ local_x
    lower_x = guide_x .- half_width .* normal_x
    lower_y = centre_y .- half_width .* normal_y
    upper_x = guide_x .+ half_width .* normal_x
    upper_y = centre_y .+ half_width .* normal_y

    zero_config = CASES["zero"].config
    horn_x = collect(range(0.0, HORN_LENGTH_MM; length=300))
    horn_half = horn_height_mm.(Ref(zero_config), horn_x) ./ 2
    body = plot(
        xlabel="x, mm", ylabel="y, mm",
        title="Maximum-delay Al channel: physical walls",
        aspect_ratio=:equal, gridalpha=0.25, legend=:topright,
        xlims=(-5, maximum_outlet + maximum_config.receiver_length_mm + 5),
        ylims=(-maximum_config.receiver_half_height_mm - 3,
               maximum_config.receiver_half_height_mm + 3),
    )
    horn_polygon_x = vcat(horn_x, reverse(horn_x))
    horn_polygon_y = vcat(horn_half, -reverse(horn_half))
    plot!(body, horn_polygon_x, horn_polygon_y; seriestype=:shape,
          color=:black, fillalpha=0.18, linewidth=2, label="collector")
    guide_polygon_x = vcat(lower_x, reverse(upper_x))
    guide_polygon_y = vcat(lower_y, reverse(upper_y))
    plot!(body, guide_polygon_x, guide_polygon_y; seriestype=:shape,
          color=:royalblue, fillalpha=0.22, linewidth=2, label="sin⁴ guide")
    receiver_x = [maximum_outlet,
                  maximum_outlet + maximum_config.receiver_length_mm,
                  maximum_outlet + maximum_config.receiver_length_mm,
                  maximum_outlet, maximum_outlet]
    receiver_y = [-maximum_config.receiver_half_height_mm,
                  -maximum_config.receiver_half_height_mm,
                  maximum_config.receiver_half_height_mm,
                  maximum_config.receiver_half_height_mm,
                  -maximum_config.receiver_half_height_mm]
    plot!(body, receiver_x, receiver_y; color=:gray35, linewidth=1.5,
          linestyle=:dash, label="Al receiver")
    scatter!(body, [maximum_outlet + maximum_config.target_distance_mm], [0.0];
             marker=:star5, markersize=7, color=:red, label="60 mm target")

    centreline_plot = plot(
        xlabel="x, mm", ylabel="y, mm",
        title="Frozen delay states (physical units)",
        aspect_ratio=:equal, gridalpha=0.25, legend=:topright,
        xlims=(HORN_LENGTH_MM - 5, HORN_LENGTH_MM + AXIAL_LENGTH_MM + 5),
        ylims=(-5, 40),
    )
    colors = Dict("zero" => :gray35, "middle" => :darkorange, "maximum" => :royalblue)
    for name in ("zero", "middle", "maximum")
        case = CASES[name]
        guide_axial = name == "zero" ? AXIAL_LENGTH_MM : case.config.rounded_axial_length_mm
        local_x = collect(range(0.0, guide_axial; length=500))
        centre_y = guide_center_y_mm.(Ref(case.config), Ref(case.variant), local_x)
        plot!(centreline_plot, HORN_LENGTH_MM .+ local_x, centre_y;
              color=colors[name], linewidth=2.5,
              label="$name centreline, ΔL=$(case.extra_path_mm) mm")
    end
    x_max = HORN_LENGTH_MM + AXIAL_LENGTH_MM
    vline!(centreline_plot, [x_max]; color=:black, linestyle=:dash, label="radiator plane")
    plot(body, centreline_plot; layout=(2, 1), size=(900, 850))
end

function run_analyze_stage()
    data = Dict(name => JLD2.load(signal_path(name)) for name in keys(CASES))
    zero = data["zero"]
    equal_path = data["equal_path"]
    rows = [state_row(name, data[name], zero, equal_path)
            for name in ("zero", "middle", "maximum", "equal_path")]
    mkpath(OUTPUT_ROOT)
    summary_path = joinpath(OUTPUT_ROOT, "aluminium_horn_ttd_summary.csv")
    write_rows(summary_path, rows)

    maximum_row = rows[3]
    passed = maximum_row.peak_over_equal_path >= 0.8 &&
             maximum_row.broadening_over_equal_path <= 1.25 &&
             maximum_row.correlation_over_equal_path >= 0.90 &&
             maximum_row.postcursor_ratio <= 0.10 &&
             maximum_row.transverse_energy_ratio <= 0.10
    verdict_path = joinpath(OUTPUT_ROOT, "verdict.txt")
    open(verdict_path, "w") do io
        println(io, passed ? "PASS: advance to a small coupled Al aperture." :
                           "STOP: do not build the Al aperture; revise the isolated channel.")
        println(io, "peak/equal_path=$(maximum_row.peak_over_equal_path)")
        println(io, "Bt/equal_path=$(maximum_row.broadening_over_equal_path)")
        println(io, "rho/equal_path=$(maximum_row.correlation_over_equal_path)")
        println(io, "postcursor=$(maximum_row.postcursor_ratio)")
        println(io, "transverse_energy=$(maximum_row.transverse_energy_ratio)")
    end

    geometry_path = joinpath(OUTPUT_ROOT, "aluminium_horn_ttd_geometry.png")
    savefig(geometry_plot(), geometry_path)

    time_us = zero["time_s"] .* 1e6
    waveform_plot = plot(
        xlabel="time, μs", ylabel="target vₓ, m/s",
        title="All-Al horn–TTD pilot: 5-cycle pulse at 60 mm",
        gridalpha=0.25,
    )
    for (name, color) in (("zero", :gray35), ("middle", :darkorange),
                          ("maximum", :royalblue), ("equal_path", :seagreen))
        plot!(waveform_plot, time_us, probe_signal(data[name], "target_axis");
              label=name, color=color, linewidth=name == "maximum" ? 2.5 : 1.8)
    end
    waveform_path = joinpath(OUTPUT_ROOT, "aluminium_horn_ttd_transient.png")
    savefig(waveform_plot, waveform_path)

    println("[+] maximum peak/equal=$(maximum_row.peak_over_equal_path), " *
            "peak/zero=$(maximum_row.peak_over_zero)")
    println("[+] Bt=$(maximum_row.broadening_over_equal_path), " *
            "rho=$(maximum_row.correlation_over_equal_path), " *
            "post=$(maximum_row.postcursor_ratio), " *
            "transverse=$(maximum_row.transverse_energy_ratio)")
    println("[+] measured delay=$(maximum_row.envelope_delay_us_over_zero) μs, " *
            "target=$(maximum_row.target_delay_us) μs, gate passed=$passed")
    println("[+] $summary_path")
    println("[+] $geometry_path")
    println("[+] $waveform_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "convert" ? run_convert_stage() :
    STAGE == "solve" ? run_solve_stage() :
    STAGE == "analyze" ? run_analyze_stage() :
    error("unknown stage: $STAGE")
end

end
