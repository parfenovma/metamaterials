module AluminiumHornTTDOutputJacobian

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_AL_HORN_TTD_OUTPUT_JACOBIAN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_output_jacobian_242khz"),
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
const ALUMINIUM_CP_M_S = 6122.102437409232
const THROAT_HEIGHT_MM = 3.5
const BASE_OUTPUT_HEIGHT_MM = 7.0
const HORN_LENGTH_MM = 65.41
const GUIDE_AXIAL_MM = 132.5
const MAXIMUM_PATH_MM = 155.53

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
using .HornPointRadiatorMesher

const WEIGHT_PATH = joinpath(
    PROJECT_ROOT, "tmp", "aluminium_horn_ttd_full_diffuser_aperture_242khz",
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
const LOW_OUTPUT_HEIGHT_MM = BASE_OUTPUT_HEIGHT_MM * EXTREMES.low / EXTREMES.weighted_mean
const HIGH_OUTPUT_HEIGHT_MM = BASE_OUTPUT_HEIGHT_MM * EXTREMES.high / EXTREMES.weighted_mean

function required_common_diffuser_length_mm()
    maximum_ratio = max(
        HIGH_OUTPUT_HEIGHT_MM / THROAT_HEIGHT_MM,
        THROAT_HEIGHT_MM / LOW_OUTPUT_HEIGHT_MM,
    )
    ALUMINIUM_CP_M_S * log(maximum_ratio) /
    (8 * LOWER_FREQUENCY_HZ * 0.05) * 1e3
end

const DIFFUSER_LENGTH_MM = required_common_diffuser_length_mm()

function physical_config(output_height_mm)
    HornPointRadiatorConfig(
        input_height_mm=7.0,
        throat_height_mm=THROAT_HEIGHT_MM,
        horn_length_mm=HORN_LENGTH_MM,
        lower_band_frequency_hz=LOWER_FREQUENCY_HZ,
        pressure_wave_speed_m_s=ALUMINIUM_CP_M_S,
        straight_guide_length_mm=MAXIMUM_PATH_MM,
        rounded_axial_length_mm=GUIDE_AXIAL_MM,
        receiver_length_mm=85.0,
        receiver_half_height_mm=45.0,
        target_distance_mm=60.0,
        diffuser_output_height_mm=Float64(output_height_mm),
        diffuser_length_mm=DIFFUSER_LENGTH_MM,
        profile_samples=360,
        arc_integration_samples=8001,
        minimum_inner_radius_mm=1.0,
    )
end

const CASES = Dict(
    "baseline_edge" => (
        variant=:collector_straight_diffuser,
        config=physical_config(BASE_OUTPUT_HEIGHT_MM),
    ),
    "physical_edge" => (
        variant=:collector_straight_diffuser,
        config=physical_config(LOW_OUTPUT_HEIGHT_MM),
    ),
    "baseline_center" => (
        variant=:collector_smooth_diffuser,
        config=physical_config(BASE_OUTPUT_HEIGHT_MM),
    ),
    "physical_center" => (
        variant=:collector_smooth_diffuser,
        config=physical_config(HIGH_OUTPUT_HEIGHT_MM),
    ),
)

if STAGE == "mesh"
    using Gmsh: gmsh
elseif STAGE in ("harmonic", "analyze")
    using JLD2
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
    include(joinpath(@__DIR__, "horn_point_radiator_harmonic_solver.jl"))
    using .SinusoidalMaterialLens
    using .HornPointRadiatorHarmonicSolver
end

mesh_path(case_name) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(case_name).msh")
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
        output_height_mm=case.config.diffuser_output_height_mm,
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

function run_analyze_stage()
    data = Dict(name => JLD2.load(harmonic_path(name)) for name in keys(CASES))
    ux = Dict(name => ComplexF64(item["target_displacement_m"][1])
              for (name, item) in data)
    power = Dict(name => abs(Float64(item["active_input_power_w_per_m"]))
                 for (name, item) in data)
    edge_ratio = abs(ux["physical_edge"] / ux["baseline_edge"])
    center_ratio = abs(ux["physical_center"] / ux["baseline_center"])
    edge_phase_deg = wrap_phase_deg(ux["physical_edge"] / ux["baseline_edge"])
    center_phase_deg = wrap_phase_deg(ux["physical_center"] / ux["baseline_center"])
    target_edge_ratio = LOW_OUTPUT_HEIGHT_MM / BASE_OUTPUT_HEIGHT_MM
    target_center_ratio = HIGH_OUTPUT_HEIGHT_MM / BASE_OUTPUT_HEIGHT_MM
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

    summary_path = joinpath(OUTPUT_ROOT, "output_jacobian_harmonic_summary.csv")
    open(summary_path, "w") do io
        println(io, "common_diffuser_length_mm,edge_output_height_mm,center_output_height_mm,target_edge_amplitude_ratio,target_center_amplitude_ratio,measured_edge_amplitude_ratio,measured_center_amplitude_ratio,extreme_ratio_relative_error,edge_phase_shift_deg,center_phase_shift_deg,edge_equal_power_efficiency_ratio,center_equal_power_efficiency_ratio,baseline_edge_power_w_per_m,physical_edge_power_w_per_m,baseline_center_power_w_per_m,physical_center_power_w_per_m,passed")
        println(io, join((
            DIFFUSER_LENGTH_MM, LOW_OUTPUT_HEIGHT_MM, HIGH_OUTPUT_HEIGHT_MM,
            target_edge_ratio, target_center_ratio, edge_ratio, center_ratio,
            extreme_ratio_error, edge_phase_deg, center_phase_deg,
            edge_efficiency_ratio, center_efficiency_ratio,
            power["baseline_edge"], power["physical_edge"],
            power["baseline_center"], power["physical_center"], passed,
        ), ','))
    end
    verdict_path = joinpath(OUTPUT_ROOT, "output_jacobian_harmonic_verdict.txt")
    open(verdict_path, "w") do io
        println(io, passed ?
            "PASS: output transition area is a usable passive amplitude regulator." :
            "STOP: output transition area does not reproduce the required weights.")
        println(io, "extreme_ratio_error=$extreme_ratio_error")
        println(io, "edge_phase_shift_deg=$edge_phase_deg")
        println(io, "center_phase_shift_deg=$center_phase_deg")
    end
    println("[+] output physical edge ratio=$edge_ratio (target=$target_edge_ratio)")
    println("[+] output physical centre ratio=$center_ratio (target=$target_center_ratio)")
    println("[+] extreme-ratio error=$(100extreme_ratio_error)%")
    println("[+] phase shifts edge/centre=$edge_phase_deg / $center_phase_deg deg")
    println("[+] passed=$passed")
    println("[+] $summary_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "harmonic" ? run_harmonic_stage() :
    STAGE == "analyze" ? run_analyze_stage() :
    error("unknown stage: $STAGE")
end

end
