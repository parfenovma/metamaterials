module AluminiumHornTTDSmallAperture

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_AL_HORN_TTD_SMALL_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_small_aperture_242khz"),
)

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "compare")
const VARIANT = argument("variant", "lens")
VARIANT in ("lens", "uniform", "corrected") ||
    error("--variant must be lens, uniform, or corrected")

const FREQUENCY_HZ = 242.0e3
const ALUMINIUM_CP_M_S = 6122.102437409232
const GUIDE_DELAY_US_PER_MM = 0.20861664549831324
const CENTERS_MM = [-57.4, -28.7, 0.0, 28.7, 57.4]
const FOCAL_DISTANCE_MM = 60.0
const AXIAL_LENGTH_MM = 132.5

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
include(joinpath(@__DIR__, "monotonic_horn_lens.jl"))
include(joinpath(@__DIR__, "horn_monotonic_array_3d_mesher.jl"))
using .MonotonicHornLens
using .HornMonotonicArray3DMesher

const DISTANCE_MM = hypot.(FOCAL_DISTANCE_MM, CENTERS_MM)
const TARGET_DELAY_US = (maximum(DISTANCE_MM) .- DISTANCE_MM) ./ ALUMINIUM_CP_M_S .* 1e3
const EXTRA_PATH_MM = TARGET_DELAY_US ./ GUIDE_DELAY_US_PER_MM
const BEND_AMPLITUDE_MM = [
    smooth_amplitude_mm(AXIAL_LENGTH_MM, extra_path) for extra_path in EXTRA_PATH_MM
]
const CORRECTED_EXTRA_PATH_MM = [
    0.0,
    18.54533608692065,
    19.981330145747584,
    18.54533608692065,
    0.0,
]
const CORRECTED_BEND_AMPLITUDE_MM = [
    smooth_amplitude_mm(AXIAL_LENGTH_MM, extra_path)
    for extra_path in CORRECTED_EXTRA_PATH_MM
]
const LENS_CONFIG = HornMonotonicArray3DConfig(
    input_height_mm=7.0,
    throat_height_mm=3.5,
    channel_depth_mm=7.0,
    channel_pitch_mm=28.7,
    horn_length_mm=65.41,
    guide_axial_length_mm=AXIAL_LENGTH_MM,
    bend_amplitude_mm=BEND_AMPLITUDE_MM,
    receiver_length_mm=85.0,
    receiver_half_width_mm=64.0,
    receiver_half_height_mm=6.0,
    focal_distance_mm=FOCAL_DISTANCE_MM,
    profile_samples=160,
)
const CONFIG_AMPLITUDES_MM = VARIANT == "lens" ? BEND_AMPLITUDE_MM :
                             VARIANT == "corrected" ? CORRECTED_BEND_AMPLITUDE_MM :
                             zeros(length(BEND_AMPLITUDE_MM))
const CONFIG = HornMonotonicArray3DConfig(
    input_height_mm=LENS_CONFIG.input_height_mm,
    throat_height_mm=LENS_CONFIG.throat_height_mm,
    channel_depth_mm=LENS_CONFIG.channel_depth_mm,
    channel_pitch_mm=LENS_CONFIG.channel_pitch_mm,
    horn_length_mm=LENS_CONFIG.horn_length_mm,
    guide_axial_length_mm=LENS_CONFIG.guide_axial_length_mm,
    bend_amplitude_mm=CONFIG_AMPLITUDES_MM,
    receiver_length_mm=LENS_CONFIG.receiver_length_mm,
    receiver_half_width_mm=LENS_CONFIG.receiver_half_width_mm,
    receiver_half_height_mm=LENS_CONFIG.receiver_half_height_mm,
    focal_distance_mm=LENS_CONFIG.focal_distance_mm,
    profile_samples=LENS_CONFIG.profile_samples,
)

if STAGE == "mesh"
    using Gmsh: gmsh
elseif STAGE in ("solve", "matrix")
    using JLD2
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
    include(joinpath(@__DIR__, "horn_monotonic_array_3d_harmonic_solver.jl"))
    using .SinusoidalMaterialLens
    using .HornMonotonicArray3DHarmonicSolver
elseif STAGE in ("compare", "corrected_compare", "matrix_compare")
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    using .ImpulseRiskAnalysis: contiguous_width
end

mesh_path(variant=VARIANT) = joinpath(OUTPUT_ROOT, "mesh_$(variant)_half.msh")
result_path(variant=VARIANT) = joinpath(OUTPUT_ROOT, "$(variant)_half_harmonic.jld2")
matrix_path(variant=VARIANT) = joinpath(OUTPUT_ROOT, "$(variant)_half_response_matrix.jld2")

function half_probe_points(config)
    centers = channel_centers_mm(config)
    selected = findall(>=(0.0), centers)
    names = String[]
    x_mm = Float64[]
    y_mm = Float64[]
    z_mm = Float64[]
    local_x_mm = config.guide_axial_length_mm - 1.0
    for (source_index, channel_index) in enumerate(selected)
        center_y_mm = centers[channel_index]
        push!(names, "preout_$source_index")
        push!(x_mm, config.horn_length_mm + local_x_mm)
        push!(y_mm, center_y_mm)
        push!(z_mm, channel_center_z_mm(config, channel_index, local_x_mm))
        push!(names, "receiver_$source_index")
        push!(x_mm, outlet_x_mm(config) + 2.0)
        push!(y_mm, center_y_mm)
        push!(z_mm, 0.0)
    end
    push!(names, "focus")
    push!(x_mm, focus_x_mm(config))
    push!(y_mm, 0.0)
    push!(z_mm, 0.0)
    (; names, x_mm, y_mm, z_mm)
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        result = build_horn_monotonic_array_3d_half_mesh(
            mesh_path();
            config=CONFIG,
            size_path_mm=1.00,
            size_focus_mm=1.80,
            size_receiver_mm=2.60,
        )
        mkpath(OUTPUT_ROOT)
        open(joinpath(OUTPUT_ROOT, "mesh_summary_$(VARIANT).csv"), "w") do io
            println(io, "variant,full_channels,half_sources,nodes,elements,outlet_x_mm,focus_x_mm")
            println(io, join((
                VARIANT,
                length(CONFIG.bend_amplitude_mm),
                result.source_count,
                result.node_count,
                result.element_count,
                outlet_x_mm(CONFIG),
                focus_x_mm(CONFIG),
            ), ','))
        end
    finally
        gmsh.finalize()
    end
end

function run_solve_stage()
    probes = half_probe_points(CONFIG)
    material = aluminium_6061()
    source_count = (length(CONFIG.bend_amplitude_mm) + 1) ÷ 2
    result = solve_horn_monotonic_array_3d_harmonic(
        mesh_path(),
        source_count,
        outlet_x_mm(CONFIG),
        probes,
        material;
        config=HornMonotonicArray3DHarmonicConfig(
            frequency_hz=FREQUENCY_HZ,
            focal_distance_mm=FOCAL_DISTANCE_MM,
            scan_x_min_mm=35.0,
            scan_x_max_mm=80.0,
            scan_y_half_width_mm=62.0,
            scan_y_min_mm=0.0,
            scan_step_mm=2.0,
            symmetry_y=true,
        ),
        compute_source_power=true,
    )
    jldsave(
        result_path();
        format_version=1,
        variant=VARIANT,
        frequency_hz=FREQUENCY_HZ,
        full_channel_count=length(CONFIG.bend_amplitude_mm),
        half_source_count=source_count,
        centers_mm=channel_centers_mm(CONFIG),
        target_delay_us=TARGET_DELAY_US,
        guide_delay_us_per_mm=GUIDE_DELAY_US_PER_MM,
        extra_path_mm=VARIANT == "lens" ? EXTRA_PATH_MM :
                      VARIANT == "corrected" ? CORRECTED_EXTRA_PATH_MM :
                      zeros(length(EXTRA_PATH_MM)),
        bend_amplitude_mm=CONFIG.bend_amplitude_mm,
        outlet_x_mm=outlet_x_mm(CONFIG),
        focal_distance_mm=FOCAL_DISTANCE_MM,
        probe_names=result.probe_names,
        probe_x_mm=result.probe_x_mm,
        probe_y_mm=result.probe_y_mm,
        probe_z_mm=result.probe_z_mm,
        probe_displacement_m=result.probe_displacement_m,
        focus_displacement_m=result.focus_displacement_m,
        scan_x_mm=result.scan_x_mm,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=result.scan_ux_m,
        scan_uy_m=result.scan_uy_m,
        scan_uz_m=result.scan_uz_m,
        active_input_power_w=result.active_input_power_w,
        reactive_input_power_var=result.reactive_input_power_var,
    )
    println("[+] $VARIANT focus |ux|=$(abs(result.focus_displacement_m[1]) * 1e9) nm")
    println("[+] $VARIANT active input power=$(result.active_input_power_w) W")
    println("[+] $(result_path())")
end

function run_matrix_stage()
    probes = half_probe_points(CONFIG)
    material = aluminium_6061()
    source_count = (length(CONFIG.bend_amplitude_mm) + 1) ÷ 2
    result = solve_horn_monotonic_array_3d_response_matrix(
        mesh_path(),
        source_count,
        outlet_x_mm(CONFIG),
        probes,
        material;
        config=HornMonotonicArray3DHarmonicConfig(
            frequency_hz=FREQUENCY_HZ,
            focal_distance_mm=FOCAL_DISTANCE_MM,
            scan_y_half_width_mm=0.0,
            scan_y_min_mm=0.0,
            scan_step_mm=2.0,
            symmetry_y=true,
        ),
    )
    jldsave(
        matrix_path();
        format_version=1,
        variant=VARIANT,
        frequency_hz=FREQUENCY_HZ,
        source_names=["centre", "inner_pair", "outer_pair"],
        focus_response_m=result.focus_response_m,
        probe_names=result.probe_names,
        probe_response_m=result.probe_response_m,
        source_normal_displacement_integral_m3=result.source_normal_displacement_integral_m3,
    )
    reconstructed_nm = abs(sum(result.focus_response_m[1, :])) * 1e9
    println("[+] $VARIANT matrix reconstructs focus |ux|=$reconstructed_nm nm")
    println("[+] $(matrix_path())")
end

function ideal_scalar_gain()
    distance_m = DISTANCE_MM .* 1e-3
    phase = exp.(-im * 2pi * FREQUENCY_HZ .* distance_m ./ ALUMINIUM_CP_M_S)
    spreading = inv.(sqrt.(distance_m))
    abs(sum(spreading) / sum(spreading .* phase))
end

function run_compare_stage(selected_variant="lens")
    lens = JLD2.load(result_path(selected_variant))
    uniform = JLD2.load(result_path("uniform"))
    lens_focus_nm = abs(lens["focus_displacement_m"][1]) * 1e9
    uniform_focus_nm = abs(uniform["focus_displacement_m"][1]) * 1e9
    pressure_gain = lens_focus_nm / uniform_focus_nm
    lens_power = abs(Float64(lens["active_input_power_w"]))
    uniform_power = abs(Float64(uniform["active_input_power_w"]))
    power_normalized_gain = pressure_gain * sqrt(uniform_power / lens_power)
    ideal_gain = ideal_scalar_gain()
    ideal_efficiency = power_normalized_gain / ideal_gain

    x_mm = lens["scan_x_mm"]
    y_positive_mm = lens["scan_y_mm"]
    lens_x_nm = abs.(lens["scan_ux_m"]) .* 1e9
    uniform_x_nm = abs.(uniform["scan_ux_m"]) .* 1e9
    y_mm = vcat(-reverse(y_positive_mm[2:end]), y_positive_mm)
    lens_map_nm = vcat(reverse(lens_x_nm[2:end, :]; dims=1), lens_x_nm)
    uniform_map_nm = vcat(reverse(uniform_x_nm[2:end, :]; dims=1), uniform_x_nm)
    target_x_index = argmin(abs.(x_mm .- FOCAL_DISTANCE_MM))
    lens_profile = lens_map_nm[:, target_x_index]
    uniform_profile = uniform_map_nm[:, target_x_index]
    lens_width = contiguous_width(y_mm, lens_profile)
    lens_map_for_peak = map(value -> isfinite(value) ? value : -Inf, lens_map_nm)
    local_y_index, local_x_index = Tuple(argmax(lens_map_for_peak))
    local_peak_x_mm = x_mm[local_x_index]
    local_peak_y_mm = y_mm[local_y_index]
    focus_vector = lens["focus_displacement_m"]
    transverse_ratio = (abs2(focus_vector[2]) + abs2(focus_vector[3])) / abs2(focus_vector[1])
    matrix_ceiling = if selected_variant == "corrected"
        matrix = JLD2.load(matrix_path("lens"))
        response = vec(matrix["focus_response_m"][1, :])
        sum(abs, response) / abs(sum(vec(
            JLD2.load(matrix_path("uniform"))["focus_response_m"][1, :],
        )))
    else
        NaN
    end
    ceiling_efficiency = selected_variant == "corrected" ?
                         pressure_gain / matrix_ceiling : NaN
    passed = if selected_variant == "corrected"
        ceiling_efficiency >= 0.90 && power_normalized_gain >= 1.0 &&
        abs(local_peak_x_mm - FOCAL_DISTANCE_MM) <= 10.0 &&
        abs(local_peak_y_mm) <= 4.0 && transverse_ratio <= 0.10
    else
        power_normalized_gain >= 1.35 && ideal_efficiency >= 0.75 &&
        abs(local_peak_x_mm - FOCAL_DISTANCE_MM) <= 10.0 &&
        abs(local_peak_y_mm) <= 4.0 && transverse_ratio <= 0.10
    end

    output_prefix = selected_variant == "corrected" ? "corrected_small_aperture" :
                    "small_aperture"
    summary_path = joinpath(OUTPUT_ROOT, "$(output_prefix)_summary.csv")
    open(summary_path, "w") do io
        println(io, "variant,lens_focus_ux_nm,uniform_focus_ux_nm,pressure_gain,lens_active_power_w,uniform_active_power_w,power_normalized_gain,ideal_scalar_gain,ideal_efficiency,matrix_ceiling,ceiling_efficiency,lens_fwhm_mm,local_peak_x_mm,local_peak_y_mm,transverse_energy_ratio,passed")
        println(io, join((
            selected_variant, lens_focus_nm, uniform_focus_nm, pressure_gain, lens_power,
            uniform_power, power_normalized_gain, ideal_gain, ideal_efficiency,
            matrix_ceiling, ceiling_efficiency, lens_width.width, local_peak_x_mm,
            local_peak_y_mm, transverse_ratio, passed,
        ), ','))
    end

    finite_lens = filter(isfinite, vec(lens_map_nm))
    finite_uniform = filter(isfinite, vec(uniform_map_nm))
    common_limit = max(maximum(finite_lens), maximum(finite_uniform))
    lens_panel = heatmap(
        x_mm, y_mm, lens_map_nm;
        xlabel="distance from radiator, mm", ylabel="y, mm",
        title="$(titlecase(selected_variant)) five-channel Al lens", color=:viridis,
        colorbar_title="|uₓ|, nm", clims=(0, common_limit), aspect_ratio=:equal,
    )
    uniform_panel = heatmap(
        x_mm, y_mm, uniform_map_nm;
        xlabel="distance from radiator, mm", ylabel="y, mm",
        title="Matched five-channel straight control", color=:viridis,
        colorbar_title="|uₓ|, nm", clims=(0, common_limit), aspect_ratio=:equal,
    )
    for panel in (lens_panel, uniform_panel)
        scatter!(panel, [FOCAL_DISTANCE_MM], [0.0]; marker=:xcross,
                 color=:white, markerstrokewidth=2, label=false)
    end
    profile_panel = plot(
        y_mm, lens_profile;
        label="$selected_variant lens", linewidth=2.5,
        xlabel="y at x=60 mm, mm", ylabel="|uₓ|, nm",
        title="Physical focus profile", gridalpha=0.25,
    )
    plot!(profile_panel, y_mm, uniform_profile;
          label="straight control", linewidth=2)
    figure_path = joinpath(OUTPUT_ROOT, "$(output_prefix)_focus.png")
    savefig(plot(lens_panel, uniform_panel, profile_panel;
                 layout=(2, 2), size=(1450, 950), margin=5Plots.mm), figure_path)

    geometry_path = joinpath(OUTPUT_ROOT, "small_aperture_geometry.csv")
    open(geometry_path, "w") do io
        println(io, "channel,y_mm,target_delay_us,extra_path_mm,bend_amplitude_mm,inner_radius_mm")
        for index in eachindex(CENTERS_MM)
            radius = minimum_inner_radius_mm(
                AXIAL_LENGTH_MM, BEND_AMPLITUDE_MM[index], LENS_CONFIG.throat_height_mm,
            )
            println(io, join((index, CENTERS_MM[index], TARGET_DELAY_US[index],
                              EXTRA_PATH_MM[index], BEND_AMPLITUDE_MM[index], radius), ','))
        end
    end
    println("[+] $selected_variant five-channel pressure gain=$pressure_gain")
    println("[+] power-normalized gain=$power_normalized_gain, " *
            "ideal=$ideal_gain, efficiency=$ideal_efficiency")
    println("[+] local peak=($local_peak_x_mm, $local_peak_y_mm) mm, " *
            "FWHM=$(lens_width.width) mm, transverse=$transverse_ratio, passed=$passed")
    println("[+] $summary_path")
    println("[+] $figure_path")
end

wrap_phase_deg(value) = rad2deg(atan(sin(angle(value)), cos(angle(value))))

function run_matrix_compare_stage()
    lens = JLD2.load(matrix_path("lens"))
    uniform = JLD2.load(matrix_path("uniform"))
    lens_focus = vec(lens["focus_response_m"][1, :])
    uniform_focus = vec(uniform["focus_response_m"][1, :])
    source_names = String.(lens["source_names"])
    current_gain = abs(sum(lens_focus)) / abs(sum(uniform_focus))
    coherence = abs(sum(lens_focus)) / sum(abs, lens_focus)
    aligned_ceiling = sum(abs, lens_focus) / abs(sum(uniform_focus))
    outer_reference = lens_focus[end]
    phase_error_deg = [wrap_phase_deg(value / outer_reference) for value in lens_focus]
    phase_slope_deg_per_mm = 360 * FREQUENCY_HZ * GUIDE_DELAY_US_PER_MM * 1e-6
    current_half_path_mm = [EXTRA_PATH_MM[3], EXTRA_PATH_MM[4], EXTRA_PATH_MM[5]]
    carrier_cycle_mm = 360 / phase_slope_deg_per_mm
    correction_mm = zeros(3)
    corrected_path_mm = copy(current_half_path_mm)
    for index in 1:2
        candidates = [phase_error_deg[index] / phase_slope_deg_per_mm +
                      cycle * carrier_cycle_mm for cycle in -2:2]
        feasible = filter(delta -> current_half_path_mm[index] + delta >= 0.0, candidates)
        isempty(feasible) && error("no non-negative path correction for $(source_names[index])")
        correction_mm[index] = feasible[argmin(abs.(feasible))]
        corrected_path_mm[index] += correction_mm[index]
    end
    corrected_amplitude_mm = [
        smooth_amplitude_mm(AXIAL_LENGTH_MM, path) for path in corrected_path_mm
    ]

    summary_path = joinpath(OUTPUT_ROOT, "small_aperture_response_matrix.csv")
    open(summary_path, "w") do io
        println(io, "source,lens_focus_amplitude_nm,lens_focus_phase_deg,uniform_focus_amplitude_nm,uniform_focus_phase_deg,phase_error_to_outer_deg,current_extra_path_mm,calculated_correction_mm,corrected_extra_path_mm,corrected_amplitude_mm")
        for index in eachindex(source_names)
            println(io, join((
                source_names[index], abs(lens_focus[index]) * 1e9,
                rad2deg(angle(lens_focus[index])), abs(uniform_focus[index]) * 1e9,
                rad2deg(angle(uniform_focus[index])), phase_error_deg[index],
                current_half_path_mm[index], correction_mm[index], corrected_path_mm[index],
                corrected_amplitude_mm[index],
            ), ','))
        end
    end
    verdict_path = joinpath(OUTPUT_ROOT, "small_aperture_matrix_verdict.txt")
    open(verdict_path, "w") do io
        println(io, "current_gain=$current_gain")
        println(io, "coherence_efficiency=$coherence")
        println(io, "phase_aligned_ceiling=$aligned_ceiling")
        println(io, "phase_slope_deg_per_mm=$phase_slope_deg_per_mm")
        println(io, "carrier_cycle_mm=$carrier_cycle_mm")
        println(io, "next=one calculated passive path correction, then stop or advance")
    end
    println("[+] matrix current gain=$current_gain, coherence=$coherence, " *
            "phase-aligned ceiling=$aligned_ceiling")
    println("[+] phase errors to outer=$(phase_error_deg) deg")
    println("[+] path corrections=$(correction_mm) mm")
    println("[+] corrected half paths=$(corrected_path_mm) mm")
    println("[+] $summary_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "solve" ? run_solve_stage() :
    STAGE == "matrix" ? run_matrix_stage() :
    STAGE == "compare" ? run_compare_stage() :
    STAGE == "corrected_compare" ? run_compare_stage("corrected") :
    STAGE == "matrix_compare" ? run_matrix_compare_stage() :
    error("unknown stage: $STAGE")
end

end
