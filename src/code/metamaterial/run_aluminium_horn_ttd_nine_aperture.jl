module AluminiumHornTTDNineAperture

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "compare")
const VARIANT = argument("variant", "lens")
VARIANT in ("lens", "uniform") || error("--variant must be lens or uniform")
const CHANNEL_COUNT = parse(Int, argument("channels", "9"))
CHANNEL_COUNT in (9, 15) || error("--channels must be 9 or 15")
const DIFFUSER = lowercase(argument("diffuser", "false")) in ("1", "true", "yes")
const WEIGHT_FILE = argument("weights-file", "15_power_optimal_weights.jld2")
const WEIGHT_LABEL = argument("weight-label", "carrier")
const OUTPUT_ROOT = get(
    ENV,
    DIFFUSER ? "METAMATERIALS_AL_HORN_TTD_FULL_DIFFUSER_OUTPUT" :
    CHANNEL_COUNT == 9 ? "METAMATERIALS_AL_HORN_TTD_NINE_OUTPUT" :
                         "METAMATERIALS_AL_HORN_TTD_FULL_OUTPUT",
    joinpath(
        PROJECT_ROOT,
        "tmp",
        DIFFUSER ? "aluminium_horn_ttd_full_diffuser_aperture_242khz" :
        CHANNEL_COUNT == 9 ? "aluminium_horn_ttd_nine_aperture_242khz" :
                             "aluminium_horn_ttd_full_aperture_242khz",
    ),
)

const FREQUENCY_HZ = parse(Float64, argument("frequency-hz", "242000"))
const CARRIER_FREQUENCY_HZ = 242.0e3
const ALUMINIUM_CP_M_S = 6122.102437409232
const GUIDE_DELAY_US_PER_MM = 0.20861664549831324
const FOCAL_DISTANCE_MM = 60.0
const AXIAL_LENGTH_MM = 132.5
const CHANNEL_PITCH_MM = 114.8 / (CHANNEL_COUNT - 1)
const CENTERS_MM = collect(-((CHANNEL_COUNT - 1) ÷ 2):((CHANNEL_COUNT - 1) ÷ 2)) .*
                   CHANNEL_PITCH_MM

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
include(joinpath(@__DIR__, "monotonic_horn_lens.jl"))
include(joinpath(@__DIR__, "horn_monotonic_array_3d_mesher.jl"))
using .MonotonicHornLens
using .HornMonotonicArray3DMesher

const DISTANCE_MM = hypot.(FOCAL_DISTANCE_MM, CENTERS_MM)
const TARGET_DELAY_US = (maximum(DISTANCE_MM) .- DISTANCE_MM) ./ ALUMINIUM_CP_M_S .* 1e3
const GROUP_DELAY_PATH_MM = TARGET_DELAY_US ./ GUIDE_DELAY_US_PER_MM

function transferred_phase_correction_mm(abs_y_mm)
    y_mm = abs(Float64(abs_y_mm))
    y_mm <= 28.7 && return 1.9456364425991912 +
        (5.607486503910028 - 1.9456364425991912) * y_mm / 28.7
    5.607486503910028 * (57.4 - y_mm) / 28.7
end

const PHASE_CORRECTION_MM = [
    isapprox(abs(y_mm), 57.4; atol=1e-12) ? 0.0 :
    transferred_phase_correction_mm(y_mm)
    for y_mm in CENTERS_MM
]
const EXTRA_PATH_MM = GROUP_DELAY_PATH_MM .+ PHASE_CORRECTION_MM
const BEND_AMPLITUDE_MM = [
    smooth_amplitude_mm(AXIAL_LENGTH_MM, path_mm) for path_mm in EXTRA_PATH_MM
]
const CONFIG = HornMonotonicArray3DConfig(
    input_height_mm=7.0,
    throat_height_mm=3.5,
    channel_depth_mm=7.0,
    channel_pitch_mm=CHANNEL_PITCH_MM,
    horn_length_mm=65.41,
    guide_axial_length_mm=AXIAL_LENGTH_MM,
    bend_amplitude_mm=VARIANT == "lens" ? BEND_AMPLITUDE_MM : zeros(CHANNEL_COUNT),
    diffuser_output_height_mm=DIFFUSER ? 7.0 : 3.5,
    diffuser_length_mm=DIFFUSER ? 65.41 : 0.0,
    receiver_length_mm=85.0,
    receiver_half_width_mm=64.0,
    receiver_half_height_mm=6.0,
    focal_distance_mm=FOCAL_DISTANCE_MM,
    profile_samples=160,
)

if STAGE == "mesh"
    using Gmsh: gmsh
elseif STAGE in ("solve", "matrix", "weighted")
    using JLD2
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
    include(joinpath(@__DIR__, "horn_monotonic_array_3d_harmonic_solver.jl"))
    using .SinusoidalMaterialLens
    using .HornMonotonicArray3DHarmonicSolver
elseif STAGE in ("compare", "matrix_compare", "weighted_compare")
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    using .ImpulseRiskAnalysis: contiguous_width
end

mesh_path(variant=VARIANT) = joinpath(OUTPUT_ROOT, "mesh_$(variant)_half.msh")
frequency_suffix(frequency_hz=FREQUENCY_HZ) =
    isapprox(frequency_hz, CARRIER_FREQUENCY_HZ; atol=1e-6) ? "" :
    "_f$(round(Int, frequency_hz))hz"
result_path(variant=VARIANT; frequency_hz=FREQUENCY_HZ) = joinpath(
    OUTPUT_ROOT, "$(variant)_half_harmonic$(frequency_suffix(frequency_hz)).jld2",
)
matrix_path(; frequency_hz=FREQUENCY_HZ) = joinpath(
    OUTPUT_ROOT, "lens_half_response_matrix$(frequency_suffix(frequency_hz)).jld2",
)
weighted_result_path(; frequency_hz=FREQUENCY_HZ) = joinpath(
    OUTPUT_ROOT,
    "lens_half_weighted" * (WEIGHT_LABEL == "carrier" ? "" : "_$(WEIGHT_LABEL)") *
    "_harmonic$(frequency_suffix(frequency_hz)).jld2",
)
weight_design_path() = joinpath(OUTPUT_ROOT, WEIGHT_FILE)

function half_probe_points(config)
    centers = channel_centers_mm(config)
    selected = findall(>=(0.0), centers)
    names = String["focus"]
    x_mm = Float64[focus_x_mm(config)]
    y_mm = Float64[0.0]
    z_mm = Float64[0.0]
    for (source_index, channel_index) in enumerate(selected)
        push!(names, "preout_$source_index")
        push!(x_mm, outlet_x_mm(config) - 1.0)
        push!(y_mm, centers[channel_index])
        push!(z_mm, channel_center_z_mm(
            config, channel_index, config.guide_axial_length_mm - 1.0,
        ))
    end
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
            println(io, join((VARIANT, CHANNEL_COUNT, result.source_count,
                              result.node_count, result.element_count,
                              outlet_x_mm(CONFIG), focus_x_mm(CONFIG)), ','))
        end
    finally
        gmsh.finalize()
    end
end

function run_solve_stage()
    source_count = (CHANNEL_COUNT + 1) ÷ 2
    result = solve_horn_monotonic_array_3d_harmonic(
        mesh_path(), source_count, outlet_x_mm(CONFIG), half_probe_points(CONFIG),
        aluminium_6061();
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
        channel_count=CHANNEL_COUNT,
        centers_mm=CENTERS_MM,
        target_delay_us=TARGET_DELAY_US,
        group_delay_path_mm=GROUP_DELAY_PATH_MM,
        phase_correction_mm=VARIANT == "lens" ? PHASE_CORRECTION_MM : zeros(CHANNEL_COUNT),
        extra_path_mm=VARIANT == "lens" ? EXTRA_PATH_MM : zeros(CHANNEL_COUNT),
        bend_amplitude_mm=CONFIG.bend_amplitude_mm,
        focus_displacement_m=result.focus_displacement_m,
        scan_x_mm=result.scan_x_mm,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=result.scan_ux_m,
        scan_uy_m=result.scan_uy_m,
        scan_uz_m=result.scan_uz_m,
        active_input_power_w=result.active_input_power_w,
        reactive_input_power_var=result.reactive_input_power_var,
    )
    println("[+] $VARIANT $CHANNEL_COUNT-channel focus |ux|=$(abs(result.focus_displacement_m[1]) * 1e9) nm")
    println("[+] $VARIANT active input power=$(result.active_input_power_w) W")
    println("[+] $(result_path())")
end

function run_matrix_stage()
    VARIANT == "lens" || error("matrix stage is only needed for --variant=lens")
    source_count = (CHANNEL_COUNT + 1) ÷ 2
    result = solve_horn_monotonic_array_3d_response_matrix(
        mesh_path("lens"), source_count, outlet_x_mm(CONFIG), half_probe_points(CONFIG),
        aluminium_6061();
        config=HornMonotonicArray3DHarmonicConfig(
            frequency_hz=FREQUENCY_HZ,
            focal_distance_mm=FOCAL_DISTANCE_MM,
            scan_y_half_width_mm=0.0,
            scan_y_min_mm=0.0,
            scan_step_mm=2.0,
            symmetry_y=true,
        ),
    )
    half_centers_mm = CENTERS_MM[CENTERS_MM .>= 0.0]
    jldsave(
        matrix_path();
        format_version=1,
        channel_count=CHANNEL_COUNT,
        half_centers_mm,
        focus_response_m=result.focus_response_m,
        source_normal_displacement_integral_m3=result.source_normal_displacement_integral_m3,
    )
    println("[+] $CHANNEL_COUNT-channel matrix reconstructs focus |ux|=" *
            "$(abs(sum(result.focus_response_m[1, :])) * 1e9) nm")
    println("[+] $(matrix_path())")
end

function run_weighted_stage()
    CHANNEL_COUNT == 15 || error("power-weight verification is defined for 15 channels")
    DIFFUSER || error("power-weight verification is defined for the output-diffuser lens")
    design = JLD2.load(weight_design_path())
    weights = Float64.(design["pressure_weights"])
    source_count = (CHANNEL_COUNT + 1) ÷ 2
    length(weights) == source_count || error("weight file does not match the half-domain")
    result = solve_horn_monotonic_array_3d_harmonic(
        mesh_path("lens"), source_count, outlet_x_mm(CONFIG), half_probe_points(CONFIG),
        aluminium_6061();
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
        source_weights=weights,
        compute_source_power=true,
    )
    JLD2.jldsave(
        weighted_result_path();
        format_version=1,
        frequency_hz=FREQUENCY_HZ,
        channel_count=CHANNEL_COUNT,
        centers_mm=CENTERS_MM,
        half_centers_mm=CENTERS_MM[CENTERS_MM .>= 0.0],
        pressure_weights=weights,
        focus_displacement_m=result.focus_displacement_m,
        scan_x_mm=result.scan_x_mm,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=result.scan_ux_m,
        scan_uy_m=result.scan_uy_m,
        scan_uz_m=result.scan_uz_m,
        active_input_power_w=result.active_input_power_w,
        reactive_input_power_var=result.reactive_input_power_var,
        source_normal_displacement_integral_m3=
            result.source_normal_displacement_integral_m3,
    )
    println("[+] weighted lens focus |ux|=$(abs(result.focus_displacement_m[1]) * 1e9) nm")
    println("[+] weighted active input power=$(result.active_input_power_w) W")
    println("[+] $(weighted_result_path())")
end

function ideal_scalar_gain()
    distance_m = DISTANCE_MM .* 1e-3
    phase = exp.(-im * 2pi * FREQUENCY_HZ .* distance_m ./ ALUMINIUM_CP_M_S)
    spreading = inv.(sqrt.(distance_m))
    abs(sum(spreading) / sum(spreading .* phase))
end

function run_compare_stage()
    lens = JLD2.load(result_path("lens"))
    uniform = JLD2.load(result_path("uniform"))
    abrupt_uniform = DIFFUSER ? JLD2.load(joinpath(
        PROJECT_ROOT,
        "tmp",
        "aluminium_horn_ttd_full_aperture_242khz",
        "uniform_half_harmonic.jld2",
    )) : uniform
    lens_focus_nm = abs(lens["focus_displacement_m"][1]) * 1e9
    uniform_focus_nm = abs(uniform["focus_displacement_m"][1]) * 1e9
    pressure_gain = lens_focus_nm / uniform_focus_nm
    lens_power = abs(Float64(lens["active_input_power_w"]))
    uniform_power = abs(Float64(uniform["active_input_power_w"]))
    power_gain = pressure_gain * sqrt(uniform_power / lens_power)
    abrupt_uniform_focus_nm = abs(abrupt_uniform["focus_displacement_m"][1]) * 1e9
    abrupt_uniform_power = abs(Float64(abrupt_uniform["active_input_power_w"]))
    pressure_gain_over_abrupt = lens_focus_nm / abrupt_uniform_focus_nm
    power_gain_over_abrupt = pressure_gain_over_abrupt *
                              sqrt(abrupt_uniform_power / lens_power)
    ideal_gain = ideal_scalar_gain()
    ideal_efficiency = power_gain / ideal_gain

    x_mm = lens["scan_x_mm"]
    y_positive_mm = lens["scan_y_mm"]
    y_mm = vcat(-reverse(y_positive_mm[2:end]), y_positive_mm)
    lens_half = abs.(lens["scan_ux_m"]) .* 1e9
    uniform_half = abs.(uniform["scan_ux_m"]) .* 1e9
    lens_map_nm = vcat(reverse(lens_half[2:end, :]; dims=1), lens_half)
    uniform_map_nm = vcat(reverse(uniform_half[2:end, :]; dims=1), uniform_half)
    target_x_index = argmin(abs.(x_mm .- FOCAL_DISTANCE_MM))
    lens_profile = lens_map_nm[:, target_x_index]
    uniform_profile = uniform_map_nm[:, target_x_index]
    width = contiguous_width(y_mm, lens_profile)
    focus_vector = lens["focus_displacement_m"]
    transverse_ratio = (abs2(focus_vector[2]) + abs2(focus_vector[3])) / abs2(focus_vector[1])
    required_gain = CHANNEL_COUNT == 15 ? 2.0 : 1.5
    gate_gain = DIFFUSER ? power_gain_over_abrupt : power_gain
    passed = gate_gain >= required_gain && power_gain >= 1.5 &&
             ideal_efficiency >= 0.50 &&
             abs(width.peak_coordinate) <= 4.0 && transverse_ratio <= 0.10

    summary_path = joinpath(OUTPUT_ROOT, "$(CHANNEL_COUNT)_aperture_summary.csv")
    open(summary_path, "w") do io
        println(io, "diffuser,lens_focus_ux_nm,matched_uniform_focus_ux_nm,matched_pressure_gain,lens_active_power_w,matched_uniform_active_power_w,matched_power_normalized_gain,abrupt_uniform_focus_ux_nm,abrupt_uniform_active_power_w,pressure_gain_over_abrupt,power_normalized_gain_over_abrupt,ideal_scalar_gain,ideal_efficiency,fwhm_mm,profile_peak_y_mm,transverse_energy_ratio,passed")
        println(io, join((DIFFUSER, lens_focus_nm, uniform_focus_nm, pressure_gain, lens_power,
                          uniform_power, power_gain, abrupt_uniform_focus_nm,
                          abrupt_uniform_power, pressure_gain_over_abrupt,
                          power_gain_over_abrupt, ideal_gain, ideal_efficiency,
                          width.width, width.peak_coordinate, transverse_ratio, passed), ','))
    end
    geometry_path = joinpath(OUTPUT_ROOT, "$(CHANNEL_COUNT)_aperture_geometry.csv")
    open(geometry_path, "w") do io
        println(io, "channel,y_mm,target_delay_us,group_delay_path_mm,phase_correction_mm,extra_path_mm,bend_amplitude_mm")
        for index in eachindex(CENTERS_MM)
            println(io, join((index, CENTERS_MM[index], TARGET_DELAY_US[index],
                              GROUP_DELAY_PATH_MM[index], PHASE_CORRECTION_MM[index],
                              EXTRA_PATH_MM[index], BEND_AMPLITUDE_MM[index]), ','))
        end
    end

    finite_lens = filter(isfinite, vec(lens_map_nm))
    finite_uniform = filter(isfinite, vec(uniform_map_nm))
    common_limit = max(maximum(finite_lens), maximum(finite_uniform))
    lens_panel = heatmap(
        x_mm, y_mm, lens_map_nm; xlabel="distance from radiator, mm", ylabel="y, mm",
        title="Response-informed $CHANNEL_COUNT-channel Al lens" *
              (DIFFUSER ? " + diffuser" : ""), color=:viridis,
        colorbar_title="|uₓ|, nm", clims=(0, common_limit), aspect_ratio=:equal,
    )
    uniform_panel = heatmap(
        x_mm, y_mm, uniform_map_nm; xlabel="distance from radiator, mm", ylabel="y, mm",
        title="Matched $CHANNEL_COUNT-channel straight control" *
              (DIFFUSER ? " + diffuser" : ""), color=:viridis,
        colorbar_title="|uₓ|, nm", clims=(0, common_limit), aspect_ratio=:equal,
    )
    for panel in (lens_panel, uniform_panel)
        scatter!(panel, [FOCAL_DISTANCE_MM], [0.0]; marker=:xcross,
                 color=:white, markerstrokewidth=2, label=false)
    end
    profile_panel = plot(
        y_mm, lens_profile; label="$CHANNEL_COUNT-channel lens", linewidth=2.5,
        xlabel="y at x=60 mm, mm", ylabel="|uₓ|, nm",
        title="Physical focus profile", gridalpha=0.25,
    )
    plot!(profile_panel, y_mm, uniform_profile; label="straight control", linewidth=2)
    figure_path = joinpath(OUTPUT_ROOT, "$(CHANNEL_COUNT)_aperture_focus.png")
    savefig(plot(lens_panel, uniform_panel, profile_panel;
                 layout=(2, 2), size=(1450, 950), margin=5Plots.mm), figure_path)
    println("[+] $CHANNEL_COUNT-channel pressure gain=$pressure_gain")
    println("[+] power-normalized gain=$power_gain, ideal=$ideal_gain, " *
            "efficiency=$ideal_efficiency, passed=$passed")
    DIFFUSER && println("[+] versus abrupt uniform: pressure=" *
                        "$pressure_gain_over_abrupt, equal-power=$power_gain_over_abrupt")
    println("[+] FWHM=$(width.width) mm, peak y=$(width.peak_coordinate) mm")
    println("[+] $summary_path")
    println("[+] $figure_path")
end

wrap_phase_deg(value) = rad2deg(atan(sin(angle(value)), cos(angle(value))))

function best_bounded_phase_correction(focus_response, current_path_mm)
    phase_slope_deg_per_mm = 360 * FREQUENCY_HZ * GUIDE_DELAY_US_PER_MM * 1e-6
    cycle_mm = 360 / phase_slope_deg_per_mm
    best = nothing
    for target_deg in range(-180.0, 180.0; length=7201)
        correction_mm = zeros(length(focus_response))
        feasible = true
        for index in eachindex(focus_response)
            error_deg = wrap_phase_deg(focus_response[index] / cis(deg2rad(target_deg)))
            candidates = [error_deg / phase_slope_deg_per_mm + cycle * cycle_mm
                          for cycle in -2:2]
            candidates = filter(candidates) do delta_mm
                0.0 <= current_path_mm[index] + delta_mm <= 23.03
            end
            if isempty(candidates)
                feasible = false
                break
            end
            correction_mm[index] = candidates[argmin(abs.(candidates))]
        end
        feasible || continue
        objective = correction_mm[1]^2 + 2sum(abs2, correction_mm[2:end])
        if isnothing(best) || objective < best.objective
            best = (; target_deg, correction_mm, objective)
        end
    end
    isnothing(best) && error("no bounded carrier correction fits the isolated-pilot path range")
    best
end

function run_matrix_compare_stage()
    matrix = JLD2.load(matrix_path())
    lens_focus = vec(matrix["focus_response_m"][1, :])
    uniform = JLD2.load(result_path("uniform"))
    uniform_focus = ComplexF64(uniform["focus_displacement_m"][1])
    current_gain = abs(sum(lens_focus)) / abs(uniform_focus)
    coherence = abs(sum(lens_focus)) / sum(abs, lens_focus)
    aligned_ceiling = sum(abs, lens_focus) / abs(uniform_focus)
    half_indices = findall(>=(0.0), CENTERS_MM)
    current_half_path_mm = EXTRA_PATH_MM[half_indices]
    correction = best_bounded_phase_correction(lens_focus, current_half_path_mm)
    corrected_path_mm = current_half_path_mm .+ correction.correction_mm
    corrected_amplitude_mm = [
        smooth_amplitude_mm(AXIAL_LENGTH_MM, path_mm) for path_mm in corrected_path_mm
    ]
    phase_error_deg = [wrap_phase_deg(value / cis(deg2rad(correction.target_deg)))
                       for value in lens_focus]
    correction_allowed = aligned_ceiling >= 2.0

    summary_path = joinpath(OUTPUT_ROOT, "$(CHANNEL_COUNT)_response_matrix.csv")
    open(summary_path, "w") do io
        println(io, "group,y_mm,focus_amplitude_nm,focus_phase_deg,target_phase_deg,phase_error_deg,current_extra_path_mm,correction_mm,corrected_extra_path_mm,corrected_amplitude_mm")
        for index in eachindex(lens_focus)
            println(io, join((index, CENTERS_MM[half_indices[index]],
                              abs(lens_focus[index]) * 1e9, rad2deg(angle(lens_focus[index])),
                              correction.target_deg, phase_error_deg[index],
                              current_half_path_mm[index], correction.correction_mm[index],
                              corrected_path_mm[index], corrected_amplitude_mm[index]), ','))
        end
    end
    verdict_path = joinpath(OUTPUT_ROOT, "$(CHANNEL_COUNT)_matrix_verdict.txt")
    open(verdict_path, "w") do io
        println(io, "current_gain=$current_gain")
        println(io, "coherence_efficiency=$coherence")
        println(io, "phase_aligned_ceiling=$aligned_ceiling")
        println(io, "correction_allowed=$correction_allowed")
        println(io, "target_phase_deg=$(correction.target_deg)")
    end
    println("[+] full matrix current gain=$current_gain, coherence=$coherence, " *
            "phase-aligned ceiling=$aligned_ceiling")
    println("[+] correction allowed=$correction_allowed, target phase=" *
            "$(correction.target_deg) deg")
    println("[+] corrections=$(correction.correction_mm) mm")
    println("[+] corrected paths=$corrected_path_mm mm")
    println("[+] $summary_path")
end

function run_weighted_compare_stage()
    CHANNEL_COUNT == 15 || error("power-weight verification is defined for 15 channels")
    DIFFUSER || error("power-weight verification is defined for the output-diffuser lens")
    weighted = JLD2.load(weighted_result_path())
    design = JLD2.load(weight_design_path())
    abrupt_uniform = JLD2.load(joinpath(
        PROJECT_ROOT,
        "tmp",
        "aluminium_horn_ttd_full_aperture_242khz",
        "uniform_half_harmonic.jld2",
    ))
    diffuser_uniform = JLD2.load(result_path("uniform"))
    focus_nm = abs(weighted["focus_displacement_m"][1]) * 1e9
    input_power_w = abs(Float64(weighted["active_input_power_w"]))
    abrupt_focus_nm = abs(abrupt_uniform["focus_displacement_m"][1]) * 1e9
    abrupt_power_w = abs(Float64(abrupt_uniform["active_input_power_w"]))
    diffuser_focus_nm = abs(diffuser_uniform["focus_displacement_m"][1]) * 1e9
    diffuser_power_w = abs(Float64(diffuser_uniform["active_input_power_w"]))
    gain_over_abrupt = focus_nm / abrupt_focus_nm * sqrt(abrupt_power_w / input_power_w)
    gain_over_diffuser = focus_nm / diffuser_focus_nm * sqrt(diffuser_power_w / input_power_w)
    predicted_focus_nm = if haskey(design, "predicted_focus_amplitude_m")
        Float64(design["predicted_focus_amplitude_m"]) * 1e9
    else
        matrix = JLD2.load(matrix_path())
        response = vec(matrix["focus_response_m"][1, :])
        abs(sum(response .* Float64.(design["pressure_weights"]))) * 1e9
    end
    prediction_error = focus_nm / predicted_focus_nm - 1

    x_mm = weighted["scan_x_mm"]
    y_positive_mm = weighted["scan_y_mm"]
    y_mm = vcat(-reverse(y_positive_mm[2:end]), y_positive_mm)
    half_nm = abs.(weighted["scan_ux_m"]) .* 1e9
    map_nm = vcat(reverse(half_nm[2:end, :]; dims=1), half_nm)
    target_x_index = argmin(abs.(x_mm .- FOCAL_DISTANCE_MM))
    profile_nm = map_nm[:, target_x_index]
    width = contiguous_width(y_mm, profile_nm)
    focus_vector = weighted["focus_displacement_m"]
    transverse_ratio = (abs2(focus_vector[2]) + abs2(focus_vector[3])) /
                       abs2(focus_vector[1])
    passed = gain_over_abrupt >= 2.0 && abs(width.peak_coordinate) <= 4.0 &&
             transverse_ratio <= 0.10 && abs(prediction_error) <= 0.01

    result_label = WEIGHT_LABEL == "carrier" ? "weighted" : "weighted_$(WEIGHT_LABEL)"
    summary_path = joinpath(OUTPUT_ROOT, "15_$(result_label)_direct_summary.csv")
    open(summary_path, "w") do io
        println(io, "weighted_focus_ux_nm,weighted_active_power_w,abrupt_uniform_focus_ux_nm,abrupt_uniform_active_power_w,equal_power_gain_over_abrupt,diffuser_uniform_focus_ux_nm,diffuser_uniform_active_power_w,equal_power_gain_over_diffuser,predicted_focus_ux_nm,prediction_relative_error,fwhm_mm,profile_peak_y_mm,transverse_energy_ratio,passed")
        println(io, join((focus_nm, input_power_w, abrupt_focus_nm, abrupt_power_w,
                          gain_over_abrupt, diffuser_focus_nm, diffuser_power_w,
                          gain_over_diffuser, predicted_focus_nm, prediction_error,
                          width.width, width.peak_coordinate, transverse_ratio, passed), ','))
    end
    field_panel = heatmap(
        x_mm, y_mm, map_nm; xlabel="distance from diffuser, mm", ylabel="y, mm",
        title="15-channel Al lens: optimized power weights", color=:viridis,
        colorbar_title="|uₓ|, nm", aspect_ratio=:equal,
    )
    scatter!(field_panel, [FOCAL_DISTANCE_MM], [0.0]; marker=:xcross,
             color=:white, markerstrokewidth=2, label=false)
    profile_panel = plot(
        y_mm, profile_nm; label="weighted lens", linewidth=2.5,
        xlabel="y at x=60 mm, mm", ylabel="|uₓ|, nm",
        title="Physical focal profile", gridalpha=0.25,
    )
    figure_path = joinpath(OUTPUT_ROOT, "15_$(result_label)_direct_focus.png")
    savefig(plot(field_panel, profile_panel; layout=(2, 1), size=(1100, 1150),
                 margin=5Plots.mm), figure_path)
    println("[+] direct equal-power gain over abrupt uniform=$gain_over_abrupt")
    println("[+] direct equal-power gain over diffuser uniform=$gain_over_diffuser")
    println("[+] matrix prediction error=$(100prediction_error)%")
    println("[+] FWHM=$(width.width) mm, peak y=$(width.peak_coordinate) mm, passed=$passed")
    println("[+] $summary_path")
    println("[+] $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "solve" ? run_solve_stage() :
    STAGE == "matrix" ? run_matrix_stage() :
    STAGE == "weighted" ? run_weighted_stage() :
    STAGE == "compare" ? run_compare_stage() :
    STAGE == "matrix_compare" ? run_matrix_compare_stage() :
    STAGE == "weighted_compare" ? run_weighted_compare_stage() :
    error("unknown stage: $STAGE")
end

end
