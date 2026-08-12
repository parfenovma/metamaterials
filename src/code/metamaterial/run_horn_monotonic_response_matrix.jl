module HornMonotonicResponseMatrix

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HALF_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_array_3d_half_symmetry")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_RESPONSE_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_response_matrix"),
)

using JLD2
ENV["GKSwstype"] = "100"
using Plots

if !isdefined(parentmodule(@__MODULE__), :HornMonotonicArray3DHalfSymmetry)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "run_horn_monotonic_array_3d_half_symmetry.jl"),
    )
end

using ..HornMonotonicArray3DHalfSymmetry
using ..HornMonotonicArray3DHarmonicSolver:
    HornMonotonicArray3DHarmonicConfig,
    solve_horn_monotonic_array_3d_harmonic,
    solve_horn_monotonic_array_3d_response_matrix
using ..HornMonotonicArray3DConvergence:
    convergence_probe_points_mm,
    average_channel_outputs
using ..HornMonotonicArray3DMesher: channel_centers_mm, outlet_x_mm
using ..SinusoidalMaterialLens: photopolymer
using ..ImpulseRiskAnalysis: contiguous_width

export energy_normalize,
       focus_optimal_weights,
       profile_metrics,
       input_power_w,
       optimize_focus_weights

const LEVEL_ID = "phase2"
const VERIFY_LEVEL_ID = "phase4"
const CARRIER_FREQUENCY_HZ = 242.0e3
const SOURCE_MULTIPLICITY = Float64[1.0; fill(2.0, 7)]

function frequency_tag(frequency_hz::Real)
    frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    tenths_khz = round(Int, Float64(frequency_hz) / 100)
    whole_khz, tenth_khz = divrem(tenths_khz, 10)
    "$(whole_khz)p$(tenth_khz)khz"
end

frequency_suffix(frequency_hz::Real) =
    isapprox(frequency_hz, CARRIER_FREQUENCY_HZ; atol=1e-6, rtol=0) ?
    "" : "_$(frequency_tag(frequency_hz))"

response_path(variant, frequency_hz::Real=CARRIER_FREQUENCY_HZ) = joinpath(
    OUTPUT_ROOT,
    "response_matrix_$(variant)_$(LEVEL_ID)$(frequency_suffix(frequency_hz)).jld2",
)
basis_path(variant, source_index, frequency_hz::Real=CARRIER_FREQUENCY_HZ) = joinpath(
    OUTPUT_ROOT,
    "basis_$(variant)_$(LEVEL_ID)_source$(source_index)$(frequency_suffix(frequency_hz)).jld2",
)
field_path(variant) = joinpath(HALF_ROOT, "field_half_$(variant)_$(LEVEL_ID).jld2")
mesh_path(variant) = joinpath(HALF_ROOT, "meshes", "mesh_half_$(variant)_$(LEVEL_ID).msh")
verify_field_path(variant) =
    joinpath(OUTPUT_ROOT, "weighted_$(variant)_$(VERIFY_LEVEL_ID).jld2")
verify_baseline_path(variant) =
    joinpath(HALF_ROOT, "field_half_$(variant)_$(VERIFY_LEVEL_ID).jld2")
verify_mesh_path(variant) =
    joinpath(HALF_ROOT, "meshes", "mesh_half_$(variant)_$(VERIFY_LEVEL_ID).msh")

function parse_options(args=ARGS)
    options = Dict{String, String}()
    for argument in args
        startswith(argument, "--") || throw(ArgumentError("unknown argument: $argument"))
        fields = split(argument[3:end], '='; limit=2)
        length(fields) == 2 || throw(ArgumentError("expected --name=value: $argument"))
        options[fields[1]] = fields[2]
    end
    stage = get(options, "stage", "analyze")
    variant = get(options, "variant", "all")
    source_index = parse(Int, get(options, "source", "0"))
    frequency_hz = parse(Float64, get(options, "frequency-khz", "242.0")) * 1e3
    frequency_hz > 0 || throw(ArgumentError("--frequency-khz must be positive"))
    stage in ("solve", "aggregate", "verify", "analyze") ||
        throw(ArgumentError("--stage must be solve, aggregate, verify, or analyze"))
    variant in ("lens", "uniform", "all") ||
        throw(ArgumentError("--variant must be lens, uniform, or all"))
    stage == "solve" && (variant == "all" || !(1 <= source_index <= 8)) &&
        throw(ArgumentError("solve requires one variant and --source=1,...,8"))
    stage == "aggregate" && variant == "all" &&
        throw(ArgumentError("aggregate requires one variant"))
    stage == "verify" && variant == "all" &&
        throw(ArgumentError("verify requires one variant"))
    stage == "analyze" && variant != "all" &&
        throw(ArgumentError("analyze consumes both variants"))
    (; stage, variant, source_index, frequency_hz)
end

function energy_normalize(weights, multiplicity=SOURCE_MULTIPLICITY; budget=sum(multiplicity))
    length(weights) == length(multiplicity) || throw(DimensionMismatch("one multiplicity per weight"))
    all(>=(0), weights) || throw(ArgumentError("source weights must be non-negative"))
    energy = sum(multiplicity .* abs2.(weights))
    energy > 0 || throw(ArgumentError("at least one source weight must be positive"))
    Float64.(weights) .* sqrt(Float64(budget) / energy)
end

focus_amplitude(focus_response, weights) = abs(sum(focus_response .* weights))

"""Global non-negative real optimum for one complex focal observable."""
function focus_optimal_weights(
    focus_response,
    multiplicity=SOURCE_MULTIPLICITY;
    phase_samples::Integer=7200,
)
    length(focus_response) == length(multiplicity) ||
        throw(DimensionMismatch("one multiplicity per source response"))
    phase_samples >= 360 || throw(ArgumentError("phase sweep is undersampled"))
    best_weights = ones(Float64, length(focus_response))
    best_amplitude = -Inf
    for phase in range(0.0, 2pi; length=phase_samples + 1)[1:end-1]
        projected = max.(real.(cis(-phase) .* focus_response), 0.0)
        all(iszero, projected) && continue
        candidate = energy_normalize(projected ./ multiplicity, multiplicity)
        amplitude = focus_amplitude(focus_response, candidate)
        if amplitude > best_amplitude
            best_amplitude = amplitude
            best_weights = candidate
        end
    end
    best_weights
end

function reconstruct_full_profile(nonnegative_coordinate, nonnegative_complex_profile)
    length(nonnegative_coordinate) == length(nonnegative_complex_profile) ||
        throw(DimensionMismatch("coordinate and profile lengths differ"))
    first(nonnegative_coordinate) == 0 ||
        throw(ArgumentError("half-domain profile must start at zero"))
    coordinate = vcat(-reverse(nonnegative_coordinate[2:end]), nonnegative_coordinate)
    complex_profile = vcat(
        reverse(nonnegative_complex_profile[2:end]),
        nonnegative_complex_profile,
    )
    (; coordinate, complex_profile, amplitude=abs.(complex_profile))
end

function sidelobe_ratio(profile)
    center_peak = argmax(profile)
    left_minimum = center_peak
    while left_minimum > firstindex(profile) + 1
        left_minimum -= 1
        profile[left_minimum] <= profile[left_minimum - 1] &&
            profile[left_minimum] <= profile[left_minimum + 1] && break
    end
    right_minimum = center_peak
    while right_minimum < lastindex(profile) - 1
        right_minimum += 1
        profile[right_minimum] <= profile[right_minimum - 1] &&
            profile[right_minimum] <= profile[right_minimum + 1] && break
    end
    outside = vcat(
        collect(firstindex(profile):(left_minimum - 1)),
        collect((right_minimum + 1):lastindex(profile)),
    )
    isempty(outside) ? 0.0 : maximum(profile[outside]) / profile[center_peak]
end

function profile_metrics(coordinate_mm, complex_profile_m)
    amplitude_nm = abs.(complex_profile_m) .* 1e9
    width = contiguous_width(coordinate_mm, amplitude_nm)
    (
        coordinate_mm=Float64.(coordinate_mm),
        amplitude_nm,
        fwhm_mm=width.width,
        peak_y_mm=width.peak_coordinate,
        peak_nm=width.peak_amplitude,
        sidelobe_amplitude_ratio=sidelobe_ratio(amplitude_nm),
    )
end

function weighted_profile(scan_y_mm, scan_ux_m, weights)
    half = scan_ux_m * weights
    full = reconstruct_full_profile(scan_y_mm, half)
    profile_metrics(full.coordinate, full.complex_profile)
end

function input_power_w(
    source_normal_displacement_integral_m3,
    weights,
    frequency_hz,
    pressure_amplitude_pa,
)
    size(source_normal_displacement_integral_m3) == (length(weights), length(weights)) ||
        throw(DimensionMismatch("source-work matrix must be square in source weights"))
    normal_displacement = source_normal_displacement_integral_m3 * weights
    velocity_integral = im * 2pi * frequency_hz .* normal_displacement
    complex_power = 0.5sum(
        (-pressure_amplitude_pa .* weights) .* conj.(velocity_integral),
    )
    (; complex_va=complex_power, active_w=real(complex_power), reactive_var=imag(complex_power))
end

function optimize_focus_weights(
    focus_response_m,
    scan_y_mm,
    scan_ux_m,
    multiplicity=SOURCE_MULTIPLICITY;
    maximum_fwhm_mm::Real=6.0,
    maximum_sidelobe_ratio::Real=0.40,
    blend_samples::Integer=1001,
)
    target = focus_optimal_weights(focus_response_m, multiplicity)
    equal = energy_normalize(ones(length(target)), multiplicity)
    best = nothing
    for alpha in range(0.0, 1.0; length=blend_samples)
        weights = energy_normalize((1 - alpha) .* equal .+ alpha .* target, multiplicity)
        metrics = weighted_profile(scan_y_mm, scan_ux_m, weights)
        feasible = metrics.fwhm_mm <= maximum_fwhm_mm &&
                   metrics.sidelobe_amplitude_ratio <= maximum_sidelobe_ratio &&
                   abs(metrics.peak_y_mm) <= 1.0
        feasible || continue
        amplitude_m = focus_amplitude(focus_response_m, weights)
        if isnothing(best) || amplitude_m > best.focus_amplitude_m
            best = (; weights, alpha, focus_amplitude_m=amplitude_m, profile=metrics)
        end
    end
    isnothing(best) && error("no focus-weight candidate satisfies FWHM/sidelobe constraints")
    merge(best, (; unconstrained_weights=target))
end

function run_solve(variant, source_index, frequency_hz=CARRIER_FREQUENCY_HZ)
    mkpath(OUTPUT_ROOT)
    config = HornMonotonicArray3DHalfSymmetry.array_config(variant)
    positive_indices = 8:15
    probes = convergence_probe_points_mm(
        config;
        channel_indices=positive_indices,
        nonnegative_half=true,
    )
    solver_config = HornMonotonicArray3DHarmonicConfig(
        frequency_hz=Float64(frequency_hz),
        scan_x_min_mm=config.focal_distance_mm,
        scan_x_max_mm=config.focal_distance_mm,
        scan_y_half_width_mm=20.0,
        scan_y_min_mm=0.0,
        scan_step_mm=1.0,
        symmetry_y=true,
    )
    source_weights = zeros(Float64, length(positive_indices))
    source_weights[source_index] = 1.0
    result = solve_horn_monotonic_array_3d_harmonic(
        mesh_path(variant),
        length(positive_indices),
        outlet_x_mm(config),
        probes,
        photopolymer();
        config=solver_config,
        source_weights,
    )
    channel_response_m = average_channel_outputs(
        result.probe_displacement_m,
        probes.channel_groups,
    )
    jldsave(
        basis_path(variant, source_index, frequency_hz);
        format_version=1,
        level_id=LEVEL_ID,
        variant,
        source_index,
        frequency_hz=solver_config.frequency_hz,
        pressure_amplitude_pa=solver_config.pressure_amplitude_pa,
        channel_response_m,
        focus_response_m=result.focus_displacement_m[1],
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=vec(result.scan_ux_m),
        scan_uy_m=vec(result.scan_uy_m),
        scan_uz_m=vec(result.scan_uz_m),
        source_normal_displacement_integral_m3=result.source_normal_displacement_integral_m3,
        complex_input_power_va=result.complex_input_power_va,
    )
    println("[+] $variant source $source_index at $(frequency_hz / 1e3) kHz focus contribution=$(abs(result.focus_displacement_m[1]) * 1e9) nm")
    println("[+] $(basis_path(variant, source_index, frequency_hz))")
end

function run_aggregate(variant, frequency_hz=CARRIER_FREQUENCY_HZ)
    config = HornMonotonicArray3DHalfSymmetry.array_config(variant)
    positive_indices = 8:15
    columns = [
        JLD2.load(basis_path(variant, index, frequency_hz))
        for index in eachindex(positive_indices)
    ]
    all(isapprox(column["frequency_hz"], frequency_hz; atol=1e-6, rtol=0) for column in columns) ||
        error("basis columns do not match the requested frequency")
    response_matrix_m = hcat(getindex.(columns, "channel_response_m")...)
    focus_response_m = ComplexF64[getindex(column, "focus_response_m") for column in columns]
    scan_ux_m = hcat(getindex.(columns, "scan_ux_m")...)
    scan_uy_m = hcat(getindex.(columns, "scan_uy_m")...)
    scan_uz_m = hcat(getindex.(columns, "scan_uz_m")...)
    source_normal_displacement_integral_m3 = hcat(
        getindex.(columns, "source_normal_displacement_integral_m3")...,
    )
    at_carrier = isapprox(frequency_hz, CARRIER_FREQUENCY_HZ; atol=1e-6, rtol=0)
    baseline_channel_error, baseline_focus_error = if at_carrier
        baseline = JLD2.load(field_path(variant))
        channel_error = maximum(abs.(
            sum(response_matrix_m; dims=2)[:] .-
            baseline["nonnegative_channel_output_ux_m"],
        )) / maximum(abs, baseline["nonnegative_channel_output_ux_m"])
        focus_error = abs(
            sum(focus_response_m) - baseline["focus_displacement_m"][1],
        ) / abs(baseline["focus_displacement_m"][1])
        channel_error <= 1e-8 || error("response matrix failed channel superposition check")
        focus_error <= 1e-8 || error("response matrix failed focus superposition check")
        (channel_error, focus_error)
    else
        (NaN, NaN)
    end
    centers_mm = channel_centers_mm(config)[positive_indices]
    jldsave(
        response_path(variant, frequency_hz);
        format_version=1,
        level_id=LEVEL_ID,
        variant,
        frequency_hz=columns[1]["frequency_hz"],
        pressure_amplitude_pa=columns[1]["pressure_amplitude_pa"],
        nonnegative_channel_indices=collect(positive_indices),
        nonnegative_channel_centers_mm=centers_mm,
        source_multiplicity=SOURCE_MULTIPLICITY,
        response_matrix_m,
        focus_response_m,
        scan_y_mm=columns[1]["scan_y_mm"],
        scan_ux_m,
        scan_uy_m,
        scan_uz_m,
        source_normal_displacement_integral_m3,
        baseline_channel_superposition_relative_error=baseline_channel_error,
        baseline_focus_superposition_relative_error=baseline_focus_error,
    )
    if at_carrier
        println("[+] $variant response matrix focus reconstruction error=$baseline_focus_error")
    else
        println("[+] $variant response matrix assembled at $(frequency_hz / 1e3) kHz")
    end
    println("[+] $(response_path(variant, frequency_hz))")
end

function run_verify(variant)
    design = JLD2.load(joinpath(OUTPUT_ROOT, "horn_aperture_precompensation.jld2"))
    weights = Float64.(design["source_weights"])
    config = HornMonotonicArray3DHalfSymmetry.array_config(variant)
    positive_indices = 8:15
    probes = convergence_probe_points_mm(
        config;
        channel_indices=positive_indices,
        nonnegative_half=true,
    )
    solver_config = HornMonotonicArray3DHarmonicConfig(
        scan_x_min_mm=config.focal_distance_mm,
        scan_x_max_mm=config.focal_distance_mm,
        scan_y_half_width_mm=12.0,
        scan_y_min_mm=0.0,
        scan_step_mm=1.0,
        symmetry_y=true,
    )
    result = solve_horn_monotonic_array_3d_harmonic(
        verify_mesh_path(variant),
        length(positive_indices),
        outlet_x_mm(config),
        probes,
        photopolymer();
        config=solver_config,
        source_weights=weights,
    )
    channel_output_ux_m = average_channel_outputs(
        result.probe_displacement_m,
        probes.channel_groups,
    )
    jldsave(
        verify_field_path(variant);
        format_version=1,
        level_id=VERIFY_LEVEL_ID,
        variant,
        frequency_hz=solver_config.frequency_hz,
        pressure_amplitude_pa=solver_config.pressure_amplitude_pa,
        source_weights=weights,
        nonnegative_channel_output_ux_m=channel_output_ux_m,
        focus_displacement_m=result.focus_displacement_m,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=vec(result.scan_ux_m),
    )
    println("[+] verified $variant weighted focus=$(abs(result.focus_displacement_m[1]) * 1e9) nm")
    println("[+] $(verify_field_path(variant))")
end

positive_power(power) = abs(power.active_w)

function run_analysis()
    mkpath(OUTPUT_ROOT)
    lens = JLD2.load(response_path("lens"))
    uniform = JLD2.load(response_path("uniform"))
    multiplicity = Float64.(lens["source_multiplicity"])
    equal = energy_normalize(ones(length(multiplicity)), multiplicity)
    optimized = optimize_focus_weights(
        lens["focus_response_m"],
        lens["scan_y_mm"],
        lens["scan_ux_m"],
        multiplicity,
    )
    corrected = optimized.weights

    lens_equal_profile = weighted_profile(lens["scan_y_mm"], lens["scan_ux_m"], equal)
    lens_corrected_profile = weighted_profile(lens["scan_y_mm"], lens["scan_ux_m"], corrected)
    uniform_equal_profile = weighted_profile(uniform["scan_y_mm"], uniform["scan_ux_m"], equal)
    uniform_corrected_profile = weighted_profile(
        uniform["scan_y_mm"],
        uniform["scan_ux_m"],
        corrected,
    )
    lens_equal_focus_m = focus_amplitude(lens["focus_response_m"], equal)
    lens_corrected_focus_m = focus_amplitude(lens["focus_response_m"], corrected)
    uniform_equal_focus_m = focus_amplitude(uniform["focus_response_m"], equal)
    uniform_corrected_focus_m = focus_amplitude(uniform["focus_response_m"], corrected)

    nominal_equal_gain = lens_equal_focus_m / uniform_equal_focus_m
    nominal_corrected_gain = lens_corrected_focus_m / uniform_corrected_focus_m
    exact_source_power_available = all(isfinite, real.(
        lens["source_normal_displacement_integral_m3"],
    )) && all(isfinite, real.(uniform["source_normal_displacement_integral_m3"]))
    powers = if exact_source_power_available
        lens_equal_power = input_power_w(
            lens["source_normal_displacement_integral_m3"], equal,
            lens["frequency_hz"], lens["pressure_amplitude_pa"],
        )
        lens_corrected_power = input_power_w(
            lens["source_normal_displacement_integral_m3"], corrected,
            lens["frequency_hz"], lens["pressure_amplitude_pa"],
        )
        uniform_equal_power = input_power_w(
            uniform["source_normal_displacement_integral_m3"], equal,
            uniform["frequency_hz"], uniform["pressure_amplitude_pa"],
        )
        uniform_corrected_power = input_power_w(
            uniform["source_normal_displacement_integral_m3"], corrected,
            uniform["frequency_hz"], uniform["pressure_amplitude_pa"],
        )
        (
            positive_power(lens_equal_power),
            positive_power(lens_corrected_power),
            positive_power(uniform_equal_power),
            positive_power(uniform_corrected_power),
        )
    else
        (1.0, 1.0, 1.0, 1.0)
    end
    energy_equal_gain = nominal_equal_gain * sqrt(powers[3] / powers[1])
    energy_corrected_gain = nominal_corrected_gain * sqrt(powers[4] / powers[2])
    gain_improvement = energy_corrected_gain / energy_equal_gain - 1
    lens_equal_power_focus_m = lens_equal_focus_m / sqrt(powers[1])
    lens_corrected_power_focus_m = lens_corrected_focus_m / sqrt(powers[2])
    lens_focus_improvement = lens_corrected_power_focus_m / lens_equal_power_focus_m - 1
    proxy_passed = gain_improvement >= 0.15 &&
                   energy_corrected_gain >= 2.2 &&
                   lens_corrected_profile.fwhm_mm <= 6.0 &&
                   lens_corrected_profile.sidelobe_amplitude_ratio <= 0.40

    verification_available = isfile(verify_field_path("lens")) &&
                             isfile(verify_field_path("uniform"))
    if verification_available
        verified_lens = JLD2.load(verify_field_path("lens"))
        verified_uniform = JLD2.load(verify_field_path("uniform"))
        verified_lens_baseline = JLD2.load(verify_baseline_path("lens"))
        verified_uniform_baseline = JLD2.load(verify_baseline_path("uniform"))
        verified_lens_focus_nm = abs(verified_lens["focus_displacement_m"][1]) * 1e9
        verified_uniform_focus_nm = abs(verified_uniform["focus_displacement_m"][1]) * 1e9
        verified_baseline_lens_focus_nm =
            abs(verified_lens_baseline["focus_displacement_m"][1]) * 1e9
        verified_baseline_uniform_focus_nm =
            abs(verified_uniform_baseline["focus_displacement_m"][1]) * 1e9
        verified_gain = verified_lens_focus_nm / verified_uniform_focus_nm
        verified_baseline_gain =
            verified_baseline_lens_focus_nm / verified_baseline_uniform_focus_nm
        verified_lens_improvement =
            verified_lens_focus_nm / verified_baseline_lens_focus_nm - 1
        verified_gain_improvement = verified_gain / verified_baseline_gain - 1
        verified_half = reconstruct_full_profile(
            verified_lens["scan_y_mm"],
            verified_lens["scan_ux_m"],
        )
        verified_profile = profile_metrics(
            verified_half.coordinate,
            verified_half.complex_profile,
        )
        verification_path = joinpath(OUTPUT_ROOT, "phase4_verification_summary.csv")
        open(verification_path, "w") do io
            println(io, "baseline_lens_focus_nm,weighted_lens_focus_nm,lens_focus_improvement,baseline_uniform_focus_nm,weighted_uniform_focus_nm,baseline_gain,weighted_gain,gain_improvement,weighted_fwhm_mm")
            println(io, join((
                verified_baseline_lens_focus_nm,
                verified_lens_focus_nm,
                verified_lens_improvement,
                verified_baseline_uniform_focus_nm,
                verified_uniform_focus_nm,
                verified_baseline_gain,
                verified_gain,
                verified_gain_improvement,
                verified_profile.fwhm_mm,
            ), ','))
        end
        println("[+] phase4 weighted lens improvement=$(100verified_lens_improvement)%")
        println("[+] phase4 matched gain: $verified_baseline_gain -> $verified_gain")
        println("[+] $verification_path")
    end

    weights_path = joinpath(OUTPUT_ROOT, "source_weights.csv")
    open(weights_path, "w") do io
        println(io, "source_group,positive_channel_index,center_y_mm,multiplicity,equal_weight,corrected_weight,unconstrained_weight")
        for index in eachindex(corrected)
            println(io, join((
                index,
                lens["nonnegative_channel_indices"][index],
                lens["nonnegative_channel_centers_mm"][index],
                multiplicity[index],
                equal[index],
                corrected[index],
                optimized.unconstrained_weights[index],
            ), ','))
        end
    end

    matrix_path = joinpath(OUTPUT_ROOT, "response_matrix.csv")
    open(matrix_path, "w") do io
        println(io, "variant,output_group,input_group,output_center_y_mm,input_center_y_mm,ux_amplitude_nm,ux_phase_deg")
        for (variant, data) in (("lens", lens), ("uniform", uniform))
            matrix = data["response_matrix_m"]
            centers = data["nonnegative_channel_centers_mm"]
            for output_group in axes(matrix, 1), input_group in axes(matrix, 2)
                value = matrix[output_group, input_group]
                println(io, join((
                    variant,
                    output_group,
                    input_group,
                    centers[output_group],
                    centers[input_group],
                    abs(value) * 1e9,
                    rad2deg(angle(value)),
                ), ','))
            end
        end
    end

    summary_path = joinpath(OUTPUT_ROOT, "precompensation_proxy_summary.csv")
    open(summary_path, "w") do io
        println(io, "baseline_lens_focus_nm,corrected_lens_focus_nm,baseline_uniform_focus_nm,corrected_uniform_focus_nm,baseline_nominal_gain,corrected_nominal_gain,baseline_equal_energy_gain,corrected_equal_energy_gain,equal_energy_gain_improvement,lens_equal_energy_focus_improvement,corrected_fwhm_mm,corrected_sidelobe_amplitude_ratio,blend_alpha,exact_source_power_available,proxy_passed")
        println(io, join((
            lens_equal_focus_m * 1e9,
            lens_corrected_focus_m * 1e9,
            uniform_equal_focus_m * 1e9,
            uniform_corrected_focus_m * 1e9,
            nominal_equal_gain,
            nominal_corrected_gain,
            energy_equal_gain,
            energy_corrected_gain,
            gain_improvement,
            lens_focus_improvement,
            lens_corrected_profile.fwhm_mm,
            lens_corrected_profile.sidelobe_amplitude_ratio,
            optimized.alpha,
            exact_source_power_available,
            proxy_passed,
        ), ','))
    end

    design_path = joinpath(OUTPUT_ROOT, "horn_aperture_precompensation.jld2")
    jldsave(
        design_path;
        format_version=1,
        frequency_hz=lens["frequency_hz"],
        source_multiplicity=multiplicity,
        source_weights=corrected,
        unconstrained_source_weights=optimized.unconstrained_weights,
        blend_alpha=optimized.alpha,
        nominal_equal_gain,
        nominal_corrected_gain,
        energy_equal_gain,
        energy_corrected_gain,
        gain_improvement,
        lens_focus_improvement,
        corrected_fwhm_mm=lens_corrected_profile.fwhm_mm,
        corrected_sidelobe_amplitude_ratio=lens_corrected_profile.sidelobe_amplitude_ratio,
        exact_source_power_available,
        proxy_passed,
    )

    weights_panel = plot(
        lens["nonnegative_channel_centers_mm"],
        corrected;
        marker=:circle,
        linewidth=2.5,
        xlabel="channel centre y, mm",
        ylabel="pressure multiplier",
        title="Symmetric source correction",
        label="corrected",
        gridalpha=0.25,
    )
    hline!(weights_panel, [1.0]; linestyle=:dash, label="equal drive")
    profile_panel = plot(
        lens_equal_profile.coordinate_mm,
        lens_equal_profile.amplitude_nm;
        linewidth=2.2,
        xlabel="y at x=35 mm, mm",
        ylabel="|ux|, nm",
        title="Lens focal profile",
        label="equal drive",
        gridalpha=0.25,
    )
    plot!(
        profile_panel,
        lens_corrected_profile.coordinate_mm,
        lens_corrected_profile.amplitude_nm;
        linewidth=2.7,
        label="corrected",
    )
    gain_panel = bar(
        ["equal", "corrected"],
        [energy_equal_gain, energy_corrected_gain];
        ylabel="matched gain at equal nominal energy",
        title="Precompensation proxy",
        label=false,
        gridalpha=0.25,
    )
    hline!(gain_panel, [2.2]; linestyle=:dash, label="2.2 target")
    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_response_matrix.png")
    savefig(plot(
        weights_panel,
        profile_panel,
        gain_panel;
        layout=(1, 3),
        size=(1700, 560),
        margin=5Plots.mm,
    ), figure_path)

    println("[+] equal-nominal-energy matched gain: $energy_equal_gain -> $energy_corrected_gain")
    println("[+] equal-nominal-energy gain improvement=$(100gain_improvement)%")
    println("[+] lens equal-nominal-energy focal amplitude improvement=$(100lens_focus_improvement)%")
    println("[+] corrected FWHM=$(lens_corrected_profile.fwhm_mm) mm")
    println("[+] corrected sidelobe=$(lens_corrected_profile.sidelobe_amplitude_ratio)")
    println("[+] proxy gate passed=$proxy_passed")
    println("[+] $summary_path")
    println("[+] $figure_path")
end

function main(args=ARGS)
    options = parse_options(args)
    if options.stage == "solve"
        run_solve(options.variant, options.source_index, options.frequency_hz)
    elseif options.stage == "aggregate"
        run_aggregate(options.variant, options.frequency_hz)
    elseif options.stage == "verify"
        run_verify(options.variant)
    else
        run_analysis()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
