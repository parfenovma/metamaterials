module HornMonotonicGroupDelay

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_GROUP_DELAY_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_group_delay"),
)
const DESIGN_PATH = joinpath(
    PROJECT_ROOT,
    "tmp",
    "horn_monotonic_aperture",
    "monotonic_aperture_impulse.jld2",
)
const FREQUENCIES_HZ = [193.8e3, 242.0e3, 290.1e3]

using JLD2
using Statistics: mean
ENV["GKSwstype"] = "100"
using Plots

if !isdefined(parentmodule(@__MODULE__), :HornMonotonicResponseMatrix)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "run_horn_monotonic_response_matrix.jl"),
    )
end

using ..HornMonotonicResponseMatrix
using ..MonotonicHornLens:
    MonotonicHornLensConfig,
    geometric_delays_s,
    smooth_centerline_mm,
    smooth_amplitude_mm,
    minimum_inner_radius_mm,
    minimum_common_axial_length_mm

export unwrap_phase,
       linear_group_delay_s,
       isotonic_nonincreasing,
       effective_group_speed_m_s,
       corrected_monotonic_paths_mm,
       analyze_group_delay

function unwrap_phase(phase_rad)
    isempty(phase_rad) && return Float64[]
    result = Float64.(phase_rad)
    for index in 2:length(result)
        result[index] += 2pi * round((result[index - 1] - result[index]) / (2pi))
    end
    result
end

function linear_group_delay_s(frequencies_hz, transfer)
    length(frequencies_hz) == length(transfer) ||
        throw(DimensionMismatch("one complex transfer value per frequency"))
    length(frequencies_hz) >= 3 || throw(ArgumentError("at least three anchors are required"))
    issorted(frequencies_hz) || throw(ArgumentError("frequencies must be sorted"))
    all(!iszero, transfer) || throw(ArgumentError("transfer contains a zero"))
    omega = 2pi .* Float64.(frequencies_hz)
    phase_rad = unwrap_phase(angle.(transfer))
    centered_omega = omega .- mean(omega)
    slope_s = sum(centered_omega .* (phase_rad .- mean(phase_rad))) /
              sum(abs2, centered_omega)
    fitted_phase_rad = mean(phase_rad) .+ slope_s .* centered_omega
    (
        delay_s=-slope_s,
        phase_rad,
        fitted_phase_rad,
        maximum_residual_deg=maximum(abs.(rad2deg.(phase_rad .- fitted_phase_rad))),
    )
end

function isotonic_nondecreasing(values, weights=ones(length(values)))
    length(values) == length(weights) || throw(DimensionMismatch("one weight per value"))
    all(>(0), weights) || throw(ArgumentError("isotonic weights must be positive"))
    block_value = Float64[]
    block_weight = Float64[]
    block_first = Int[]
    block_last = Int[]
    for index in eachindex(values)
        push!(block_value, Float64(values[index]))
        push!(block_weight, Float64(weights[index]))
        push!(block_first, index)
        push!(block_last, index)
        while length(block_value) >= 2 && block_value[end - 1] > block_value[end]
            merged_weight = block_weight[end - 1] + block_weight[end]
            merged_value = (
                block_weight[end - 1] * block_value[end - 1] +
                block_weight[end] * block_value[end]
            ) / merged_weight
            block_value[end - 1] = merged_value
            block_weight[end - 1] = merged_weight
            block_last[end - 1] = block_last[end]
            pop!(block_value)
            pop!(block_weight)
            pop!(block_first)
            pop!(block_last)
        end
    end
    result = zeros(Float64, length(values))
    for block in eachindex(block_value)
        result[block_first[block]:block_last[block]] .= block_value[block]
    end
    result
end

isotonic_nonincreasing(values, weights=ones(length(values))) =
    -isotonic_nondecreasing(-Float64.(values), weights)

function effective_group_speed_m_s(extra_path_mm, relative_delay_s)
    length(extra_path_mm) == length(relative_delay_s) ||
        throw(DimensionMismatch("one delay per path length"))
    extra_path_m = Float64.(extra_path_mm) .* 1e-3
    delay_s = Float64.(relative_delay_s)
    slope_s_m = sum(extra_path_m .* delay_s) / sum(abs2, extra_path_m)
    slope_s_m > 0 || error("measured delay has the wrong sign for the extra path")
    inv(slope_s_m)
end

function corrected_monotonic_paths_mm(
    old_extra_path_mm,
    measured_delay_s,
    target_delay_s,
    group_speed_m_s,
    ;
    reliable=trues(length(old_extra_path_mm)),
)
    length(old_extra_path_mm) == length(measured_delay_s) == length(target_delay_s) ||
        throw(DimensionMismatch("path and delay arrays must have the same length"))
    group_speed_m_s > 0 || throw(ArgumentError("group speed must be positive"))
    length(reliable) == length(old_extra_path_mm) ||
        throw(DimensionMismatch("one reliability flag per path"))
    any(reliable) || throw(ArgumentError("at least one channel must be reliable"))
    target_delay_s = Float64.(target_delay_s)
    raw_mm = fill(NaN, length(old_extra_path_mm))
    raw_mm[reliable] .= Float64.(old_extra_path_mm)[reliable] .+
                        group_speed_m_s .*
                        (target_delay_s[reliable] .- measured_delay_s[reliable]) .* 1e3
    for index in eachindex(raw_mm)
        reliable[index] && continue
        left = findlast(reliable[begin:(index - 1)])
        right_local = findfirst(reliable[(index + 1):end])
        right = isnothing(right_local) ? nothing : index + right_local
        if isnothing(left)
            raw_mm[index] = raw_mm[right] +
                            group_speed_m_s *
                            (target_delay_s[index] - target_delay_s[right]) * 1e3
        elseif isnothing(right)
            raw_mm[index] = raw_mm[left] +
                            group_speed_m_s *
                            (target_delay_s[index] - target_delay_s[left]) * 1e3
        else
            fraction = (target_delay_s[index] - target_delay_s[left]) /
                       (target_delay_s[right] - target_delay_s[left])
            raw_mm[index] = raw_mm[left] + fraction * (raw_mm[right] - raw_mm[left])
        end
    end
    raw_mm = max.(raw_mm, 0.0)
    raw_mm[end] = 0.0
    corrected_mm = isotonic_nonincreasing(raw_mm)
    (; raw_mm, corrected_mm)
end

reconstruct_symmetric(nonnegative) = vcat(reverse(nonnegative[2:end]), nonnegative)

function response_data(frequency_hz)
    lens = JLD2.load(HornMonotonicResponseMatrix.response_path("lens", frequency_hz))
    uniform = JLD2.load(HornMonotonicResponseMatrix.response_path("uniform", frequency_hz))
    (; lens, uniform)
end

function analyze_group_delay(; frequencies_hz=FREQUENCIES_HZ)
    mkpath(OUTPUT_ROOT)
    frequencies_hz = Float64.(frequencies_hz)
    issorted(frequencies_hz) || throw(ArgumentError("frequency anchors must be sorted"))
    data = response_data.(frequencies_hz)
    channel_count = size(first(data).lens["response_matrix_m"], 1)
    device_transfer = Matrix{ComplexF64}(undef, length(frequencies_hz), channel_count)
    for frequency_index in eachindex(frequencies_hz)
        lens_matrix = data[frequency_index].lens["response_matrix_m"]
        uniform_matrix = data[frequency_index].uniform["response_matrix_m"]
        size(lens_matrix) == size(uniform_matrix) == (channel_count, channel_count) ||
            error("response matrices have inconsistent dimensions")
        for channel_index in 1:channel_count
            device_transfer[frequency_index, channel_index] =
                lens_matrix[channel_index, channel_index] /
                uniform_matrix[channel_index, channel_index]
        end
    end

    delay_fit = [
        linear_group_delay_s(frequencies_hz, device_transfer[:, channel])
        for channel in 1:channel_count
    ]
    absolute_delay_s = getproperty.(delay_fit, :delay_s)
    measured_delay_s = absolute_delay_s .- absolute_delay_s[end]
    maximum_phase_fit_residual_deg = maximum(getproperty.(delay_fit, :maximum_residual_deg))
    minimum_anchor_amplitude_ratio = vec(minimum(abs.(device_transfer); dims=1))
    delay_reliable = getproperty.(delay_fit, :maximum_residual_deg) .<= 30.0 .&&
                     minimum_anchor_amplitude_ratio .>= 0.50
    delay_reliable[end] || error("zero-extra-path reference channel failed the delay gate")

    design = JLD2.load(DESIGN_PATH)
    positive_indices = 8:15
    centers_mm = Float64.(design["centers_mm"])[positive_indices]
    old_extra_path_mm = Float64.(design["extra_path_mm"])[positive_indices]
    target_delay_s = geometric_delays_s(MonotonicHornLensConfig())[positive_indices]
    group_speed_m_s = effective_group_speed_m_s(
        old_extra_path_mm[delay_reliable],
        measured_delay_s[delay_reliable],
    )
    correction = corrected_monotonic_paths_mm(
        old_extra_path_mm,
        measured_delay_s,
        target_delay_s,
        group_speed_m_s,
        reliable=delay_reliable,
    )
    corrected_extra_path_mm = correction.corrected_mm
    predicted_corrected_delay_s = target_delay_s .+
        (corrected_extra_path_mm .- correction.raw_mm) .* 1e-3 ./ group_speed_m_s

    old_axial_length_mm = Float64(design["common_axial_length_mm"])
    required_axial_length_mm = minimum_common_axial_length_mm(maximum(corrected_extra_path_mm))
    corrected_axial_length_mm = max(old_axial_length_mm, required_axial_length_mm)
    corrected_amplitude_mm = [
        smooth_amplitude_mm(corrected_axial_length_mm, extra)
        for extra in corrected_extra_path_mm
    ]
    corrected_inner_radius_mm = [
        minimum_inner_radius_mm(corrected_axial_length_mm, amplitude, 1.6)
        for amplitude in corrected_amplitude_mm
    ]
    monotonic_gate = all(diff(corrected_extra_path_mm) .<= 1e-10) &&
                     corrected_extra_path_mm[end] == 0.0
    curvature_gate = minimum(corrected_inner_radius_mm) >= 3.93 - 1e-8
    monotonic_gate || error("corrected path is not monotonic")
    curvature_gate || error("corrected path violates the inner-radius gate")

    old_error_us = (measured_delay_s .- target_delay_s) .* 1e6
    corrected_error_us = (predicted_corrected_delay_s .- target_delay_s) .* 1e6
    rms_error_before_us = sqrt(mean(abs2, old_error_us[delay_reliable]))
    rms_error_after_us = sqrt(mean(abs2, corrected_error_us))
    carrier_index = argmin(abs.(frequencies_hz .- 242.0e3))
    carrier_transfer = device_transfer[carrier_index, :]
    carrier_frequency_hz = frequencies_hz[carrier_index]
    pitch_mm = 4.8
    throat_width_mm = 1.6
    receiver_speed_m_s = MonotonicHornLensConfig().receiver_speed_m_s
    carrier_pressure_wavelength_mm = receiver_speed_m_s / carrier_frequency_hz * 1e3

    channels_path = joinpath(OUTPUT_ROOT, "group_delay_path_correction.csv")
    open(channels_path, "w") do io
        println(io, "source_group,center_y_mm,delay_reliable,target_delay_us,measured_delay_us,delay_error_us,old_extra_path_mm,raw_corrected_extra_path_mm,corrected_extra_path_mm,predicted_corrected_delay_us,predicted_corrected_error_us,corrected_bend_amplitude_mm,minimum_inner_radius_mm,minimum_anchor_amplitude_ratio,carrier_amplitude_ratio,carrier_phase_deg,center_y_over_pitch,target_delay_cycles,measured_delay_cycles,old_extra_path_over_lambda_p,corrected_extra_path_over_lambda_p,bend_amplitude_over_pitch,inner_radius_over_throat_width")
        for index in 1:channel_count
            println(io, join((
                index,
                centers_mm[index],
                delay_reliable[index],
                target_delay_s[index] * 1e6,
                measured_delay_s[index] * 1e6,
                old_error_us[index],
                old_extra_path_mm[index],
                correction.raw_mm[index],
                corrected_extra_path_mm[index],
                predicted_corrected_delay_s[index] * 1e6,
                corrected_error_us[index],
                corrected_amplitude_mm[index],
                corrected_inner_radius_mm[index],
                minimum_anchor_amplitude_ratio[index],
                abs(carrier_transfer[index]),
                rad2deg(angle(carrier_transfer[index])),
                centers_mm[index] / pitch_mm,
                carrier_frequency_hz * target_delay_s[index],
                carrier_frequency_hz * measured_delay_s[index],
                old_extra_path_mm[index] / carrier_pressure_wavelength_mm,
                corrected_extra_path_mm[index] / carrier_pressure_wavelength_mm,
                corrected_amplitude_mm[index] / pitch_mm,
                corrected_inner_radius_mm[index] / throat_width_mm,
            ), ','))
        end
    end

    anchors_path = joinpath(OUTPUT_ROOT, "frequency_anchor_transfer.csv")
    open(anchors_path, "w") do io
        println(io, "frequency_khz,source_group,center_y_mm,amplitude_ratio,wrapped_phase_deg,unwrapped_phase_deg,phase_fit_residual_deg,frequency_over_f0,center_y_over_pitch")
        for channel in 1:channel_count
            fit = delay_fit[channel]
            for frequency_index in eachindex(frequencies_hz)
                println(io, join((
                    frequencies_hz[frequency_index] / 1e3,
                    channel,
                    centers_mm[channel],
                    abs(device_transfer[frequency_index, channel]),
                    rad2deg(angle(device_transfer[frequency_index, channel])),
                    rad2deg(fit.phase_rad[frequency_index]),
                    rad2deg(fit.phase_rad[frequency_index] - fit.fitted_phase_rad[frequency_index]),
                    frequencies_hz[frequency_index] / carrier_frequency_hz,
                    centers_mm[channel] / pitch_mm,
                ), ','))
            end
        end
    end

    summary_path = joinpath(OUTPUT_ROOT, "group_delay_summary.csv")
    open(summary_path, "w") do io
        println(io, "lower_frequency_khz,carrier_frequency_khz,upper_frequency_khz,reliable_channel_count,maximum_target_delay_us,maximum_reliable_measured_delay_us,effective_group_speed_m_s,reliable_rms_delay_error_before_us,predicted_rms_delay_error_after_us,maximum_phase_fit_residual_deg,maximum_reliable_phase_fit_residual_deg,old_maximum_extra_path_mm,corrected_maximum_extra_path_mm,old_common_axial_length_mm,corrected_common_axial_length_mm,minimum_inner_radius_mm,minimum_anchor_amplitude_ratio,minimum_carrier_amplitude_ratio,monotonic_gate,curvature_gate,lower_frequency_over_f0,upper_frequency_over_f0,maximum_target_delay_cycles,effective_group_speed_over_receiver_speed,old_maximum_extra_path_over_lambda_p,corrected_maximum_extra_path_over_lambda_p,corrected_axial_length_over_lambda_p,minimum_inner_radius_over_throat_width")
        println(io, join((
            frequencies_hz[1] / 1e3,
            frequencies_hz[carrier_index] / 1e3,
            frequencies_hz[end] / 1e3,
            count(delay_reliable),
            maximum(target_delay_s) * 1e6,
            maximum(measured_delay_s[delay_reliable]) * 1e6,
            group_speed_m_s,
            rms_error_before_us,
            rms_error_after_us,
            maximum_phase_fit_residual_deg,
            maximum(getproperty.(delay_fit, :maximum_residual_deg)[delay_reliable]),
            maximum(old_extra_path_mm),
            maximum(corrected_extra_path_mm),
            old_axial_length_mm,
            corrected_axial_length_mm,
            minimum(corrected_inner_radius_mm),
            minimum(minimum_anchor_amplitude_ratio),
            minimum(abs.(carrier_transfer)),
            monotonic_gate,
            curvature_gate,
            frequencies_hz[1] / carrier_frequency_hz,
            frequencies_hz[end] / carrier_frequency_hz,
            carrier_frequency_hz * maximum(target_delay_s),
            group_speed_m_s / receiver_speed_m_s,
            maximum(old_extra_path_mm) / carrier_pressure_wavelength_mm,
            maximum(corrected_extra_path_mm) / carrier_pressure_wavelength_mm,
            corrected_axial_length_mm / carrier_pressure_wavelength_mm,
            minimum(corrected_inner_radius_mm) / throat_width_mm,
        ), ','))
    end

    design_path = joinpath(OUTPUT_ROOT, "corrected_monotonic_lens_v1.jld2")
    full_centers_mm = vcat(-reverse(centers_mm[2:end]), centers_mm)
    full_target_delay_s = reconstruct_symmetric(target_delay_s)
    full_measured_delay_s = reconstruct_symmetric(measured_delay_s)
    full_predicted_corrected_delay_s = reconstruct_symmetric(predicted_corrected_delay_s)
    full_old_extra_path_mm = reconstruct_symmetric(old_extra_path_mm)
    full_raw_corrected_extra_path_mm = reconstruct_symmetric(correction.raw_mm)
    full_corrected_extra_path_mm = reconstruct_symmetric(corrected_extra_path_mm)
    full_corrected_amplitude_mm = reconstruct_symmetric(corrected_amplitude_mm)
    full_corrected_inner_radius_mm = reconstruct_symmetric(corrected_inner_radius_mm)
    full_delay_reliable = reconstruct_symmetric(delay_reliable)
    full_minimum_anchor_amplitude_ratio = reconstruct_symmetric(
        minimum_anchor_amplitude_ratio,
    )
    jldsave(
        design_path;
        format_version=1,
        frequencies_hz,
        centers_mm=full_centers_mm,
        target_delay_s=full_target_delay_s,
        measured_delay_s=full_measured_delay_s,
        predicted_corrected_delay_s=full_predicted_corrected_delay_s,
        old_extra_path_mm=full_old_extra_path_mm,
        raw_corrected_extra_path_mm=full_raw_corrected_extra_path_mm,
        corrected_extra_path_mm=full_corrected_extra_path_mm,
        corrected_path_length_mm=corrected_axial_length_mm .+ full_corrected_extra_path_mm,
        corrected_bend_amplitude_mm=full_corrected_amplitude_mm,
        corrected_inner_radius_mm=full_corrected_inner_radius_mm,
        old_common_axial_length_mm=old_axial_length_mm,
        corrected_common_axial_length_mm=corrected_axial_length_mm,
        effective_group_speed_m_s=group_speed_m_s,
        rms_delay_error_before_us=rms_error_before_us,
        rms_delay_error_after_us=rms_error_after_us,
        maximum_phase_fit_residual_deg,
        delay_reliable=full_delay_reliable,
        minimum_anchor_amplitude_ratio=full_minimum_anchor_amplitude_ratio,
        device_transfer,
        monotonic_gate,
        curvature_gate,
        # Compatibility aliases for the existing monotonic-array mesh drivers.
        common_axial_length_mm=corrected_axial_length_mm,
        extra_path_mm=full_corrected_extra_path_mm,
        path_length_mm=corrected_axial_length_mm .+ full_corrected_extra_path_mm,
        bend_amplitude_mm=full_corrected_amplitude_mm,
        inner_radius_mm=full_corrected_inner_radius_mm,
    )

    geometry_path = joinpath(OUTPUT_ROOT, "corrected_monotonic_lens_v1_geometry.csv")
    open(geometry_path, "w") do io
        println(io, "element_index,center_y_mm,target_delay_us,extra_path_mm,path_length_mm,bend_amplitude_mm,minimum_inner_radius_mm,delay_reliable,center_y_over_pitch,target_delay_cycles,extra_path_over_lambda_p,path_length_over_lambda_p,bend_amplitude_over_pitch,inner_radius_over_throat_width")
        for index in eachindex(full_centers_mm)
            println(io, join((
                index,
                full_centers_mm[index],
                full_target_delay_s[index] * 1e6,
                full_corrected_extra_path_mm[index],
                corrected_axial_length_mm + full_corrected_extra_path_mm[index],
                full_corrected_amplitude_mm[index],
                full_corrected_inner_radius_mm[index],
                full_delay_reliable[index],
                full_centers_mm[index] / pitch_mm,
                carrier_frequency_hz * full_target_delay_s[index],
                full_corrected_extra_path_mm[index] / carrier_pressure_wavelength_mm,
                (corrected_axial_length_mm + full_corrected_extra_path_mm[index]) /
                carrier_pressure_wavelength_mm,
                full_corrected_amplitude_mm[index] / pitch_mm,
                full_corrected_inner_radius_mm[index] / throat_width_mm,
            ), ','))
        end
    end

    delay_panel = plot(
        centers_mm,
        target_delay_s .* 1e6;
        marker=:circle,
        linewidth=2.5,
        xlabel="channel centre y, mm",
        ylabel="relative group delay, us",
        title="Measured delay and one correction",
        label="target",
        gridalpha=0.25,
    )
    plot!(
        delay_panel,
        centers_mm[delay_reliable],
        measured_delay_s[delay_reliable] .* 1e6;
        marker=:diamond,
        label="measured, reliable",
    )
    scatter!(
        delay_panel,
        centers_mm[.!delay_reliable],
        measured_delay_s[.!delay_reliable] .* 1e6;
        marker=:x,
        markersize=8,
        markerstrokewidth=3,
        label="rejected by phase-linearity gate",
    )
    plot!(
        delay_panel,
        centers_mm,
        predicted_corrected_delay_s .* 1e6;
        marker=:square,
        label="predicted corrected",
    )

    path_panel = plot(
        centers_mm,
        old_extra_path_mm;
        marker=:circle,
        linewidth=2.5,
        xlabel="channel centre y, mm",
        ylabel="extra path, mm",
        title="Monotonic unwrapped path",
        label="v0",
        gridalpha=0.25,
    )
    plot!(path_panel, centers_mm, corrected_extra_path_mm; marker=:square, label="corrected v1")

    phase_panel = plot(
        xlabel="frequency, kHz",
        ylabel="unwrapped device phase, deg",
        title="Lens / straight diagonal phase",
        gridalpha=0.25,
    )
    for channel in 1:channel_count
        plot!(
            phase_panel,
            frequencies_hz ./ 1e3,
            rad2deg.(delay_fit[channel].phase_rad);
            marker=:circle,
            linewidth=1.7,
            label="y=$(round(centers_mm[channel]; digits=1)) mm",
        )
    end

    amplitude_panel = plot(
        xlabel="channel centre y, mm",
        ylabel="|lens / straight|",
        title="Diagonal displacement ratio",
        gridalpha=0.25,
    )
    for frequency_index in eachindex(frequencies_hz)
        plot!(
            amplitude_panel,
            centers_mm,
            abs.(device_transfer[frequency_index, :]);
            marker=:circle,
            linewidth=1.9,
            label="$(frequencies_hz[frequency_index] / 1e3) kHz",
        )
    end

    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_group_delay.png")
    savefig(plot(
        delay_panel,
        path_panel,
        phase_panel,
        amplitude_panel;
        layout=(2, 2),
        size=(1500, 1050),
        margin=5Plots.mm,
    ), figure_path)

    geometry_figure = plot(
        xlabel="axial x, mm",
        ylabel="aperture y, mm",
        zlabel="bend z, mm",
        title="Corrected monotonic lens v1 — guide centrelines",
        camera=(35, 28),
        legend=false,
        gridalpha=0.2,
        size=(1400, 900),
    )
    axial_x_mm = range(0.0, corrected_axial_length_mm; length=401)
    for index in eachindex(full_centers_mm)
        plot!(
            geometry_figure,
            axial_x_mm,
            fill(full_centers_mm[index], length(axial_x_mm)),
            smooth_centerline_mm.(
                axial_x_mm,
                corrected_axial_length_mm,
                full_corrected_amplitude_mm[index],
            );
            linewidth=2.4,
            line_z=fill(full_corrected_extra_path_mm[index], length(axial_x_mm)),
            colorbar=false,
        )
    end
    geometry_figure_path = joinpath(
        OUTPUT_ROOT,
        "corrected_monotonic_lens_v1_geometry.png",
    )
    savefig(geometry_figure, geometry_figure_path)

    println("[+] maximum measured delay=$(maximum(measured_delay_s) * 1e6) us")
    println("[+] effective guide group speed=$group_speed_m_s m/s")
    println("[+] delay RMS error: $rms_error_before_us -> $rms_error_after_us us")
    println("[+] extra path max: $(maximum(old_extra_path_mm)) -> $(maximum(corrected_extra_path_mm)) mm")
    println("[+] common axial length: $old_axial_length_mm -> $corrected_axial_length_mm mm")
    println("[+] minimum inner radius=$(minimum(corrected_inner_radius_mm)) mm")
    println("[+] $summary_path")
    println("[+] $figure_path")
    println("[+] $geometry_figure_path")
    (;
        summary_path,
        channels_path,
        anchors_path,
        design_path,
        geometry_path,
        figure_path,
        geometry_figure_path,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_group_delay()
end

end
