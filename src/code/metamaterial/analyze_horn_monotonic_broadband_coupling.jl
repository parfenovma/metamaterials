module HornMonotonicBroadbandCoupling

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RESPONSE_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_response_matrix")
const V1_GATE_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_v1_carrier")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_COUPLING_OUTPUT",
    V1_GATE_ROOT,
)
const FREQUENCIES_HZ = [193.8e3, 242.0e3, 290.1e3]

using JLD2
using LinearAlgebra: norm
using Statistics: mean
ENV["GKSwstype"] = "100"
using Plots

if !isdefined(parentmodule(@__MODULE__), :HornMonotonicGroupDelay)
    Base.include(
        parentmodule(@__MODULE__),
        joinpath(@__DIR__, "run_horn_monotonic_group_delay.jl"),
    )
end
using ..HornMonotonicGroupDelay: linear_group_delay_s

export column_crosstalk_ratio, analyze_broadband_coupling

function column_crosstalk_ratio(matrix, column_index)
    diagonal = abs(matrix[column_index, column_index])
    diagonal > 0 || return Inf
    off_diagonal_energy = sum(abs2, matrix[:, column_index]) - diagonal^2
    sqrt(max(off_diagonal_energy, 0.0)) / diagonal
end

frequency_suffix(frequency_hz) = frequency_hz == 242.0e3 ? "" :
    frequency_hz == 193.8e3 ? "_193p8khz" : "_290p1khz"

function load_response(variant, frequency_hz)
    JLD2.load(joinpath(
        RESPONSE_ROOT,
        "response_matrix_$(variant)_phase2$(frequency_suffix(frequency_hz)).jld2",
    ))
end

function analyze_broadband_coupling()
    mkpath(OUTPUT_ROOT)
    lens = load_response.(Ref("lens"), FREQUENCIES_HZ)
    uniform = load_response.(Ref("uniform"), FREQUENCIES_HZ)
    centers_mm = Float64.(lens[2]["nonnegative_channel_centers_mm"])
    channel_count = length(centers_mm)

    diagonal_delay_us = zeros(channel_count)
    focus_delay_us = zeros(channel_count)
    diagonal_residual_deg = zeros(channel_count)
    focus_residual_deg = zeros(channel_count)
    minimum_diagonal_ratio = zeros(channel_count)
    minimum_focus_ratio = zeros(channel_count)
    carrier_lens_crosstalk = zeros(channel_count)
    carrier_uniform_crosstalk = zeros(channel_count)

    for channel in 1:channel_count
        diagonal_transfer = ComplexF64[
            lens[index]["response_matrix_m"][channel, channel] /
            uniform[index]["response_matrix_m"][channel, channel]
            for index in eachindex(FREQUENCIES_HZ)
        ]
        focus_transfer = ComplexF64[
            lens[index]["focus_response_m"][channel] /
            uniform[index]["focus_response_m"][channel]
            for index in eachindex(FREQUENCIES_HZ)
        ]
        diagonal_fit = linear_group_delay_s(FREQUENCIES_HZ, diagonal_transfer)
        focus_fit = linear_group_delay_s(FREQUENCIES_HZ, focus_transfer)
        diagonal_delay_us[channel] = diagonal_fit.delay_s * 1e6
        focus_delay_us[channel] = focus_fit.delay_s * 1e6
        diagonal_residual_deg[channel] = diagonal_fit.maximum_residual_deg
        focus_residual_deg[channel] = focus_fit.maximum_residual_deg
        minimum_diagonal_ratio[channel] = minimum(abs.(diagonal_transfer))
        minimum_focus_ratio[channel] = minimum(abs.(focus_transfer))
        carrier_lens_crosstalk[channel] = column_crosstalk_ratio(
            lens[2]["response_matrix_m"],
            channel,
        )
        carrier_uniform_crosstalk[channel] = column_crosstalk_ratio(
            uniform[2]["response_matrix_m"],
            channel,
        )
    end

    reliable = diagonal_residual_deg .<= 30.0 .&& minimum_diagonal_ratio .>= 0.50
    reliable_focus = focus_residual_deg .<= 35.0 .&& minimum_focus_ratio .>= 0.45
    joint_reliable = reliable .&& reliable_focus
    delay_observable_rms_difference_us = sqrt(mean(abs2,
        diagonal_delay_us[joint_reliable] .- focus_delay_us[joint_reliable],
    ))
    maximum_lens_crosstalk_ratio = maximum(carrier_lens_crosstalk)
    maximum_uniform_crosstalk_ratio = maximum(carrier_uniform_crosstalk)
    crosstalk_dominant = maximum_lens_crosstalk_ratio >= 0.20

    v1_gate = JLD2.load(joinpath(V1_GATE_ROOT, "v1_carrier_gate.jld2"))
    phase_delay_compatible = v1_gate["maximum_phase_error_deg"] <= 10.0 &&
                             v1_gate["focus_change"] > 0
    dispersion_intercept_conflict = !phase_delay_compatible && !crosstalk_dominant

    channels_path = joinpath(OUTPUT_ROOT, "broadband_coupling_channels.csv")
    open(channels_path, "w") do io
        println(io, "source_group,center_y_mm,diagonal_delay_us,focus_delay_us,delay_difference_us,diagonal_phase_residual_deg,focus_phase_residual_deg,minimum_diagonal_amplitude_ratio,minimum_focus_amplitude_ratio,carrier_lens_crosstalk_ratio,carrier_uniform_crosstalk_ratio,diagonal_delay_reliable,focus_delay_reliable,center_y_over_pitch,diagonal_delay_cycles,focus_delay_cycles")
        for channel in 1:channel_count
            println(io, join((
                channel,
                centers_mm[channel],
                diagonal_delay_us[channel],
                focus_delay_us[channel],
                diagonal_delay_us[channel] - focus_delay_us[channel],
                diagonal_residual_deg[channel],
                focus_residual_deg[channel],
                minimum_diagonal_ratio[channel],
                minimum_focus_ratio[channel],
                carrier_lens_crosstalk[channel],
                carrier_uniform_crosstalk[channel],
                reliable[channel],
                reliable_focus[channel],
                centers_mm[channel] / 4.8,
                diagonal_delay_us[channel] * 1e-6 * 242.0e3,
                focus_delay_us[channel] * 1e-6 * 242.0e3,
            ), ','))
        end
    end

    summary_path = joinpath(OUTPUT_ROOT, "broadband_coupling_summary.csv")
    open(summary_path, "w") do io
        println(io, "reliable_diagonal_channel_count,reliable_focus_channel_count,joint_reliable_channel_count,delay_observable_rms_difference_us,maximum_lens_crosstalk_ratio,maximum_uniform_crosstalk_ratio,crosstalk_dominant,v1_focus_change,v1_maximum_phase_error_deg,phase_delay_compatible,dispersion_intercept_conflict")
        println(io, join((
            count(reliable),
            count(reliable_focus),
            count(joint_reliable),
            delay_observable_rms_difference_us,
            maximum_lens_crosstalk_ratio,
            maximum_uniform_crosstalk_ratio,
            crosstalk_dominant,
            v1_gate["focus_change"],
            v1_gate["maximum_phase_error_deg"],
            phase_delay_compatible,
            dispersion_intercept_conflict,
        ), ','))
    end

    delay_panel = plot(
        centers_mm,
        diagonal_delay_us;
        marker=:square,
        linewidth=2.5,
        xlabel="channel centre y, mm",
        ylabel="excess group delay, us",
        title="Two delay observables",
        label="diagonal preout",
        gridalpha=0.25,
    )
    plot!(delay_panel, centers_mm, focus_delay_us;
          marker=:circle, linewidth=2.2, label="source contribution at focus")

    residual_panel = plot(
        centers_mm,
        diagonal_residual_deg;
        marker=:square,
        linewidth=2.5,
        xlabel="channel centre y, mm",
        ylabel="three-anchor residual, deg",
        title="Phase linearity",
        label="diagonal preout",
        gridalpha=0.25,
    )
    plot!(residual_panel, centers_mm, focus_residual_deg;
          marker=:circle, linewidth=2.2, label="focus contribution")
    hline!(residual_panel, [30.0]; linestyle=:dash, color=:black, label="30 deg gate")

    crosstalk_panel = plot(
        centers_mm,
        carrier_lens_crosstalk;
        marker=:square,
        linewidth=2.5,
        xlabel="channel centre y, mm",
        ylabel="off-diagonal L2 / diagonal",
        title="Carrier response-matrix cross-talk",
        label="lens",
        gridalpha=0.25,
    )
    plot!(crosstalk_panel, centers_mm, carrier_uniform_crosstalk;
          marker=:circle, linewidth=2.2, label="matched straight")
    hline!(crosstalk_panel, [0.20]; linestyle=:dash, color=:black, label="dominance gate")

    figure_path = joinpath(OUTPUT_ROOT, "horn_monotonic_broadband_coupling.png")
    savefig(plot(
        delay_panel,
        residual_panel,
        crosstalk_panel;
        layout=(1, 3),
        size=(1800, 600),
        margin=5Plots.mm,
    ), figure_path)

    println("[+] max lens carrier cross-talk=$maximum_lens_crosstalk_ratio")
    println("[+] diagonal/focus delay RMS difference=$delay_observable_rms_difference_us us")
    println("[+] cross-talk dominant=$crosstalk_dominant")
    println("[+] phase-delay compatible=$phase_delay_compatible")
    println("[+] dispersion/intercept conflict=$dispersion_intercept_conflict")
    println("[+] $summary_path")
    (; summary_path, channels_path, figure_path, dispersion_intercept_conflict)
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_broadband_coupling()
end

end
