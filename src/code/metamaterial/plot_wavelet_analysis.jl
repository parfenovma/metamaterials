ENV["GKSwstype"] = "100"

using JLD2
using Plots

include(joinpath(@__DIR__, "wavelet_analysis.jl"))
using .WaveletAnalysis

function normalized_power_db(result::CWTResult, normalization)
    10.0 .* log10.(wavelet_power(result) ./ normalization .+ eps(Float64))
end

function energy_masked_ridge_us(result::CWTResult, ridge_s, energy_floor_relative)
    power = wavelet_power(result)
    energy = [
        sum(view(power, frequency_index, :) .* view(
            result.cone_of_influence_valid,
            frequency_index,
            :,
        ))
        for frequency_index in eachindex(result.frequency_hz)
    ]
    relative_energy = maximum(energy) > 0.0 ? energy ./ maximum(energy) : zero.(energy)
    ridge_us = copy(ridge_s) .* 1.0e6
    ridge_us[relative_energy .< energy_floor_relative] .= NaN
    ridge_us
end

function plot_wavelet_analysis(
    result_path::AbstractString;
    output_dir::AbstractString=joinpath(dirname(result_path), "figures"),
    dynamic_range_db::Real=50.0,
)
    saved = JLD2.load(result_path)
    reference = saved["reference_cwt"]
    sample = saved["sample_cwt"]
    comparison = saved["comparison"]
    config = saved["morlet_config"]
    metadata = saved["sample_metadata"]
    carrier_hz = metadata.carrier_frequency_hz
    time_us = reference.time_s .* 1.0e6
    frequency_khz = reference.frequency_hz ./ 1.0e3
    reference_power = wavelet_power(reference)
    normalization = maximum(reference_power)
    reference_db = normalized_power_db(reference, normalization)
    sample_db = normalized_power_db(sample, normalization)
    color_limits = (-Float64(dynamic_range_db), 0.0)
    interior_frequency_ticks = collect(
        ceil(first(frequency_khz) / 200.0) * 200.0:200.0:
        floor(last(frequency_khz) / 200.0) * 200.0,
    )
    frequency_ticks_khz = unique(sort(vcat(
        first(frequency_khz),
        interior_frequency_ticks,
        last(frequency_khz),
    )))

    omega0 = Float64(config["omega0"])
    coi_factor = Float64(config["cone_of_influence_factor"])
    energy_floor = Float64(config["reference_energy_floor_relative"])
    scales_s = morlet_scale_s.(reference.frequency_hz, omega0)
    left_coi_us = (first(reference.time_s) .+ coi_factor .* scales_s) .* 1.0e6
    right_coi_us = (last(reference.time_s) .- coi_factor .* scales_s) .* 1.0e6

    common_heatmap = (
        xlabel="Frequency, kHz",
        ylabel="Time, μs",
        color=:viridis,
        clims=color_limits,
        xlims=(first(frequency_khz), last(frequency_khz)),
        ylims=(first(time_us), last(time_us)),
        xticks=frequency_ticks_khz,
        framestyle=:box,
        grid=false,
        guidefontsize=10,
        tickfontsize=8,
        titlefontsize=11,
        left_margin=8Plots.mm,
        right_margin=5Plots.mm,
        bottom_margin=6Plots.mm,
        legend=:topright,
    )

    reference_ridge_us = energy_masked_ridge_us(
        reference,
        comparison.reference_peak_s,
        energy_floor,
    )
    sample_ridge_us = energy_masked_ridge_us(
        sample,
        comparison.sample_peak_s,
        energy_floor,
    )

    reference_plot = heatmap(
        frequency_khz,
        time_us,
        permutedims(reference_db);
        common_heatmap...,
        title="Rectangular reference, dB",
        colorbar_title="dB",
    )
    plot!(
        reference_plot,
        frequency_khz,
        reference_ridge_us;
        label="peak ridge",
        color=:white,
        linewidth=1.8,
    )
    plot!(reference_plot, frequency_khz, left_coi_us; label="COI", color=:white, linestyle=:dash)
    plot!(reference_plot, frequency_khz, right_coi_us; label="", color=:white, linestyle=:dash)
    vline!(reference_plot, [carrier_hz / 1.0e3]; label="carrier", color=:white, linestyle=:dot)

    sample_plot = heatmap(
        frequency_khz,
        time_us,
        permutedims(sample_db);
        common_heatmap...,
        title="Sample $(metadata.profile_type), dB relative to reference maximum",
        colorbar_title="dB",
    )
    plot!(
        sample_plot,
        frequency_khz,
        sample_ridge_us;
        label="peak ridge",
        color=:white,
        linewidth=1.8,
    )
    plot!(sample_plot, frequency_khz, left_coi_us; label="COI", color=:white, linestyle=:dash)
    plot!(sample_plot, frequency_khz, right_coi_us; label="", color=:white, linestyle=:dash)
    vline!(sample_plot, [carrier_hz / 1.0e3]; label="carrier", color=:white, linestyle=:dot)

    valid = comparison.valid
    centroid_delay_us = copy(comparison.centroid_delay_s) .* 1.0e6
    correlation_delay_us = copy(comparison.correlation_delay_s) .* 1.0e6
    centroid_delay_us[.!valid] .= NaN
    correlation_delay_us[.!valid] .= NaN
    delay_plot = plot(
        frequency_khz,
        centroid_delay_us;
        label="energy centroid",
        marker=:circle,
        linewidth=2,
        xlabel="Frequency, kHz",
        ylabel="Additional delay, μs",
        title="Wavelet delays; carrier=$(carrier_hz / 1e3) kHz",
        xticks=frequency_ticks_khz,
        framestyle=:box,
        gridalpha=0.25,
        guidefontsize=10,
        tickfontsize=8,
        titlefontsize=11,
        left_margin=8Plots.mm,
        right_margin=5Plots.mm,
        bottom_margin=6Plots.mm,
    )
    plot!(
        delay_plot,
        frequency_khz,
        correlation_delay_us;
        label="envelope correlation",
        marker=:square,
        linewidth=2,
    )
    hline!(delay_plot, [0.0]; label="", color=:gray50, linestyle=:dash)
    invalid = findall(.!comparison.valid)
    isempty(invalid) || scatter!(
        delay_plot,
        frequency_khz[invalid],
        comparison.centroid_delay_s[invalid] .* 1.0e6;
        label="low confidence",
        marker=:xcross,
        markersize=3,
        markeralpha=0.45,
        color=:red,
    )

    figure = plot(
        reference_plot,
        sample_plot,
        delay_plot;
        layout=(3, 1),
        size=(1200, 1500),
    )
    mkpath(output_dir)
    output_path = joinpath(output_dir, "scalogram_$(splitext(basename(result_path))[1]).png")
    savefig(figure, output_path)
    output_path
end

function main(args=ARGS)
    isempty(args) && error("usage: julia plot_wavelet_analysis.jl RESULT.jld2 [OUTPUT_DIR]")
    output_dir = length(args) >= 2 ? args[2] : joinpath(dirname(args[1]), "figures")
    println("[+] Wavelet figure: ", plot_wavelet_analysis(args[1]; output_dir))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
