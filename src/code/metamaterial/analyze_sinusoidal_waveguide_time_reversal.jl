module AnalyzeSinusoidalWaveguideTimeReversal

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Statistics

include(joinpath(@__DIR__, "operator_time_reversal_audit.jl"))
using .OperatorTimeReversalAudit

export analyze_sinusoidal_waveguide_time_reversal

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PACKAGE_ROOT = joinpath(
    PROJECT_ROOT, "results", "sinusoidal_material_transfer_p6_242khz",
)
const INPUT_PATH = joinpath(
    PACKAGE_ROOT, "response_matrix", "base",
    "aluminium_strip_group_response_matrix.jld2",
)
const OUTPUT_ROOT = joinpath(PACKAGE_ROOT, "time_reversal_carrier")
const FREQUENCY_HZ = 242.0e3
const PRESSURE_SPEED_M_S = 6122.102437409232
const PRESSURE_WAVELENGTH_MM = 1.0e3 * PRESSURE_SPEED_M_S / FREQUENCY_HZ
const SOURCE_STRIP_COUNTS = [1; fill(2, 7)]
const DYNAMIC_RANGES_DB = (20.0, 30.0, 40.0, 60.0)

function csv_value(value)
    value isa AbstractFloat && return isfinite(value) ? string(value) :
        value == Inf ? "Inf" : value == -Inf ? "-Inf" : "NaN"
    value isa Symbol && return String(value)
    value isa Bool && return value ? "true" : "false"
    string(value)
end

function write_csv(path, rows)
    isempty(rows) && return
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function load_matrices()
    isfile(INPUT_PATH) || error("missing frozen response matrix: $INPUT_PATH")
    data = JLD2.load(INPUT_PATH)
    data["format_version"] == 1 || error("unsupported response-matrix format")
    centers_mm = Float64.(data["centers_mm"])
    matrices = Dict(
        :sinusoidal => ComplexF64.(data["selected"][:preout]),
        :straight => ComplexF64.(data["uniform"][:preout]),
    )
    all(size(matrix) == (8, 8) for matrix in values(matrices)) ||
        throw(DimensionMismatch("expected two frozen 8 x 8 matrices"))
    length(centers_mm) == 8 || throw(DimensionMismatch(
        "expected eight saved receiver coordinates",
    ))
    matrices, centers_mm
end

function normalized_matrices(matrices, normalization)
    if normalization == :as_solved
        Dict(role => copy(matrix) for (role, matrix) in matrices)
    elseif normalization == :nominal_area
        scale = sqrt.(SOURCE_STRIP_COUNTS)'
        Dict(role => matrix ./ scale for (role, matrix) in matrices)
    else
        throw(ArgumentError("unknown source normalization: $normalization"))
    end
end

function coherence_metrics(kernel)
    coherence = coherence_kernel(kernel)
    off_diagonal = [
        coherence[row, column]
        for row in axes(coherence, 1), column in axes(coherence, 2)
        if row != column
    ]
    (
        mean_offdiagonal_coherence=mean(off_diagonal),
        rms_offdiagonal_coherence=sqrt(mean(abs2, off_diagonal)),
        maximum_offdiagonal_coherence=maximum(off_diagonal),
    )
end

function analyze_case(matrix, coordinates_lambda)
    singular = singular_summary(matrix, DYNAMIC_RANGES_DB)
    kernel = time_reversal_kernel(matrix)
    coherence = coherence_metrics(kernel)
    spacing = only(unique(round.(diff(coordinates_lambda); digits=12)))
    focuses = [
        focus_metrics(kernel, focal_bin; bin_spacing=spacing)
        for focal_bin in eachindex(coordinates_lambda)
    ]
    (; matrix, singular, kernel, coherence, focuses)
end

function summary_row(geometry, normalization, result)
    focuses = result.focuses
    singular = result.singular
    (
        geometry,
        source_normalization=normalization,
        stable_rank=singular.stable_rank,
        entropy_rank=singular.entropy_rank,
        condition_number=cond(result.matrix),
        modes_20db=singular.observable_modes[20.0],
        modes_30db=singular.observable_modes[30.0],
        modes_40db=singular.observable_modes[40.0],
        modes_60db=singular.observable_modes[60.0],
        median_sampled_fwhm_lambda=median(getproperty.(focuses, :fwhm)),
        median_focus_contrast_db=median(getproperty.(focuses, :focus_contrast_db)),
        median_peak_outside_fwhm_db=median(getproperty.(focuses, :peak_outside_fwhm_db)),
        all_reduced_targets_are_global_maxima=all(
            getproperty.(focuses, :target_is_global_maximum),
        ),
        mean_offdiagonal_coherence=result.coherence.mean_offdiagonal_coherence,
        rms_offdiagonal_coherence=result.coherence.rms_offdiagonal_coherence,
        maximum_offdiagonal_coherence=result.coherence.maximum_offdiagonal_coherence,
    )
end

function relative_rows(summary_rows)
    rows = NamedTuple[]
    for normalization in (:as_solved, :nominal_area)
        lens = only(filter(row ->
            row.geometry == :sinusoidal && row.source_normalization == normalization,
            summary_rows,
        ))
        straight = only(filter(row ->
            row.geometry == :straight && row.source_normalization == normalization,
            summary_rows,
        ))
        push!(rows, (
            source_normalization=normalization,
            stable_rank_ratio=lens.stable_rank / straight.stable_rank,
            entropy_rank_ratio=lens.entropy_rank / straight.entropy_rank,
            condition_number_ratio=lens.condition_number / straight.condition_number,
            focus_contrast_change_db=
                lens.median_focus_contrast_db - straight.median_focus_contrast_db,
            peak_outside_fwhm_change_db=
                lens.median_peak_outside_fwhm_db - straight.median_peak_outside_fwhm_db,
            mean_offdiagonal_coherence_change=
                lens.mean_offdiagonal_coherence - straight.mean_offdiagonal_coherence,
        ))
    end
    rows
end

function plot_singular_values(results)
    figure = plot(
        xlabel="even-sector mode index", ylabel="relative singular value / dB",
        title="242 kHz: sinusoidal versus straight waveguides",
        size=(900, 560), legend=:bottomleft, ylim=(-22, 1), xticks=1:8,
    )
    for (geometry, label, color) in (
        (:sinusoidal, "sinusoidal", :darkorange),
        (:straight, "straight", :steelblue),
    )
        values = results[(:nominal_area, geometry)].singular.singular_value_db
        plot!(figure, eachindex(values), values;
            marker=:circle, linewidth=2.5, label, color)
    end
    hline!(figure, [-20.0]; color=:black, linestyle=:dash,
        label="20 dB dynamic range")
    savefig(figure, joinpath(OUTPUT_ROOT, "singular_spectrum.png"))
end

function plot_coherence(results, coordinates_mm)
    panels = Any[]
    for (geometry, title) in ((:sinusoidal, "sinusoidal"), (:straight, "straight"))
        coherence = coherence_kernel(results[(:nominal_area, geometry)].kernel)
        push!(panels, heatmap(
            coordinates_mm, coordinates_mm, coherence;
            xlabel="target |y| / mm", ylabel="observation |y| / mm",
            title, aspect_ratio=:equal, color=:viridis, clims=(0, 1),
            colorbar_title="|C|",
        ))
    end
    savefig(plot(
        panels...; layout=(1, 2), size=(1150, 500),
        left_margin=6Plots.mm, bottom_margin=5Plots.mm,
    ),
        joinpath(OUTPUT_ROOT, "time_reversal_coherence.png"))
end

function plot_focus_profiles(results, coordinates_mm)
    target_bins = (1, 5)
    panels = Any[]
    for target_bin in target_bins
        panel = plot(
            xlabel="observation |y| / mm", ylabel="normalized intensity / dB",
            title="target |y| = $(coordinates_mm[target_bin]) mm",
            legend=:bottomleft, ylim=(-45, 2),
        )
        for (geometry, label, color) in (
            (:sinusoidal, "sinusoidal", :darkorange),
            (:straight, "straight", :steelblue),
        )
            profile = normalized_focus_intensity(
                results[(:nominal_area, geometry)].kernel, target_bin,
            )
            plot!(panel, coordinates_mm, 10 .* log10.(max.(profile, 10.0^(-4.5)));
                marker=:circle, linewidth=2.5, label, color)
        end
        vline!(panel, [coordinates_mm[target_bin]];
            color=:black, linestyle=:dot, label="target")
        push!(panels, panel)
    end
    savefig(plot(
        panels...; layout=(1, 2), size=(1150, 480),
        left_margin=6Plots.mm, bottom_margin=5Plots.mm,
    ),
        joinpath(OUTPUT_ROOT, "time_reversal_profiles.png"))
end

function write_readme(summary_rows, comparison_rows, coordinates_mm)
    summaries = Dict(
        (row.source_normalization, row.geometry) => row for row in summary_rows
    )
    nominal_lens = summaries[(:nominal_area, :sinusoidal)]
    nominal_straight = summaries[(:nominal_area, :straight)]
    nominal_comparison = only(filter(row ->
        row.source_normalization == :nominal_area, comparison_rows,
    ))
    robust = all(
        row.stable_rank_ratio < 1 && row.entropy_rank_ratio < 1 &&
            row.condition_number_ratio > 1
        for row in comparison_rows
    )
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# Sinusoidal aluminium waveguides: carrier time-reversal audit")
        println(io)
        println(io, "## Result")
        println(io)
        println(io, "At `242 kHz`, the selected sinusoidal lens does **not** add spatial channels relative to the matched straight strips. Both matrices retain all eight reflection-even modes within 20 dB, while the sinusoidal spectrum is less balanced.")
        println(io)
        println(io, "Under the nominal source-area normalization:")
        println(io)
        println(io, "| geometry | stable rank | entropy rank | condition number | modes @20 dB | median contrast / dB | median sidelobe outside FWHM / dB |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|")
        for geometry in (:sinusoidal, :straight)
            row = summaries[(:nominal_area, geometry)]
            println(io, "| $(geometry) | $(round(row.stable_rank; digits=3)) | $(round(row.entropy_rank; digits=3)) | $(round(row.condition_number; digits=3)) | $(row.modes_20db) | $(round(row.median_focus_contrast_db; digits=2)) | $(round(row.median_peak_outside_fwhm_db; digits=2)) |")
        end
        println(io)
        println(io, "The sinusoidal stable rank is $(round(100nominal_comparison.stable_rank_ratio; digits=1))% of the straight value, its condition number is `$(round(nominal_comparison.condition_number_ratio; digits=2))x` larger, and its median time-reversal contrast changes by $(round(nominal_comparison.focus_contrast_change_db; digits=2)) dB. The sign of all three comparisons is $(robust ? "the same" : "not the same") under the original as-solved traction convention.")
        println(io)
        println(io, "This does not contradict the existing `1.633x` carrier focus gain under the fixed uniform drive. That design optimized one coherent focal combination; it did not optimize the full singular spectrum or general-purpose channel conditioning. In book terminology, it redirects energy into one useful wavefront but does not increase the carrier information capacity of this even-sector basis.")
        println(io)
        println(io, "## What the kernel means")
        println(io)
        println(io, "`K = H H^dagger` is computed on the plane five millimetres after the strips. Its columns are virtual phase-conjugate refocusing profiles over `|y|`. Because every source is a mirror-symmetric strip pair, an off-axis target represents two equal physical foci at `+y` and `-y`; signed lateral localization has not been demonstrated.")
        println(io)
        println(io, "The sampled half-aperture spacing is `$(round(coordinates_mm[2] - coordinates_mm[1]; digits=3)) mm = $(round((coordinates_mm[2] - coordinates_mm[1]) / PRESSURE_WAVELENGTH_MM; digits=3)) lambdaP`. Every sampled half-maximum width is one bin, so the present matrix cannot compare sub-bin resolution. Contrast and sidelobe metrics remain resolved.")
        println(io)
        println(io, "## Scope and next physical calculation")
        println(io)
        println(io, "This is a one-frequency result. A genuine broadband time-reversal test needs the same eight source-group responses and eight output samples on the frozen 17-frequency grid, plus complex source work for exact power normalization. Odd source groups are also required before claiming signed lateral localization. Until those data exist, no conclusion about temporal compression, Q, ringing, or broadband super-resolution is made.")
    end
end

function analyze_sinusoidal_waveguide_time_reversal()
    mkpath(OUTPUT_ROOT)
    raw_matrices, coordinates_mm = load_matrices()
    coordinates_lambda = coordinates_mm ./ PRESSURE_WAVELENGTH_MM
    results = Dict{Tuple{Symbol,Symbol},Any}()
    summary_rows = NamedTuple[]
    singular_rows = NamedTuple[]
    focus_rows = NamedTuple[]
    profile_rows = NamedTuple[]

    for normalization in (:as_solved, :nominal_area)
        matrices = normalized_matrices(raw_matrices, normalization)
        for geometry in (:sinusoidal, :straight)
            result = analyze_case(matrices[geometry], coordinates_lambda)
            results[(normalization, geometry)] = result
            push!(summary_rows, summary_row(geometry, normalization, result))
            for mode_index in eachindex(result.singular.singular_values)
                push!(singular_rows, (
                    geometry, source_normalization=normalization, mode_index,
                    singular_value=result.singular.singular_values[mode_index],
                    relative_db=result.singular.singular_value_db[mode_index],
                ))
            end
            for (focal_bin, metrics) in enumerate(result.focuses)
                push!(focus_rows, (
                    geometry, source_normalization=normalization, focal_bin,
                    focal_abs_y_mm=coordinates_mm[focal_bin],
                    peak_bin=metrics.peak_bin,
                    target_is_global_maximum=metrics.target_is_global_maximum,
                    sampled_fwhm_bins=metrics.fwhm_bins,
                    sampled_fwhm_lambda=metrics.fwhm,
                    focus_contrast_db=metrics.focus_contrast_db,
                    peak_outside_fwhm_db=metrics.peak_outside_fwhm_db,
                ))
                for observation_bin in eachindex(metrics.profile)
                    push!(profile_rows, (
                        geometry, source_normalization=normalization, focal_bin,
                        observation_bin,
                        observation_abs_y_mm=coordinates_mm[observation_bin],
                        normalized_intensity=metrics.profile[observation_bin],
                    ))
                end
            end
        end
    end
    comparison_rows = relative_rows(summary_rows)

    write_csv(joinpath(OUTPUT_ROOT, "architecture_summary.csv"), summary_rows)
    write_csv(joinpath(OUTPUT_ROOT, "relative_comparison.csv"), comparison_rows)
    write_csv(joinpath(OUTPUT_ROOT, "singular_values.csv"), singular_rows)
    write_csv(joinpath(OUTPUT_ROOT, "focus_metrics.csv"), focus_rows)
    write_csv(joinpath(OUTPUT_ROOT, "focus_profiles.csv"), profile_rows)
    plot_singular_values(results)
    plot_coherence(results, coordinates_mm)
    plot_focus_profiles(results, coordinates_mm)
    write_readme(summary_rows, comparison_rows, coordinates_mm)
    jldsave(
        joinpath(OUTPUT_ROOT, "sinusoidal_waveguide_time_reversal.jld2");
        frequency_hz=FREQUENCY_HZ,
        pressure_wavelength_mm=PRESSURE_WAVELENGTH_MM,
        coordinates_mm, coordinates_lambda, source_strip_counts=SOURCE_STRIP_COUNTS,
        results, summary_rows, comparison_rows,
    )
    for row in comparison_rows
        println(
            "[+] $(row.source_normalization): stable-rank ratio=" *
            "$(round(row.stable_rank_ratio; digits=3)), condition ratio=" *
            "$(round(row.condition_number_ratio; digits=3)), contrast delta=" *
            "$(round(row.focus_contrast_change_db; digits=3)) dB",
        )
    end
    summary_rows
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    analyze_sinusoidal_waveguide_time_reversal()
end

end
