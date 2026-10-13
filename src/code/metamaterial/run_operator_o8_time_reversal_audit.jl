module RunOperatorO8TimeReversalAudit

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Statistics

include(joinpath(@__DIR__, "operator_time_reversal_audit.jl"))
using .OperatorTimeReversalAudit

export run_o8_time_reversal_audit

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const INPUT_ROOT = joinpath(PROJECT_ROOT, "results", "operator_space", "o8_plane_strain_v0")
const INPUT_BUNDLE = joinpath(INPUT_ROOT, "pzt_physical_bundle.jld2")
const METADATA_BUNDLE = joinpath(INPUT_ROOT, "o8_plane_strain_boundaries.jld2")
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o8_time_reversal_audit_v0",
)
const ARCHITECTURES = ((id=:q00, states=[1]), (id=:q4, states=collect(1:4)))
const STATE_LABELS = ("q00", "q10", "q01", "q11")
const DYNAMIC_RANGES_DB = (20.0, 30.0, 40.0, 60.0)
const APERTURE_WAVELENGTHS = 8.0

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

function load_inputs()
    isfile(INPUT_BUNDLE) || error("missing O8 PZT bundle: $INPUT_BUNDLE")
    isfile(METADATA_BUNDLE) || error("missing O8 metadata bundle: $METADATA_BUNDLE")
    physical = JLD2.load(INPUT_BUNDLE)
    physical["passed"] || error("O8 PZT physical branch did not pass")
    metadata = JLD2.load(METADATA_BUNDLE)
    outputs = ComplexF64.(physical["common_scale"] .* physical["outputs"])
    frequencies_hz = Float64.(metadata["frequencies_hz"])
    frequency_ratios = Float64.(metadata["frequency_ratios"])
    size(outputs) == (17, 3, 4, 7) || throw(DimensionMismatch(
        "expected frozen 17 x 3 x 4 x 7 O8 transfer, got $(size(outputs))",
    ))
    outputs, frequencies_hz, frequency_ratios, physical["common_scale"]
end

function aperture_coordinates(output_count)
    spacing = APERTURE_WAVELENGTHS / output_count
    coordinates = collect(range(
        -APERTURE_WAVELENGTHS / 2 + spacing / 2,
        APERTURE_WAVELENGTHS / 2 - spacing / 2;
        length=output_count,
    ))
    coordinates, spacing
end

function architecture_rows(operators, summaries, focus_rows, novelty)
    rows = NamedTuple[]
    for architecture in ARCHITECTURES
        id = architecture.id
        selected = filter(row -> row.architecture == id, focus_rows)
        summary = summaries[id]
        push!(rows, (
            architecture=id,
            physical_states=length(architecture.states),
            operator_columns=size(operators[id], 2),
            stable_rank=summary.stable_rank,
            entropy_rank=summary.entropy_rank,
            modes_20db=summary.observable_modes[20.0],
            modes_30db=summary.observable_modes[30.0],
            modes_40db=summary.observable_modes[40.0],
            modes_60db=summary.observable_modes[60.0],
            median_fwhm_lambda=median(getproperty.(selected, :fwhm_lambda)),
            median_focus_contrast_db=median(getproperty.(selected, :focus_contrast_db)),
            median_peak_outside_fwhm_db=median(getproperty.(selected, :peak_outside_fwhm_db)),
            all_targets_are_global_maxima=all(getproperty.(selected, :target_is_global_maximum)),
            added_state_novelty_20db=id == :q4 ? novelty[20.0] : NaN,
            added_state_novelty_30db=id == :q4 ? novelty[30.0] : NaN,
            added_state_novelty_40db=id == :q4 ? novelty[40.0] : NaN,
        ))
    end
    rows
end

function plot_singular_spectra(summaries)
    figure = plot(
        xlabel="spatial mode index", ylabel="relative singular value / dB",
        title="O8 PZT: observable spatial modes", size=(900, 560),
        legend=:topright, ylim=(-80, 2), xticks=1:17,
    )
    for (id, label, color) in ((:q00, "q00 only", :steelblue), (:q4, "four states", :darkorange))
        values = summaries[id].singular_value_db
        plot!(figure, eachindex(values), values; marker=:circle, linewidth=2.5,
            label, color)
    end
    for level in (20, 40, 60)
        hline!(figure, [-level]; color=:gray, linestyle=:dash, linewidth=1,
            label=level == 20 ? "dynamic-range levels" : "")
    end
    savefig(figure, joinpath(OUTPUT_ROOT, "singular_spectrum.png"))
end

function plot_center_focus(kernels, coordinates)
    focal_bin = cld(length(coordinates), 2)
    figure = plot(
        xlabel="output coordinate / lambdaP", ylabel="normalized intensity / dB",
        title="O8 virtual time reversal: centre focus", size=(900, 560),
        legend=:bottomleft, ylim=(-60, 3),
    )
    for (id, label, color) in ((:q00, "q00 only", :steelblue), (:q4, "four states", :darkorange))
        profile = normalized_focus_intensity(kernels[id], focal_bin)
        profile_db = 10 .* log10.(max.(profile, 1.0e-6))
        plot!(figure, coordinates, profile_db; marker=:circle, linewidth=2.5,
            label, color)
    end
    hline!(figure, [-3.0103]; color=:black, linestyle=:dash, label="half maximum")
    vline!(figure, [coordinates[focal_bin]]; color=:black, linestyle=:dot,
        label="target")
    savefig(figure, joinpath(OUTPUT_ROOT, "center_time_reversal_focus.png"))
end

function plot_coherence_kernels(kernels, coordinates)
    panels = Any[]
    for (id, title) in ((:q00, "q00 only"), (:q4, "four states"))
        push!(panels, heatmap(
            coordinates, coordinates, coherence_kernel(kernels[id]);
            xlabel="target / lambdaP", ylabel="observation / lambdaP",
            title, aspect_ratio=:equal, color=:viridis, clims=(0, 1),
            colorbar_title="|C|",
        ))
    end
    figure = plot(panels...; layout=(1, 2), size=(1150, 500))
    savefig(figure, joinpath(OUTPUT_ROOT, "time_reversal_coherence.png"))
end

function write_readme(
    architecture_summary, weights, frequencies_hz, novelty, common_scale,
)
    rows = Dict(row.architecture => row for row in architecture_summary)
    q00 = rows[:q00]
    q4 = rows[:q4]
    mode_delta = q4.modes_30db - q00.modes_30db
    width_change = 100 * (q4.median_fwhm_lambda / q00.median_fwhm_lambda - 1)
    interpretation = mode_delta > 0 ?
        "Q4 exposes $mode_delta additional spatial mode(s) at 30 dB." :
        mode_delta == 0 ?
            "Q4 does not increase the 30 dB mode count; any detector gain must come from conditioning or defect-specific coding rather than a larger integer mode count." :
            "Q4 loses $(-mode_delta) observable mode(s) at 30 dB under fixed total energy."
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O8 book-guided time-reversal and modal audit")
        println(io)
        println(io, "This is a post-hoc diagnostic of the frozen, physically passing O8 PZT/RL transfer. It does not tune the head and does not open a new detector split.")
        println(io)
        println(io, "The pulse-weighted operator is")
        println(io)
        println(io, "```math")
        println(io, "A_Q=[w_f H_q(f)]_{q,f}/\\sqrt{Q}, \\qquad K_Q=A_QA_Q^\\dagger.")
        println(io, "```")
        println(io)
        println(io, "The `1/sqrt(Q)` factor keeps total incident energy fixed. The seven weights are the normalized amplitude spectrum of the project-standard five-cycle Hann burst. No state-wise or frequency-wise transfer normalization is used; the inherited common O8 scale is `$(common_scale)`.")
        println(io)
        println(io, "| architecture | stable rank | entropy rank | modes @20 dB | @30 dB | @40 dB | median FWHM / lambdaP | median focus contrast / dB |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
        for row in architecture_summary
            println(io, "| $(row.architecture) | $(round(row.stable_rank; digits=3)) | $(round(row.entropy_rank; digits=3)) | $(row.modes_20db) | $(row.modes_30db) | $(row.modes_40db) | $(round(row.median_fwhm_lambda; digits=3)) | $(round(row.median_focus_contrast_db; digits=2)) |")
        end
        println(io)
        println(io, interpretation)
        println(io, " The Q4 median FWHM change relative to q00 is $(round(width_change; digits=1))%.")
        println(io, " Additional-state energy outside the q00 observable subspace is $(round(100novelty[20.0]; digits=2))% / $(round(100novelty[30.0]; digits=2))% / $(round(100novelty[40.0]; digits=2))% at 20/30/40 dB.")
        println(io)
        println(io, "## Scope boundary")
        println(io)
        println(io, "This kernel tests spatial focusing and mode observability in the already lossy O8 transfer. Seven frequency anchors are insufficient to estimate modal Q, decay time, or ringing. Those require a dense-band solve or a measured impulse response. Dynamic-range thresholds here are sensitivity diagnostics, not a replacement for the existing stochastic O8 detector replay.")
        println(io)
        println(io, "Frequencies (kHz): `$(join(round.(frequencies_hz ./ 1e3; digits=3), ", "))`.")
        println(io, "Pulse weights: `$(join(round.(weights; digits=6), ", "))`.")
    end
end

function run_o8_time_reversal_audit()
    mkpath(OUTPUT_ROOT)
    transfers, frequencies_hz, frequency_ratios, common_scale = load_inputs()
    weights = normalized_pulse_weights(frequency_ratios; cycles=5.0)
    coordinates, spacing = aperture_coordinates(size(transfers, 1))

    operators = Dict(
        architecture.id => weighted_operator(transfers, architecture.states, weights)
        for architecture in ARCHITECTURES
    )
    kernels = Dict(id => time_reversal_kernel(operator) for (id, operator) in operators)
    summaries = Dict(
        id => singular_summary(operator, DYNAMIC_RANGES_DB)
        for (id, operator) in operators
    )
    added_operator = weighted_operator(transfers, 2:4, weights)
    novelty = Dict(
        level => subspace_novelty(operators[:q00], added_operator, level)
        for level in (20.0, 30.0, 40.0)
    )

    singular_rows = NamedTuple[]
    for architecture in ARCHITECTURES
        summary = summaries[architecture.id]
        for mode_index in eachindex(summary.singular_values)
            push!(singular_rows, (
                architecture=architecture.id, mode_index,
                singular_value=summary.singular_values[mode_index],
                relative_db=summary.singular_value_db[mode_index],
            ))
        end
    end

    focus_rows = NamedTuple[]
    focus_profile_rows = NamedTuple[]
    for architecture in ARCHITECTURES, focal_bin in axes(transfers, 1)
        metrics = focus_metrics(kernels[architecture.id], focal_bin; bin_spacing=spacing)
        push!(focus_rows, (
            architecture=architecture.id, focal_bin,
            focal_coordinate_lambda=coordinates[focal_bin],
            peak_bin=metrics.peak_bin,
            target_is_global_maximum=metrics.target_is_global_maximum,
            fwhm_bins=metrics.fwhm_bins, fwhm_lambda=metrics.fwhm,
            focus_contrast_db=metrics.focus_contrast_db,
            peak_outside_fwhm_db=metrics.peak_outside_fwhm_db,
        ))
        for observation_bin in eachindex(metrics.profile)
            push!(focus_profile_rows, (
                architecture=architecture.id, focal_bin, observation_bin,
                observation_coordinate_lambda=coordinates[observation_bin],
                normalized_intensity=metrics.profile[observation_bin],
            ))
        end
    end
    architecture_summary = architecture_rows(
        operators, summaries, focus_rows, novelty,
    )
    novelty_rows = [(
        dynamic_range_db=level,
        added_state_energy_outside_q00_subspace=novelty[level],
    ) for level in sort!(collect(keys(novelty)))]
    pulse_rows = [(
        frequency_index=index, frequency_hz=frequencies_hz[index],
        frequency_ratio=frequency_ratios[index], pulse_amplitude_weight=weights[index],
    ) for index in eachindex(weights)]

    write_csv(joinpath(OUTPUT_ROOT, "architecture_summary.csv"), architecture_summary)
    write_csv(joinpath(OUTPUT_ROOT, "singular_values.csv"), singular_rows)
    write_csv(joinpath(OUTPUT_ROOT, "focus_metrics.csv"), focus_rows)
    write_csv(joinpath(OUTPUT_ROOT, "focus_profiles.csv"), focus_profile_rows)
    write_csv(joinpath(OUTPUT_ROOT, "subspace_novelty.csv"), novelty_rows)
    write_csv(joinpath(OUTPUT_ROOT, "pulse_weights.csv"), pulse_rows)
    plot_singular_spectra(summaries)
    plot_center_focus(kernels, coordinates)
    plot_coherence_kernels(kernels, coordinates)
    write_readme(architecture_summary, weights, frequencies_hz, novelty, common_scale)
    jldsave(
        joinpath(OUTPUT_ROOT, "o8_time_reversal_audit.jld2");
        transfers, frequencies_hz, frequency_ratios, weights, coordinates,
        operators, kernels, summaries, novelty, architecture_summary, focus_rows,
    )
    for row in architecture_summary
        println(
            "[+] $(row.architecture): stable_rank=$(round(row.stable_rank; digits=3)), " *
            "modes30=$(row.modes_30db), median_FWHM=$(round(row.median_fwhm_lambda; digits=3)) lambdaP",
        )
    end
    architecture_summary
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_o8_time_reversal_audit()
end

end
