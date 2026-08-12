module RunSinusoidalMaterialLens

ENV["GKSwstype"] = "100"

using JLD2
using Plots

include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
using .SinusoidalMaterialLens

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_242khz"),
)
const PULSE_CYCLES = parse(
    Float64,
    get(ENV, "METAMATERIALS_MATERIAL_LENS_PULSE_CYCLES", "4.0"),
)

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function configuration_summary(candidate)
    (
        configuration=String(candidate.configuration),
        lens_material=candidate.lens_material.name,
        output_material=candidate.output_material.name,
        matching_material=isnothing(candidate.matching_material) ? "none" : candidate.matching_material.name,
        matching_thickness_mm=candidate.matching_thickness_mm,
        cell_length_mm=first(candidate.cells).length_mm,
        harmonic_focus_amplitude=candidate.harmonic_focus_amplitude,
        pulse_peak_amplitude=candidate.pulse_peak_amplitude,
        pulse_peak_time_us=candidate.pulse_peak_time_s * 1e6,
        uniform_pulse_peak_amplitude=candidate.uniform_pulse_peak_amplitude,
        focus_gain_over_uniform=candidate.pulse_peak_amplitude / candidate.uniform_pulse_peak_amplitude,
    )
end

function save_pulse_plot(path, candidates, config)
    figure = plot(
        xlabel="time, microseconds",
        ylabel="reduced longitudinal amplitude",
        title="$(config.pulse_cycles)-cycle pulse at the target focus",
        size=(1050, 650),
        gridalpha=0.25,
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    for candidate in candidates
        response = pulse_response(
            candidate.cells,
            candidate.lens_material,
            candidate.output_material,
            config;
            matching_material=candidate.matching_material,
            matching_thickness_mm=candidate.matching_thickness_mm,
        )
        plot!(
            figure,
            response.time_s .* 1e6,
            response.focused;
            linewidth=2,
            label=replace(String(candidate.configuration), '_' => ' '),
        )
    end
    savefig(figure, path)
end

function save_geometry_plot(path, candidates)
    panels = Plots.Plot[]
    for candidate in candidates
        rows = design_rows(candidate)
        gaps = getproperty.(rows, :minimum_gap_mm)
        periods = getproperty.(rows, :periods)
        centers = getproperty.(rows, :center_y_mm)
        panel = scatter(
            centers,
            gaps;
            marker_z=periods,
            color=:viridis,
            markersize=9,
            xlabel="element center y, mm",
            ylabel="minimum gap, mm",
            colorbar_title="periods",
            title=replace(String(candidate.configuration), '_' => ' '),
            label=false,
            ylim=(1.0, 7.4),
            gridalpha=0.25,
        )
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(length(panels), 1), size=(1000, 950)), path)
end

function run_material_lens_search()
    mkpath(OUTPUT_ROOT)
    config = MaterialLensConfig(pulse_cycles=PULSE_CYCLES)
    matching = geometric_matching_material(photopolymer(), aluminium_6061())
    quarter_wave = quarter_wave_thickness_mm(matching, config.center_frequency_hz)
    thicknesses = sort(unique(vcat(collect(2.0:0.5:6.0), quarter_wave)))

    candidates = LensDesignCandidate[]
    for configuration in (:polymer_direct, :polymer_matched, :aluminium)
        println("=== Optimizing $(configuration) ===")
        candidate = optimize_configuration(
            configuration;
            config,
            matching_thicknesses_mm=thicknesses,
            shortlist_count=16,
        )
        push!(candidates, candidate)
        rows = design_rows(candidate; slot_width_mm=config.slot_width_mm)
        path = joinpath(OUTPUT_ROOT, "selected_$(configuration).csv")
        write_csv(path, rows)
        println("[+] $(configuration): peak=$(candidate.pulse_peak_amplitude), gain=$(candidate.pulse_peak_amplitude / candidate.uniform_pulse_peak_amplitude)")
        println("[+] $path")
    end

    summaries = configuration_summary.(candidates)
    summary_path = joinpath(OUTPUT_ROOT, "material_lens_summary.csv")
    pulse_path = joinpath(OUTPUT_ROOT, "material_lens_pulses.png")
    geometry_path = joinpath(OUTPUT_ROOT, "material_lens_geometries.png")
    data_path = joinpath(OUTPUT_ROOT, "material_lens_designs.jld2")
    write_csv(summary_path, summaries)
    save_pulse_plot(pulse_path, candidates, config)
    save_geometry_plot(geometry_path, candidates)
    jldsave(
        data_path;
        format_version=1,
        config,
        candidates,
        matching_quarter_wave_thickness_mm=quarter_wave,
    )
    println("[+] Summary: $summary_path")
    println("[+] Pulses: $pulse_path")
    println("[+] Designs: $data_path")
    candidates
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_material_lens_search()
end

end
