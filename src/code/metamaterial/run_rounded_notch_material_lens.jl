module RunRoundedNotchMaterialLens

ENV["GKSwstype"] = "100"

using JLD2
using Plots

include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
using .SinusoidalMaterialLens
include(joinpath(@__DIR__, "material_lens_mesher.jl"))
using .MaterialLensMesher

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_ROUNDED_NOTCH_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "rounded_notch_material_lens_242khz"),
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

function summary_row(candidate)
    (
        configuration=String(candidate.configuration),
        lens_material=candidate.lens_material.name,
        output_material=candidate.output_material.name,
        matching_thickness_mm=candidate.matching_thickness_mm,
        cell_length_mm=first(candidate.cells).length_mm,
        harmonic_focus_amplitude=candidate.harmonic_focus_amplitude,
        pulse_peak_amplitude=candidate.pulse_peak_amplitude,
        pulse_peak_time_us=1e6 * candidate.pulse_peak_time_s,
        uniform_pulse_peak_amplitude=candidate.uniform_pulse_peak_amplitude,
        focus_gain_over_uniform=candidate.pulse_peak_amplitude /
                                candidate.uniform_pulse_peak_amplitude,
    )
end

function save_pulses(path, candidates, config)
    figure = plot(
        xlabel="time, microseconds",
        ylabel="reduced longitudinal amplitude",
        title="Rounded-notch lens: four-cycle pulse at target",
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
            1e6 .* response.time_s,
            response.focused;
            linewidth=2,
            label=replace(String(candidate.configuration), '_' => ' '),
        )
    end
    savefig(figure, path)
end

function representative_cell(candidate)
    profiled = filter(cell -> cell.minimum_gap_mm < cell.height_mm, candidate.cells)
    isempty(profiled) ? first(candidate.cells) : profiled[argmin(getproperty.(profiled, :minimum_gap_mm))]
end

function save_cell_shapes(path, candidates)
    panels = Plots.Plot[]
    for candidate in candidates
        cell = representative_cell(candidate)
        points = MaterialLensMesher.cell_boundary_points(cell, cell.height_mm / 2)
        shape = Shape(first.(points), last.(points))
        panel = plot(
            shape;
            fillcolor=:steelblue,
            fillalpha=0.25,
            linecolor=:navy,
            linewidth=2,
            label=false,
            xlabel="x, mm",
            ylabel="y, mm",
            title="$(replace(String(candidate.configuration), '_' => ' ')); " *
                  "n=$(cell.notch_count), w=$(cell.notch_width_mm), g=$(cell.minimum_gap_mm) mm",
            aspect_ratio=:equal,
            xlims=(-0.5, cell.length_mm + 0.5),
            ylims=(-0.4, cell.height_mm + 0.4),
            gridalpha=0.2,
        )
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(length(panels), 1), size=(1050, 850)), path)
end

function save_aperture(path, candidates)
    panels = Plots.Plot[]
    for candidate in candidates
        rows = rounded_design_rows(candidate)
        panel = scatter(
            getproperty.(rows, :center_y_mm),
            getproperty.(rows, :minimum_gap_mm);
            marker_z=getproperty.(rows, :notch_count),
            color=:viridis,
            markersize=8,
            colorbar_title="notches",
            xlabel="element center y, mm",
            ylabel="tip ligament, mm",
            title=replace(String(candidate.configuration), '_' => ' '),
            ylim=(1.5, 7.4),
            label=false,
            gridalpha=0.25,
        )
        push!(panels, panel)
    end
    savefig(plot(panels...; layout=(3, 1), size=(1000, 950)), path)
end

function run_search()
    mkpath(OUTPUT_ROOT)
    config = MaterialLensConfig(transfer_segments=150)
    candidates = RoundedLensDesignCandidate[]
    for configuration in (:polymer_direct, :polymer_matched, :aluminium)
        println("=== Optimizing rounded notches: $configuration ===")
        candidate = optimize_rounded_configuration(
            configuration;
            config,
            lengths_mm=[16.0, 20.0, 24.0, 28.0, 32.0, 36.0],
            gaps_mm=[2.0, 3.0, 4.0, 5.0],
            notch_counts=[3, 4, 5],
            notch_widths_mm=[1.2, 1.6, 2.0],
            matching_thicknesses_mm=[1.0, 1.4, 1.8, 2.2],
            shortlist_count=14,
        )
        push!(candidates, candidate)
        selection_path = joinpath(OUTPUT_ROOT, "selected_$(configuration).csv")
        write_csv(selection_path, rounded_design_rows(candidate; slot_width_mm=config.slot_width_mm))
        println("[+] $configuration: peak=$(candidate.pulse_peak_amplitude), " *
                "gain=$(candidate.pulse_peak_amplitude / candidate.uniform_pulse_peak_amplitude)")
        println("[+] $selection_path")
    end

    summary_path = joinpath(OUTPUT_ROOT, "rounded_notch_summary.csv")
    pulse_path = joinpath(OUTPUT_ROOT, "rounded_notch_pulses.png")
    cell_path = joinpath(OUTPUT_ROOT, "rounded_notch_cells.png")
    aperture_path = joinpath(OUTPUT_ROOT, "rounded_notch_apertures.png")
    data_path = joinpath(OUTPUT_ROOT, "rounded_notch_designs.jld2")
    write_csv(summary_path, summary_row.(candidates))
    save_pulses(pulse_path, candidates, config)
    save_cell_shapes(cell_path, candidates)
    save_aperture(aperture_path, candidates)
    jldsave(data_path; format_version=1, config, candidates)
    println("[+] Summary: $summary_path")
    println("[+] Cell shapes: $cell_path")
    candidates
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_search()
end

end
