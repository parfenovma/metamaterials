module PathPackingBound

ENV["GKSwstype"] = "100"

using Plots
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_PATH_PACKING_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "path_packing_bound"),
)
const STRIP_WIDTH_MM = 3.2
const PITCH_MM = 3.8
const CENTERLINE_AMPLITUDE_MM = (PITCH_MM - STRIP_WIDTH_MM) / 2
const PRINTED_MIN_RADIUS_MM = 0.35
const TARGETS_MM = (
    (name="group_matched_trim", extra_path_mm=6.20),
    (name="phase_matched_path", extra_path_mm=10.42),
)

function arc_ratio(slope_amplitude; samples=20000)
    angles = ((0:(samples - 1)) .+ 0.5) .* (2pi / samples)
    sum(sqrt.(1 .+ slope_amplitude^2 .* cos.(angles).^2)) / samples
end

function best_sinusoidal_meander(length_mm, center_radius_mm)
    amplitude = CENTERLINE_AMPLITUDE_MM
    maximum_periods = floor(Int, length_mm / (2pi * sqrt(amplitude * center_radius_mm)))
    maximum_periods <= 0 && return (periods=0, extra_path_mm=0.0, arc_length_mm=length_mm)
    slope_amplitude = 2pi * maximum_periods * amplitude / length_mm
    arc_length = length_mm * arc_ratio(slope_amplitude)
    (periods=maximum_periods, extra_path_mm=arc_length - length_mm, arc_length_mm=arc_length)
end

function minimum_length_for_target(target_extra_mm, center_radius_mm)
    for length_mm in 5.0:0.1:350.0
        result = best_sinusoidal_meander(length_mm, center_radius_mm)
        result.extra_path_mm >= target_extra_mm && return (
            length_mm,
            result.periods,
            result.extra_path_mm,
        )
    end
    (NaN, 0, NaN)
end

function write_rows(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function run()
    mkpath(OUTPUT_ROOT)
    radius_cases = (
        (name="optimistic_centerline", center_radius_mm=PRINTED_MIN_RADIUS_MM),
        (
            name="inner_boundary_printable",
            center_radius_mm=PRINTED_MIN_RADIUS_MM + STRIP_WIDTH_MM / 2,
        ),
    )
    summary_rows = [begin
        length_mm, periods, achieved = minimum_length_for_target(
            target.extra_path_mm,
            radius.center_radius_mm,
        )
        (
            radius_case=radius.name,
            center_radius_mm=radius.center_radius_mm,
            target=target.name,
            target_extra_path_mm=target.extra_path_mm,
            minimum_cell_length_mm=length_mm,
            sinusoidal_periods=periods,
            achieved_extra_path_mm=achieved,
        )
    end for radius in radius_cases for target in TARGETS_MM]
    write_rows(joinpath(OUTPUT_ROOT, "path_packing_summary.csv"), summary_rows)

    length_grid = collect(5.0:0.5:300.0)
    sweep_rows = [begin
        result = best_sinusoidal_meander(length_mm, radius.center_radius_mm)
        (
            radius_case=radius.name,
            center_radius_mm=radius.center_radius_mm,
            cell_length_mm=length_mm,
            periods=result.periods,
            extra_path_mm=result.extra_path_mm,
        )
    end for radius in radius_cases for length_mm in length_grid]
    write_rows(joinpath(OUTPUT_ROOT, "path_packing_sweep.csv"), sweep_rows)

    panel = plot(
        xlabel="cell length, mm",
        ylabel="maximum extra centerline path, mm",
        title="Optimistic in-plane sinusoidal packing bound",
    )
    for radius in radius_cases
        selected = filter(row -> row.radius_case == radius.name, sweep_rows)
        plot!(
            panel,
            getproperty.(selected, :cell_length_mm),
            getproperty.(selected, :extra_path_mm);
            linewidth=2,
            label=@sprintf("%s, R=%.2f mm", radius.name, radius.center_radius_mm),
        )
    end
    for target in TARGETS_MM
        hline!(panel, [target.extra_path_mm]; linestyle=:dash, label=target.name)
    end
    savefig(panel, joinpath(OUTPUT_ROOT, "path_packing_bound.png"))

    for row in summary_rows
        @printf(
            "[+] %s / %s: L >= %.1f mm, N=%d, extra=%.2f mm\n",
            row.radius_case,
            row.target,
            row.minimum_cell_length_mm,
            row.sinusoidal_periods,
            row.achieved_extra_path_mm,
        )
    end
    println("[+] $OUTPUT_ROOT")
    (; summary_rows, sweep_rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
