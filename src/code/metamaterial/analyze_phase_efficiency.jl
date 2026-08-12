module PhaseEfficiencyAnalysis

using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const LIBRARY_ROOT = joinpath(PROJECT_ROOT, "tmp", "target_frequency_library")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "phase_efficiency")
const FREQUENCIES_KHZ = [242, 244]
const TARGET_PHASES_DEG = collect(0.0:5.0:355.0)
const TOLERANCES_DEG = [5.0, 10.0, 15.0, 20.0]

function read_csv_rows(path::AbstractString)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1) if !isempty(strip(line))]
end

circular_error_deg(a, b) = abs(mod(a - b + 180.0, 360.0) - 180.0)

function library_rows(frequency_khz::Integer)
    path = joinpath(LIBRARY_ROOT, "library_$(frequency_khz)khz.csv")
    [
        (
            case_id=row["case_id"],
            periods=parse(Int, row["periods"]),
            length_mm=parse(Float64, row["length_mm"]),
            gap_mm=parse(Float64, row["gap_mm"]),
            amplitude=parse(Float64, row["amplitude"]),
            phase_deg=mod(rad2deg(parse(Float64, row["phase_rad"])), 360.0),
            valid=lowercase(row["spectral_valid"]) == "true",
        )
        for row in read_csv_rows(path)
    ]
end

function best_for_phase(rows, target_phase_deg, tolerance_deg)
    candidates = filter(rows) do row
        row.valid && circular_error_deg(row.phase_deg, target_phase_deg) <= tolerance_deg
    end
    isempty(candidates) && return nothing
    candidates[argmax(getproperty.(candidates, :amplitude))]
end

function efficiency_rows(frequency_khz::Integer)
    library = library_rows(frequency_khz)
    [
        let best = best_for_phase(library, target, tolerance)
            (
                frequency_khz,
                target_phase_deg=target,
                tolerance_deg=tolerance,
                available=!isnothing(best),
                max_amplitude=isnothing(best) ? NaN : best.amplitude,
                best_phase_deg=isnothing(best) ? NaN : best.phase_deg,
                phase_error_deg=isnothing(best) ? NaN : circular_error_deg(best.phase_deg, target),
                case_id=isnothing(best) ? "" : best.case_id,
                periods=isnothing(best) ? 0 : best.periods,
                length_mm=isnothing(best) ? NaN : best.length_mm,
                gap_mm=isnothing(best) ? NaN : best.gap_mm,
            )
        end
        for tolerance in TOLERANCES_DEG
        for target in TARGET_PHASES_DEG
    ]
end

function pareto_candidates(rows, target_phase_deg)
    candidates = [
        merge(row, (phase_error_deg=circular_error_deg(row.phase_deg, target_phase_deg),))
        for row in rows if row.valid
    ]
    filter(candidates) do candidate
        !any(candidates) do other
            other.phase_error_deg <= candidate.phase_error_deg &&
            other.amplitude >= candidate.amplitude &&
            (other.phase_error_deg < candidate.phase_error_deg || other.amplitude > candidate.amplitude)
        end
    end |> values -> sort(values; by=value -> value.phase_error_deg)
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function save_efficiency_plot(rows, frequency_khz)
    library = library_rows(frequency_khz)
    figure = scatter(
        getproperty.(library, :phase_deg),
        getproperty.(library, :amplitude);
        xlabel="Phase modulo 2pi, degrees",
        ylabel="|H|",
        title="Best transmission available near each phase, $(frequency_khz) kHz",
        color=:gray,
        alpha=0.45,
        markersize=5,
        label="library elements",
        xlims=(0, 360),
        xticks=0:60:360,
        gridalpha=0.25,
        size=(1000, 650),
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    for (tolerance, marker) in zip(TOLERANCES_DEG, (:circle, :diamond, :square, :utriangle))
        current = filter(row -> row.tolerance_deg == tolerance && row.available, rows)
        plot!(
            figure,
            getproperty.(current, :target_phase_deg),
            getproperty.(current, :max_amplitude);
            linewidth=2,
            marker,
            markersize=3,
            label="±$(Int(tolerance))° tolerance",
        )
    end
    path = joinpath(OUTPUT_ROOT, "phase_efficiency_$(frequency_khz)khz.png")
    savefig(figure, path)
    path
end

function analyze()
    mkpath(OUTPUT_ROOT)
    all_rows = NamedTuple[]
    for frequency_khz in FREQUENCIES_KHZ
        rows = efficiency_rows(frequency_khz)
        append!(all_rows, rows)
        println("[+] $(save_efficiency_plot(rows, frequency_khz))")
    end
    path = joinpath(OUTPUT_ROOT, "phase_efficiency.csv")
    write_csv(path, all_rows)

    required_path = joinpath(PROJECT_ROOT, "tmp", "lens_prototype_242khz_9elements", "selected_elements.csv")
    required = read_csv_rows(required_path)
    unique_targets = sort(unique(mod(rad2deg(parse(Float64, row["target_phase_rad"])), 360.0) for row in required))
    library = library_rows(242)
    pareto_rows = NamedTuple[]
    for target in unique_targets
        for candidate in pareto_candidates(library, target)
            push!(pareto_rows, (
                target_phase_deg=target,
                candidate.case_id,
                candidate.periods,
                candidate.length_mm,
                candidate.gap_mm,
                candidate.phase_deg,
                candidate.phase_error_deg,
                candidate.amplitude,
            ))
        end
    end
    pareto_path = joinpath(OUTPUT_ROOT, "lens_phase_pareto_242khz.csv")
    write_csv(pareto_path, pareto_rows)
    println("[+] $path")
    println("[+] $pareto_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze()
end

end
