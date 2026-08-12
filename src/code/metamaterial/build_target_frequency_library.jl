module TargetFrequencyLibrary

using JLD2
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "target_frequency_library")
const TARGET_FREQUENCIES_HZ = [242.0e3, 244.0e3]
const MINIMUM_WORKING_AMPLITUDE = 0.1

function read_csv_rows(path::AbstractString)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1) if !isempty(strip(line))]
end

number_label(value::Real) = replace(string(Float64(value)), "." => "p")

function characteristic_path(root, case_id, amplitude_mm, periods)
    joinpath(
        root,
        "characteristics",
        case_id,
        "characteristics_data_sin_A_$(Float64(amplitude_mm))_N_$(periods)_STG_F_220.0.jld2",
    )
end

function candidate(row, root)
    periods = parse(Int, row["periods"])
    length_mm = parse(Float64, row["length_mm"])
    amplitude_mm = parse(Float64, row["amplitude_mm"])
    gap_mm = haskey(row, "gap_mm") ?
             parse(Float64, row["gap_mm"]) :
             parse(Float64, row["minimum_gap_mm"])
    case_id = row["case_id"]
    (
        case_id=case_id,
        periods=periods,
        length_mm=length_mm,
        gap_mm=gap_mm,
        amplitude_mm=amplitude_mm,
        characteristic_path=characteristic_path(root, case_id, amplitude_mm, periods),
    )
end

geometry_key(item) = (item.periods, item.length_mm, item.gap_mm)

function collect_candidates()
    gap_root = joinpath(PROJECT_ROOT, "tmp", "gap_pilot_220khz")
    period_root = joinpath(PROJECT_ROOT, "tmp", "period_length_pilot_220khz")
    combined = Dict{Tuple{Int, Float64, Float64}, NamedTuple}()

    for row in read_csv_rows(joinpath(gap_root, "gap_pilot.csv"))
        item = candidate(row, gap_root)
        combined[geometry_key(item)] = item
    end
    for row in read_csv_rows(joinpath(period_root, "period_length_pilot.csv"))
        item = candidate(row, period_root)
        get!(combined, geometry_key(item), item)
    end
    sort(collect(values(combined)); by=item -> geometry_key(item))
end

function nearest_point(path::AbstractString, frequency_hz::Real)
    isfile(path) || error("characteristics not found: $path")
    transfer = JLD2.load(path)["calibrated_transmission"]
    index = argmin(abs.(transfer.frequency_hz .- frequency_hz))
    (
        sampled_frequency_hz=transfer.frequency_hz[index],
        transfer=ComplexF64(transfer.transfer[index]),
        amplitude=transfer.amplitude[index],
        phase_rad=transfer.phase_rad[index],
        spectral_valid=transfer.valid[index],
    )
end

function library_rows(candidates, frequency_hz)
    [
        let point = nearest_point(item.characteristic_path, frequency_hz)
            (
                case_id=item.case_id,
                profile_family="sinusoidal",
                periods=item.periods,
                length_mm=item.length_mm,
                gap_mm=item.gap_mm,
                amplitude_mm=item.amplitude_mm,
                requested_frequency_hz=frequency_hz,
                sampled_frequency_hz=point.sampled_frequency_hz,
                H_real=real(point.transfer),
                H_imag=imag(point.transfer),
                amplitude=point.amplitude,
                phase_rad=point.phase_rad,
                spectral_valid=point.spectral_valid,
                valid=point.spectral_valid && point.amplitude >= MINIMUM_WORKING_AMPLITUDE,
                source_characteristics=item.characteristic_path,
            )
        end
        for item in candidates
    ]
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

function phase_coverage(rows, threshold)
    phases = sort(collect(
        mod(row.phase_rad, 2.0 * pi)
        for row in rows
        if row.spectral_valid && row.amplitude >= threshold
    ))
    isempty(phases) && return (count=0, covered_arc_rad=0.0, missing_arc_rad=2.0 * pi)
    length(phases) == 1 && return (count=1, covered_arc_rad=0.0, missing_arc_rad=2.0 * pi)
    circular_gaps = vcat(diff(phases), first(phases) + 2.0 * pi - last(phases))
    missing = maximum(circular_gaps)
    (count=length(phases), covered_arc_rad=2.0 * pi - missing, missing_arc_rad=missing)
end

function save_complex_plot(rows, frequency_hz)
    valid_rows = filter(row -> row.valid, rows)
    limit = 1.15 * maximum(getproperty.(valid_rows, :amplitude))
    figure = scatter(
        getproperty.(valid_rows, :H_real),
        getproperty.(valid_rows, :H_imag);
        marker_z=getproperty.(valid_rows, :amplitude),
        color=:viridis,
        colorbar_title="|H|",
        xlabel="Re H",
        ylabel="Im H",
        title="Target-frequency element library, $(frequency_hz / 1e3) kHz",
        aspect_ratio=:equal,
        xlims=(-limit, limit),
        ylims=(-limit, limit),
        markersize=8,
        legend=false,
        gridalpha=0.25,
        size=(850, 720),
    )
    hline!(figure, [0.0]; color=:gray, linestyle=:dot)
    vline!(figure, [0.0]; color=:gray, linestyle=:dot)
    path = joinpath(OUTPUT_ROOT, "library_$(Int(round(frequency_hz / 1e3)))khz_complex_H.png")
    savefig(figure, path)
    path
end

function build_library()
    mkpath(OUTPUT_ROOT)
    candidates = collect_candidates()
    summary = NamedTuple[]
    for frequency_hz in TARGET_FREQUENCIES_HZ
        rows = library_rows(candidates, frequency_hz)
        label = Int(round(frequency_hz / 1e3))
        csv_path = joinpath(OUTPUT_ROOT, "library_$(label)khz.csv")
        write_csv(csv_path, rows)
        figure_path = save_complex_plot(rows, frequency_hz)
        for threshold in (0.1, 0.2, 0.3, 0.4, 0.5)
            coverage = phase_coverage(rows, threshold)
            push!(summary, (
                frequency_hz,
                amplitude_threshold=threshold,
                element_count=coverage.count,
                covered_arc_rad=coverage.covered_arc_rad,
                covered_arc_deg=rad2deg(coverage.covered_arc_rad),
                missing_arc_rad=coverage.missing_arc_rad,
            ))
        end
        best = rows[argmax(getproperty.(rows, :amplitude))]
        println("[+] $csv_path")
        println("[+] $figure_path")
        println("[+] Best at $(label) kHz: $(best.case_id), |H|=$(best.amplitude), phase=$(best.phase_rad) rad")
    end
    summary_path = joinpath(OUTPUT_ROOT, "phase_coverage_summary.csv")
    write_csv(summary_path, summary)
    println("[+] $summary_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    build_library()
end

end
