module IdealHuygensAperture

ENV["GKSwstype"] = "100"

using Plots

include(joinpath(@__DIR__, "lens_design.jl"))
using .LensDesign

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_IDEAL_APERTURE_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "ideal_huygens_aperture_242khz"),
)
const FREQUENCY_HZ = parse(
    Float64,
    get(ENV, "METAMATERIALS_IDEAL_APERTURE_FREQUENCY_HZ", "242000"),
)
const LONGITUDINAL_SPEED_M_PER_S = parse(
    Float64,
    get(ENV, "METAMATERIALS_IDEAL_APERTURE_CP_M_PER_S", "2340"),
)

function nearest_odd_count(aperture_mm, pitch_mm, slot_width_mm)
    raw = max(3, round(Int, (aperture_mm + slot_width_mm) / pitch_mm))
    isodd(raw) ? raw : raw + (raw * pitch_mm < aperture_mm ? 1 : -1)
end

function aperture_width_mm(config::LensConfig)
    pitch = config.element_width_mm + config.slot_width_mm
    config.element_count * pitch - config.slot_width_mm
end

function make_config(aperture_mm, focal_distance_mm, pitch_mm, slot_width_mm)
    pitch_mm > slot_width_mm || throw(ArgumentError("pitch must exceed slot width"))
    LensConfig(
        frequency_hz=FREQUENCY_HZ,
        longitudinal_speed_m_per_s=LONGITUDINAL_SPEED_M_PER_S,
        element_count=nearest_odd_count(aperture_mm, pitch_mm, slot_width_mm),
        element_width_mm=pitch_mm - slot_width_mm,
        slot_width_mm=slot_width_mm,
        focal_distance_mm=focal_distance_mm,
        aperture_quadrature_points=9,
    )
end

function candidate_configs()
    result = LensConfig[]
    for aperture_mm in (74.0, 110.0, 125.0, 140.0)
        focal_distances = aperture_mm == 74.0 ? (35.0, 60.0) : (60.0,)
        for focal_distance_mm in focal_distances
            for (pitch_mm, slot_width_mm) in ((8.2, 1.2), (4.8, 0.6))
                push!(result, make_config(
                    aperture_mm,
                    focal_distance_mm,
                    pitch_mm,
                    slot_width_mm,
                ))
            end
        end
    end
    result
end

function ideal_selection(config::LensConfig; transmission_amplitude=1.0)
    centers = lens_centers_mm(config)
    k = 2pi / (1000.0 * config.longitudinal_speed_m_per_s / config.frequency_hz)
    entries = [
        LibraryEntry(
            "ideal_$(index)",
            0,
            0.0,
            0.0,
            0.0,
            transmission_amplitude * cis(
                -k * (hypot(config.focal_distance_mm, center) - config.focal_distance_mm),
            ),
        )
        for (index, center) in enumerate(centers)
    ]
    LensSelection(
        centers,
        entries,
        field_at_from_entries(config.focal_distance_mm, 0.0, centers, entries, config),
    )
end

function field_at_from_entries(x_mm, y_mm, centers, entries, config)
    sum(
        entry.transfer * LensDesign.element_kernel(center, x_mm, y_mm, config)
        for (center, entry) in zip(centers, entries)
    )
end

function contiguous_width(coordinate, amplitude)
    peak_index = argmax(amplitude)
    threshold = amplitude[peak_index] / sqrt(2.0)
    left = peak_index
    right = peak_index
    while left > firstindex(amplitude) && amplitude[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(amplitude) && amplitude[right + 1] >= threshold
        right += 1
    end
    (width=coordinate[right] - coordinate[left], start=coordinate[left], stop=coordinate[right])
end

function aperture_metrics(config::LensConfig)
    selection = ideal_selection(config)
    axial_x = collect(range(max(2.0, 0.25 * config.focal_distance_mm),
                            1.9 * config.focal_distance_mm; step=0.25))
    axial = abs.([field_at(x, 0.0, selection; config) for x in axial_x])
    peak_index = argmax(axial)
    peak_x = axial_x[peak_index]
    transverse_y = collect(range(-aperture_width_mm(config) / 2,
                                 aperture_width_mm(config) / 2; step=0.25))
    transverse = abs.([field_at(peak_x, y, selection; config) for y in transverse_y])
    dof = contiguous_width(axial_x, axial)
    fwhm = contiguous_width(transverse_y, transverse)
    pitch = config.element_width_mm + config.slot_width_mm
    (
        requested_frequency_hz=config.frequency_hz,
        wavelength_mm=1000.0 * config.longitudinal_speed_m_per_s / config.frequency_hz,
        element_count=config.element_count,
        pitch_mm=pitch,
        slot_width_mm=config.slot_width_mm,
        aperture_mm=aperture_width_mm(config),
        focal_distance_mm=config.focal_distance_mm,
        numerical_aperture=(aperture_width_mm(config) / 2) /
                           hypot(config.focal_distance_mm, aperture_width_mm(config) / 2),
        target_gain=abs(field_at(config.focal_distance_mm, 0.0, selection; config)),
        axial_peak_gain=axial[peak_index],
        axial_peak_x_mm=peak_x,
        dof_minus_3db_mm=dof.width,
        transverse_fwhm_mm=fwhm.width,
        required_phase_span_rad=required_phase_span(config),
    )
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

function qualifying_index(rows)
    valid = findall(row -> row.target_gain >= 2.8 && row.dof_minus_3db_mm <= 30.0, rows)
    isempty(valid) && return argmax(getproperty.(rows, :target_gain))
    sort(valid; by=index -> (
        rows[index].aperture_mm,
        rows[index].focal_distance_mm,
        -rows[index].target_gain,
    ))[1]
end

function save_summary_plot(path, rows, selected_index)
    colors = [row.pitch_mm < 6.0 ? :royalblue : :darkorange for row in rows]
    markers = [row.focal_distance_mm < 50.0 ? :diamond : :circle for row in rows]
    gain = scatter(
        getproperty.(rows, :aperture_mm),
        getproperty.(rows, :target_gain);
        markercolor=colors,
        markershape=markers,
        markersize=7,
        label=false,
        xlabel="actual aperture D, mm",
        ylabel="ideal gain at target",
        title="Ideal phase-only Huygens aperture",
        gridalpha=0.25,
    )
    hline!(gain, [2.8]; color=:gray, linestyle=:dash, label="G=2.8")
    scatter!(gain, [rows[selected_index].aperture_mm], [rows[selected_index].target_gain];
             marker=:star5, markersize=13, color=:red, label="selected")
    dof = scatter(
        getproperty.(rows, :aperture_mm),
        getproperty.(rows, :dof_minus_3db_mm);
        markercolor=colors,
        markershape=markers,
        markersize=7,
        label=false,
        xlabel="actual aperture D, mm",
        ylabel="axial DOF -3 dB, mm",
        gridalpha=0.25,
    )
    hline!(dof, [30.0]; color=:gray, linestyle=:dash, label="30 mm")
    scatter!(dof, [rows[selected_index].aperture_mm], [rows[selected_index].dof_minus_3db_mm];
             marker=:star5, markersize=13, color=:red, label="selected")
    savefig(plot(gain, dof; layout=(2, 1), size=(1050, 900), margin=5Plots.mm), path)
end

function save_selected_artifacts(output_root, config, row)
    selection = ideal_selection(config)
    x_mm = collect(range(2.0, max(100.0, 1.7 * config.focal_distance_mm); step=0.5))
    y_limit = max(45.0, 0.65 * aperture_width_mm(config))
    y_mm = collect(range(-y_limit, y_limit; step=0.5))
    field = field_grid(x_mm, y_mm, selection; config)
    amplitude = abs.(field)
    figure = heatmap(
        x_mm,
        y_mm,
        amplitude;
        xlabel="distance after lens x, mm",
        ylabel="transverse coordinate y, mm",
        title="Ideal broadband phase-only lens: |u| / |u_inc|",
        color=:viridis,
        colorbar_title="amplitude",
        aspect_ratio=:equal,
        size=(1200, 700),
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    scatter!(figure, [config.focal_distance_mm], [0.0]; marker=:xcross,
             markersize=9, markerstrokewidth=3, color=:white, label="target")
    savefig(figure, joinpath(output_root, "selected_ideal_field.png"))

    phases = angle.(getproperty.(selection.entries, :transfer))
    phase_rows = [(
        position_index=index,
        center_y_mm=center,
        required_phase_rad=phase,
        required_phase_deg=rad2deg(phase),
    ) for (index, (center, phase)) in enumerate(zip(selection.centers_mm, phases))]
    write_csv(joinpath(output_root, "selected_required_phases.csv"), phase_rows)

    phase_plot = plot(
        selection.centers_mm,
        phases;
        marker=:circle,
        linewidth=2,
        xlabel="element center y, mm",
        ylabel="required phase, rad",
        title="Physical-cell targets for selected aperture",
        label=false,
        gridalpha=0.25,
    )
    savefig(phase_plot, joinpath(output_root, "selected_required_phases.png"))

    color_limit = maximum(abs, field)
    animation = @animate for temporal_phase in range(0.0, 2pi; length=25)[1:end-1]
        instantaneous = real.(field .* cis(-temporal_phase))
        frame = heatmap(
            x_mm,
            y_mm,
            instantaneous;
            xlabel="distance after lens x, mm",
            ylabel="transverse coordinate y, mm",
            title="Ideal broadband lens at 242 kHz",
            color=:balance,
            clim=(-color_limit, color_limit),
            colorbar_title="u / |u_inc|",
            aspect_ratio=:equal,
            size=(1000, 600),
        )
        scatter!(frame, [config.focal_distance_mm], [0.0]; marker=:xcross,
                 markersize=8, markerstrokewidth=2, color=:black, label=false)
        frame
    end
    gif(animation, joinpath(output_root, "selected_ideal_wave.gif"); fps=12, show_msg=false)
    println("Selected aperture: $row")
end

function main()
    configs = candidate_configs()
    rows = aperture_metrics.(configs)
    selected_index = qualifying_index(rows)
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "ideal_aperture_summary.csv"), rows)
    save_summary_plot(joinpath(OUTPUT_ROOT, "ideal_aperture_summary.png"), rows, selected_index)
    save_selected_artifacts(OUTPUT_ROOT, configs[selected_index], rows[selected_index])
    foreach(println, rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
