ENV["GKSwstype"] = "100"

using JLD2
using Plots

include(joinpath(@__DIR__, "lens_design.jl"))
using .LensDesign

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const LENS_FREQUENCY_HZ = parse(Float64, get(ENV, "METAMATERIALS_LENS_FREQUENCY_HZ", "220000"))
const LENS_ELEMENT_COUNT = parse(Int, get(ENV, "METAMATERIALS_LENS_ELEMENT_COUNT", "7"))
const LENS_FOCAL_DISTANCE_MM = parse(Float64, get(ENV, "METAMATERIALS_LENS_FOCAL_DISTANCE_MM", "60"))
const LENS_MINIMUM_AMPLITUDE = parse(Float64, get(ENV, "METAMATERIALS_LENS_MINIMUM_AMPLITUDE", "0"))
const LENS_OPTIMIZER = get(ENV, "METAMATERIALS_LENS_OPTIMIZER", "amplitude")
const LENS_MINIMUM_FOCUS_AMPLITUDE = parse(
    Float64,
    get(ENV, "METAMATERIALS_LENS_MINIMUM_FOCUS_AMPLITUDE", "1.0"),
)
const LIBRARY_PATH = get(
    ENV,
    "METAMATERIALS_LENS_LIBRARY",
    joinpath(PROJECT_ROOT, "tmp", "gap_pilot_220khz", "gap_pilot.csv"),
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_LENS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "lens_prototype_$(Int(round(LENS_FREQUENCY_HZ / 1e3)))khz"),
)

function write_selection_csv(path, selection, phase_values)
    open(path, "w") do io
        println(io, "position_index,center_y_mm,case_id,periods,length_mm,gap_mm,amplitude_mm,H_real,H_imag,H_amplitude,H_phase_rad,target_phase_rad,phase_error_rad")
        for (index, (center, entry, target, error)) in enumerate(zip(
            selection.centers_mm,
            selection.entries,
            phase_values.target,
            phase_values.errors,
        ))
            println(io, join((
                index,
                center,
                entry.case_id,
                entry.periods,
                entry.length_mm,
                entry.gap_mm,
                entry.amplitude_mm,
                real(entry.transfer),
                imag(entry.transfer),
                abs(entry.transfer),
                angle(entry.transfer),
                target,
                error,
            ), ','))
        end
    end
end

function find_local_peak(amplitude, x_mm, y_mm; x_range=(35.0, 85.0), y_limit=25.0)
    best_amplitude = -Inf
    best_index = CartesianIndex(1, 1)
    for index in CartesianIndices(amplitude)
        y = y_mm[index[1]]
        x = x_mm[index[2]]
        if x_range[1] <= x <= x_range[2] && abs(y) <= y_limit && amplitude[index] > best_amplitude
            best_amplitude = amplitude[index]
            best_index = index
        end
    end
    (
        amplitude=best_amplitude,
        x_mm=x_mm[best_index[2]],
        y_mm=y_mm[best_index[1]],
        index=best_index,
    )
end

function axial_depth_of_focus(amplitude, x_mm, y_mm)
    center_y_index = argmin(abs.(y_mm))
    axial = amplitude[center_y_index, :]
    peak_index = argmax(axial)
    level = axial[peak_index] / sqrt(2.0)
    first_index = peak_index
    last_index = peak_index
    while first_index > firstindex(axial) && axial[first_index - 1] >= level
        first_index -= 1
    end
    while last_index < lastindex(axial) && axial[last_index + 1] >= level
        last_index += 1
    end
    (
        depth_mm=x_mm[last_index] - x_mm[first_index],
        start_mm=x_mm[first_index],
        stop_mm=x_mm[last_index],
        level_amplitude=level,
    )
end

function save_selection_plot(path, selection, phase_values)
    centers = selection.centers_mm
    amplitudes = abs.(getproperty.(selection.entries, :transfer))
    phase_plot = plot(
        centers,
        phase_values.target;
        marker=:circle,
        linewidth=2,
        label="required",
        xlabel="Element center y, mm",
        ylabel="phase, rad",
        title="Selected phase profile",
        gridalpha=0.25,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    plot!(
        phase_plot,
        centers,
        phase_values.selected;
        marker=:diamond,
        linewidth=2,
        label="selected H",
    )
    amplitude_plot = bar(
        string.(round.(centers; digits=1)),
        amplitudes;
        xlabel="Element center y, mm",
        ylabel="|H|",
        title="Selected transmission amplitudes",
        label=false,
        ylim=(0.0, max(0.32, 1.12 * maximum(amplitudes))),
        gridalpha=0.25,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    figure = plot(phase_plot, amplitude_plot; layout=(2, 1), size=(1000, 900))
    savefig(figure, path)
end

function save_field_plot(path, amplitude, x_mm, y_mm, config, peak)
    figure = heatmap(
        x_mm,
        y_mm,
        amplitude;
        xlabel="Distance after lens x, mm",
        ylabel="Transverse coordinate y, mm",
        title="Independent-element prediction: |u(x,y)| / |u_inc|",
        color=:viridis,
        colorbar_title="amplitude",
        aspect_ratio=:equal,
        size=(1150, 650),
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    scatter!(
        figure,
        [config.focal_distance_mm],
        [0.0];
        marker=:xcross,
        markersize=9,
        markerstrokewidth=3,
        color=:white,
        label="target",
    )
    scatter!(
        figure,
        [peak.x_mm],
        [peak.y_mm];
        marker=:circle,
        markersize=6,
        color=:red,
        label="local maximum",
    )
    savefig(figure, path)
end

function save_profile_plot(path, amplitude, reference_amplitude, x_mm, y_mm, config, peak)
    center_y_index = argmin(abs.(y_mm))
    peak_x_index = peak.index[2]
    axial = plot(
        x_mm,
        amplitude[center_y_index, :];
        linewidth=2,
        label="lens",
        xlabel="x, mm",
        ylabel="on-axis amplitude",
        title="Axial field",
        gridalpha=0.25,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    plot!(axial, x_mm, reference_amplitude[center_y_index, :]; linewidth=2, label="uniform aperture")
    vline!(axial, [config.focal_distance_mm]; linestyle=:dash, color=:gray, label="target")

    transverse = plot(
        y_mm,
        amplitude[:, peak_x_index];
        linewidth=2,
        label=false,
        xlabel="y, mm",
        ylabel="amplitude",
        title="Transverse profile at x=$(round(peak.x_mm; digits=1)) mm",
        gridalpha=0.25,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    figure = plot(axial, transverse; layout=(2, 1), size=(1000, 900))
    savefig(figure, path)
end

function save_animation(path, field, x_mm, y_mm, config)
    color_limit = maximum(abs, field)
    phases = range(0.0, 2.0 * pi; length=33)[1:end-1]
    animation = @animate for phase in phases
        instantaneous = real.(field .* cis(-phase))
        figure = heatmap(
            x_mm,
            y_mm,
            instantaneous;
            xlabel="x after lens, mm",
            ylabel="y, mm",
            title="Harmonic field at $(round(config.frequency_hz / 1e3; digits=1)) kHz; phase=$(round(phase; digits=2)) rad",
            color=:balance,
            clim=(-color_limit, color_limit),
            colorbar_title="u / |u_inc|",
            aspect_ratio=:equal,
            size=(950, 560),
            left_margin=7Plots.mm,
            bottom_margin=7Plots.mm,
        )
        scatter!(
            figure,
            [config.focal_distance_mm],
            [0.0];
            marker=:xcross,
            markersize=8,
            markerstrokewidth=2,
            color=:black,
            label=false,
        )
        figure
    end
    gif(animation, path; fps=12, show_msg=false)
end

function run_lens_prototype()
    config = LensConfig(
        frequency_hz=LENS_FREQUENCY_HZ,
        element_count=LENS_ELEMENT_COUNT,
        focal_distance_mm=LENS_FOCAL_DISTANCE_MM,
    )
    entries = filter(
        entry -> abs(entry.transfer) >= LENS_MINIMUM_AMPLITUDE,
        load_gap_library(LIBRARY_PATH),
    )
    isempty(entries) && error("no library entries pass the requested amplitude threshold")
    objective = LensObjectiveConfig(
        axial_offsets_mm=vcat(collect(-30.0:2.0:-10.0), collect(10.0:2.0:40.0)),
        transverse_offsets_mm=[-24.0, -16.0, -12.0, 12.0, 16.0, 24.0],
        axial_weight=1.0,
        sidelobe_weight=0.2,
        minimum_focus_amplitude=LENS_MINIMUM_FOCUS_AMPLITUDE,
    )
    selection = if LENS_OPTIMIZER == "amplitude"
        select_lens(entries; config)
    elseif LENS_OPTIMIZER == "contrast"
        select_lens_contrast(entries; config, objective)
    else
        error("unknown lens optimizer: $LENS_OPTIMIZER")
    end
    phase_values = phase_design_values(selection; config)

    x_mm = collect(range(2.0, 100.0; length=197))
    y_mm = collect(range(-45.0, 45.0; length=181))
    field = field_grid(x_mm, y_mm, selection; config)
    reference_field = reference_field_grid(x_mm, y_mm; config)
    amplitude = abs.(field)
    reference_amplitude = abs.(reference_field)
    peak = find_local_peak(amplitude, x_mm, y_mm)
    depth_of_focus = axial_depth_of_focus(amplitude, x_mm, y_mm)

    target_field = field_at(config.focal_distance_mm, 0.0, selection; config)
    target_y_index = argmin(abs.(y_mm))
    target_x_index = argmin(abs.(x_mm .- config.focal_distance_mm))
    target_reference = reference_field[target_y_index, target_x_index]
    rms_phase_error = sqrt(sum(abs2, phase_values.errors) / length(phase_values.errors))

    mkpath(OUTPUT_ROOT)
    selection_csv = joinpath(OUTPUT_ROOT, "selected_elements.csv")
    selection_png = joinpath(OUTPUT_ROOT, "lens_selection.png")
    field_png = joinpath(OUTPUT_ROOT, "lens_field_amplitude.png")
    profiles_png = joinpath(OUTPUT_ROOT, "lens_field_profiles.png")
    animation_gif = joinpath(OUTPUT_ROOT, "lens_wave.gif")
    result_jld2 = joinpath(OUTPUT_ROOT, "lens_prototype.jld2")

    write_selection_csv(selection_csv, selection, phase_values)
    save_selection_plot(selection_png, selection, phase_values)
    save_field_plot(field_png, amplitude, x_mm, y_mm, config, peak)
    save_profile_plot(profiles_png, amplitude, reference_amplitude, x_mm, y_mm, config, peak)
    save_animation(animation_gif, field, x_mm, y_mm, config)
    JLD2.jldsave(
        result_jld2;
        format_version=1,
        config,
        selection,
        phase_values,
        x_mm,
        y_mm,
        field,
        reference_field,
        target_field,
        target_reference,
        peak,
        rms_phase_error,
        depth_of_focus,
        optimizer=LENS_OPTIMIZER,
        objective,
    )

    println("[+] Selected elements: $selection_csv")
    println("[+] Target amplitude / incident: $(abs(target_field))")
    println("[+] Target amplitude / uniform aperture: $(abs(target_field) / abs(target_reference))")
    println("[+] RMS phase error: $rms_phase_error rad")
    println("[+] Local peak: amplitude=$(peak.amplitude), x=$(peak.x_mm) mm, y=$(peak.y_mm) mm")
    println("[+] Axial -3 dB depth: $(depth_of_focus.depth_mm) mm ($(depth_of_focus.start_mm)..$(depth_of_focus.stop_mm) mm)")
    println("[+] Field figure: $field_png")
    println("[+] Animation: $animation_gif")
    result_jld2
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_lens_prototype()
end
