module BalancedHuygensPilot

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_BALANCED_HUYGENS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "balanced_huygens_pilot_242khz"),
)
const DAMPING_SCALE = parse(
    Float64,
    get(ENV, "METAMATERIALS_BALANCED_HUYGENS_DAMPING_SCALE", "0.25"),
)
const SERIES_HEIGHTS_MM = (0.45, 0.65, 0.85, 1.05)
const SHUNT_LENGTHS_MM = (2.4, 3.0, 3.6, 4.2)
const FREQUENCIES_HZ = collect(226.0e3:2.0e3:258.0e3)
const TARGET_FREQUENCY_HZ = 242.0e3
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "balanced_huygens_mesher.jl"))
using .BalancedHuygensMesher

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) && startswith(REQUESTED_STAGE, "solve-")
    include(joinpath(@__DIR__, "harmonic_solver.jl"))
    using .HarmonicElasticity
    using JLD2
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
    include(joinpath(@__DIR__, "lens_design.jl"))
    using .LensDesign
end

number_label(value) = replace(@sprintf("%.2f", Float64(value)), "." => "p")
series_id(height) = "series_hs$(number_label(height))"
shunt_id(length_mm) = "shunt_lm$(number_label(length_mm))"
balanced_id(height, length_mm) =
    "balanced_hs$(number_label(height))_lm$(number_label(length_mm))"

config(height, length_mm) = BalancedHuygensConfig(
    series_height_mm=Float64(height),
    shunt_mass_length_mm=Float64(length_mm),
)

function balanced_cases()
    [
        (
            id=balanced_id(height, length_mm),
            series_height_mm=Float64(height),
            shunt_mass_length_mm=Float64(length_mm),
            config=config(height, length_mm),
        )
        for height in SERIES_HEIGHTS_MM
        for length_mm in SHUNT_LENGTHS_MM
    ]
end

function model_specs()
    specs = [(id="reference", variant=:r0, config=config(0.85, 3.0))]
    append!(specs, [
        (id=series_id(height), variant=:series, config=config(height, 3.0))
        for height in SERIES_HEIGHTS_MM
    ])
    append!(specs, [
        (id=shunt_id(length_mm), variant=:shunt, config=config(0.85, length_mm))
        for length_mm in SHUNT_LENGTHS_MM
    ])
    append!(specs, [
        (id=case.id, variant=:balanced, config=case.config)
        for case in balanced_cases()
    ])
    specs
end

mesh_path(id) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(id).msh")
model_path(id) = joinpath(OUTPUT_ROOT, "models", "model_$(id).json")
response_path(id) = joinpath(OUTPUT_ROOT, "harmonic", "response_$(id).jld2")

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for spec in model_specs()
            build_balanced_huygens_mesh(
                mesh_path(spec.id);
                config=spec.config,
                variant=spec.variant,
                size_min_mm=0.06,
                size_max_mm=0.28,
            )
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    for spec in model_specs()
        convert_mesh(mesh_path(spec.id); output_dir=joinpath(OUTPUT_ROOT, "models"))
    end
end

harmonic_config() = HarmonicConfig(
    rayleigh_alpha=DAMPING_SCALE * 79560.0,
    rayleigh_beta=DAMPING_SCALE * 2.5e-9,
)

function save_response(id)
    points = solve_harmonic_sweep(model_path(id), FREQUENCIES_HZ; config=harmonic_config())
    path = response_path(id)
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        frequency_hz=getproperty.(points, :frequency_hz),
        right_displacement=getproperty.(points, :right_displacement),
        source_power_w_per_m=getproperty.(points, :source_power_w_per_m),
        left_absorbed_power_w_per_m=getproperty.(points, :left_absorbed_power_w_per_m),
        right_absorbed_power_w_per_m=getproperty.(points, :right_absorbed_power_w_per_m),
        internal_dissipated_power_w_per_m=getproperty.(points, :internal_dissipated_power_w_per_m),
    )
    println("[+] $path")
end

function run_solve_stage(index_text)
    index = parse(Int, index_text)
    specs = model_specs()
    1 <= index <= length(specs) || error("invalid model index")
    save_response(specs[index].id)
end

function run_solve_parallel_stage()
    stages = ["solve-$index" for index in eachindex(model_specs())]
    semaphore = Base.Semaphore(min(4, length(stages), Sys.CPU_THREADS))
    @sync for stage in stages
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command(stage))
            finally
                Base.release(semaphore)
            end
        end
    end
end

circular_error(a, b) = mod(a - b + pi, 2pi) - pi

function unwrap_phase(values)
    raw = angle.(values)
    result = Float64[first(raw)]
    for index in 2:length(raw)
        push!(result, last(result) + mod(raw[index] - raw[index - 1] + pi, 2pi) - pi)
    end
    result
end

function contiguous_bandwidth(amplitude, center_index; threshold=0.8)
    amplitude[center_index] >= threshold || return 0.0
    left = center_index
    right = center_index
    while left > firstindex(amplitude) && amplitude[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(amplitude) && amplitude[right + 1] >= threshold
        right += 1
    end
    FREQUENCIES_HZ[right] - FREQUENCIES_HZ[left]
end

function phase_coverage(phases)
    length(phases) <= 1 && return 0.0
    sorted = sort(mod.(phases, 2pi))
    gaps = [diff(sorted); first(sorted) + 2pi - last(sorted)]
    2pi - maximum(gaps)
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

function load_response(id, reference)
    saved = JLD2.load(response_path(id))
    response = saved["right_displacement"] ./ reference["right_displacement"]
    power_ratio = saved["right_absorbed_power_w_per_m"] ./
                  reference["right_absorbed_power_w_per_m"]
    balance = saved["source_power_w_per_m"] .-
              saved["left_absorbed_power_w_per_m"] .-
              saved["right_absorbed_power_w_per_m"] .-
              saved["internal_dissipated_power_w_per_m"]
    (
        response,
        amplitude=abs.(response),
        phase=unwrap_phase(response),
        power_ratio,
        balance,
        source_power=saved["source_power_w_per_m"],
    )
end

function lens_config()
    LensConfig(
        frequency_hz=TARGET_FREQUENCY_HZ,
        longitudinal_speed_m_per_s=2340.0,
        element_count=15,
        element_width_mm=4.2,
        slot_width_mm=0.6,
        focal_distance_mm=35.0,
        aperture_quadrature_points=9,
    )
end

function select_reduced_lens(rows)
    config = lens_config()
    centers = lens_centers_mm(config)
    absolute_centers = sort(unique(abs.(centers)))
    k = 2pi / (1000.0 * config.longitudinal_speed_m_per_s / config.frequency_hz)
    base_targets = [
        -k * (hypot(config.focal_distance_mm, center) - config.focal_distance_mm)
        for center in absolute_centers
    ]
    eligible = findall(row -> row.target_amplitude >= 0.45, rows)
    isempty(eligible) && (eligible = collect(eachindex(rows)))
    best = nothing
    for offset in range(-pi, pi; length=721)
        indices = Int[]
        for target in base_targets
            scores = [
                circular_error(rows[index].target_phase_rad, target + offset)^2 +
                0.35 * max(0.0, 0.9 - rows[index].target_amplitude)^2
                for index in eligible
            ]
            push!(indices, eligible[argmin(scores)])
        end
        entries = [
            begin
                group = findfirst(==(abs(center)), absolute_centers)
                row = rows[indices[group]]
                LibraryEntry(
                    row.case_id,
                    0,
                    row.series_height_mm,
                    row.shunt_mass_length_mm,
                    0.0,
                    complex(row.target_H_real, row.target_H_imag),
                )
            end
            for center in centers
        ]
        selection = LensSelection(centers, entries, 0.0 + 0.0im)
        focus = field_at(config.focal_distance_mm, 0.0, selection; config)
        errors = [
            circular_error(rows[indices[group]].target_phase_rad, base_targets[group] + offset)
            for group in eachindex(absolute_centers)
        ]
        candidate = (
            gain=abs(focus),
            offset,
            indices=copy(indices),
            selection=LensSelection(centers, entries, focus),
            rms_phase_error=sqrt(sum(abs2, errors) / length(errors)),
            minimum_amplitude=minimum(rows[index].target_amplitude for index in indices),
            errors,
            targets=base_targets .+ offset,
        )
        if isnothing(best) || candidate.gain > best.gain
            best = candidate
        end
    end
    best, absolute_centers
end

function save_analysis_plot(path, balanced_rows, balanced_spectra, series_points, shunt_points, selected_ids)
    complex_plane = scatter(
        real.(getproperty.(balanced_rows, :target_transfer)),
        imag.(getproperty.(balanced_rows, :target_transfer));
        marker_z=getproperty.(balanced_rows, :target_amplitude),
        color=:viridis,
        colorbar_title="|H|",
        xlabel="Re H",
        ylabel="Im H",
        title="Balanced cells at 242 kHz, s=$(DAMPING_SCALE)",
        aspect_ratio=:equal,
        label="balanced",
        gridalpha=0.25,
    )
    scatter!(complex_plane, real.(series_points), imag.(series_points);
             marker=:square, color=:royalblue, label="series only")
    scatter!(complex_plane, real.(shunt_points), imag.(shunt_points);
             marker=:utriangle, color=:darkorange, label="shunt only")
    circle = range(0, 2pi; length=241)
    plot!(complex_plane, cos.(circle), sin.(circle); color=:gray, linestyle=:dash, label=false)

    spectra_plot = plot(
        xlabel="frequency, kHz",
        ylabel="|H_cal|",
        title="Balanced-cell spectra",
        gridalpha=0.25,
        legend=false,
    )
    for (row, amplitude) in zip(balanced_rows, balanced_spectra)
        selected = row.case_id in selected_ids
        plot!(spectra_plot, FREQUENCIES_HZ ./ 1e3, amplitude;
              color=selected ? :firebrick : :gray70,
              alpha=selected ? 0.95 : 0.4,
              linewidth=selected ? 2.1 : 0.9,
              label=false)
    end
    hline!(spectra_plot, [0.8]; color=:black, linestyle=:dash, label=false)

    target_degrees = [-143.2716, -24.4423, 83.0705, 176.8089,
                      -105.9786, -48.1282, -12.1972, 0.0]
    phase_amplitude = scatter(
        rad2deg.(getproperty.(balanced_rows, :target_phase_rad)),
        getproperty.(balanced_rows, :target_amplitude);
        marker_z=getproperty.(balanced_rows, :bandwidth_0p8_hz) ./ 1e3,
        color=:plasma,
        colorbar_title="BW |H|>=0.8, kHz",
        xlabel="phase at 242 kHz, deg",
        ylabel="|H_cal|",
        title="Phase--amplitude coverage",
        label=false,
        gridalpha=0.25,
    )
    vline!(phase_amplitude, target_degrees; color=:gray70, alpha=0.45,
           linestyle=:dot, label=false)
    hline!(phase_amplitude, [0.8]; color=:black, linestyle=:dash, label=false)

    robust_plot = scatter(
        rad2deg.(getproperty.(balanced_rows, :target_phase_rad)),
        getproperty.(balanced_rows, :minimum_amplitude_234_250khz);
        marker_z=getproperty.(balanced_rows, :target_amplitude),
        color=:viridis,
        colorbar_title="|H(242)|",
        xlabel="phase at 242 kHz, deg",
        ylabel="min |H|, 234--250 kHz",
        title="Useful-band robustness",
        label=false,
        gridalpha=0.25,
    )
    hline!(robust_plot, [0.8]; color=:black, linestyle=:dash, label=false)
    savefig(plot(complex_plane, spectra_plot, phase_amplitude, robust_plot;
                 layout=(2, 2), size=(1300, 1000), margin=5Plots.mm), path)
end

function save_lens_artifacts(best, rows, absolute_centers)
    config = lens_config()
    x_mm = collect(2.0:0.5:90.0)
    y_mm = collect(-45.0:0.5:45.0)
    amplitude = abs.(field_grid(x_mm, y_mm, best.selection; config))
    heat = heatmap(
        x_mm,
        y_mm,
        amplitude;
        xlabel="distance after lens x, mm",
        ylabel="transverse coordinate y, mm",
        title="Reduced lens from balanced physical library",
        color=:viridis,
        colorbar_title="|u| / |u_inc|",
        aspect_ratio=:equal,
        size=(1150, 650),
    )
    scatter!(heat, [config.focal_distance_mm], [0.0]; marker=:xcross,
             markersize=9, markerstrokewidth=3, color=:white, label="target")
    savefig(heat, joinpath(OUTPUT_ROOT, "balanced_huygens_lens_field.png"))

    selected = [
        begin
            row = rows[best.indices[group]]
            (
                absolute_center_y_mm=absolute_centers[group],
                case_id=row.case_id,
                target_phase_rad=best.targets[group],
                selected_phase_rad=row.target_phase_rad,
                phase_error_rad=best.errors[group],
                amplitude=row.target_amplitude,
                power_ratio=row.target_power_ratio,
            )
        end
        for group in eachindex(absolute_centers)
    ]
    write_csv(joinpath(OUTPUT_ROOT, "selected_balanced_cells.csv"), selected)
end

function run_analysis_stage()
    reference = JLD2.load(response_path("reference"))
    target_index = findfirst(==(TARGET_FREQUENCY_HZ), FREQUENCIES_HZ)
    band_indices = findall(frequency -> 234.0e3 <= frequency <= 250.0e3, FREQUENCIES_HZ)
    balanced_rows = NamedTuple[]
    balanced_spectra = Vector{Float64}[]
    response_rows = NamedTuple[]
    ablation_rows = NamedTuple[]

    series_points = ComplexF64[]
    for height in SERIES_HEIGHTS_MM
        loaded = load_response(series_id(height), reference)
        push!(series_points, loaded.response[target_index])
        push!(ablation_rows, (
            variant="series",
            parameter_mm=height,
            H_real=real(loaded.response[target_index]),
            H_imag=imag(loaded.response[target_index]),
            amplitude=loaded.amplitude[target_index],
            phase_rad=angle(loaded.response[target_index]),
            power_ratio=loaded.power_ratio[target_index],
        ))
    end
    shunt_points = ComplexF64[]
    for length_mm in SHUNT_LENGTHS_MM
        loaded = load_response(shunt_id(length_mm), reference)
        push!(shunt_points, loaded.response[target_index])
        push!(ablation_rows, (
            variant="shunt",
            parameter_mm=length_mm,
            H_real=real(loaded.response[target_index]),
            H_imag=imag(loaded.response[target_index]),
            amplitude=loaded.amplitude[target_index],
            phase_rad=angle(loaded.response[target_index]),
            power_ratio=loaded.power_ratio[target_index],
        ))
    end

    for case in balanced_cases()
        loaded = load_response(case.id, reference)
        push!(balanced_spectra, loaded.amplitude)
        push!(balanced_rows, (
            case_id=case.id,
            series_height_mm=case.series_height_mm,
            shunt_mass_length_mm=case.shunt_mass_length_mm,
            target_transfer=loaded.response[target_index],
            target_H_real=real(loaded.response[target_index]),
            target_H_imag=imag(loaded.response[target_index]),
            target_amplitude=loaded.amplitude[target_index],
            target_phase_rad=angle(loaded.response[target_index]),
            target_power_ratio=loaded.power_ratio[target_index],
            minimum_amplitude_234_250khz=minimum(loaded.amplitude[band_indices]),
            bandwidth_0p8_hz=contiguous_bandwidth(loaded.amplitude, target_index),
            max_relative_power_balance_error=maximum(abs.(loaded.balance) ./
                max.(abs.(loaded.source_power), eps())),
        ))
        for index in eachindex(FREQUENCIES_HZ)
            push!(response_rows, (
                case_id=case.id,
                frequency_hz=FREQUENCIES_HZ[index],
                amplitude=loaded.amplitude[index],
                phase_rad=loaded.phase[index],
                power_ratio=loaded.power_ratio[index],
            ))
        end
    end

    best, absolute_centers = select_reduced_lens(balanced_rows)
    selected_ids = Set(balanced_rows[index].case_id for index in best.indices)
    mkpath(OUTPUT_ROOT)
    serializable_rows = [Base.structdiff(row, NamedTuple{(:target_transfer,)}) for row in balanced_rows]
    write_csv(joinpath(OUTPUT_ROOT, "balanced_huygens_library.csv"), serializable_rows)
    write_csv(joinpath(OUTPUT_ROOT, "balanced_huygens_response.csv"), response_rows)
    write_csv(joinpath(OUTPUT_ROOT, "balanced_huygens_ablation_242khz.csv"), ablation_rows)
    save_analysis_plot(
        joinpath(OUTPUT_ROOT, "balanced_huygens_library.png"),
        balanced_rows,
        balanced_spectra,
        series_points,
        shunt_points,
        selected_ids,
    )
    save_lens_artifacts(best, balanced_rows, absolute_centers)

    accepted_phases = [row.target_phase_rad for row in balanced_rows if row.target_amplitude >= 0.8]
    summary = (
        damping_scale=DAMPING_SCALE,
        balanced_case_count=length(balanced_rows),
        maximum_target_amplitude=maximum(getproperty.(balanced_rows, :target_amplitude)),
        phase_coverage_amplitude_0p8_rad=phase_coverage(accepted_phases),
        reduced_lens_gain=best.gain,
        reduced_lens_rms_phase_error_rad=best.rms_phase_error,
        reduced_lens_minimum_amplitude=best.minimum_amplitude,
    )
    write_csv(joinpath(OUTPUT_ROOT, "balanced_huygens_summary.csv"), [summary])
    println(summary)
    println("Selected case ids: $(join(sort(collect(selected_ids)), ", "))")
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_solve_parallel_stage()
    elseif startswith(stage, "solve-")
        run_solve_stage(stage[7:end])
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown balanced-Huygens stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Balanced Huygens stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_balanced_huygens_pilot.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
