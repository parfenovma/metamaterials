module HuygensPairPilot

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HUYGENS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "huygens_pair_pilot_242khz"),
)
const DAMPING_SCALE = parse(Float64, get(ENV, "METAMATERIALS_HUYGENS_DAMPING_SCALE", "0.25"))
const MESH_SIZE_MIN_MM = parse(Float64, get(ENV, "METAMATERIALS_HUYGENS_MESH_MIN_MM", "0.06"))
const MESH_SIZE_MAX_MM = parse(Float64, get(ENV, "METAMATERIALS_HUYGENS_MESH_MAX_MM", "0.28"))
const FREQUENCIES_HZ = collect(
    parse(Float64, get(ENV, "METAMATERIALS_HUYGENS_START_HZ", "220000")):
    parse(Float64, get(ENV, "METAMATERIALS_HUYGENS_STEP_HZ", "2000")):
    parse(Float64, get(ENV, "METAMATERIALS_HUYGENS_STOP_HZ", "264000")),
)
const TARGET_FREQUENCY_HZ = 242.0e3
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "huygens_pair_mesher.jl"))
using .HuygensPairMesher

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

struct PairCase
    id::String
    config::HuygensPairConfig
end

number_label(value) = replace(@sprintf("%.1f", Float64(value)), "." => "p")

function pair_cases()
    cases = PairCase[]
    for front_length in (2.4, 3.0, 3.6, 4.2), rear_length in (2.4, 3.0, 3.6, 4.2)
        push!(cases, PairCase(
            "lf$(number_label(front_length))_lr$(number_label(rear_length))_nf0p6_nr0p6",
            HuygensPairConfig(
                front_mass_length_mm=front_length,
                rear_mass_length_mm=rear_length,
            ),
        ))
    end
    for (front_neck, rear_neck) in ((0.4, 0.4), (0.4, 0.8), (0.8, 0.4), (0.8, 0.8))
        push!(cases, PairCase(
            "lf3p0_lr3p0_nf$(number_label(front_neck))_nr$(number_label(rear_neck))",
            HuygensPairConfig(
                front_neck_width_mm=front_neck,
                rear_neck_width_mm=rear_neck,
            ),
        ))
    end
    cases
end

mesh_path(id) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(id).msh")
model_path(id) = joinpath(OUTPUT_ROOT, "models", "model_$(id).json")
response_path(id) = joinpath(OUTPUT_ROOT, "harmonic", "response_$(id).jld2")

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    cases = pair_cases()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        build_huygens_pair_mesh(
            mesh_path("reference");
            config=first(cases).config,
            variant=:r0,
            size_min_mm=MESH_SIZE_MIN_MM,
            size_max_mm=MESH_SIZE_MAX_MM,
        )
        for case in cases
            build_huygens_pair_mesh(
                mesh_path(case.id);
                config=case.config,
                variant=:pair,
                size_min_mm=MESH_SIZE_MIN_MM,
                size_max_mm=MESH_SIZE_MAX_MM,
            )
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    convert_mesh(mesh_path("reference"); output_dir=joinpath(OUTPUT_ROOT, "models"))
    for case in pair_cases()
        convert_mesh(mesh_path(case.id); output_dir=joinpath(OUTPUT_ROOT, "models"))
    end
end

function harmonic_config()
    HarmonicConfig(
        rayleigh_alpha=DAMPING_SCALE * 79560.0,
        rayleigh_beta=DAMPING_SCALE * 2.5e-9,
    )
end

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

function run_solve_stage(specification)
    if specification == "reference"
        save_response("reference")
        return
    end
    index = parse(Int, specification)
    cases = pair_cases()
    1 <= index <= length(cases) || error("invalid case index")
    save_response(cases[index].id)
end

function run_solve_parallel_stage()
    stages = ["solve-reference"; ["solve-$index" for index in eachindex(pair_cases())]]
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

function unwrap_phase(values)
    raw = angle.(values)
    result = Float64[first(raw)]
    for index in 2:length(raw)
        push!(result, last(result) + mod(raw[index] - raw[index - 1] + pi, 2pi) - pi)
    end
    result
end

circular_error(a, b) = mod(a - b + pi, 2pi) - pi

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
    isempty(phases) && return 0.0
    length(phases) == 1 && return 0.0
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
    eligible = findall(row -> row.target_amplitude >= 0.55, rows)
    isempty(eligible) && (eligible = eachindex(rows))
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
                LibraryEntry(row.case_id, 0, row.front_mass_length_mm,
                             row.rear_mass_length_mm, 0.0,
                             complex(row.target_H_real, row.target_H_imag))
            end
            for center in centers
        ]
        selection = LensSelection(centers, entries, 0.0 + 0.0im)
        focus = field_at(config.focal_distance_mm, 0.0, selection; config)
        phase_errors = [
            circular_error(rows[indices[group]].target_phase_rad,
                           base_targets[group] + offset)
            for group in eachindex(absolute_centers)
        ]
        candidate = (
            gain=abs(focus),
            offset=offset,
            indices=copy(indices),
            selection=LensSelection(centers, entries, focus),
            rms_phase_error=sqrt(sum(abs2, phase_errors) / length(phase_errors)),
            minimum_amplitude=minimum(rows[index].target_amplitude for index in indices),
            phase_errors=phase_errors,
            targets=base_targets .+ offset,
        )
        if isnothing(best) || candidate.gain > best.gain
            best = candidate
        end
    end
    best, absolute_centers
end

function save_library_plot(path, rows, spectra, selected_ids)
    target_phases_deg = [-143.2716, -24.4423, 83.0705, 176.8089,
                         -105.9786, -48.1282, -12.1972, 0.0]
    complex_plane = scatter(
        getproperty.(rows, :target_H_real),
        getproperty.(rows, :target_H_imag);
        marker_z=getproperty.(rows, :target_amplitude),
        color=:viridis,
        colorbar_title="|H|",
        xlabel="Re H",
        ylabel="Im H",
        title="Physical H(242 kHz), s=$(DAMPING_SCALE)",
        aspect_ratio=:equal,
        label=false,
        gridalpha=0.25,
    )
    circle = range(0, 2pi; length=241)
    plot!(complex_plane, cos.(circle), sin.(circle); color=:gray, linestyle=:dash, label=false)

    spectra_plot = plot(
        xlabel="frequency, kHz",
        ylabel="|H_cal|",
        title="Broadband response of pair cells",
        gridalpha=0.25,
        legend=false,
    )
    for (row, amplitude) in zip(rows, spectra)
        selected = row.case_id in selected_ids
        plot!(spectra_plot, FREQUENCIES_HZ ./ 1e3, amplitude;
              color=selected ? :firebrick : :gray70,
              alpha=selected ? 0.9 : 0.35,
              linewidth=selected ? 2.0 : 0.8,
              label=false)
    end
    hline!(spectra_plot, [0.8]; color=:black, linestyle=:dash, label=false)

    phase_amplitude = scatter(
        rad2deg.(getproperty.(rows, :target_phase_rad)),
        getproperty.(rows, :target_amplitude);
        marker_z=getproperty.(rows, :bandwidth_0p8_hz) ./ 1e3,
        color=:plasma,
        colorbar_title="BW |H|>=0.8, kHz",
        xlabel="phase at 242 kHz, deg",
        ylabel="|H_cal|",
        title="Phase--amplitude library",
        label=false,
        gridalpha=0.25,
    )
    vline!(phase_amplitude, target_phases_deg; color=:gray70, alpha=0.45,
           linestyle=:dot, label=false)
    hline!(phase_amplitude, [0.8]; color=:black, linestyle=:dash, label=false)

    band_plot = scatter(
        rad2deg.(getproperty.(rows, :target_phase_rad)),
        getproperty.(rows, :minimum_amplitude_234_250khz);
        marker_z=getproperty.(rows, :target_amplitude),
        color=:viridis,
        colorbar_title="|H(242)|",
        xlabel="phase at 242 kHz, deg",
        ylabel="min |H|, 234--250 kHz",
        title="Useful-band robustness",
        label=false,
        gridalpha=0.25,
    )
    hline!(band_plot, [0.8]; color=:black, linestyle=:dash, label=false)
    savefig(plot(complex_plane, spectra_plot, phase_amplitude, band_plot;
                 layout=(2, 2), size=(1300, 1000), margin=5Plots.mm), path)
end

function save_lens_artifacts(output_root, best, rows, absolute_centers)
    config = lens_config()
    x_mm = collect(2.0:0.5:90.0)
    y_mm = collect(-45.0:0.5:45.0)
    field = field_grid(x_mm, y_mm, best.selection; config)
    amplitude = abs.(field)
    heat = heatmap(
        x_mm,
        y_mm,
        amplitude;
        xlabel="distance after lens x, mm",
        ylabel="transverse coordinate y, mm",
        title="Reduced lens from first physical pair library",
        color=:viridis,
        colorbar_title="|u| / |u_inc|",
        aspect_ratio=:equal,
        size=(1150, 650),
    )
    scatter!(heat, [config.focal_distance_mm], [0.0]; marker=:xcross,
             markersize=9, markerstrokewidth=3, color=:white, label="target")
    savefig(heat, joinpath(output_root, "huygens_pair_lens_field.png"))

    selection_rows = [
        begin
            row = rows[best.indices[group]]
            (
                absolute_center_y_mm=absolute_centers[group],
                case_id=row.case_id,
                target_phase_rad=best.targets[group],
                selected_phase_rad=row.target_phase_rad,
                phase_error_rad=best.phase_errors[group],
                amplitude=row.target_amplitude,
                power_ratio=row.target_power_ratio,
            )
        end
        for group in eachindex(absolute_centers)
    ]
    write_csv(joinpath(output_root, "selected_physical_cells.csv"), selection_rows)
end

function run_analysis_stage()
    reference = JLD2.load(response_path("reference"))
    target_index = findfirst(==(TARGET_FREQUENCY_HZ), FREQUENCIES_HZ)
    isnothing(target_index) && error("frequency grid must include 242 kHz")
    band_indices = findall(frequency -> 234.0e3 <= frequency <= 250.0e3, FREQUENCIES_HZ)
    rows = NamedTuple[]
    spectra = Vector{Float64}[]
    response_rows = NamedTuple[]
    for case in pair_cases()
        saved = JLD2.load(response_path(case.id))
        response = saved["right_displacement"] ./ reference["right_displacement"]
        amplitude = abs.(response)
        phase = unwrap_phase(response)
        power_ratio = saved["right_absorbed_power_w_per_m"] ./
                      reference["right_absorbed_power_w_per_m"]
        balance = saved["source_power_w_per_m"] .-
                  saved["left_absorbed_power_w_per_m"] .-
                  saved["right_absorbed_power_w_per_m"] .-
                  saved["internal_dissipated_power_w_per_m"]
        push!(spectra, amplitude)
        push!(rows, (
            case_id=case.id,
            front_mass_length_mm=case.config.front_mass_length_mm,
            rear_mass_length_mm=case.config.rear_mass_length_mm,
            front_neck_width_mm=case.config.front_neck_width_mm,
            rear_neck_width_mm=case.config.rear_neck_width_mm,
            target_H_real=real(response[target_index]),
            target_H_imag=imag(response[target_index]),
            target_amplitude=amplitude[target_index],
            target_phase_rad=angle(response[target_index]),
            target_power_ratio=power_ratio[target_index],
            minimum_amplitude_234_250khz=minimum(amplitude[band_indices]),
            bandwidth_0p8_hz=contiguous_bandwidth(amplitude, target_index),
            max_relative_power_balance_error=maximum(abs.(balance) ./
                max.(abs.(saved["source_power_w_per_m"]), eps())),
        ))
        for index in eachindex(FREQUENCIES_HZ)
            push!(response_rows, (
                case_id=case.id,
                frequency_hz=FREQUENCIES_HZ[index],
                amplitude=amplitude[index],
                phase_rad=phase[index],
                power_ratio=power_ratio[index],
            ))
        end
    end

    best, absolute_centers = select_reduced_lens(rows)
    selected_ids = Set(rows[index].case_id for index in best.indices)
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "huygens_pair_library.csv"), rows)
    write_csv(joinpath(OUTPUT_ROOT, "huygens_pair_response.csv"), response_rows)
    save_library_plot(joinpath(OUTPUT_ROOT, "huygens_pair_library.png"), rows, spectra, selected_ids)
    save_lens_artifacts(OUTPUT_ROOT, best, rows, absolute_centers)
    accepted_phases = [row.target_phase_rad for row in rows if row.target_amplitude >= 0.8]
    summary = (
        damping_scale=DAMPING_SCALE,
        case_count=length(rows),
        maximum_target_amplitude=maximum(getproperty.(rows, :target_amplitude)),
        minimum_target_amplitude=minimum(getproperty.(rows, :target_amplitude)),
        phase_coverage_amplitude_0p8_rad=phase_coverage(accepted_phases),
        reduced_lens_gain=best.gain,
        reduced_lens_rms_phase_error_rad=best.rms_phase_error,
        reduced_lens_minimum_amplitude=best.minimum_amplitude,
    )
    write_csv(joinpath(OUTPUT_ROOT, "huygens_pair_summary.csv"), [summary])
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
        error("unknown Huygens-pair stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Huygens pair stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_huygens_pair_pilot.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
