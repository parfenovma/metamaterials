module CombinedKMKOpenCell

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_KMK_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "combined_kmk_open_cell"),
)
const MAX_WORKERS = parse(Int, get(ENV, "METAMATERIALS_KMK_WORKERS", "3"))
const FRAME_COUNT = parse(Int, get(ENV, "METAMATERIALS_KMK_FRAMES", "24"))
const FRAME_RATE = parse(Int, get(ENV, "METAMATERIALS_KMK_FPS", "12"))
const VARIANTS = (:solid, :k_only, :m_only, :combined)
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "combined_kmk_mesher.jl"))
using .CombinedKMKMesher

is_solver_stage(stage) = stage in ("coarse", "refine", "validate") ||
                         startswith(stage, "coarse-") || startswith(stage, "refine-")

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) &&
       (is_solver_stage(REQUESTED_STAGE) || REQUESTED_STAGE in ("analyze", "field"))
    include(joinpath(@__DIR__, "modal_harmonic_solver.jl"))
    include(joinpath(@__DIR__, "modal_projection.jl"))
    using .ModalHarmonicElasticity
    using .ElasticModalProjection
    using Gridap
    using JLD2
end

if REQUESTED_STAGE in ("analyze", "field")
    ENV["GKSwstype"] = "100"
    using Plots
end

if REQUESTED_STAGE == "field"
    using FFMPEG
end

const CONFIG = CombinedKMKConfig(
    lead_length_mm=parse(Float64, get(ENV, "METAMATERIALS_KMK_LEAD_MM", "12.0")),
)

mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")
model_path(variant) = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
point_path(stage, variant, frequency_hz) = joinpath(
    OUTPUT_ROOT,
    String(stage),
    "$(variant)_$(round(Int, frequency_hz))hz.jld2",
)

coarse_frequencies_hz() = Float64[
    220.0e3, 225.0e3, 230.0e3, 235.0e3, 240.0e3, 242.0e3,
    244.0e3, 246.0e3, 248.0e3, 250.0e3, 252.0e3,
]
coarse_specs() = [(variant=:combined, frequency_hz=f) for f in coarse_frequencies_hz()]

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            build_combined_kmk_mesh(
                mesh_path(variant);
                config=CONFIG,
                variant,
            )
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    for variant in VARIANTS
        convert_mesh(mesh_path(variant); output_dir=dirname(model_path(variant)))
    end
end

function project_all(state, tag)
    right_modes = filter(
        mode -> mode.kind == :propagating && mode.direction == :right,
        state.modes,
    )
    [
        begin
            left_candidates = filter(
                mode -> mode.kind == :propagating &&
                        mode.direction == :left &&
                        mode.parity == right_mode.parity,
                state.modes,
            )
            left_mode = first(sort(left_candidates; by=mode -> abs(
                mode.wavenumber_per_m + right_mode.wavenumber_per_m,
            )))
            projection = project_boundary_mode_pair(
                state,
                tag,
                right_mode,
                left_mode;
                quadrature_degree=6,
            )
            (; right_mode, left_mode, projection...)
        end
        for right_mode in right_modes
    ]
end

function solve_scattering(variant, frequency_hz; return_state=false)
    harmonic_config = ModalHarmonicElasticity.HarmonicElasticity.HarmonicConfig(
        rayleigh_alpha=0.0,
        rayleigh_beta=0.0,
        element_order=2,
        quadrature_degree=4,
    )
    state = solve_symmetry_reduced_modal_state(
        model_path(variant),
        frequency_hz;
        config=harmonic_config,
        port_height_m=CONFIG.stiffness.height_mm * 1e-3,
        port_element_count=80,
    )
    left = project_all(state, "Source")
    right = project_all(state, "Microphone")
    input_index = findfirst(pair -> pair.right_mode.parity == :symmetric, left)
    input_index === nothing && error("symmetric input mode was not found")
    input_left = left[input_index]
    input_right = right[input_index]
    incident = abs2(input_left.right_amplitude)
    reflection = abs2(input_left.left_amplitude) / incident
    transmission = abs2(input_right.right_amplitude) / incident
    right_incoming = abs2(input_right.left_amplitude) / incident
    conversion = sum(
        abs2(left_pair.left_amplitude) + abs2(right_pair.right_amplitude)
        for (index, (left_pair, right_pair)) in enumerate(zip(left, right))
        if index != input_index
    ) / incident
    transmission_amplitude = input_right.right_amplitude / input_left.right_amplitude
    length_m = total_length_mm(CONFIG) * 1e-3
    wavenumber = real(state.right_mode.wavenumber_per_m)
    deembedded = transmission_amplitude * exp(im * wavenumber * length_m)
    row = (
        variant,
        frequency_hz=Float64(frequency_hz),
        incident_power=incident,
        T00=transmission,
        R00=reflection,
        C0=conversion,
        right_incoming,
        amplitude=abs(transmission_amplitude),
        phase_rad=angle(deembedded),
        balance=transmission + reflection + conversion,
    )
    return_state ? (; row, state) : row
end

function solve_and_save(stage, spec)
    row = solve_scattering(spec.variant, spec.frequency_hz)
    path = point_path(stage, spec.variant, spec.frequency_hz)
    mkpath(dirname(path))
    JLD2.jldsave(path; row)
    println(
        "[+] $stage $(spec.variant) ", round(spec.frequency_hz / 1e3; digits=2),
        " kHz: T=", round(row.T00; digits=4),
        ", R=", round(row.R00; digits=4),
        ", phase=", round(row.phase_rad; digits=3),
    )
end

function load_stage_rows(stage, specs)
    [JLD2.load(point_path(stage, spec.variant, spec.frequency_hz))["row"] for spec in specs]
end

function feature_frequency_hz()
    rows = load_stage_rows(:coarse, coarse_specs())
    rows[argmax(row.T00 for row in rows)].frequency_hz
end

function refine_frequencies_hz()
    center = feature_frequency_hz()
    lower = max(220.0e3, center - 8.0e3)
    # Above approximately 252.8 kHz two more y-even port branches propagate;
    # the current rank-one boundary is intentionally not used there.
    upper = min(246.0e3, center + 8.0e3)
    sort(unique(vcat(
        collect(lower:1.0e3:upper),
        242.0e3,
    )))
end

function refine_specs()
    frequencies = refine_frequencies_hz()
    specs = NamedTuple[]
    for variant in (:combined, :m_only, :k_only), frequency_hz in frequencies
        push!(specs, (; variant, frequency_hz))
    end
    push!(specs, (variant=:solid, frequency_hz=242.0e3))
    specs
end

function run_parallel_specs(stage, specs)
    indices = [
        index for index in eachindex(specs)
        if !isfile(point_path(Symbol(stage), specs[index].variant, specs[index].frequency_hz))
    ]
    isempty(indices) && return
    for batch_start in 1:MAX_WORKERS:length(indices)
        batch = indices[batch_start:min(batch_start + MAX_WORKERS - 1, end)]
        @sync for index in batch
            @async run(child_command("$stage-$index"))
        end
    end
end

function run_validation_stage()
    frequencies = (242.0e3, 245.0e3, 246.0e3)
    rows = [solve_scattering(:combined, frequency_hz) for frequency_hz in frequencies]
    path = joinpath(OUTPUT_ROOT, "port_length_validation.csv")
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "lead_mm,frequency_hz,T00,R00,C0,balance,deembedded_phase_rad")
        for row in rows
            println(
                io,
                join((
                    CONFIG.lead_length_mm,
                    row.frequency_hz,
                    row.T00,
                    row.R00,
                    row.C0,
                    row.balance,
                    row.phase_rad,
                ), ','),
            )
        end
    end
    for row in rows
        println(
            "[+] lead=", CONFIG.lead_length_mm,
            " mm, f=", row.frequency_hz / 1e3,
            " kHz: T=", round(row.T00; digits=5),
            ", R=", round(row.R00; digits=5),
            ", balance=", round(row.balance; digits=5),
            ", phase=", round(row.phase_rad; digits=4),
        )
    end
    println("[+] $path")
end

function unwrap_phase(values)
    result = Float64[first(values)]
    for raw in Iterators.drop(values, 1)
        push!(result, last(result) + mod(raw - last(result) + pi, 2pi) - pi)
    end
    result
end

function add_phase_and_delay(rows)
    sorted = sort(rows; by=row -> row.frequency_hz)
    phase = unwrap_phase(getproperty.(sorted, :phase_rad))
    omega = 2pi .* getproperty.(sorted, :frequency_hz)
    delay = similar(phase)
    delay[1] = -(phase[2] - phase[1]) / (omega[2] - omega[1])
    for index in 2:(length(phase) - 1)
        delay[index] = -(phase[index + 1] - phase[index - 1]) /
                       (omega[index + 1] - omega[index - 1])
    end
    delay[end] = -(phase[end] - phase[end - 1]) / (omega[end] - omega[end - 1])
    [merge(row, (phase_rad=phase[index], group_delay_s=delay[index]))
     for (index, row) in enumerate(sorted)]
end

function analyzed_rows()
    raw = load_stage_rows(:refine, refine_specs())
    rows = NamedTuple[]
    for variant in (:combined, :m_only, :k_only)
        append!(rows, add_phase_and_delay(filter(row -> row.variant == variant, raw)))
    end
    append!(rows, filter(row -> row.variant == :solid, raw))
    rows
end

function write_csv(rows)
    path = joinpath(OUTPUT_ROOT, "combined_kmk_modal_scattering.csv")
    open(path, "w") do io
        println(io, "variant,frequency_hz,T00,R00,C0,right_incoming,amplitude,deembedded_phase_rad,group_delay_s,power_sum")
        for row in rows
            println(io, join((
                row.variant,
                row.frequency_hz,
                row.T00,
                row.R00,
                row.C0,
                row.right_incoming,
                row.amplitude,
                row.phase_rad,
                hasproperty(row, :group_delay_s) ? row.group_delay_s : NaN,
                row.balance,
            ), ','))
        end
    end
    println("[+] $path")
end

function variant_rows(rows, variant)
    sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
end

function save_summary_figure(rows)
    combined = variant_rows(rows, :combined)
    frequencies = getproperty.(combined, :frequency_hz) ./ 1e3
    power_panel = plot(
        frequencies,
        getproperty.(combined, :T00);
        linewidth=2,
        marker=:circle,
        label="T00",
        xlabel="frequency, kHz",
        ylabel="power fraction",
        title="Absolute KMK scattering",
        gridalpha=0.25,
    )
    plot!(power_panel, frequencies, getproperty.(combined, :R00);
          linewidth=2, marker=:circle, label="R00")
    plot!(power_panel, frequencies, getproperty.(combined, :C0);
          linewidth=2, marker=:circle, label="C0")
    vline!(power_panel, [242.0]; color=:black, linestyle=:dash, label="242 kHz")

    ablation_panel = plot(
        xlabel="frequency, kHz",
        ylabel="T00",
        title="Ablation of the open cell",
        gridalpha=0.25,
        legend=:bottomright,
        ylims=(0, 1.05),
    )
    styles = Dict(:combined => (:royalblue, :circle), :m_only => (:darkorange, :diamond),
                  :k_only => (:seagreen, :square))
    for variant in (:combined, :m_only, :k_only)
        current = variant_rows(rows, variant)
        color, marker = styles[variant]
        plot!(
            ablation_panel,
            getproperty.(current, :frequency_hz) ./ 1e3,
            getproperty.(current, :T00);
            color,
            marker,
            linewidth=2,
            label=String(variant),
        )
    end
    vline!(ablation_panel, [242.0]; color=:black, linestyle=:dash, label="242 kHz")

    phase_panel = plot(
        xlabel="frequency, kHz",
        ylabel="de-embedded phase, rad",
        title="Transmission phase",
        gridalpha=0.25,
        legend=:topleft,
    )
    for variant in (:combined, :m_only, :k_only)
        current = variant_rows(rows, variant)
        color, marker = styles[variant]
        plot!(
            phase_panel,
            getproperty.(current, :frequency_hz) ./ 1e3,
            getproperty.(current, :phase_rad);
            color,
            marker,
            linewidth=2,
            label=String(variant),
        )
    end
    vline!(phase_panel, [242.0]; color=:black, linestyle=:dash, label="242 kHz")

    delay_values = [
        row.T00 >= 0.05 && abs(row.balance - 1) <= 0.02 ?
        row.group_delay_s * 1e6 : NaN
        for row in combined
    ]
    delay_panel = plot(
        frequencies,
        delay_values;
        linewidth=2,
        marker=:circle,
        label="KMK (T00 >= 0.05)",
        xlabel="frequency, kHz",
        ylabel="excess group delay, μs",
        title="Phase-slope diagnostic",
        gridalpha=0.25,
    )
    vline!(delay_panel, [242.0]; color=:black, linestyle=:dash, label="242 kHz")
    figure = plot(
        power_panel,
        ablation_panel,
        phase_panel,
        delay_panel;
        layout=(2, 2),
        size=(1500, 950),
        margin=5Plots.mm,
        plot_title="Open K/2–M–K/2 cell — lossless modal S-parameters",
    )
    path = joinpath(OUTPUT_ROOT, "combined_kmk_modal_scattering.png")
    savefig(figure, path)
    println("[+] $path")
end

function run_analysis_stage()
    rows = analyzed_rows()
    write_csv(rows)
    save_summary_figure(rows)
    println("\nAt 242 kHz:")
    for variant in (:solid, :k_only, :m_only, :combined)
        current = filter(
            row -> row.variant == variant && row.frequency_hz == 242.0e3,
            rows,
        )
        isempty(current) && continue
        row = only(current)
        println(
            "  $variant: T00=", round(row.T00; digits=4),
            ", R00=", round(row.R00; digits=4),
            ", C0=", round(row.C0; digits=3),
            ", phase=", round(row.phase_rad; digits=3), " rad",
        )
    end
    combined = variant_rows(rows, :combined)
    feature = combined[argmax(row.T00 for row in combined)]
    println(
        "Best combined transmission point: ", feature.frequency_hz / 1e3,
        " kHz, T00=", round(feature.T00; digits=4),
        ", R00=", round(feature.R00; digits=4),
    )
end

function visualization_field(state)
    data = only(Gridap.Visualization.visualization_data(
        state.domain,
        "KMK_field";
        order=1,
        cellfields=Dict("displacement" => state.displacement),
    ))
    coordinates = collect(Gridap.Geometry.get_node_coordinates(data.grid))
    displacement = collect(data.nodaldata["displacement"])
    (
        x_mm=Float64[point[1] * 1e3 for point in coordinates],
        y_mm=Float64[point[2] * 1e3 for point in coordinates],
        ux=ComplexF64[value[1] for value in displacement],
        uy=ComplexF64[value[2] for value in displacement],
    )
end

function field_panel(field, phase; zoom=false)
    instantaneous_x = real.(field.ux .* exp(im * phase))
    instantaneous_y = real.(field.uy .* exp(im * phase))
    normalization = maximum(sqrt.(abs2.(field.ux) .+ abs2.(field.uy)))
    scale_mm = zoom ? 0.18 : 0.10
    x = field.x_mm .+ scale_mm .* instantaneous_x ./ normalization
    y = field.y_mm .+ scale_mm .* instantaneous_y ./ normalization
    cell_start = CONFIG.lead_length_mm
    cell_end = cell_start + cell_length_mm(CONFIG)
    scatter(
        x,
        y;
        marker_z=instantaneous_x ./ normalization,
        color=:balance,
        clims=(-1, 1),
        colorbar=false,
        markersize=zoom ? 2.0 : 1.2,
        markerstrokewidth=0,
        label=false,
        aspect_ratio=:equal,
        xlims=zoom ? (cell_start - 0.8, cell_end + 0.8) : (-0.5, total_length_mm(CONFIG) + 0.5),
        ylims=(-2.6, 2.6),
        xlabel="x, mm",
        ylabel="y, mm",
        title=zoom ? "KMK cell (deformation rescaled)" : "Full modal-port model",
        grid=false,
        framestyle=:box,
    )
end

function encode_animation(frame_dir)
    pattern = joinpath(frame_dir, "frame_%03d.png")
    mp4_path = joinpath(OUTPUT_ROOT, "combined_kmk_field.mp4")
    gif_path = joinpath(OUTPUT_ROOT, "combined_kmk_field.gif")
    FFMPEG.ffmpeg_exe(Cmd(String[
        "-y", "-framerate", string(FRAME_RATE), "-i", pattern,
        "-c:v", "libx264", "-pix_fmt", "yuv420p", mp4_path,
    ]))
    FFMPEG.ffmpeg_exe(Cmd(String[
        "-y", "-framerate", string(FRAME_RATE), "-i", pattern,
        "-vf", "fps=$(FRAME_RATE),scale=1200:-1:flags=lanczos", "-loop", "0", gif_path,
    ]))
    println("[+] $mp4_path")
    println("[+] $gif_path")
end

function run_field_stage()
    rows = analyzed_rows()
    combined = variant_rows(rows, :combined)
    selected = combined[argmax(row.T00 for row in combined)]
    solved = solve_scattering(:combined, selected.frequency_hz; return_state=true)
    field = visualization_field(solved.state)
    frame_dir = joinpath(OUTPUT_ROOT, "field_frames")
    mkpath(frame_dir)
    for frame_index in 1:FRAME_COUNT
        phase = 2pi * (frame_index - 1) / FRAME_COUNT
        figure = plot(
            field_panel(field, phase; zoom=false),
            field_panel(field, phase; zoom=true);
            layout=(2, 1),
            size=(1500, 800),
            margin=4Plots.mm,
            plot_title=string(
                "Open KMK field at ", selected.frequency_hz / 1e3,
                " kHz; T00=", round(selected.T00; digits=3),
                ", R00=", round(selected.R00; digits=3),
            ),
        )
        savefig(figure, joinpath(frame_dir, "frame_$(lpad(frame_index, 3, '0')).png"))
    end
    encode_animation(frame_dir)
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "coarse"
        run_parallel_specs("coarse", coarse_specs())
    elseif startswith(stage, "coarse-")
        index = parse(Int, stage[8:end])
        solve_and_save(:coarse, coarse_specs()[index])
    elseif stage == "refine"
        run_parallel_specs("refine", refine_specs())
    elseif startswith(stage, "refine-")
        index = parse(Int, stage[8:end])
        solve_and_save(:refine, refine_specs()[index])
    elseif stage == "analyze"
        run_analysis_stage()
    elseif stage == "validate"
        run_validation_stage()
    elseif stage == "field"
        run_field_stage()
    else
        error("unknown KMK stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "coarse", "refine", "analyze", "field")
            println("\n=== Open KMK stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_combined_kmk_open_cell.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
