module SideMassModalPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MODAL_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "side_mass_modal_pilot"),
)
const TARGET_FREQUENCY_HZ = parse(
    Float64,
    get(ENV, "METAMATERIALS_MODAL_TARGET_HZ", "242000"),
)
const MODE_COUNT = parse(Int, get(ENV, "METAMATERIALS_MODAL_COUNT", "12"))
const FRAME_COUNT = parse(Int, get(ENV, "METAMATERIALS_MODAL_FRAMES", "24"))
const FRAME_RATE = parse(Int, get(ENV, "METAMATERIALS_MODAL_FPS", "12"))
const VARIANTS = (:b, :bd)
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "side_mass_mesher.jl"))
using .SideMassMesher

const CONFIG = SideMassConfig(bright_neck_width_mm=0.7)

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif REQUESTED_STAGE in ("solve-b", "solve-bd", "analyze")
    include(joinpath(@__DIR__, "modal_solver.jl"))
    using .ConservativeElasticModes
    using JLD2
end

if REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using FFMPEG
    using Plots
end

mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")
model_path(variant) = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
result_path(variant) = joinpath(OUTPUT_ROOT, "modes_$(variant).jld2")
free_result_path(variant) = joinpath(OUTPUT_ROOT, "modes_$(variant)_free_ends.jld2")
summary_path() = joinpath(OUTPUT_ROOT, "modal_summary.csv")

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            build_side_mass_mesh(
                mesh_path(variant);
                config=CONFIG,
                variant,
                size_min_mm=0.10,
                size_max_mm=0.45,
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

function solve_and_save(variant, path; clamp_tags, target_frequency_hz)
    modes = solve_conservative_modes(
        model_path(variant);
        config=ModeConfig(
            target_frequency_hz=target_frequency_hz,
            mode_count=MODE_COUNT,
            tolerance=1.0e-9,
            clamp_tags=clamp_tags,
        ),
    )
    metrics = [side_mass_mode_metrics(mode, CONFIG; variant) for mode in modes]
    mkpath(OUTPUT_ROOT)
    JLD2.jldsave(path; modes, metrics, variant=String(variant))
    println("[+] $path")
    for (index, mode) in enumerate(modes)
        metric = metrics[index]
        println(
            "mode $index: ",
            round(mode.frequency_hz / 1e3; digits=3),
            " kHz, parity=", round(mode.parity_y; digits=3),
            ", bright=", round(metric.bright_fraction; digits=3),
            ", dark=", round(metric.dark_fraction; digits=3),
            ", residual=", mode.relative_residual,
        )
    end
end

function run_solve_stage(variant)
    solve_and_save(
        variant,
        result_path(variant);
        clamp_tags=["Source", "Microphone"],
        target_frequency_hz=TARGET_FREQUENCY_HZ,
    )
    # For B, centre the free-end search on the localized cluster found by the
    # clamped diagnostic; for BD the target resonance itself is the useful test.
    free_target_hz = variant == :b ? 217.0e3 : TARGET_FREQUENCY_HZ
    solve_and_save(
        variant,
        free_result_path(variant);
        clamp_tags=String[],
        target_frequency_hz=free_target_hz,
    )
end

function run_parallel_solve_stage()
    @sync for variant in VARIANTS
        @async run(child_command("solve-$(variant)"))
    end
end

function write_summary(data, free_data)
    path = summary_path()
    open(path, "w") do io
        println(io, "variant,boundary_condition,mode_index,frequency_hz,parity_y,longitudinal_fraction,bright_fraction,dark_fraction,side_mass_fraction,relative_residual")
        for (boundary, current_data) in (("clamped", data), ("free", free_data))
            for variant in VARIANTS
                modes = current_data[variant]["modes"]
                metrics = current_data[variant]["metrics"]
                for (index, mode) in enumerate(modes)
                    metric = metrics[index]
                    println(io, join((
                        uppercase(String(variant)),
                        boundary,
                        index,
                        mode.frequency_hz,
                        mode.parity_y,
                        mode.longitudinal_fraction,
                        metric.bright_fraction,
                        metric.dark_fraction,
                        metric.side_mass_fraction,
                        mode.relative_residual,
                    ), ','))
                end
            end
        end
    end
    println("[+] $path")
end

function save_boundary_sensitivity(data, free_data)
    path = joinpath(OUTPUT_ROOT, "boundary_sensitivity.csv")
    open(path, "w") do io
        println(io, "variant,boundary_condition,frequency_hz,parity_y,bright_fraction,dark_fraction,side_mass_fraction")
        for variant in VARIANTS
            metric_name = variant == :b ? :bright_fraction : :side_mass_fraction
            for (boundary, current_data) in (("clamped", data), ("free", free_data))
                modes = current_data[variant]["modes"]
                metrics = current_data[variant]["metrics"]
                index = select_mode(modes, metrics; parity_sign=1, metric_name)
                mode = modes[index]
                metric = metrics[index]
                println(io, join((
                    uppercase(String(variant)),
                    boundary,
                    mode.frequency_hz,
                    mode.parity_y,
                    metric.bright_fraction,
                    metric.dark_fraction,
                    metric.side_mass_fraction,
                ), ','))
            end
        end
    end
    println("[+] $path")
end

function select_mode(modes, metrics; parity_sign=1, metric_name=:side_mass_fraction)
    candidates = findall(eachindex(modes)) do index
        parity_sign * modes[index].parity_y >= 0.8
    end
    isempty(candidates) && error("no mode with requested parity")
    candidates[argmax(getproperty(metrics[index], metric_name) for index in candidates)]
end

function selected_triplet(data)
    modes_b = data[:b]["modes"]
    metrics_b = data[:b]["metrics"]
    modes_bd = data[:bd]["modes"]
    metrics_bd = data[:bd]["metrics"]
    index_b = select_mode(modes_b, metrics_b; parity_sign=1, metric_name=:bright_fraction)
    index_bd_even = select_mode(modes_bd, metrics_bd; parity_sign=1)
    index_bd_odd = select_mode(modes_bd, metrics_bd; parity_sign=-1)
    [
        (label="B: most bright-localized even mode", variant=:b,
         mode=modes_b[index_b], metric=metrics_b[index_b]),
        (label="BD: most localized even mode", variant=:bd,
         mode=modes_bd[index_bd_even], metric=metrics_bd[index_bd_even]),
        (label="BD: most localized odd partner", variant=:bd,
         mode=modes_bd[index_bd_odd], metric=metrics_bd[index_bd_odd]),
    ]
end

function spectrum_panel(data)
    panel = plot(
        xlabel="frequency, kHz",
        ylabel="side-mass kinetic fraction",
        title="Clamped-end conservative modes near 242 kHz",
        xlims=(200, 275),
        ylims=(0, 0.48),
        gridalpha=0.25,
        legend=:topleft,
    )
    styles = Dict(:b => (:darkorange, :circle), :bd => (:royalblue, :diamond))
    for variant in VARIANTS
        modes = data[variant]["modes"]
        metrics = data[variant]["metrics"]
        color, marker = styles[variant]
        even = [mode.parity_y >= 0 for mode in modes]
        scatter!(
            panel,
            getproperty.(modes, :frequency_hz) ./ 1e3,
            getproperty.(metrics, :side_mass_fraction);
            marker,
            markercolor=[flag ? color : :transparent for flag in even],
            markerstrokecolor=color,
            markersize=7,
            label="$(uppercase(String(variant))) (filled: y-even)",
        )
    end
    vline!(panel, [TARGET_FREQUENCY_HZ / 1e3]; color=:black, linestyle=:dash, label="target")
    panel
end

function mode_panel(selected, phase)
    mode = selected.mode
    x_mm = mode.node_x_m .* 1e3
    y_mm = mode.node_y_m .* 1e3
    magnitude = sqrt.(mode.displacement_x .^ 2 .+ mode.displacement_y .^ 2)
    normalization = maximum(magnitude)
    ux = mode.displacement_x ./ normalization
    uy = mode.displacement_y ./ normalization
    oscillation = cos(phase)
    deformation_mm = 0.45
    title = string(
        selected.label,
        "\n",
        round(mode.frequency_hz / 1e3; digits=3),
        " kHz; parity=", round(mode.parity_y; digits=2),
        "; B=", round(selected.metric.bright_fraction; digits=2),
        "; D=", round(selected.metric.dark_fraction; digits=2),
    )
    scatter(
        x_mm .+ deformation_mm .* oscillation .* ux,
        y_mm .+ deformation_mm .* oscillation .* uy;
        marker_z=oscillation .* ux,
        color=:balance,
        clims=(-1, 1),
        colorbar=false,
        markersize=1.7,
        markerstrokewidth=0,
        label=false,
        aspect_ratio=:equal,
        xlims=(-0.5, 30.5),
        ylims=(-0.7, 7.3),
        xlabel="x, mm",
        ylabel="y, mm",
        title,
        grid=false,
        framestyle=:box,
    )
end

function save_static_figure(data, selected)
    panels = [spectrum_panel(data); [mode_panel(item, 0.0) for item in selected]]
    figure = plot(
        panels...;
        layout=(4, 1),
        size=(1450, 1450),
        margin=4Plots.mm,
        plot_title="Modal pilot — modes are M-normalized; deformation is rescaled for visibility",
    )
    path = joinpath(OUTPUT_ROOT, "modal_pilot.png")
    savefig(figure, path)
    println("[+] $path")
end

function encode_animation(frame_dir)
    input_pattern = joinpath(frame_dir, "frame_%03d.png")
    mp4_path = joinpath(OUTPUT_ROOT, "selected_eigenmodes.mp4")
    gif_path = joinpath(OUTPUT_ROOT, "selected_eigenmodes.gif")
    FFMPEG.ffmpeg_exe(Cmd(String[
        "-y", "-framerate", string(FRAME_RATE), "-i", input_pattern,
        "-c:v", "libx264", "-pix_fmt", "yuv420p", mp4_path,
    ]))
    FFMPEG.ffmpeg_exe(Cmd(String[
        "-y", "-framerate", string(FRAME_RATE), "-i", input_pattern,
        "-vf", "fps=$(FRAME_RATE),scale=1200:-1:flags=lanczos", "-loop", "0", gif_path,
    ]))
    println("[+] $mp4_path")
    println("[+] $gif_path")
end

function save_animation(selected)
    frame_dir = joinpath(OUTPUT_ROOT, "animation_frames")
    mkpath(frame_dir)
    for frame_index in 1:FRAME_COUNT
        phase = 2pi * (frame_index - 1) / FRAME_COUNT
        panels = [mode_panel(item, phase) for item in selected]
        figure = plot(
            panels...;
            layout=(3, 1),
            size=(1450, 1000),
            margin=4Plots.mm,
            plot_title="Standing eigenmodes; displayed deformation is independently rescaled",
        )
        path = joinpath(frame_dir, "frame_$(lpad(frame_index, 3, '0')).png")
        savefig(figure, path)
        println("[frame $frame_index/$FRAME_COUNT] $path")
    end
    encode_animation(frame_dir)
end

function run_analysis_stage()
    data = Dict(variant => JLD2.load(result_path(variant)) for variant in VARIANTS)
    free_data = Dict(variant => JLD2.load(free_result_path(variant)) for variant in VARIANTS)
    write_summary(data, free_data)
    save_boundary_sensitivity(data, free_data)
    selected = selected_triplet(data)
    save_static_figure(data, selected)
    save_animation(selected)
    println("\nSelected modes:")
    for item in selected
        println(
            "  ", item.label, ": ", round(item.mode.frequency_hz / 1e3; digits=3),
            " kHz, parity=", round(item.mode.parity_y; digits=3),
            ", B=", round(item.metric.bright_fraction; digits=3),
            ", D=", round(item.metric.dark_fraction; digits=3),
        )
    end
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_parallel_solve_stage()
    elseif stage == "solve-b"
        run_solve_stage(:b)
    elseif stage == "solve-bd"
        run_solve_stage(:bd)
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown modal-pilot stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Modal pilot stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_side_mass_modal_pilot.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
