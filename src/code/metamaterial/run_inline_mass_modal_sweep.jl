module InlineMassModalSweep

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_INLINE_MASS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "inline_mass_modal_sweep"),
)
const TARGET_FREQUENCY_HZ = parse(
    Float64,
    get(ENV, "METAMATERIALS_INLINE_MASS_TARGET_HZ", "250000"),
)
const MODE_COUNT = parse(Int, get(ENV, "METAMATERIALS_INLINE_MASS_MODE_COUNT", "8"))
const MAX_SOLVERS = parse(Int, get(ENV, "METAMATERIALS_INLINE_MASS_WORKERS", "4"))
const FRAME_COUNT = parse(Int, get(ENV, "METAMATERIALS_INLINE_MASS_FRAMES", "24"))
const FRAME_RATE = parse(Int, get(ENV, "METAMATERIALS_INLINE_MASS_FPS", "12"))
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

starts_with_solve(stage) = startswith(stage, "solve-")

include(joinpath(@__DIR__, "inline_mass_component_mesher.jl"))
using .InlineMassComponentMesher

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) &&
       (starts_with_solve(REQUESTED_STAGE) || REQUESTED_STAGE == "analyze")
    include(joinpath(@__DIR__, "modal_solver.jl"))
    using .ConservativeElasticModes
    using JLD2
end

if REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using FFMPEG
    using Plots
end

Base.@kwdef struct SweepCase
    id::String
    family::Symbol
    parameter_mm::Float64
    config::InlineMassConfig
end

function unique_cases(cases)
    seen = Set{String}()
    [case for case in cases if case.id ∉ seen && (push!(seen, case.id); true)]
end

function case_id(prefix, value)
    encoded = replace(string(round(Float64(value); digits=2)), "." => "p")
    "$(prefix)_$(encoded)mm"
end

function sweep_cases()
    cases = SweepCase[]
    for height in (0.35, 0.45, 0.55, 0.65, 0.75, 0.85)
        push!(cases, SweepCase(
            id=case_id("neck_h", height),
            family=:neck_height,
            parameter_mm=height,
            config=InlineMassConfig(neck_height_mm=height),
        ))
    end
    for length in (1.4, 1.7, 2.0, 2.3, 2.6, 2.9)
        push!(cases, SweepCase(
            id=case_id("mass_l", length),
            family=:mass_length,
            parameter_mm=length,
            config=InlineMassConfig(mass_length_mm=length),
        ))
    end
    for length in (0.35, 0.45, 0.55, 0.70, 0.85)
        push!(cases, SweepCase(
            id=case_id("neck_l", length),
            family=:neck_length,
            parameter_mm=length,
            config=InlineMassConfig(neck_length_mm=length),
        ))
    end
    # The one-factor pilot deliberately starts from a conservative printable
    # geometry. These combined cases follow the observed f ~ sqrt(K/M) trend
    # toward the target without changing the component topology.
    for mass_length in (1.0, 1.2, 1.4), neck_height in (0.75, 0.85, 0.95)
        mass_code = replace(string(round(mass_length; digits=2)), "." => "p")
        neck_code = replace(string(round(neck_height; digits=2)), "." => "p")
        push!(cases, SweepCase(
            id="combined_m$(mass_code)_h$(neck_code)",
            family=:combined,
            parameter_mm=mass_length,
            config=InlineMassConfig(
                mass_length_mm=mass_length,
                neck_length_mm=0.35,
                neck_height_mm=neck_height,
            ),
        ))
    end
    for mass_height in (1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.4, 2.5, 2.6, 2.7)
        push!(cases, SweepCase(
            id=case_id("target_mass_h", mass_height),
            family=:mass_height_target,
            parameter_mm=mass_height,
            config=InlineMassConfig(
                mass_length_mm=1.4,
                mass_height_mm=mass_height,
                neck_length_mm=0.35,
                neck_height_mm=0.95,
            ),
        ))
    end
    unique_cases(cases)
end

mesh_path(case) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(case.id).msh")
model_path(case) = joinpath(OUTPUT_ROOT, "models", "model_$(case.id).json")
result_path(case) = joinpath(OUTPUT_ROOT, "modes", "$(case.id).jld2")

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for case in sweep_cases()
            build_inline_mass_component_mesh(
                mesh_path(case);
                config=case.config,
                size_min_mm=0.045,
                size_max_mm=0.15,
            )
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    for case in sweep_cases()
        convert_mesh(mesh_path(case); output_dir=dirname(model_path(case)))
    end
end

function solve_case(index)
    cases = sweep_cases()
    1 <= index <= length(cases) || error("case index $index is out of range")
    case = cases[index]
    modes = solve_conservative_modes(
        model_path(case);
        config=ModeConfig(
            element_order=2,
            quadrature_degree=4,
            target_frequency_hz=TARGET_FREQUENCY_HZ,
            mode_count=MODE_COUNT,
            tolerance=1.0e-9,
            clamp_tags=["FixedInterface"],
        ),
    )
    metrics = [inline_mass_mode_metrics(mode, case.config) for mode in modes]
    selected_index = select_inline_mass_mode(modes, metrics)
    mkpath(dirname(result_path(case)))
    JLD2.jldsave(
        result_path(case);
        modes,
        metrics,
        selected_index,
        case_id=case.id,
        family=String(case.family),
        parameter_mm=case.parameter_mm,
        config=case.config,
    )
    selected = modes[selected_index]
    metric = metrics[selected_index]
    println(
        "[+] $(case.id): ", round(selected.frequency_hz / 1e3; digits=3),
        " kHz, mass=", round(metric.mass_fraction; digits=3),
        ", coherence=", round(metric.translation_coherence; digits=3),
        ", longitudinal=", round(selected.longitudinal_fraction; digits=3),
        ", parity_y=", round(selected.parity_y; digits=3),
    )
end

function run_parallel_solve_stage()
    indices = collect(eachindex(sweep_cases()))
    for batch_start in 1:MAX_SOLVERS:length(indices)
        batch = indices[batch_start:min(batch_start + MAX_SOLVERS - 1, end)]
        @sync for index in batch
            @async run(child_command("solve-$index"))
        end
    end
end

function loaded_rows()
    rows = NamedTuple[]
    for case in sweep_cases()
        data = JLD2.load(result_path(case))
        index = data["selected_index"]
        mode = data["modes"][index]
        metric = data["metrics"][index]
        push!(rows, (
            case,
            mode,
            metric,
            frequency_hz=mode.frequency_hz,
            detuning_hz=mode.frequency_hz - 242.0e3,
        ))
    end
    rows
end

function write_summary(rows)
    path = joinpath(OUTPUT_ROOT, "inline_mass_modal_sweep.csv")
    open(path, "w") do io
        println(io, "case_id,family,parameter_mm,mass_length_mm,mass_height_mm,neck_length_mm,neck_height_mm,frequency_hz,detuning_from_242_hz,mass_fraction,translation_coherence,longitudinal_fraction,parity_y,relative_residual,score")
        for row in rows
            case = row.case
            mode = row.mode
            metric = row.metric
            println(io, join((
                case.id,
                case.family,
                case.parameter_mm,
                case.config.mass_length_mm,
                case.config.mass_height_mm,
                case.config.neck_length_mm,
                case.config.neck_height_mm,
                row.frequency_hz,
                row.detuning_hz,
                metric.mass_fraction,
                metric.translation_coherence,
                mode.longitudinal_fraction,
                mode.parity_y,
                mode.relative_residual,
                metric.score,
            ), ','))
        end
    end
    println("[+] $path")
end

function best_row(rows)
    acceptable = filter(rows) do row
        row.frequency_hz > 242.0e3 &&
            row.metric.translation_coherence > 0.8 &&
            row.mode.longitudinal_fraction > 0.7
    end
    isempty(acceptable) && (acceptable = rows)
    acceptable[argmin(abs(row.frequency_hz - TARGET_FREQUENCY_HZ) for row in acceptable)]
end

function tuning_panel(rows, family, title, xlabel)
    selected = sort(filter(row -> row.case.family == family, rows); by=row -> row.case.parameter_mm)
    plot(
        getproperty.(getproperty.(selected, :case), :parameter_mm),
        getproperty.(selected, :frequency_hz) ./ 1e3;
        marker=:circle,
        linewidth=2,
        xlabel,
        ylabel="selected f_M, kHz",
        title,
        label=false,
        gridalpha=0.25,
    )
end

function mode_panel(row, phase)
    mode = row.mode
    metric = row.metric
    magnitude = sqrt.(mode.displacement_x .^ 2 .+ mode.displacement_y .^ 2)
    normalization = maximum(magnitude)
    oscillation = cos(phase)
    deformation_mm = 0.35
    x_mm = mode.node_x_m .* 1e3 .+
           deformation_mm .* oscillation .* mode.displacement_x ./ normalization
    y_mm = mode.node_y_m .* 1e3 .+
           deformation_mm .* oscillation .* mode.displacement_y ./ normalization
    color = oscillation .* mode.displacement_x ./ normalization
    total_length = component_length_mm(row.case.config)
    scatter(
        x_mm,
        y_mm;
        marker_z=color,
        color=:balance,
        clims=(-1, 1),
        colorbar=false,
        markersize=2.0,
        markerstrokewidth=0,
        label=false,
        aspect_ratio=:equal,
        xlims=(-0.5, total_length + 0.5),
        ylims=(-2.2, 2.2),
        xlabel="x, mm",
        ylabel="y, mm",
        title=string(
            row.case.id, ": ", round(mode.frequency_hz / 1e3; digits=2),
            " kHz | mass fraction=", round(metric.mass_fraction; digits=2),
            " | rigid-x=", round(metric.translation_coherence; digits=2),
        ),
        grid=false,
        framestyle=:box,
    )
end

function save_static_figure(rows, best)
    panels = [
        tuning_panel(rows, :neck_height, "Stiffness control", "neck height, mm"),
        tuning_panel(rows, :mass_length, "Mass control", "mass length, mm"),
        tuning_panel(rows, :neck_length, "Compliance control", "neck length, mm"),
        tuning_panel(rows, :mass_height_target, "Target refinement", "mass height, mm"),
    ]
    for panel in panels
        hline!(panel, [242.0]; color=:black, linestyle=:dash, label="242 kHz")
        hline!(panel, [TARGET_FREQUENCY_HZ / 1e3]; color=:darkorange, linestyle=:dot,
               label="tuning target")
    end
    push!(panels, mode_panel(best, 0.0))
    figure = plot(
        panels...;
        layout=(2, 3),
        size=(1650, 900),
        margin=5Plots.mm,
        plot_title="Fixed-interface inline mass M — quadratic FEM",
    )
    path = joinpath(OUTPUT_ROOT, "inline_mass_modal_sweep.png")
    savefig(figure, path)
    println("[+] $path")
end

function encode_animation(frame_dir)
    input_pattern = joinpath(frame_dir, "frame_%03d.png")
    mp4_path = joinpath(OUTPUT_ROOT, "inline_mass_mode.mp4")
    gif_path = joinpath(OUTPUT_ROOT, "inline_mass_mode.gif")
    FFMPEG.ffmpeg_exe(Cmd(String[
        "-y", "-framerate", string(FRAME_RATE), "-i", input_pattern,
        "-c:v", "libx264", "-pix_fmt", "yuv420p", mp4_path,
    ]))
    FFMPEG.ffmpeg_exe(Cmd(String[
        "-y", "-framerate", string(FRAME_RATE), "-i", input_pattern,
        "-vf", "fps=$(FRAME_RATE),scale=1000:-1:flags=lanczos", "-loop", "0", gif_path,
    ]))
    println("[+] $mp4_path")
    println("[+] $gif_path")
end

function save_animation(best)
    frame_dir = joinpath(OUTPUT_ROOT, "animation_frames")
    mkpath(frame_dir)
    for frame_index in 1:FRAME_COUNT
        phase = 2pi * (frame_index - 1) / FRAME_COUNT
        figure = plot(
            mode_panel(best, phase);
            size=(1100, 500),
            margin=5Plots.mm,
            plot_title="Longitudinal standing mode of fixed-interface M (deformation rescaled)",
        )
        path = joinpath(frame_dir, "frame_$(lpad(frame_index, 3, '0')).png")
        savefig(figure, path)
    end
    encode_animation(frame_dir)
end

function run_analysis_stage()
    rows = loaded_rows()
    write_summary(rows)
    best = best_row(rows)
    save_static_figure(rows, best)
    save_animation(best)
    println("\nSelected tuned M:")
    println("  case: $(best.case.id)")
    println("  frequency: $(round(best.frequency_hz / 1e3; digits=3)) kHz")
    println("  detuning from 242 kHz: $(round(best.detuning_hz / 1e3; digits=3)) kHz")
    println("  mass fraction: $(round(best.metric.mass_fraction; digits=4))")
    println("  rigid-x coherence: $(round(best.metric.translation_coherence; digits=4))")
    println("  longitudinal fraction: $(round(best.mode.longitudinal_fraction; digits=4))")
    println("  parity_y: $(round(best.mode.parity_y; digits=4))")
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_parallel_solve_stage()
    elseif starts_with_solve(stage)
        solve_case(parse(Int, stage[7:end]))
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown inline-mass stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Inline mass stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_inline_mass_modal_sweep.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
