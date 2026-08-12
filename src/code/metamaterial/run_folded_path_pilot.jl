module FoldedPathPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_FOLDED_PATH_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "folded_path_pilot"),
)
const MAX_WORKERS = parse(Int, get(ENV, "METAMATERIALS_FOLDED_PATH_WORKERS", "2"))
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "folded_path_mesher.jl"))
using .FoldedPathMesher

const VARIANTS = FOLDED_PATH_VARIANTS
const CONFIG = FoldedPathConfig()

if REQUESTED_STAGE in ("mesh", "render")
    using Gmsh: gmsh
end
if REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) &&
       (REQUESTED_STAGE in ("solve", "analyze") || startswith(REQUESTED_STAGE, "solve-worker-"))
    include(joinpath(@__DIR__, "transient_model_solver.jl"))
    using .TransientModelElasticity
end
if REQUESTED_STAGE in ("render", "analyze")
    ENV["GKSwstype"] = "100"
    using Plots
end
if REQUESTED_STAGE == "analyze"
    using JLD2
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    include(joinpath(@__DIR__, "spectral_analysis.jl"))
    using .ImpulseRiskAnalysis
    using .SpectralAnalysis
end

mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")
model_path(variant) = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
signal_path(variant) = joinpath(OUTPUT_ROOT, "signals", "$(variant).jld2")

function child_command(stage)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            build_folded_path_mesh(mesh_path(variant); config=CONFIG, variant)
        end
    finally
        gmsh.finalize()
    end
end

function triangle_segments_mm(path)
    gmsh.clear()
    gmsh.open(path)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    node_index = Dict(Int(tag) => index for (index, tag) in enumerate(node_tags))
    x = coordinates[1:3:end] .* 1e3
    y = coordinates[2:3:end] .* 1e3
    element_types, _, element_nodes = gmsh.model.mesh.getElements(2)
    x_segments = Float64[]
    y_segments = Float64[]
    for (element_type, connectivity) in zip(element_types, element_nodes)
        name, _, _, node_count, _, _ = gmsh.model.mesh.getElementProperties(element_type)
        startswith(name, "Triangle") || continue
        for offset in 1:node_count:length(connectivity)
            indices = getindex.(Ref(node_index), Int.(connectivity[offset:(offset + 2)]))
            append!(x_segments, (x[indices[1]], x[indices[2]], x[indices[3]], x[indices[1]], NaN))
            append!(y_segments, (y[indices[1]], y[indices[2]], y[indices[3]], y[indices[1]], NaN))
        end
    end
    x_segments, y_segments
end

function run_render_stage()
    labels = Dict(
        :straight_device => "straight device-length reference",
        :straight_unfolded => "straight equal-path reference",
        :sharp_fold => "sharp mirror-like fold proxy",
        :rounded_fold => "rounded guided fold",
    )
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    panels = Any[]
    try
        for variant in VARIANTS
            x, y = triangle_segments_mm(mesh_path(variant))
            push!(panels, plot(
                x,
                y;
                color=:steelblue,
                linewidth=0.22,
                label=false,
                aspect_ratio=:equal,
                xlabel="x, mm",
                ylabel="y, mm",
                title=labels[variant],
                grid=false,
            ))
        end
    finally
        gmsh.finalize()
    end
    figure = plot(
        panels...;
        layout=(2, 2),
        size=(1500, 900),
        margin=4Plots.mm,
        plot_title="Global-fold transient pilot — equal unfolded path 30.42 mm",
    )
    path = joinpath(OUTPUT_ROOT, "folded_path_mesh_gallery.png")
    mkpath(dirname(path))
    savefig(figure, path)
    println("[+] $path")
end

function run_convert_stage()
    for variant in VARIANTS
        convert_mesh(mesh_path(variant); output_dir=dirname(model_path(variant)))
    end
end

function solve_variant(variant)
    transient_config = TransientModelConfig(
        frequency_hz=242.0e3,
        pulse_cycles=5.0,
        final_time_s=80.0e-6,
        samples_per_period=30,
        element_order=1,
        quadrature_degree=2,
    )
    run_transient_model(
        model_path(variant),
        signal_path(variant);
        id=variant,
        material=TransientMaterialConfig(),
        config=transient_config,
    )
end

function run_worker(worker_index)
    for index in worker_index:MAX_WORKERS:length(VARIANTS)
        variant = VARIANTS[index]
        isfile(signal_path(variant)) && continue
        solve_variant(variant)
    end
end

function run_solve_stage()
    any(!isfile(signal_path(variant)) for variant in VARIANTS) || return
    @sync for worker_index in 1:MAX_WORKERS
        @async run(child_command("solve-worker-$worker_index"))
    end
end

load_signal(variant) = JLD2.load(signal_path(variant))

function variant_metrics(variant, signal, reference)
    time_s = signal["time_s"]
    reference_time_s = reference["time_s"]
    length(time_s) == length(reference_time_s) || error("transient sample counts differ")
    maximum(abs.(time_s .- reference_time_s)) < 1e-15 || error("transient time grids differ")
    normal = signal["right_normal_velocity_m_per_s"]
    tangent = signal["right_tangent_velocity_m_per_s"]
    reference_normal = reference["right_normal_velocity_m_per_s"]
    dt_s = time_s[2] - time_s[1]
    pulse = pulse_metrics(normal, reference_normal, reference_normal, dt_s)
    normal_energy = sum(abs2, normal) * dt_s
    reference_energy = sum(abs2, reference_normal) * dt_s
    tangent_energy = sum(abs2, tangent) * dt_s
    envelope = analytic_envelope(normal)
    reference_envelope = analytic_envelope(reference_normal)
    (
        variant,
        peak_ratio=pulse.gain_peak,
        normal_energy_ratio=normal_energy / reference_energy,
        broadening_ratio=pulse.broadening_ratio,
        pulse_correlation=pulse.pulse_correlation,
        postcursor_ratio=pulse.postcursor_ratio,
        envelope_delay_us=envelope_delay(time_s, reference_normal, normal) * 1e6,
        envelope_peak_shift_us=(time_s[argmax(envelope)] - time_s[argmax(reference_envelope)]) * 1e6,
        tangent_to_normal_peak=maximum(abs, tangent) / maximum(abs, normal),
        tangent_to_normal_energy=tangent_energy / normal_energy,
        primary_gate=pulse.gain_peak >= 0.8,
        pulse_gate=pulse.broadening_ratio <= 1.25 &&
                   pulse.pulse_correlation >= 0.90 &&
                   pulse.postcursor_ratio <= 0.10,
    )
end

function write_summary(rows)
    path = joinpath(OUTPUT_ROOT, "folded_path_summary.csv")
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
    println("[+] $path")
end

function save_result_figure(rows, signals)
    colors = Dict(
        :straight_device => :black,
        :straight_unfolded => :gray40,
        :sharp_fold => :firebrick,
        :rounded_fold => :royalblue,
    )
    labels = Dict(
        :straight_device => "straight device",
        :straight_unfolded => "straight equal path",
        :sharp_fold => "sharp fold",
        :rounded_fold => "rounded fold",
    )
    waveform_panel = plot(
        xlabel="time, μs", ylabel="mean output vx, m/s",
        title="Useful output component", gridalpha=0.25,
    )
    envelope_panel = plot(
        xlabel="time, μs", ylabel="envelope / equal-path peak",
        title="Analytic envelopes", gridalpha=0.25,
    )
    reference_peak = maximum(analytic_envelope(
        signals[:straight_unfolded]["right_normal_velocity_m_per_s"],
    ))
    for variant in VARIANTS
        data = signals[variant]
        time_us = data["time_s"] .* 1e6
        waveform = data["right_normal_velocity_m_per_s"]
        plot!(waveform_panel, time_us, waveform;
              color=colors[variant], linewidth=2, label=labels[variant], xlims=(15, 75))
        plot!(envelope_panel, time_us, analytic_envelope(waveform) ./ reference_peak;
              color=colors[variant], linewidth=2, label=labels[variant], xlims=(15, 75))
    end

    x = collect(1:length(VARIANTS))
    tick_labels = [replace(labels[variant], " " => "\n") for variant in VARIANTS]
    amplitude_panel = plot(
        x,
        getproperty.(rows, :peak_ratio);
        marker=:circle,
        linewidth=2,
        label="peak ratio",
        xticks=(x, tick_labels),
        ylabel="relative to equal path",
        title="Amplitude and output energy",
        gridalpha=0.25,
    )
    plot!(amplitude_panel, x, getproperty.(rows, :normal_energy_ratio);
          marker=:square, linewidth=2, label="normal energy ratio")
    hline!(amplitude_panel, [0.8]; color=:gray, linestyle=:dash, label="primary gate")

    fidelity_panel = plot(
        x,
        getproperty.(rows, :broadening_ratio);
        marker=:circle,
        linewidth=2,
        label="B_t",
        xticks=(x, tick_labels),
        title="Pulse fidelity and converted output",
        gridalpha=0.25,
    )
    plot!(fidelity_panel, x, getproperty.(rows, :pulse_correlation);
          marker=:square, linewidth=2, label="correlation")
    plot!(fidelity_panel, x, getproperty.(rows, :postcursor_ratio);
          marker=:diamond, linewidth=2, label="postcursor")
    plot!(fidelity_panel, x, getproperty.(rows, :tangent_to_normal_energy);
          marker=:utriangle, linewidth=2, label="tangent/normal energy")

    figure = plot(
        waveform_panel,
        envelope_panel,
        amplitude_panel,
        fidelity_panel;
        layout=(2, 2),
        size=(1500, 920),
        margin=5Plots.mm,
        plot_title="Five-cycle elastic global-fold screening at 242 kHz",
    )
    path = joinpath(OUTPUT_ROOT, "folded_path_transient.png")
    savefig(figure, path)
    println("[+] $path")
end

function run_analysis_stage()
    signals = Dict(variant => load_signal(variant) for variant in VARIANTS)
    reference = signals[:straight_unfolded]
    rows = [variant_metrics(variant, signals[variant], reference) for variant in VARIANTS]
    write_summary(rows)
    save_result_figure(rows, signals)
    for row in rows
        println(
            "[+] $(row.variant): peak=", round(row.peak_ratio; digits=3),
            ", energy=", round(row.normal_energy_ratio; digits=3),
            ", Bt=", round(row.broadening_ratio; digits=3),
            ", rho=", round(row.pulse_correlation; digits=3),
            ", post=", round(row.postcursor_ratio; digits=3),
            ", tangent E=", round(row.tangent_to_normal_energy; digits=3),
        )
    end
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "render"
        run_render_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_solve_stage()
    elseif startswith(stage, "solve-worker-")
        worker_index = parse(Int, split(stage, "solve-worker-"; limit=2)[2])
        1 <= worker_index <= MAX_WORKERS || error("invalid worker index")
        run_worker(worker_index)
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown folded-path stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "render", "convert", "solve", "analyze")
            println("\n=== Folded-path stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_folded_path_pilot.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
