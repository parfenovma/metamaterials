module HornPointRadiatorPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_POINT_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_point_radiator_pilot"),
)
const MAX_WORKERS = parse(Int, get(ENV, "METAMATERIALS_HORN_POINT_WORKERS", "2"))
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
using .HornPointRadiatorMesher

const CONFIG = HornPointRadiatorConfig()
const DEVICE_CONFIG = HornPointRadiatorConfig(straight_guide_length_mm=26.0)
const GENTLE_CONFIG = HornPointRadiatorConfig(rounded_axial_length_mm=26.0)
const VARIANTS = HORN_POINT_VARIANTS

case_config(variant) = variant in (:collector_gentle, :collector_smooth) ? GENTLE_CONFIG :
                       variant == :collector_device ? DEVICE_CONFIG : CONFIG

if REQUESTED_STAGE in ("mesh", "render")
    using Gmsh: gmsh
end
if REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) &&
       (REQUESTED_STAGE in ("solve", "analyze") ||
        startswith(REQUESTED_STAGE, "solve-worker-"))
    include(joinpath(@__DIR__, "transient_model_solver.jl"))
    include(joinpath(@__DIR__, "horn_point_radiator_transient_solver.jl"))
    using .TransientModelElasticity
    using .HornPointRadiatorTransientSolver
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
            build_horn_point_radiator_mesh(
                mesh_path(variant);
                config=case_config(variant),
                variant,
            )
        end
    finally
        gmsh.finalize()
    end
    println("[+] epsilon_ad=$(adiabatic_parameter(CONFIG))")
    println("[+] rounded guide amplitude=$(rounded_guide_amplitude_mm(CONFIG)) mm")
    println("[+] rounded guide arc=$(rounded_guide_length_mm(CONFIG, rounded_guide_amplitude_mm(CONFIG))) mm")
end

function triangle_segments_mm(path)
    gmsh.clear()
    gmsh.open(path)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    node_index = Dict(Int(tag) => index for (index, tag) in enumerate(node_tags))
    x = coordinates[1:3:end] .* 1e3
    y = coordinates[2:3:end] .* 1e3
    element_types, _, element_nodes = gmsh.model.mesh.getElements(2)
    x_segments, y_segments = Float64[], Float64[]
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
        :abrupt_straight => "abrupt 3.2 → 1.6 mm control",
        :collector_device => "collector + 26 mm zero-delay guide",
        :collector_straight => "log-cosine collector + straight guide",
        :collector_rounded => "collector + equal-path rounded guide",
        :collector_gentle => "collector + wavelength-scale gentle guide",
        :collector_smooth => "collector + zero-curvature sin⁴ guide",
    )
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    panels = Any[]
    try
        for variant in VARIANTS
            x, y = triangle_segments_mm(mesh_path(variant))
            push!(panels, plot(
                x, y;
                color=:steelblue, linewidth=0.18, label=false,
                aspect_ratio=:equal, xlabel="x, mm", ylabel="y, mm",
                title=labels[variant], grid=false,
            ))
        end
    finally
        gmsh.finalize()
    end
    path = joinpath(OUTPUT_ROOT, "horn_point_mesh_gallery.png")
    savefig(plot(
        panels...;
        layout=(3, 2), size=(1450, 1050), margin=4Plots.mm,
        plot_title="Collector → narrow guide → point radiator",
    ), path)
    println("[+] $path")
end

function run_convert_stage()
    for variant in VARIANTS
        convert_mesh(mesh_path(variant); output_dir=dirname(model_path(variant)))
    end
end

function solve_variant(variant)
    geometry = case_config(variant)
    probes = probe_points_mm(geometry, variant)
    run_horn_point_radiator_transient(
        model_path(variant),
        signal_path(variant);
        id=variant,
        probe_names=probes.names,
        probe_x_mm=probes.x_mm,
        probe_y_mm=probes.y_mm,
        outlet_x_mm=outlet_x_mm(geometry, variant),
        material=TransientMaterialConfig(),
        config=TransientModelConfig(
            frequency_hz=242.0e3,
            pulse_cycles=5.0,
            final_time_s=85.0e-6,
            samples_per_period=30,
            element_order=1,
            quadrature_degree=2,
        ),
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
    @sync for worker_index in 1:min(MAX_WORKERS, length(VARIANTS))
        @async run(child_command("solve-worker-$worker_index"))
    end
end

function probe_index(data, name)
    index = findfirst(==(name), data["probe_names"])
    isnothing(index) && error("probe $name is absent")
    index
end

probe_signal(data, name, component=:x) = begin
    key = component == :x ? "probe_velocity_x_m_per_s" : "probe_velocity_y_m_per_s"
    vec(data[key][probe_index(data, name), :])
end

function variant_metrics(variant, data, reference, abrupt)
    time_s = data["time_s"]
    reference_signal = probe_signal(reference, "target_axis")
    target = probe_signal(data, "target_axis")
    target_y = probe_signal(data, "target_axis", :y)
    dt_s = time_s[2] - time_s[1]
    pulse = pulse_metrics(target, reference_signal, reference_signal, dt_s)
    target_energy = sum(abs2, target) * dt_s
    reference_energy = sum(abs2, reference_signal) * dt_s
    transverse_energy = sum(abs2, target_y) * dt_s
    throat = probe_signal(data, "throat")
    abrupt_throat = probe_signal(abrupt, "throat")
    source = data["source_normal_velocity_m_per_s"]
    (
        variant=String(variant),
        target_peak_m_per_s=maximum(analytic_envelope(target)),
        peak_ratio_to_collector_straight=pulse.gain_peak,
        target_energy_ratio=target_energy / reference_energy,
        broadening_ratio=pulse.broadening_ratio,
        pulse_correlation=pulse.pulse_correlation,
        postcursor_ratio=pulse.postcursor_ratio,
        envelope_delay_us=envelope_delay(time_s, reference_signal, target) * 1e6,
        transverse_to_longitudinal_energy=transverse_energy / target_energy,
        throat_peak_m_per_s=maximum(analytic_envelope(throat)),
        throat_gain_over_abrupt=maximum(analytic_envelope(throat)) /
                                  maximum(analytic_envelope(abrupt_throat)),
        throat_gain_over_source=maximum(analytic_envelope(throat)) /
                                maximum(analytic_envelope(source)),
        primary_gate=pulse.gain_peak >= 0.8,
        pulse_gate=pulse.broadening_ratio <= 1.25 &&
                   pulse.pulse_correlation >= 0.90 &&
                   pulse.postcursor_ratio <= 0.10,
    )
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function transverse_profile(data)
    indices = findall(name -> startswith(name, "target_y_"), data["probe_names"])
    y_mm = data["probe_y_mm"][indices]
    peaks = [
        maximum(analytic_envelope(vec(data["probe_velocity_x_m_per_s"][index, :])))
        for index in indices
    ]
    order = sortperm(y_mm)
    y_mm[order], peaks[order]
end

function save_figure(rows, signals)
    colors = Dict(
        :abrupt_straight => :gray35,
        :collector_device => :darkorange,
        :collector_straight => :royalblue,
        :collector_rounded => :firebrick,
        :collector_gentle => :darkgreen,
        :collector_smooth => :purple,
    )
    labels = Dict(
        :abrupt_straight => "abrupt control",
        :collector_device => "collector + zero delay",
        :collector_straight => "collector + straight",
        :collector_rounded => "collector + rounded delay",
        :collector_gentle => "collector + gentle delay",
        :collector_smooth => "collector + sin⁴ delay",
    )
    target_panel = plot(
        xlabel="time, μs", ylabel="target vx, m/s",
        title="Radiated P-like pulse at 35 mm", gridalpha=0.25,
    )
    throat_panel = plot(
        xlabel="time, μs", ylabel="throat |analytic vx|, m/s",
        title="Local throat amplitude (diagnostic only)", gridalpha=0.25,
    )
    profile_panel = plot(
        xlabel="transverse y, mm", ylabel="peak |vx| / on-axis straight",
        title="Single-radiator transverse field at 35 mm", gridalpha=0.25,
    )
    straight_peak = maximum(analytic_envelope(probe_signal(signals[:collector_straight], "target_axis")))
    for variant in VARIANTS
        data = signals[variant]
        time_us = data["time_s"] .* 1e6
        target = probe_signal(data, "target_axis")
        throat = probe_signal(data, "throat")
        y_mm, profile = transverse_profile(data)
        plot!(target_panel, time_us, target; color=colors[variant], linewidth=2,
              label=labels[variant], xlims=(20, 80))
        plot!(throat_panel, time_us, analytic_envelope(throat); color=colors[variant],
              linewidth=2, label=labels[variant], xlims=(0, 55))
        plot!(profile_panel, y_mm, profile ./ straight_peak; color=colors[variant],
              linewidth=2, marker=:circle, markersize=3, label=labels[variant])
    end

    x = collect(eachindex(rows))
    ticks = (x, replace.(getproperty.(rows, :variant), "_" => "\n"))
    gate_panel = bar(
        x .- 0.18, getproperty.(rows, :peak_ratio_to_collector_straight);
        bar_width=0.34, label="target peak", xticks=ticks,
        ylabel="relative to collector + straight", title="End-to-end gate",
        gridalpha=0.25,
    )
    bar!(gate_panel, x .+ 0.18, getproperty.(rows, :target_energy_ratio);
         bar_width=0.34, label="target energy")
    hline!(gate_panel, [0.8]; color=:black, linestyle=:dash, label="peak gate")

    path = joinpath(OUTPUT_ROOT, "horn_point_transient.png")
    savefig(plot(
        target_panel, throat_panel, profile_panel, gate_panel;
        layout=(2, 2), size=(1500, 950), margin=5Plots.mm,
        plot_title="Five-cycle horn-fed point-radiator pilot at 242 kHz",
    ), path)
    println("[+] $path")
end

function run_analysis_stage()
    signals = Dict(variant => JLD2.load(signal_path(variant)) for variant in VARIANTS)
    reference = signals[:collector_straight]
    abrupt = signals[:abrupt_straight]
    rows = [variant_metrics(variant, signals[variant], reference, abrupt) for variant in VARIANTS]
    summary_path = joinpath(OUTPUT_ROOT, "horn_point_summary.csv")
    write_csv(summary_path, rows)
    save_figure(rows, signals)
    println("[+] $summary_path")
    for row in rows
        println(
            "[+] $(row.variant): target peak ratio=", round(row.peak_ratio_to_collector_straight; digits=3),
            ", energy=", round(row.target_energy_ratio; digits=3),
            ", throat/abrupt=", round(row.throat_gain_over_abrupt; digits=3),
            ", Bt=", round(row.broadening_ratio; digits=3),
            ", rho=", round(row.pulse_correlation; digits=3),
            ", transverse E=", round(row.transverse_to_longitudinal_energy; digits=3),
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
        worker = parse(Int, split(stage, "solve-worker-"; limit=2)[2])
        run_worker(worker)
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown horn-point stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "render", "convert", "solve", "analyze")
            println("\n=== Horn-point stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_horn_point_radiator_pilot.jl [--stage=...] ")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
