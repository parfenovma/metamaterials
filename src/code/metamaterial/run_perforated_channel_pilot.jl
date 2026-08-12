module PerforatedChannelPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_PERFORATED_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "perforated_channel_pilot"),
)
const MAX_WORKERS = parse(Int, get(ENV, "METAMATERIALS_PERFORATED_WORKERS", "3"))
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "perforated_channel_mesher.jl"))
using .PerforatedChannelMesher

const VARIANTS = PERFORATED_VARIANTS
const CONFIG = PerforatedChannelConfig()
const FREQUENCIES_HZ = collect(162.0e3:20.0e3:322.0e3)

is_solver_stage(stage) = stage == "lossless" || startswith(stage, "lossless-worker-")

if REQUESTED_STAGE in ("mesh", "render")
    using Gmsh: gmsh
end
if REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) && (is_solver_stage(REQUESTED_STAGE) || REQUESTED_STAGE == "analyze")
    include(joinpath(@__DIR__, "modal_harmonic_solver.jl"))
    include(joinpath(@__DIR__, "modal_projection.jl"))
    using .ModalHarmonicElasticity
    using .ElasticModalProjection
    using JLD2
end
if REQUESTED_STAGE in ("render", "analyze")
    ENV["GKSwstype"] = "100"
    using Plots
end

mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")
model_path(variant) = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
point_path(variant, frequency_hz) = joinpath(
    OUTPUT_ROOT,
    "lossless",
    "$(variant)_$(round(Int, frequency_hz))hz.jld2",
)

specs() = [
    (; variant, frequency_hz)
    for variant in VARIANTS
    for frequency_hz in FREQUENCIES_HZ
]

function child_command(stage)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            build_perforated_channel_mesh(mesh_path(variant); config=CONFIG, variant)
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
            corner_tags = connectivity[offset:(offset + 2)]
            indices = getindex.(Ref(node_index), Int.(corner_tags))
            append!(x_segments, (x[indices[1]], x[indices[2]], x[indices[3]], x[indices[1]], NaN))
            append!(y_segments, (y[indices[1]], y[indices[2]], y[indices[3]], y[indices[1]], NaN))
        end
    end
    x_segments, y_segments
end

function run_render_stage()
    labels = Dict(
        :solid => "solid reference",
        :round => "circular holes",
        :ellipse_longitudinal => "longitudinal ellipses",
        :ellipse_transverse => "transverse ellipses",
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
                linewidth=0.25,
                label=false,
                aspect_ratio=:equal,
                xlims=(CONFIG.lead_length_mm - 0.5, CONFIG.lead_length_mm + CONFIG.active_length_mm + 0.5),
                ylims=(-CONFIG.port_height_mm / 2 - 0.15, CONFIG.port_height_mm / 2 + 0.15),
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
        layout=(4, 1),
        size=(1350, 950),
        margin=3Plots.mm,
        plot_title="Equal-area perforated-channel pilot: active region",
    )
    path = joinpath(OUTPUT_ROOT, "perforated_mesh_gallery.png")
    mkpath(dirname(path))
    savefig(figure, path)
    println("[+] $path")
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

function solve_scattering(variant, frequency_hz)
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
        port_height_m=CONFIG.port_height_mm * 1e-3,
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
    wavenumber = real(state.right_mode.wavenumber_per_m)
    deembedded = transmission_amplitude * cis(wavenumber * total_length_mm(CONFIG) * 1e-3)
    (
        variant,
        frequency_hz=Float64(frequency_hz),
        T00=transmission,
        R00=reflection,
        C0=conversion,
        right_incoming,
        transmission_amplitude=ComplexF64(transmission_amplitude),
        deembedded_amplitude=abs(deembedded),
        deembedded_phase_rad=angle(deembedded),
        propagating_mode_count=length(left),
        power_sum=transmission + reflection + conversion,
        variational_balance_residual=state.metrics.balance_residual_w_per_m / incident,
    )
end

function solve_and_save(spec)
    row = solve_scattering(spec.variant, spec.frequency_hz)
    path = point_path(spec.variant, spec.frequency_hz)
    mkpath(dirname(path))
    JLD2.jldsave(path; row)
    println(
        "[+] $(spec.variant) ", round(spec.frequency_hz / 1e3; digits=1),
        " kHz: T=", round(row.T00; digits=4),
        ", R=", round(row.R00; digits=4),
        ", sum=", round(row.power_sum; digits=4),
    )
end

function run_worker(worker_index)
    for index in worker_index:MAX_WORKERS:length(specs())
        spec = specs()[index]
        isfile(point_path(spec.variant, spec.frequency_hz)) && continue
        solve_and_save(spec)
    end
end

function run_parallel_stage()
    any(!isfile(point_path(spec.variant, spec.frequency_hz)) for spec in specs()) || return
    @sync for worker_index in 1:MAX_WORKERS
        @async run(child_command("lossless-worker-$worker_index"))
    end
end

function load_rows()
    [JLD2.load(point_path(spec.variant, spec.frequency_hz))["row"] for spec in specs()]
end

function unwrap_phase(values)
    result = Float64[first(values)]
    for raw in Iterators.drop(values, 1)
        push!(result, last(result) + mod(raw - last(result) + pi, 2pi) - pi)
    end
    result
end

function phase_delay_rows(rows, variant)
    selected = sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
    references = Dict(row.frequency_hz => row for row in rows if row.variant == :solid)
    relative_transfer = [
        row.transmission_amplitude / references[row.frequency_hz].transmission_amplitude
        for row in selected
    ]
    phase = unwrap_phase(angle.(relative_transfer))
    omega = 2pi .* getproperty.(selected, :frequency_hz)
    delay = similar(phase)
    delay[1] = -(phase[2] - phase[1]) / (omega[2] - omega[1])
    for index in 2:(length(phase) - 1)
        delay[index] = -(phase[index + 1] - phase[index - 1]) /
                       (omega[index + 1] - omega[index - 1])
    end
    delay[end] = -(phase[end] - phase[end - 1]) / (omega[end] - omega[end - 1])
    trusted_delay = [
        begin
            stencil = index == 1 ? (1:2) :
                      index == length(selected) ? ((length(selected) - 1):length(selected)) :
                      ((index - 1):(index + 1))
            all(selected[neighbor].T00 >= 0.5 for neighbor in stencil) ? delay[index] : NaN
        end
        for index in eachindex(selected)
    ]
    [
        merge(row, (
            relative_amplitude=abs(relative_transfer[index]),
            relative_power=abs2(relative_transfer[index]),
            relative_phase_rad=phase[index],
            excess_group_delay_s=delay[index],
            trusted_excess_group_delay_s=trusted_delay[index],
        ))
        for (index, row) in enumerate(selected)
    ]
end

function linear_interpolate(x, values, query)
    query <= first(x) && return first(values)
    query >= last(x) && return last(values)
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function write_csv(rows)
    path = joinpath(OUTPUT_ROOT, "perforated_lossless.csv")
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
    println("[+] $path")
end

function save_figure(rows)
    transmission_panel = plot(
        xlabel="frequency, kHz", ylabel="T00", title="Absolute transmission",
        ylims=(0, 1.05), gridalpha=0.25,
    )
    reflection_panel = plot(
        xlabel="frequency, kHz", ylabel="R00", title="Absolute reflection",
        ylims=(0, 1.05), gridalpha=0.25,
    )
    phase_panel = plot(
        xlabel="frequency, kHz", ylabel="phase, rad",
        title="Transmission phase relative to solid", gridalpha=0.25,
    )
    delay_panel = plot(
        xlabel="frequency, kHz", ylabel="excess delay, μs",
        title="Coarse group-delay diagnostic", gridalpha=0.25,
    )
    styles = Dict(
        :solid => (:black, :circle),
        :round => (:royalblue, :square),
        :ellipse_longitudinal => (:darkorange, :diamond),
        :ellipse_transverse => (:seagreen, :utriangle),
    )
    labels = Dict(
        :solid => "solid",
        :round => "round",
        :ellipse_longitudinal => "ellipse ∥ x",
        :ellipse_transverse => "ellipse along y",
    )
    for variant in VARIANTS
        current = sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
        color, marker = styles[variant]
        frequencies = getproperty.(current, :frequency_hz) ./ 1e3
        plot!(transmission_panel, frequencies, getproperty.(current, :T00);
              color, marker, linewidth=2, label=labels[variant])
        plot!(reflection_panel, frequencies, getproperty.(current, :R00);
              color, marker, linewidth=2, label=labels[variant])
        plot!(phase_panel, frequencies, getproperty.(current, :relative_phase_rad);
              color, marker, linewidth=2, label=labels[variant])
        plot!(delay_panel, frequencies, getproperty.(current, :trusted_excess_group_delay_s) .* 1e6;
              color, marker, linewidth=2, label=labels[variant])
    end
    for panel in (transmission_panel, reflection_panel, phase_panel, delay_panel)
        vline!(panel, [242.0]; color=:gray, linestyle=:dash, label=false)
    end
    figure = plot(
        transmission_panel,
        reflection_panel,
        phase_panel,
        delay_panel;
        layout=(2, 2),
        size=(1450, 900),
        margin=5Plots.mm,
        plot_title="Equal-area through-hole pilot — straight 3.2 mm channel",
    )
    path = joinpath(OUTPUT_ROOT, "perforated_lossless.png")
    savefig(figure, path)
    println("[+] $path")
end

function write_summary(rows)
    path = joinpath(OUTPUT_ROOT, "perforated_summary.csv")
    open(path, "w") do io
        println(io, "variant,void_area_mm2,min_ligament_mm,min_T00,max_R00,raw_phase_span_rad,trusted_delay_at_242_us,max_trusted_abs_delay_us,max_power_error")
        for variant in VARIANTS
            current = sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
            ligaments = minimum_ligaments_mm(CONFIG, variant)
            minimum_ligament = variant == :solid ? Inf : minimum(values(ligaments))
            frequencies = getproperty.(current, :frequency_hz)
            delay_us = getproperty.(current, :trusted_excess_group_delay_s) .* 1e6
            finite_delay_us = filter(isfinite, delay_us)
            metrics = (
                void_area_mm2(CONFIG, variant),
                minimum_ligament,
                minimum(getproperty.(current, :T00)),
                maximum(getproperty.(current, :R00)),
                maximum(getproperty.(current, :relative_phase_rad)) - minimum(getproperty.(current, :relative_phase_rad)),
                linear_interpolate(frequencies, delay_us, 242.0e3),
                isempty(finite_delay_us) ? NaN : maximum(abs.(finite_delay_us)),
                maximum(abs.(getproperty.(current, :power_sum) .- 1)),
            )
            println(io, join((variant, metrics...), ','))
            println(
                "[+] $variant: min T=", round(metrics[3]; digits=4),
                ", max R=", round(metrics[4]; digits=4),
                ", delay@242=", round(metrics[6]; digits=3), " μs",
                ", phase span=", round(metrics[5]; digits=3), " rad",
            )
        end
    end
    println("[+] $path")
end

function run_analysis_stage()
    raw = load_rows()
    rows = reduce(vcat, (phase_delay_rows(raw, variant) for variant in VARIANTS))
    write_csv(rows)
    write_summary(rows)
    save_figure(rows)
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "render"
        run_render_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "lossless"
        run_parallel_stage()
    elseif startswith(stage, "lossless-worker-")
        worker_index = parse(Int, split(stage, "lossless-worker-"; limit=2)[2])
        1 <= worker_index <= MAX_WORKERS || error("invalid worker index")
        run_worker(worker_index)
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown perforated-channel stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "render", "convert", "lossless", "analyze")
            println("\n=== Perforated-channel stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_perforated_channel_pilot.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
