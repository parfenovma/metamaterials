module SideMassAblation

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_SIDE_MASS_ABLATION_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "side_mass_ablation"),
)
const FREQUENCIES_HZ = collect(
    parse(Float64, get(ENV, "METAMATERIALS_ABLATION_FREQUENCY_START_HZ", "235000")):
    parse(Float64, get(ENV, "METAMATERIALS_ABLATION_FREQUENCY_STEP_HZ", "250")):
    parse(Float64, get(ENV, "METAMATERIALS_ABLATION_FREQUENCY_STOP_HZ", "250000")),
)
const DAMPING_SCALE = parse(
    Float64,
    get(ENV, "METAMATERIALS_ABLATION_DAMPING_SCALE", "0"),
)
const VARIANTS = (:r0, :b, :bd)

const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "side_mass_mesher.jl"))
using .SideMassMesher

const CONFIG = SideMassConfig(bright_neck_width_mm=0.7)

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) &&
       (startswith(REQUESTED_STAGE, "solve-") || REQUESTED_STAGE == "fields")
    include(joinpath(@__DIR__, "harmonic_solver.jl"))
    using .HarmonicElasticity
    using JLD2
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
end

variant_name(variant) = uppercase(String(variant))
mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")
model_path(variant) = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
response_path(variant) = joinpath(OUTPUT_ROOT, "harmonic", "response_$(variant).jld2")
field_path(label, variant) = joinpath(OUTPUT_ROOT, "fields", "field_$(label)_$(variant).jld2")
selected_path() = joinpath(OUTPUT_ROOT, "selected_frequencies.csv")

function harmonic_config()
    HarmonicConfig(
        rayleigh_alpha=DAMPING_SCALE * 79560.0,
        rayleigh_beta=DAMPING_SCALE * 2.5e-9,
    )
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

function save_response(variant)
    points = solve_harmonic_sweep(
        model_path(variant),
        FREQUENCIES_HZ;
        config=harmonic_config(),
    )
    path = response_path(variant)
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        frequency_hz=getproperty.(points, :frequency_hz),
        left_displacement=getproperty.(points, :left_displacement),
        right_displacement=getproperty.(points, :right_displacement),
        left_traction_pa=getproperty.(points, :left_traction_pa),
        right_traction_pa=getproperty.(points, :right_traction_pa),
    )
    println("[+] $path")
end

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_solve_stage()
    jobs = min(
        parse(Int, get(ENV, "METAMATERIALS_GEOMETRY_JOBS", "3")),
        length(VARIANTS),
        Sys.CPU_THREADS,
    )
    semaphore = Base.Semaphore(jobs)
    @sync for variant in VARIANTS
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command("solve-$(variant)"))
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

function calibrated_response(data, reference)
    data["right_displacement"] ./ reference["right_displacement"]
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

function select_transparency_frequency(amplitude_b, amplitude_bd)
    contrast = amplitude_bd .- amplitude_b
    reliable = findall(index -> amplitude_bd[index] >= 0.05, eachindex(amplitude_bd))
    isempty(reliable) && error("BD response has no reliable points")
    reliable[argmax(contrast[reliable])]
end

function select_narrow_dip(amplitude_bd, peak_index)
    candidates = collect((peak_index + 1):length(amplitude_bd))
    isempty(candidates) && error("no frequencies above the transparency candidate")
    candidates[argmin(amplitude_bd[candidates])]
end

function run_analysis_stage()
    data = Dict(variant => JLD2.load(response_path(variant)) for variant in VARIANTS)
    reference = data[:r0]
    responses = Dict(
        variant => calibrated_response(data[variant], reference)
        for variant in VARIANTS
    )
    phases = Dict(variant => unwrap_phase(responses[variant]) for variant in VARIANTS)
    rows = [
        (
            variant=variant_name(variant),
            frequency_hz=frequency_hz,
            H_real=real(responses[variant][index]),
            H_imag=imag(responses[variant][index]),
            amplitude=abs(responses[variant][index]),
            phase_rad=phases[variant][index],
        )
        for variant in VARIANTS
        for (index, frequency_hz) in enumerate(FREQUENCIES_HZ)
    ]
    mkpath(OUTPUT_ROOT)
    csv_path = joinpath(OUTPUT_ROOT, "ablation_response.csv")
    write_csv(csv_path, rows)

    amplitude_b = abs.(responses[:b])
    amplitude_bd = abs.(responses[:bd])
    selected_index = select_transparency_frequency(amplitude_b, amplitude_bd)
    selected_frequency_hz = FREQUENCIES_HZ[selected_index]
    dip_index = select_narrow_dip(amplitude_bd, selected_index)
    dip_frequency_hz = FREQUENCIES_HZ[dip_index]
    open(selected_path(), "w") do io
        println(io, "label,frequency_hz")
        println(io, "transparency,$selected_frequency_hz")
        println(io, "narrow_dip,$dip_frequency_hz")
    end

    amplitude_plot = plot(
        xlabel="frequency, kHz",
        ylabel="calibrated |H|",
        title="R0 / B / BD ablation, damping scale=$(DAMPING_SCALE)",
        gridalpha=0.25,
        legend=:outertopright,
    )
    phase_plot = plot(
        xlabel="frequency, kHz",
        ylabel="unwrapped phase, rad",
        gridalpha=0.25,
        legend=:outertopright,
    )
    colors = Dict(:r0 => :gray, :b => :darkorange, :bd => :royalblue)
    for variant in VARIANTS
        plot!(
            amplitude_plot,
            FREQUENCIES_HZ ./ 1e3,
            abs.(responses[variant]);
            linewidth=2.5,
            color=colors[variant],
            label=variant_name(variant),
        )
        plot!(
            phase_plot,
            FREQUENCIES_HZ ./ 1e3,
            phases[variant];
            linewidth=2.5,
            color=colors[variant],
            label=variant_name(variant),
        )
    end
    vline!(amplitude_plot, [selected_frequency_hz / 1e3]; linestyle=:dash, color=:black, label=false)
    vline!(phase_plot, [selected_frequency_hz / 1e3]; linestyle=:dash, color=:black, label=false)
    vline!(amplitude_plot, [dip_frequency_hz / 1e3]; linestyle=:dot, color=:black, label=false)
    vline!(phase_plot, [dip_frequency_hz / 1e3]; linestyle=:dot, color=:black, label=false)
    figure = plot(amplitude_plot, phase_plot; layout=(2, 1), size=(1100, 900), margin=5Plots.mm)
    figure_path = joinpath(OUTPUT_ROOT, "ablation_spectra.png")
    savefig(figure, figure_path)
    println("[+] $csv_path")
    println("[+] $figure_path")
    println("Selected transparency candidate: $(selected_frequency_hz / 1e3) kHz")
    println("  |H_B|  = $(amplitude_b[selected_index])")
    println("  |H_BD| = $(amplitude_bd[selected_index])")
    println("Selected narrow dip: $(dip_frequency_hz / 1e3) kHz")
    println("  |H_BD| = $(amplitude_bd[dip_index])")
end

function save_field(label, variant, frequency_hz)
    field = solve_harmonic_field(
        model_path(variant),
        frequency_hz;
        config=harmonic_config(),
    )
    path = field_path(label, variant)
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        label=String(label),
        variant=String(variant),
        frequency_hz=field.point.frequency_hz,
        node_x_m=field.node_x_m,
        node_y_m=field.node_y_m,
        cell_node_ids=field.cell_node_ids,
        displacement_x_m=field.displacement_x_m,
        displacement_y_m=field.displacement_y_m,
        right_displacement=field.point.right_displacement,
    )
    println("[+] $path")
end

function selected_frequencies()
    lines = readlines(selected_path())
    [
        let columns = split(line, ',')
            (label=Symbol(columns[1]), frequency_hz=parse(Float64, columns[2]))
        end
        for line in Iterators.drop(lines, 1)
        if !isempty(strip(line))
    ]
end

function run_fields_stage()
    isfile(selected_path()) || error("run analyze before fields")
    for selected in selected_frequencies()
        for variant in VARIANTS
            save_field(selected.label, variant, selected.frequency_hz)
        end
    end
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_solve_stage()
    elseif startswith(stage, "solve-")
        variant = Symbol(stage[7:end])
        variant in VARIANTS || error("unknown ablation variant: $variant")
        save_response(variant)
    elseif stage == "analyze"
        run_analysis_stage()
    elseif stage == "fields"
        run_fields_stage()
    else
        error("unknown ablation stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze", "fields")
            println("\n=== Ablation stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_side_mass_ablation.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
