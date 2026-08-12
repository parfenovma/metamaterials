module FanoDarkMassTuning

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_FANO_TUNING_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "fano_dark_mass_tuning"),
)
const DARK_LENGTHS_MM = parse.(
    Float64,
    split(get(ENV, "METAMATERIALS_FANO_DARK_LENGTHS_MM", "2.94,3.00,3.06"), ','),
)
const DARK_BRIDGE_HEIGHT_MM = parse(
    Float64,
    get(ENV, "METAMATERIALS_FANO_DARK_BRIDGE_HEIGHT_MM", "0.5"),
)
const BRIGHT_NECK_WIDTH_MM = parse(
    Float64,
    get(ENV, "METAMATERIALS_FANO_BRIGHT_NECK_WIDTH_MM", "0.7"),
)
const DAMPING_SCALE = parse(
    Float64,
    get(ENV, "METAMATERIALS_FANO_DAMPING_SCALE", "0"),
)
const MESH_SIZE_MIN_MM = parse(
    Float64,
    get(ENV, "METAMATERIALS_FANO_MESH_SIZE_MIN_MM", "0.10"),
)
const MESH_SIZE_MAX_MM = parse(
    Float64,
    get(ENV, "METAMATERIALS_FANO_MESH_SIZE_MAX_MM", "0.45"),
)
const FREQUENCIES_HZ = collect(
    parse(Float64, get(ENV, "METAMATERIALS_FANO_START_HZ", "238000")):
    parse(Float64, get(ENV, "METAMATERIALS_FANO_STEP_HZ", "250")):
    parse(Float64, get(ENV, "METAMATERIALS_FANO_STOP_HZ", "247000")),
)
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "side_mass_mesher.jl"))
using .SideMassMesher
include(joinpath(@__DIR__, "side_mass_component_mesher.jl"))
using .SideMassComponentMesher

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) && startswith(REQUESTED_STAGE, "eigen-")
    include(joinpath(@__DIR__, "modal_solver.jl"))
    using .ConservativeElasticModes
    using JLD2
elseif !isnothing(REQUESTED_STAGE) && startswith(REQUESTED_STAGE, "harmonic-")
    include(joinpath(@__DIR__, "harmonic_solver.jl"))
    using .HarmonicElasticity
    using JLD2
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
end

label(value) = replace(@sprintf("%.2f", Float64(value)), "." => "p")
config(length_mm) = SideMassConfig(
    bright_neck_width_mm=BRIGHT_NECK_WIDTH_MM,
    dark_mass_length_mm=Float64(length_mm),
    dark_bridge_height_mm=DARK_BRIDGE_HEIGHT_MM,
)
case_id(length_mm) = "ld_$(label(length_mm))"
mesh_path(kind, id) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(kind)_$(id).msh")
model_path(kind, id) = joinpath(OUTPUT_ROOT, "models", "model_$(kind)_$(id).json")
eigen_path(id) = joinpath(OUTPUT_ROOT, "eigen", "eigen_$(id).jld2")
harmonic_path(id) = joinpath(OUTPUT_ROOT, "harmonic", "response_$(id).jld2")

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        baseline = config(3.0)
        build_side_mass_mesh(
            mesh_path("r0", "reference");
            config=baseline,
            variant=:r0,
            size_min_mm=MESH_SIZE_MIN_MM,
            size_max_mm=MESH_SIZE_MAX_MM,
        )
        build_bright_component_mesh(mesh_path("component", "bright"); config=baseline)
        for length_mm in DARK_LENGTHS_MM
            id = case_id(length_mm)
            current = config(length_mm)
            build_side_mass_mesh(
                mesh_path("bd", id);
                config=current,
                variant=:bd,
                size_min_mm=MESH_SIZE_MIN_MM,
                size_max_mm=MESH_SIZE_MAX_MM,
            )
            build_dark_component_mesh(mesh_path("component", id); config=current)
        end
    finally
        gmsh.finalize()
    end
end

function convert_one(kind, id)
    convert_mesh(mesh_path(kind, id); output_dir=joinpath(OUTPUT_ROOT, "models"))
end

function run_convert_stage()
    convert_one("r0", "reference")
    convert_one("component", "bright")
    for length_mm in DARK_LENGTHS_MM
        id = case_id(length_mm)
        convert_one("bd", id)
        convert_one("component", id)
    end
end

function closest_mode(modes, target_hz)
    argmin(abs(mode.frequency_hz - target_hz) for mode in modes)
end

function localized_even_mode(modes, metrics)
    candidates = findall(index ->
        modes[index].parity_y >= 0.8 &&
        230.0e3 <= modes[index].frequency_hz <= 250.0e3,
        eachindex(modes),
    )
    isempty(candidates) && error("no localized even full-cell mode in the target window")
    candidates[argmax(metrics[index].side_mass_fraction for index in candidates)]
end

function solve_component(model, target_hz)
    modes = solve_conservative_modes(
        model;
        config=ModeConfig(
            target_frequency_hz=target_hz,
            mode_count=8,
            tolerance=1.0e-10,
            clamp_tags=["FixedInterface"],
        ),
    )
    modes[closest_mode(modes, target_hz)]
end

function run_eigen_stage(index)
    1 <= index <= length(DARK_LENGTHS_MM) || error("invalid case index")
    length_mm = DARK_LENGTHS_MM[index]
    id = case_id(length_mm)
    current = config(length_mm)
    bright_mode = solve_component(model_path("component", "bright"), 246.0e3)
    dark_mode = solve_component(model_path("component", id), 246.0e3)
    full_modes = solve_conservative_modes(
        model_path("bd", id);
        config=ModeConfig(
            target_frequency_hz=242.0e3,
            mode_count=12,
            tolerance=1.0e-9,
            clamp_tags=String[],
        ),
    )
    full_metrics = [side_mass_mode_metrics(mode, current; variant=:bd) for mode in full_modes]
    selected_index = localized_even_mode(full_modes, full_metrics)
    full_mode = full_modes[selected_index]
    full_metric = full_metrics[selected_index]
    path = eigen_path(id)
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        dark_length_mm=length_mm,
        bright_frequency_hz=bright_mode.frequency_hz,
        dark_frequency_hz=dark_mode.frequency_hz,
        full_frequency_hz=full_mode.frequency_hz,
        full_parity_y=full_mode.parity_y,
        full_bright_fraction=full_metric.bright_fraction,
        full_dark_fraction=full_metric.dark_fraction,
        full_residual=full_mode.relative_residual,
    )
    println("[+] $path")
end

function run_eigen_parallel_stage()
    @sync for index in eachindex(DARK_LENGTHS_MM)
        @async run(child_command("eigen-$index"))
    end
end

function harmonic_config()
    HarmonicConfig(
        rayleigh_alpha=DAMPING_SCALE * 79560.0,
        rayleigh_beta=DAMPING_SCALE * 2.5e-9,
    )
end

function save_harmonic(path, model)
    points = solve_harmonic_sweep(model, FREQUENCIES_HZ; config=harmonic_config())
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        frequency_hz=getproperty.(points, :frequency_hz),
        right_displacement=getproperty.(points, :right_displacement),
        right_traction_pa=getproperty.(points, :right_traction_pa),
        source_power_w_per_m=getproperty.(points, :source_power_w_per_m),
        left_absorbed_power_w_per_m=getproperty.(points, :left_absorbed_power_w_per_m),
        right_absorbed_power_w_per_m=getproperty.(points, :right_absorbed_power_w_per_m),
        internal_dissipated_power_w_per_m=getproperty.(points, :internal_dissipated_power_w_per_m),
    )
    println("[+] $path")
end

function run_harmonic_stage(id)
    if id == "reference"
        save_harmonic(harmonic_path(id), model_path("r0", "reference"))
        return
    end
    index = parse(Int, id)
    1 <= index <= length(DARK_LENGTHS_MM) || error("invalid harmonic case index")
    case = case_id(DARK_LENGTHS_MM[index])
    save_harmonic(harmonic_path(case), model_path("bd", case))
end

function run_harmonic_parallel_stage()
    stages = ["harmonic-reference"; ["harmonic-$index" for index in eachindex(DARK_LENGTHS_MM)]]
    @sync for stage in stages
        @async run(child_command(stage))
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

function characteristic_indices(amplitude, center_hz)
    candidates = findall(abs.(FREQUENCIES_HZ .- center_hz) .<= 3.0e3)
    isempty(candidates) && error("no harmonic points around eigenfrequency")
    peak = candidates[argmax(amplitude[candidates])]
    dip = candidates[argmin(amplitude[candidates])]
    peak, dip
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

function run_analysis_stage()
    reference = JLD2.load(harmonic_path("reference"))
    rows = NamedTuple[]
    summary = NamedTuple[]
    amplitude_plot = plot(
        xlabel="frequency, kHz",
        ylabel="reference-calibrated |H|",
        title="Fano tuning: tb=$(BRIGHT_NECK_WIDTH_MM), td=$(DARK_BRIDGE_HEIGHT_MM) mm, s=$(DAMPING_SCALE)",
        gridalpha=0.25,
        legend=:outertopright,
    )
    phase_plot = plot(
        xlabel="frequency, kHz",
        ylabel="unwrapped phase, rad",
        gridalpha=0.25,
        legend=:outertopright,
    )
    power_plot = plot(
        xlabel="frequency, kHz",
        ylabel="Pout / Pout,R0",
        gridalpha=0.25,
        legend=:outertopright,
    )
    colors = [:darkorange, :royalblue, :seagreen]
    for (case_index, length_mm) in enumerate(DARK_LENGTHS_MM)
        id = case_id(length_mm)
        harmonic = JLD2.load(harmonic_path(id))
        eigen = JLD2.load(eigen_path(id))
        response = harmonic["right_displacement"] ./ reference["right_displacement"]
        amplitude = abs.(response)
        phase = unwrap_phase(response)
        transmitted_power = harmonic["right_absorbed_power_w_per_m"] ./
                            reference["right_absorbed_power_w_per_m"]
        power_balance_error = harmonic["source_power_w_per_m"] .-
                              harmonic["left_absorbed_power_w_per_m"] .-
                              harmonic["right_absorbed_power_w_per_m"] .-
                              harmonic["internal_dissipated_power_w_per_m"]
        peak_index, dip_index = characteristic_indices(amplitude, eigen["full_frequency_hz"])
        push!(summary, (
            dark_length_mm=length_mm,
            dark_bridge_height_mm=DARK_BRIDGE_HEIGHT_MM,
            bright_neck_width_mm=BRIGHT_NECK_WIDTH_MM,
            damping_scale=DAMPING_SCALE,
            mesh_size_min_mm=MESH_SIZE_MIN_MM,
            mesh_size_max_mm=MESH_SIZE_MAX_MM,
            bright_component_frequency_hz=eigen["bright_frequency_hz"],
            dark_component_frequency_hz=eigen["dark_frequency_hz"],
            full_even_mode_frequency_hz=eigen["full_frequency_hz"],
            peak_frequency_hz=FREQUENCIES_HZ[peak_index],
            peak_amplitude=amplitude[peak_index],
            peak_power_ratio=transmitted_power[peak_index],
            dip_frequency_hz=FREQUENCIES_HZ[dip_index],
            dip_amplitude=amplitude[dip_index],
            dip_power_ratio=transmitted_power[dip_index],
            peak_to_dip_contrast=amplitude[peak_index] / amplitude[dip_index],
            max_relative_power_balance_error=maximum(abs.(power_balance_error) ./
                max.(abs.(harmonic["source_power_w_per_m"]), eps())),
        ))
        for frequency_index in eachindex(FREQUENCIES_HZ)
            push!(rows, (
                dark_length_mm=length_mm,
                frequency_hz=FREQUENCIES_HZ[frequency_index],
                H_real=real(response[frequency_index]),
                H_imag=imag(response[frequency_index]),
                amplitude=amplitude[frequency_index],
                phase_rad=phase[frequency_index],
                transmitted_power_ratio=transmitted_power[frequency_index],
                source_power_w_per_m=harmonic["source_power_w_per_m"][frequency_index],
                left_absorbed_power_w_per_m=harmonic["left_absorbed_power_w_per_m"][frequency_index],
                right_absorbed_power_w_per_m=harmonic["right_absorbed_power_w_per_m"][frequency_index],
                power_balance_error_w_per_m=power_balance_error[frequency_index],
            ))
        end
        label_text = "ld=$(length_mm) mm"
        plot!(
            amplitude_plot,
            FREQUENCIES_HZ ./ 1e3,
            amplitude;
            linewidth=2.3,
            color=colors[case_index],
            label=label_text,
        )
        plot!(
            power_plot,
            FREQUENCIES_HZ ./ 1e3,
            transmitted_power;
            linewidth=2.3,
            color=colors[case_index],
            label=label_text,
        )
        plot!(
            phase_plot,
            FREQUENCIES_HZ ./ 1e3,
            phase;
            linewidth=2.3,
            color=colors[case_index],
            label=label_text,
        )
        scatter!(
            amplitude_plot,
            [FREQUENCIES_HZ[peak_index] / 1e3, FREQUENCIES_HZ[dip_index] / 1e3],
            [amplitude[peak_index], amplitude[dip_index]];
            color=colors[case_index],
            markersize=5,
            label=false,
        )
        vline!(
            amplitude_plot,
            [eigen["full_frequency_hz"] / 1e3];
            color=colors[case_index],
            alpha=0.4,
            linestyle=:dot,
            label=false,
        )
    end
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "fano_tuning_response.csv"), rows)
    write_csv(joinpath(OUTPUT_ROOT, "fano_tuning_summary.csv"), summary)
    figure = plot(
        amplitude_plot,
        power_plot,
        phase_plot;
        layout=(3, 1),
        size=(1150, 1200),
        margin=5Plots.mm,
    )
    figure_path = joinpath(OUTPUT_ROOT, "fano_dark_mass_tuning.png")
    savefig(figure, figure_path)
    println("[+] $figure_path")
    foreach(println, summary)
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "eigen"
        run_eigen_parallel_stage()
    elseif startswith(stage, "eigen-")
        run_eigen_stage(parse(Int, stage[7:end]))
    elseif stage == "harmonic"
        run_harmonic_parallel_stage()
    elseif startswith(stage, "harmonic-")
        run_harmonic_stage(stage[10:end])
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown Fano tuning stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "eigen", "harmonic", "analyze")
            println("\n=== Fano tuning stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_fano_dark_mass_tuning.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
