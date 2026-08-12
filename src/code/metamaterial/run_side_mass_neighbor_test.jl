module SideMassNeighborTest

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(ENV, "METAMATERIALS_NEIGHBOR_OUTPUT", joinpath(PROJECT_ROOT, "tmp", "side_mass_neighbor_test"))
const ISOLATED_ROOT = get(ENV, "METAMATERIALS_NEIGHBOR_ISOLATED_ROOT", joinpath(PROJECT_ROOT, "tmp", "side_mass_refined_lossless"))
const FREQUENCIES_HZ = collect(
    parse(Float64, get(ENV, "METAMATERIALS_NEIGHBOR_FREQUENCY_START_HZ", "235000")):
    parse(Float64, get(ENV, "METAMATERIALS_NEIGHBOR_FREQUENCY_STEP_HZ", "250")):
    parse(Float64, get(ENV, "METAMATERIALS_NEIGHBOR_FREQUENCY_STOP_HZ", "250000")),
)
const DAMPING_SCALE = parse(Float64, get(ENV, "METAMATERIALS_NEIGHBOR_DAMPING_SCALE", "0"))
const PITCH_MM = 8.2

const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "side_mass_mesher.jl"))
using .SideMassMesher
const CELL_CONFIG = SideMassConfig(bright_neck_width_mm=0.7)

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif REQUESTED_STAGE in ("solve-reference", "solve-array")
    include(joinpath(@__DIR__, "harmonic_solver.jl"))
    using .HarmonicElasticity
    using JLD2
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
end

mesh_path(kind) = joinpath(OUTPUT_ROOT, "meshes", kind, "mesh_$kind.msh")
model_path(kind) = joinpath(OUTPUT_ROOT, "models", kind, "model_$kind.json")
result_path(kind) = joinpath(OUTPUT_ROOT, "harmonic", kind, "response.jld2")

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        isfile(mesh_path("reference")) || build_side_mass_array_mesh(
            mesh_path("reference");
            config=CELL_CONFIG,
            element_count=2,
            pitch_mm=PITCH_MM,
            resonators=false,
        )
        isfile(mesh_path("dimer")) || build_side_mass_array_mesh(
            mesh_path("dimer");
            config=CELL_CONFIG,
            element_count=2,
            pitch_mm=PITCH_MM,
            resonators=true,
        )
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    for kind in ("reference", "dimer")
        isfile(model_path(kind)) || convert_mesh(mesh_path(kind); output_dir=dirname(model_path(kind)))
    end
end

function save_response(kind)
    config = HarmonicConfig(
        rayleigh_alpha=DAMPING_SCALE * 79560.0,
        rayleigh_beta=DAMPING_SCALE * 2.5e-9,
    )
    points = solve_harmonic_sweep(model_path(kind), FREQUENCIES_HZ; config)
    path = result_path(kind)
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        frequency_hz=getproperty.(points, :frequency_hz),
        right_displacement=getproperty.(points, :right_displacement),
        left_displacement=getproperty.(points, :left_displacement),
    )
    println("[+] $path")
end

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_solve_stage()
    pending = [kind for kind in ("reference", "dimer") if !isfile(result_path(kind))]
    @sync for kind in pending
        @async run(child_command(kind == "reference" ? "solve-reference" : "solve-array"))
    end
end

function read_csv_rows(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1) if !isempty(strip(line))]
end

circular_error(a, b) = mod(a - b + pi, 2pi) - pi

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
    reference = JLD2.load(result_path("reference"))
    dimer = JLD2.load(result_path("dimer"))
    transfer = dimer["right_displacement"] ./ reference["right_displacement"]
    isolated_all = read_csv_rows(joinpath(ISOLATED_ROOT, "side_mass_response.csv"))
    isolated = Dict(
        parse(Float64, row["frequency_hz"]) => row
        for row in isolated_all
        if parse(Float64, row["bright_neck_width_mm"]) == 0.7
    )
    rows = [
        let isolated_row = isolated[frequency_hz],
            isolated_amplitude = parse(Float64, isolated_row["amplitude"]),
            isolated_phase = parse(Float64, isolated_row["phase_rad"]),
            dimer_amplitude = abs(transfer[index]),
            dimer_phase = angle(transfer[index])
            (
                frequency_hz,
                isolated_amplitude,
                dimer_amplitude,
                relative_amplitude_change=(dimer_amplitude - isolated_amplitude) / isolated_amplitude,
                isolated_phase_rad=isolated_phase,
                dimer_phase_rad=dimer_phase,
                phase_change_rad=circular_error(dimer_phase, isolated_phase),
            )
        end
        for (index, frequency_hz) in enumerate(FREQUENCIES_HZ)
    ]
    mkpath(OUTPUT_ROOT)
    csv_path = joinpath(OUTPUT_ROOT, "neighbor_comparison.csv")
    write_csv(csv_path, rows)

    amplitude_plot = plot(
        getproperty.(rows, :frequency_hz) ./ 1e3,
        getproperty.(rows, :isolated_amplitude);
        linewidth=2,
        label="isolated cell",
        xlabel="frequency, kHz",
        ylabel="calibrated |H|",
        title="Neighbour interaction through common end plates",
        gridalpha=0.25,
    )
    plot!(amplitude_plot, getproperty.(rows, :frequency_hz) ./ 1e3, getproperty.(rows, :dimer_amplitude); linewidth=2, label="two-cell dimer")
    error_plot = plot(
        getproperty.(rows, :frequency_hz) ./ 1e3,
        100 .* getproperty.(rows, :relative_amplitude_change);
        linewidth=2,
        label="amplitude",
        xlabel="frequency, kHz",
        ylabel="change, % / phase ×100",
        gridalpha=0.25,
    )
    plot!(error_plot, getproperty.(rows, :frequency_hz) ./ 1e3, 100 .* getproperty.(rows, :phase_change_rad); linewidth=2, label="phase, rad ×100")
    hline!(error_plot, [5, -5]; color=:gray, linestyle=:dash, label=false)
    figure = plot(amplitude_plot, error_plot; layout=(2, 1), size=(1000, 850), margin=5Plots.mm)
    figure_path = joinpath(OUTPUT_ROOT, "neighbor_comparison.png")
    savefig(figure, figure_path)
    target = rows[argmin(abs.(getproperty.(rows, :frequency_hz) .- 242e3))]
    println("[+] $csv_path")
    println("[+] $figure_path")
    println("At 242 kHz: $target")
    println("Maximum amplitude change: $(maximum(abs, getproperty.(rows, :relative_amplitude_change)))")
    println("Maximum phase change: $(maximum(abs, getproperty.(rows, :phase_change_rad))) rad")
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_solve_stage()
    elseif stage == "solve-reference"
        save_response("reference")
    elseif stage == "solve-array"
        save_response("dimer")
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown neighbour-test stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_side_mass_neighbor_test.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
