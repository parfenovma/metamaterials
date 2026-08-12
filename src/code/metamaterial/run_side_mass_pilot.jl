module SideMassPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(ENV, "METAMATERIALS_SIDE_MASS_OUTPUT", joinpath(PROJECT_ROOT, "tmp", "side_mass_pilot"))
parse_values(name, default) = parse.(Float64, split(get(ENV, name, default), ','))
const FREQUENCIES_HZ = collect(
    parse(Float64, get(ENV, "METAMATERIALS_SIDE_MASS_FREQUENCY_START_HZ", "180000")):
    parse(Float64, get(ENV, "METAMATERIALS_SIDE_MASS_FREQUENCY_STEP_HZ", "2000")):
    parse(Float64, get(ENV, "METAMATERIALS_SIDE_MASS_FREQUENCY_STOP_HZ", "320000")),
)
const BRIGHT_NECK_WIDTHS_MM = parse_values("METAMATERIALS_SIDE_MASS_NECK_WIDTHS_MM", "0.4,0.6,0.8")
const DAMPING_SCALE = parse(Float64, get(ENV, "METAMATERIALS_SIDE_MASS_DAMPING_SCALE", "0"))

const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "side_mass_mesher.jl"))
using .SideMassMesher

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif REQUESTED_STAGE == "solve-reference" ||
       (!isnothing(REQUESTED_STAGE) && startswith(REQUESTED_STAGE, "solve-one-"))
    include(joinpath(@__DIR__, "harmonic_solver.jl"))
    using .HarmonicElasticity
    using JLD2
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
end

struct SideMassCase
    id::String
    neck_width_mm::Float64
    config::SideMassConfig
end

label(value) = replace(string(Float64(value)), "." => "p")

function side_mass_case(neck_width_mm)
    width = Float64(neck_width_mm)
    SideMassCase(
        "side_mass_Tb_$(label(width))",
        width,
        SideMassConfig(bright_neck_width_mm=width),
    )
end

side_mass_cases() = side_mass_case.(BRIGHT_NECK_WIDTHS_MM)
reference_config() = SideMassConfig()
case_dir(kind, case) = joinpath(OUTPUT_ROOT, kind, case.id)
reference_dir(kind) = joinpath(OUTPUT_ROOT, kind, "reference")
mesh_path(case) = joinpath(case_dir("meshes", case), "mesh_cell.msh")
model_path(case) = joinpath(case_dir("models", case), "model_cell.json")
result_path(case) = joinpath(case_dir("harmonic", case), "response.jld2")
reference_mesh_path() = joinpath(reference_dir("meshes"), "mesh_reference.msh")
reference_model_path() = joinpath(reference_dir("models"), "model_reference.json")
reference_result_path() = joinpath(reference_dir("harmonic"), "response.jld2")

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        isfile(reference_mesh_path()) || build_side_mass_mesh(
            reference_mesh_path();
            config=reference_config(),
            resonators=false,
        )
        for case in side_mass_cases()
            isfile(mesh_path(case)) || build_side_mass_mesh(mesh_path(case); config=case.config)
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    isfile(reference_model_path()) || convert_mesh(reference_mesh_path(); output_dir=reference_dir("models"))
    for case in side_mass_cases()
        isfile(model_path(case)) || convert_mesh(mesh_path(case); output_dir=case_dir("models", case))
    end
end

function harmonic_config()
    HarmonicConfig(
        rayleigh_alpha=DAMPING_SCALE * 79560.0,
        rayleigh_beta=DAMPING_SCALE * 2.5e-9,
    )
end

function save_response(path, model)
    points = solve_harmonic_sweep(model, FREQUENCIES_HZ; config=harmonic_config())
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

run_reference_stage() = isfile(reference_result_path()) ? println("[=] Existing reference") : save_response(reference_result_path(), reference_model_path())

function run_single_stage(index)
    cases = side_mass_cases()
    1 <= index <= length(cases) || error("invalid side-mass case index: $index")
    case = cases[index]
    isfile(result_path(case)) ? println("[=] Existing response: $(case.id)") : save_response(result_path(case), model_path(case))
end

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_solve_stage()
    run(child_command("solve-reference"))
    cases = side_mass_cases()
    pending = findall(case -> !isfile(result_path(case)), cases)
    isempty(pending) && return println("[=] All side-mass responses exist")
    jobs = min(parse(Int, get(ENV, "METAMATERIALS_GEOMETRY_JOBS", "4")), length(pending), Sys.CPU_THREADS)
    semaphore = Base.Semaphore(jobs)
    println("=== Side-mass harmonic pilot: $(length(pending)) cases, $jobs processes ===")
    @sync for index in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command("solve-one-$index"))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function unwrap_phase(phases)
    result = Float64[first(phases)]
    previous_raw = first(phases)
    for raw in Iterators.drop(phases, 1)
        push!(result, last(result) + mod(raw - previous_raw + pi, 2pi) - pi)
        previous_raw = raw
    end
    result
end

function response_rows(case, reference)
    sample = JLD2.load(result_path(case))
    transfer = sample["right_displacement"] ./ reference["right_displacement"]
    phases = unwrap_phase(angle.(transfer))
    omega = 2pi .* FREQUENCIES_HZ
    delays = similar(phases)
    delays[1] = -(phases[2] - phases[1]) / (omega[2] - omega[1])
    for index in 2:(length(phases) - 1)
        delays[index] = -(phases[index + 1] - phases[index - 1]) / (omega[index + 1] - omega[index - 1])
    end
    delays[end] = -(phases[end] - phases[end - 1]) / (omega[end] - omega[end - 1])
    [
        (
            case_id=case.id,
            bright_neck_width_mm=case.neck_width_mm,
            frequency_hz=frequency_hz,
            H_real=real(transfer[index]),
            H_imag=imag(transfer[index]),
            amplitude=abs(transfer[index]),
            phase_rad=phases[index],
            group_delay_s=delays[index],
        )
        for (index, frequency_hz) in enumerate(FREQUENCIES_HZ)
    ]
end

function extrema(rows)
    result = NamedTuple[]
    amplitude = getproperty.(rows, :amplitude)
    for index in 2:(length(rows) - 1)
        kind = if amplitude[index] > amplitude[index - 1] && amplitude[index] >= amplitude[index + 1]
            "peak"
        elseif amplitude[index] < amplitude[index - 1] && amplitude[index] <= amplitude[index + 1]
            "dip"
        else
            continue
        end
        push!(result, (
            kind,
            frequency_hz=rows[index].frequency_hz,
            amplitude=amplitude[index],
            group_delay_s=rows[index].group_delay_s,
        ))
    end
    result
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

function save_geometry_plot()
    figure = plot(
        aspect_ratio=:equal,
        xlabel="x, mm",
        ylabel="y, mm",
        title="Side-mass bright--dark cell",
        legend=false,
        xlims=(-0.5, 30.5),
        ylims=(-0.5, 7.5),
        size=(1100, 360),
    )
    for (x, y, width, height) in SideMassMesher.solid_rectangles(SideMassConfig())
        plot!(figure, Shape([x, x + width, x + width, x], [y, y, y + height, y + height]); color=:steelblue, linecolor=:steelblue)
    end
    path = joinpath(OUTPUT_ROOT, "side_mass_geometry.png")
    savefig(figure, path)
    path
end

function save_spectra(rows)
    amplitude_plot = plot(
        xlabel="frequency, kHz",
        ylabel="calibrated |H|",
        title="Side-mass cell, damping scale=$(DAMPING_SCALE)",
        gridalpha=0.25,
        legend=:outertopright,
    )
    delay_plot = plot(
        xlabel="frequency, kHz",
        ylabel="phase-derived delay, us",
        gridalpha=0.25,
        legend=:outertopright,
    )
    for case in side_mass_cases()
        current = filter(row -> row.case_id == case.id, rows)
        plot!(amplitude_plot, getproperty.(current, :frequency_hz) ./ 1e3, getproperty.(current, :amplitude); linewidth=2, label="tb=$(case.neck_width_mm) mm")
        reliable_delay = [row.amplitude >= 0.05 ? row.group_delay_s * 1e6 : NaN for row in current]
        plot!(delay_plot, getproperty.(current, :frequency_hz) ./ 1e3, reliable_delay; linewidth=2, label="tb=$(case.neck_width_mm) mm")
    end
    vline!(amplitude_plot, [242.0]; color=:gray, linestyle=:dash, label=false)
    vline!(delay_plot, [242.0]; color=:gray, linestyle=:dash, label=false)
    figure = plot(amplitude_plot, delay_plot; layout=(2, 1), size=(1050, 920), margin=5Plots.mm)
    path = joinpath(OUTPUT_ROOT, "side_mass_spectra.png")
    savefig(figure, path)
    path
end

function run_analysis_stage()
    reference = JLD2.load(reference_result_path())
    rows = reduce(vcat, response_rows(case, reference) for case in side_mass_cases())
    extrema_rows = NamedTuple[]
    for case in side_mass_cases()
        current = filter(row -> row.case_id == case.id, rows)
        for point in extrema(current)
            push!(extrema_rows, (
                case_id=case.id,
                bright_neck_width_mm=case.neck_width_mm,
                kind=point.kind,
                frequency_hz=point.frequency_hz,
                amplitude=point.amplitude,
                group_delay_s=point.group_delay_s,
            ))
        end
    end
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "side_mass_response.csv"), rows)
    isempty(extrema_rows) || write_csv(joinpath(OUTPUT_ROOT, "side_mass_extrema.csv"), extrema_rows)
    println("[+] $(save_geometry_plot())")
    println("[+] $(save_spectra(rows))")
    foreach(println, extrema_rows)
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "solve"
        run_solve_stage()
    elseif stage == "solve-reference"
        run_reference_stage()
    elseif startswith(stage, "solve-one-")
        run_single_stage(parse(Int, split(stage, '-')[3]))
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown side-mass stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "solve", "analyze")
            println("\n=== Stage: $stage ===")
            run(child_command(stage))
        end
    elseif args == ["--list"]
        foreach(case -> println(case.id), side_mass_cases())
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_side_mass_pilot.jl [--list | --stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
