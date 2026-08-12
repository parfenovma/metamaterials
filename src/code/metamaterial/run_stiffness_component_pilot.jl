module StiffnessComponentPilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_STIFFNESS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "stiffness_component_pilot"),
)
const TARGET_FREQUENCY_HZ = 242.0e3
const WORKING_BAND_HZ = (222.0e3, 262.0e3)
const MAX_WORKERS = parse(Int, get(ENV, "METAMATERIALS_STIFFNESS_WORKERS", "4"))
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "stiffness_component_mesher.jl"))
using .StiffnessComponentMesher

is_solver_stage(stage) = startswith(stage, "eigen-") || startswith(stage, "response-")

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) &&
       (is_solver_stage(REQUESTED_STAGE) || REQUESTED_STAGE == "analyze")
    include(joinpath(@__DIR__, "modal_solver.jl"))
    include(joinpath(@__DIR__, "dynamic_stiffness_solver.jl"))
    using .ConservativeElasticModes
    using .DynamicStiffnessSolver
    using JLD2
end

if REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using Plots
end

Base.@kwdef struct PilotCase
    id::String
    variant::Symbol
    config::StiffnessConfig
end

function pilot_cases()
    cases = [PilotCase(id="solid", variant=:solid, config=StiffnessConfig())]
    for ligament in (0.40, 0.50, 0.60, 0.70, 0.75)
        code = replace(string(round(ligament; digits=2)), "." => "p")
        push!(cases, PilotCase(
            id="ligament_$(code)mm",
            variant=:slotted,
            config=StiffnessConfig(ligament_height_mm=ligament),
        ))
    end
    cases
end

frequencies_hz() = sort(unique(vcat(0.0, collect(200.0e3:5.0e3:280.0e3), TARGET_FREQUENCY_HZ)))

mesh_path(case) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(case.id).msh")
model_path(case) = joinpath(OUTPUT_ROOT, "models", "model_$(case.id).json")
eigen_path(case) = joinpath(OUTPUT_ROOT, "eigen_$(case.id).jld2")
response_path(case) = joinpath(OUTPUT_ROOT, "response_$(case.id).jld2")

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for case in pilot_cases()
            build_stiffness_component_mesh(
                mesh_path(case);
                config=case.config,
                variant=case.variant,
                size_min_mm=0.045,
                size_max_mm=0.15,
            )
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    for case in pilot_cases()
        convert_mesh(mesh_path(case); output_dir=dirname(model_path(case)))
    end
end

function solve_eigen_case(index)
    cases = pilot_cases()
    1 <= index <= length(cases) || error("case index out of range")
    case = cases[index]
    modes = solve_conservative_modes(
        model_path(case);
        config=ModeConfig(
            element_order=2,
            quadrature_degree=4,
            target_frequency_hz=TARGET_FREQUENCY_HZ,
            mode_count=12,
            tolerance=1.0e-9,
            clamp_tags=["FixedInterface"],
        ),
    )
    mkpath(OUTPUT_ROOT)
    JLD2.jldsave(eigen_path(case); modes, case_id=case.id)
    symmetric = filter(mode -> mode.parity_y >= 0.8, modes)
    nearest = modes[argmin(abs(mode.frequency_hz - TARGET_FREQUENCY_HZ) for mode in modes)]
    nearest_symmetric = symmetric[
        argmin(abs(mode.frequency_hz - TARGET_FREQUENCY_HZ) for mode in symmetric)
    ]
    println(
        "[+] $(case.id): nearest=", round(nearest.frequency_hz / 1e3; digits=3),
        " kHz (parity=", round(nearest.parity_y; digits=3),
        "), symmetric=", round(nearest_symmetric.frequency_hz / 1e3; digits=3), " kHz",
    )
end

function solve_response_case(index)
    cases = pilot_cases()
    1 <= index <= length(cases) || error("case index out of range")
    case = cases[index]
    points = solve_dynamic_stiffness_sweep(
        model_path(case),
        frequencies_hz();
        component_length_m=case.config.length_mm * 1e-3,
        config=DynamicStiffnessConfig(),
    )
    mkpath(OUTPUT_ROOT)
    JLD2.jldsave(response_path(case); points, case_id=case.id)
    static_point = only(filter(point -> point.frequency_hz == 0, points))
    target_point = only(filter(point -> point.frequency_hz == TARGET_FREQUENCY_HZ, points))
    println(
        "[+] $(case.id): K0=", real(static_point.stiffness_n_per_m2),
        ", K242=", real(target_point.stiffness_n_per_m2),
    )
end

function run_parallel_stage(prefix)
    indices = collect(eachindex(pilot_cases()))
    for batch_start in 1:MAX_WORKERS:length(indices)
        batch = indices[batch_start:min(batch_start + MAX_WORKERS - 1, end)]
        @sync for index in batch
            @async run(child_command("$prefix-$index"))
        end
    end
end

function closest_mode(modes; symmetric=false)
    candidates = symmetric ? filter(mode -> mode.parity_y >= 0.8, modes) : modes
    isempty(candidates) && error("no mode satisfying requested symmetry")
    candidates[argmin(abs(mode.frequency_hz - TARGET_FREQUENCY_HZ) for mode in candidates)]
end

function loaded_rows()
    cases = pilot_cases()
    response_data = Dict(case.id => JLD2.load(response_path(case))["points"] for case in cases)
    eigen_data = Dict(case.id => JLD2.load(eigen_path(case))["modes"] for case in cases)
    solid_points = response_data["solid"]
    map(cases) do case
        points = response_data[case.id]
        static_point = only(filter(point -> point.frequency_hz == 0, points))
        target_point = only(filter(point -> point.frequency_hz == TARGET_FREQUENCY_HZ, points))
        solid_static = only(filter(point -> point.frequency_hz == 0, solid_points))
        solid_target = only(filter(point -> point.frequency_hz == TARGET_FREQUENCY_HZ, solid_points))
        modes = eigen_data[case.id]
        nearest = closest_mode(modes)
        nearest_symmetric = closest_mode(modes; symmetric=true)
        band_points = filter(
            point -> first(WORKING_BAND_HZ) <= point.frequency_hz <= last(WORKING_BAND_HZ),
            points,
        )
        band_stiffness = real.(getproperty.(band_points, :stiffness_n_per_m2))
        (
            case,
            points,
            modes,
            static_stiffness_ratio=real(static_point.stiffness_n_per_m2 /
                                        solid_static.stiffness_n_per_m2),
            target_stiffness_ratio=real(target_point.stiffness_n_per_m2 /
                                        solid_target.stiffness_n_per_m2),
            target_compliance_ratio=real(target_point.compliance_m2_per_n /
                                         solid_target.compliance_m2_per_n),
            band_min_stiffness=minimum(band_stiffness),
            band_max_stiffness=maximum(band_stiffness),
            band_stiffness_variation=maximum(band_stiffness) / minimum(band_stiffness),
            nearest,
            nearest_symmetric,
        )
    end
end

function write_summary(rows)
    path = joinpath(OUTPUT_ROOT, "stiffness_component_summary.csv")
    open(path, "w") do io
        println(io, "case_id,variant,ligament_height_mm,slot_height_mm,area_fraction,K0_over_solid,K242_over_solid,C242_over_solid,band_min_stiffness_n_per_m2,band_max_stiffness_n_per_m2,band_stiffness_variation,nonresonant_band_gate,nearest_mode_hz,nearest_mode_parity_y,nearest_symmetric_mode_hz,symmetric_distance_from_242_hz")
        for row in rows
            config = row.case.config
            println(io, join((
                row.case.id,
                row.case.variant,
                config.ligament_height_mm,
                row.case.variant == :solid ? 0.0 : slot_height_mm(config),
                row.case.variant == :solid ? 1.0 : ligament_area_fraction(config),
                row.static_stiffness_ratio,
                row.target_stiffness_ratio,
                row.target_compliance_ratio,
                row.band_min_stiffness,
                row.band_max_stiffness,
                row.band_stiffness_variation,
                row.band_min_stiffness > 0 && row.band_stiffness_variation <= 2.0,
                row.nearest.frequency_hz,
                row.nearest.parity_y,
                row.nearest_symmetric.frequency_hz,
                abs(row.nearest_symmetric.frequency_hz - TARGET_FREQUENCY_HZ),
            ), ','))
        end
    end
    println("[+] $path")
end

function compliance_panel(rows)
    panel = plot(
        xlabel="frequency, kHz",
        ylabel="Re(C_K / C_solid)",
        title="Dynamic axial compliance",
        gridalpha=0.25,
        legend=:topleft,
        ylims=(-2, 12),
    )
    solid = only(filter(row -> row.case.variant == :solid, rows))
    solid_lookup = Dict(point.frequency_hz => point for point in solid.points)
    for row in filter(row -> row.case.variant == :slotted, rows)
        frequencies = getproperty.(row.points, :frequency_hz)
        ratio = [
            real(point.compliance_m2_per_n /
                 solid_lookup[point.frequency_hz].compliance_m2_per_n)
            for point in row.points
        ]
        plot!(
            panel,
            frequencies ./ 1e3,
            ratio;
            linewidth=2,
            label="t=$(row.case.config.ligament_height_mm) mm",
        )
    end
    vspan!(panel, collect(WORKING_BAND_HZ) ./ 1e3; color=:gray, alpha=0.10, label="working band")
    vline!(panel, [TARGET_FREQUENCY_HZ / 1e3]; color=:black, linestyle=:dash, label="242 kHz")
    panel
end

function static_panel(rows)
    slotted = sort(filter(row -> row.case.variant == :slotted, rows);
                   by=row -> row.case.config.ligament_height_mm)
    plot(
        [row.case.config.ligament_height_mm for row in slotted],
        getproperty.(slotted, :static_stiffness_ratio);
        marker=:circle,
        linewidth=2,
        xlabel="ligament height, mm",
        ylabel="K(0) / K_solid(0)",
        title="Quasi-static stiffness control",
        label=false,
        gridalpha=0.25,
    )
end

function spectrum_panel(rows)
    panel = plot(
        xlabel="ligament height, mm",
        ylabel="fixed-interface mode, kHz",
        title="Component poles near working band",
        gridalpha=0.25,
        legend=:topleft,
    )
    for row in filter(row -> row.case.variant == :slotted, rows)
        x = fill(row.case.config.ligament_height_mm, length(row.modes))
        parity = getproperty.(row.modes, :parity_y)
        scatter!(
            panel,
            x,
            getproperty.(row.modes, :frequency_hz) ./ 1e3;
            markercolor=[value >= 0.8 ? :royalblue : :transparent for value in parity],
            markerstrokecolor=:royalblue,
            markersize=5,
            label=false,
        )
    end
    hspan!(panel, collect(WORKING_BAND_HZ) ./ 1e3; color=:gray, alpha=0.10, label="working band")
    hline!(panel, [TARGET_FREQUENCY_HZ / 1e3]; color=:black, linestyle=:dash, label="242 kHz")
    panel
end

function geometry_panel(config)
    panel = plot(
        aspect_ratio=:equal,
        xlims=(-0.2, config.length_mm + 0.2),
        ylims=(-config.height_mm / 2 - 0.2, config.height_mm / 2 + 0.2),
        xlabel="x, mm",
        ylabel="y, mm",
        title="Fixed topology: three rounded slots",
        legend=false,
        grid=false,
        framestyle=:box,
    )
    plot!(panel, Shape(
        [0.0, config.length_mm, config.length_mm, 0.0],
        [-config.height_mm / 2, -config.height_mm / 2,
         config.height_mm / 2, config.height_mm / 2],
    ); color=:lightgray, linecolor=:black)
    for (x, y, length, height) in slot_boxes(config)
        radius = height / 2
        theta_left = range(pi / 2, 3pi / 2; length=40)
        theta_right = range(-pi / 2, pi / 2; length=40)
        xs = vcat(
            x + radius .+ radius .* cos.(theta_left),
            x + length - radius .+ radius .* cos.(theta_right),
        )
        ys = vcat(
            y + radius .+ radius .* sin.(theta_left),
            y + radius .+ radius .* sin.(theta_right),
        )
        plot!(panel, Shape(xs, ys); color=:white, linecolor=:black)
    end
    panel
end

function save_figure(rows)
    candidate = only(filter(row -> row.case.id == "ligament_0p4mm", rows))
    figure = plot(
        geometry_panel(candidate.case.config),
        static_panel(rows),
        compliance_panel(rows),
        spectrum_panel(rows);
        layout=(2, 2),
        size=(1500, 950),
        margin=5Plots.mm,
        plot_title="Non-resonant stiffness component K/2 — quadratic FEM",
    )
    path = joinpath(OUTPUT_ROOT, "stiffness_component_pilot.png")
    savefig(figure, path)
    println("[+] $path")
end

function run_analysis_stage()
    rows = loaded_rows()
    write_summary(rows)
    save_figure(rows)
    println("\nSlotted K candidates at 242 kHz:")
    for row in filter(row -> row.case.variant == :slotted, rows)
        println(
            "  t=", row.case.config.ligament_height_mm,
            " mm: K0/Ksolid=", round(row.static_stiffness_ratio; digits=3),
            ", C242/Csolid=", round(row.target_compliance_ratio; digits=3),
            ", band variation=", round(row.band_stiffness_variation; digits=2),
            ", gate=", row.band_min_stiffness > 0 && row.band_stiffness_variation <= 2.0,
            ", nearest symmetric mode=",
            round(row.nearest_symmetric.frequency_hz / 1e3; digits=2), " kHz",
        )
    end
end


function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "eigen"
        run_parallel_stage("eigen")
    elseif stage == "response"
        run_parallel_stage("response")
    elseif startswith(stage, "eigen-")
        solve_eigen_case(parse(Int, stage[7:end]))
    elseif startswith(stage, "response-")
        solve_response_case(parse(Int, stage[10:end]))
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown stiffness-pilot stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "eigen", "response", "analyze")
            println("\n=== Stiffness component stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_stiffness_component_pilot.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
