module HornNeighbor3DGate

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_NEIGHBOR_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_neighbor_3d_242khz"),
)

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "report")
const CASE_NAME = argument("case", nothing)

using Gmsh: gmsh
include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
include(joinpath(@__DIR__, "horn_neighbor_3d_mesher.jl"))
using .HornNeighbor3DMesher

if STAGE == "solve"
    using JLD2
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
    include(joinpath(@__DIR__, "horn_neighbor_3d_harmonic_solver.jl"))
    using .SinusoidalMaterialLens
    using .HornNeighbor3DHarmonicSolver
elseif STAGE == "report"
    using JLD2
end

const CONFIG = HornNeighbor3DConfig()
const LONG_CONFIG = HornNeighbor3DConfig(
    device_guide_length_mm=38.158058,
    gentle_guide_axial_length_mm=38.158058,
    gentle_guide_path_length_mm=42.578058,
)
const CASES = Dict(
    "single_device" => (states=[:device], config=CONFIG),
    "single_gentle" => (states=[:gentle], config=CONFIG),
    "single_smooth" => (states=[:smooth], config=CONFIG),
    "DDD" => (states=[:device, :device, :device], config=CONFIG),
    "DGD" => (states=[:device, :gentle, :device], config=CONFIG),
    "DSD" => (states=[:device, :smooth, :device], config=CONFIG),
    "long_single_device" => (states=[:device], config=LONG_CONFIG),
    "long_single_smooth" => (states=[:smooth], config=LONG_CONFIG),
    "long_DDD" => (states=[:device, :device, :device], config=LONG_CONFIG),
    "long_DSD" => (states=[:device, :smooth, :device], config=LONG_CONFIG),
    "long_DSS" => (states=[:device, :smooth, :smooth], config=LONG_CONFIG),
    "long_SDD" => (states=[:smooth, :device, :device], config=LONG_CONFIG),
    "long_SDS" => (states=[:smooth, :device, :smooth], config=LONG_CONFIG),
    "long_SSS" => (states=[:smooth, :smooth, :smooth], config=LONG_CONFIG),
)

mesh_path(case_name) = joinpath(OUTPUT_ROOT, "mesh_$(case_name).msh")
result_path(case_name) = joinpath(OUTPUT_ROOT, "$(case_name)_harmonic.jld2")

function selected_case()
    isnothing(CASE_NAME) && error("--case is required for stage $STAGE")
    haskey(CASES, CASE_NAME) || error("unknown case: $CASE_NAME")
    CASE_NAME, CASES[CASE_NAME]
end

function run_mesh_stage()
    case_name, case = selected_case()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        build_horn_neighbor_3d_mesh(
            mesh_path(case_name);
            config=case.config,
            states=case.states,
            size_path_mm=0.65,
            size_receiver_mm=1.40,
        )
    finally
        gmsh.finalize()
    end
end

function run_solve_stage()
    case_name, case = selected_case()
    probes = probe_points_mm(case.config, case.states)
    if case_name in ("long_DSS", "long_SDS")
        selected = findall(startswith("preout_"), probes.names)
        probes = (
            names=probes.names[selected],
            x_mm=probes.x_mm[selected],
            y_mm=probes.y_mm[selected],
            z_mm=probes.z_mm[selected],
        )
    end
    result = solve_horn_neighbor_3d_harmonic(
        mesh_path(case_name),
        length(case.states),
        probes,
        photopolymer();
        config=HornNeighbor3DHarmonicConfig(),
    )
    jldsave(
        result_path(case_name);
        format_version=1,
        case_name,
        states=String.(case.states),
        frequency_hz=242.0e3,
        probe_names=result.probe_names,
        probe_x_mm=result.probe_x_mm,
        probe_y_mm=result.probe_y_mm,
        probe_z_mm=result.probe_z_mm,
        probe_displacement_m=result.probe_displacement_m,
    )
    println("[+] $(result_path(case_name))")
end

function probe(data, name)
    index = findfirst(==(name), data["probe_names"])
    isnothing(index) && error("probe $name was not found")
    vec(data["probe_displacement_m"][index, :])
end

phase_error_deg(value) = rad2deg(atan(sin(angle(value)), cos(angle(value))))

function transfer_row(name, value)
    (name=name, amplitude=abs(value), phase_deg=rad2deg(angle(value)))
end

function run_report_stage()
    data = Dict(case_name => JLD2.load(result_path(case_name)) for case_name in keys(CASES))
    rows = NamedTuple[]
    verdicts = NamedTuple[]
    comparisons = (
        (state_name="gentle", singleton="single_gentle", array="DGD",
         device="single_device", control="DDD"),
        (state_name="smooth", singleton="single_smooth", array="DSD",
         device="single_device", control="DDD"),
        (state_name="long_smooth", singleton="long_single_smooth", array="long_DSD",
         device="long_single_device", control="long_DDD"),
    )
    for comparison in comparisons
        state_name = comparison.state_name
        array_name = comparison.array
        isolated = probe(data[comparison.singleton], "preout_1")[1] /
                   probe(data[comparison.device], "preout_1")[1]
        neighbor = probe(data[array_name], "preout_2")[1] /
                   probe(data[comparison.control], "preout_2")[1]
        perturbation = neighbor / isolated
        receiver_isolated = probe(data[comparison.singleton], "receiver_1")[1] /
                            probe(data[comparison.device], "receiver_1")[1]
        receiver_neighbor = probe(data[array_name], "receiver_2")[1] /
                            probe(data[comparison.control], "receiver_2")[1]
        outer_left = probe(data[array_name], "preout_1")[1] /
                     probe(data[comparison.control], "preout_1")[1]
        outer_right = probe(data[array_name], "preout_3")[1] /
                      probe(data[comparison.control], "preout_3")[1]
        center_receiver = probe(data[array_name], "receiver_2")
        transverse_energy_ratio = (abs2(center_receiver[2]) + abs2(center_receiver[3])) /
                                  abs2(center_receiver[1])
        isolated_gate = abs(isolated) >= 0.75
        coupling_amplitude_gate = abs(abs(perturbation) - 1) <= 0.15
        coupling_phase_gate = abs(phase_error_deg(perturbation)) <= 15.0
        outer_gate = all(value -> abs(abs(value) - 1) <= 0.15 && abs(phase_error_deg(value)) <= 15.0,
                         (outer_left, outer_right))
        purity_gate = transverse_energy_ratio <= 0.10
        passed = isolated_gate && coupling_amplitude_gate && coupling_phase_gate && outer_gate && purity_gate
        append!(rows, [
            transfer_row("isolated_$(state_name)_over_device_preout", isolated),
            transfer_row("neighbor_$(array_name)_over_DDD_center_preout", neighbor),
            transfer_row("$(state_name)_neighbor_perturbation_over_isolated", perturbation),
            transfer_row("isolated_$(state_name)_over_device_receiver", receiver_isolated),
            transfer_row("neighbor_$(array_name)_over_DDD_center_receiver", receiver_neighbor),
            transfer_row("$(array_name)_over_DDD_left_device_preout", outer_left),
            transfer_row("$(array_name)_over_DDD_right_device_preout", outer_right),
            (name="$(array_name)_center_transverse_energy_over_longitudinal",
             amplitude=transverse_energy_ratio, phase_deg=NaN),
        ])
        push!(verdicts, (; state_name, array_name, isolated_gate, coupling_amplitude_gate,
                          coupling_phase_gate, outer_gate, purity_gate, passed,
                          transverse_energy_ratio))
    end
    mkpath(OUTPUT_ROOT)
    summary_path = joinpath(OUTPUT_ROOT, "horn_neighbor_3d_summary.csv")
    open(summary_path, "w") do io
        println(io, "metric,amplitude,phase_deg")
        for row in rows
            println(io, "$(row.name),$(row.amplitude),$(row.phase_deg)")
        end
    end
    verdict_path = joinpath(OUTPUT_ROOT, "horn_neighbor_3d_verdict.txt")
    open(verdict_path, "w") do io
        for verdict in verdicts
            println(io, "[$(verdict.state_name)]")
            for key in (:passed, :isolated_gate, :coupling_amplitude_gate,
                        :coupling_phase_gate, :outer_gate, :purity_gate,
                        :transverse_energy_ratio)
                println(io, "$key=$(getproperty(verdict, key))")
            end
        end
    end
    for row in rows
        println("[+] $(row.name): $(row.amplitude) ∠ $(row.phase_deg) deg")
    end
    for verdict in verdicts
        println("[+] $(verdict.state_name) 3D neighbour gate passed=$(verdict.passed)")
    end
    println("[+] $summary_path")
    println("[+] $verdict_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "solve" ? run_solve_stage() :
    STAGE == "report" ? run_report_stage() :
    error("unknown stage: $STAGE")
end

end
