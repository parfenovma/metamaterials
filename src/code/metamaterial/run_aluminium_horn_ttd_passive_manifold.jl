module RunAluminiumHornTTDPassiveManifold

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_manifold")
const DESIGN_PATH = let provisional = joinpath(
        PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_feed", "passive_feed_design.jld2",
    )
    isfile(provisional) ? provisional : joinpath(
        PROJECT_ROOT, "results", "aluminium_horn_ttd_p6_242khz", "passive_feed",
        "design", "passive_feed_design.jld2",
    )
end

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "solve")
const FREQUENCY_HZ = parse(Float64, argument("frequency-hz", "242000"))
const ELEMENT_ORDER = parse(Int, argument("element-order", "1"))
const EVANESCENT_MODE_COUNT = parse(Int, argument("evanescent-modes", "4"))

using JLD2
using LinearAlgebra
include(joinpath(@__DIR__, "passive_feed_manifold_mesher.jl"))
using .PassiveFeedManifoldMesher

const CONFIG = PassiveFeedManifoldConfig()
const OUTPUT_TAGS = ["Output$(lpad(index, 2, '0'))" for index in 1:CONFIG.channel_count]

if STAGE == "mesh"
    using Gmsh: gmsh
elseif STAGE == "solve"
    include(joinpath(@__DIR__, "multimode_port_harmonic_solver.jl"))
    using .MultimodePortHarmonicSolver
end

mesh_path() = joinpath(OUTPUT_ROOT, "passive_feed_manifold.msh")
result_path(frequency_hz=FREQUENCY_HZ) = joinpath(
    OUTPUT_ROOT, "passive_feed_manifold_f$(round(Int, frequency_hz))hz.jld2",
)

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        result = build_passive_feed_manifold_mesh(mesh_path(); config=CONFIG)
        open(joinpath(OUTPUT_ROOT, "passive_feed_manifold_mesh_summary.csv"), "w") do io
            println(io, "nodes,elements,input_width_mm,aperture_width_mm,channel_count,channel_width_mm,channel_gap_mm,transformer_length_mm,required_transformer_length_mm,input_straight_length_mm,output_straight_length_mm")
            println(io, join((
                result.node_count, result.element_count, CONFIG.input_width_mm,
                aperture_width_mm(CONFIG), CONFIG.channel_count,
                CONFIG.channel_width_mm, CONFIG.channel_gap_mm,
                CONFIG.transformer_length_mm, required_transformer_length_mm(CONFIG),
                CONFIG.input_straight_length_mm, CONFIG.output_straight_length_mm,
            ), ','))
        end
    finally
        gmsh.finalize()
    end
end

function fundamental_amplitude(port_amplitudes)
    symmetric = filter(item -> item.parity == :symmetric, port_amplitudes)
    length(symmetric) == 1 || error(
        "expected one propagating symmetric mode, found $(length(symmetric))",
    )
    only(symmetric).right_amplitude
end

wrap_phase_deg(value) = rad2deg(atan(sin(angle(value)), cos(angle(value))))

function run_solve_stage()
    result = solve_multimode_port_harmonic(
        mesh_path(), OUTPUT_TAGS;
        config=MultimodePortConfig(
            frequency_hz=FREQUENCY_HZ,
            element_order=ELEMENT_ORDER,
            quadrature_degree=2ELEMENT_ORDER,
            port_element_count=100,
            evanescent_modes_per_port=EVANESCENT_MODE_COUNT,
        ),
    )
    fundamental = ComplexF64[fundamental_amplitude(port) for port in result.output_amplitudes]
    centre = fundamental[(length(fundamental) + 1) ÷ 2]
    phase_deg = wrap_phase_deg.(fundamental ./ centre)
    useful_power = abs2.(fundamental)
    target = Float64.(JLD2.load(DESIGN_PATH, "full_pressure_weights"))
    target ./= norm(target)
    normalized_fundamental = fundamental / norm(fundamental)
    target_complex_coherence = abs(dot(target, normalized_fundamental))
    total_transmitted = sum(result.transmitted_power_w_per_m)
    total_useful = sum(useful_power)
    balance = (
        result.reflected_power_w_per_m + total_transmitted
    ) / result.incident_power_w_per_m
    useful_transmission = total_useful / result.incident_power_w_per_m
    all_mode_transmission = total_transmitted / result.incident_power_w_per_m
    reflection = result.reflected_power_w_per_m / result.incident_power_w_per_m
    max_phase_error = maximum(abs, phase_deg)

    parity = [String(item.parity) for item in first(result.output_amplitudes)]
    right_amplitudes = reduce(hcat, [
        ComplexF64[item.right_amplitude for item in port]
        for port in result.output_amplitudes
    ])
    left_amplitudes = reduce(hcat, [
        ComplexF64[item.left_amplitude for item in port]
        for port in result.output_amplitudes
    ])
    wavenumber_per_m = ComplexF64[
        item.wavenumber_per_m for item in first(result.output_amplitudes)
    ]
    mkpath(OUTPUT_ROOT)
    JLD2.jldsave(
        result_path();
        format_version=1,
        frequency_hz=FREQUENCY_HZ,
        element_order=ELEMENT_ORDER,
        evanescent_modes_per_port=EVANESCENT_MODE_COUNT,
        incident_power_w_per_m=result.incident_power_w_per_m,
        reflected_power_w_per_m=result.reflected_power_w_per_m,
        transmitted_power_w_per_m=result.transmitted_power_w_per_m,
        incoming_output_power_w_per_m=result.incoming_output_power_w_per_m,
        fundamental_amplitude=fundamental,
        fundamental_phase_relative_center_deg=phase_deg,
        fundamental_power_w_per_m=useful_power,
        target_pressure_weights=target,
        target_complex_coherence,
        useful_power_transmission=useful_transmission,
        all_mode_power_transmission=all_mode_transmission,
        reflection_power_fraction=reflection,
        power_balance_ratio=balance,
        maximum_fundamental_phase_error_deg=max_phase_error,
        propagating_mode_parity=parity,
        propagating_mode_wavenumber_per_m=wavenumber_per_m,
        output_right_amplitudes=right_amplitudes,
        output_left_amplitudes=left_amplitudes,
        port_basis_sizes=result.port_basis_sizes,
        port_gram_condition_numbers=result.port_gram_condition_numbers,
    )
    open(joinpath(OUTPUT_ROOT, "passive_feed_manifold_f$(round(Int, FREQUENCY_HZ))hz.csv"), "w") do io
        println(io, "channel_index,center_y_mm,target_weight,fundamental_amplitude,fundamental_phase_relative_center_deg,fundamental_power_w_per_m,total_transmitted_power_w_per_m")
        for index in eachindex(fundamental)
            println(io, join((
                index, output_center_y_mm(CONFIG, index), target[index],
                abs(fundamental[index]), phase_deg[index], useful_power[index],
                result.transmitted_power_w_per_m[index],
            ), ','))
        end
    end
    println("[+] manifold R/T_all/T_useful=$(reflection) / $(all_mode_transmission) / $(useful_transmission)")
    println("[+] target complex coherence=$target_complex_coherence")
    println("[+] fundamental phase error=$max_phase_error deg")
    println("[+] power balance=$balance")
    println("[+] $(result_path())")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "solve" ? run_solve_stage() :
    error("unknown stage: $STAGE")
end

end
