module RunAluminiumHornTTDPassiveSplitterModal

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const VARIANT = argument("variant", "calibrated_beat")
const FREQUENCY_HZ = parse(Float64, argument("frequency-hz", "242000"))
const EVANESCENT_MODE_COUNT = parse(Int, argument("evanescent-modes", "6"))
const ROOT_NAMES = Dict(
    "straight" => "aluminium_horn_ttd_passive_straight_control",
    "equal_long" => "aluminium_horn_ttd_passive_equal_long_splitter",
    "calibrated_beat" => "aluminium_horn_ttd_passive_calibrated_beat_splitter",
)
haskey(ROOT_NAMES, VARIANT) || error("unknown modal splitter variant: $VARIANT")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", ROOT_NAMES[VARIANT])
const OUTPUT_TAGS = VARIANT == "straight" ? ["Output"] : ["OutputSmall", "OutputLarge"]

using JLD2
include(joinpath(@__DIR__, "multimode_port_harmonic_solver.jl"))
using .MultimodePortHarmonicSolver

mesh_path() = joinpath(OUTPUT_ROOT, "passive_feed_splitter.msh")
result_path() = joinpath(
    OUTPUT_ROOT,
    "passive_feed_splitter_modal_f$(round(Int, FREQUENCY_HZ))hz.jld2",
)

function fundamental_amplitude(port_amplitudes)
    symmetric = filter(item -> item.parity == :symmetric, port_amplitudes)
    length(symmetric) == 1 || error("expected exactly one propagating symmetric mode")
    only(symmetric).right_amplitude
end

function run()
    result = solve_multimode_port_harmonic(
        mesh_path(), OUTPUT_TAGS;
        config=MultimodePortConfig(
            frequency_hz=FREQUENCY_HZ,
            element_order=2,
            quadrature_degree=4,
            port_element_count=100,
            evanescent_modes_per_port=EVANESCENT_MODE_COUNT,
        ),
    )
    fundamental = ComplexF64[
        fundamental_amplitude(port) for port in result.output_amplitudes
    ]
    useful_power = abs2.(fundamental)
    all_mode_transmission = sum(result.transmitted_power_w_per_m) /
                            result.incident_power_w_per_m
    useful_transmission = sum(useful_power) / result.incident_power_w_per_m
    reflection = result.reflected_power_w_per_m / result.incident_power_w_per_m
    power_balance = (
        result.reflected_power_w_per_m + sum(result.transmitted_power_w_per_m)
    ) / result.incident_power_w_per_m
    amplitude_ratio = length(fundamental) == 2 ?
        abs(fundamental[1] / fundamental[2]) : 1.0
    relative_phase_deg = length(fundamental) == 2 ?
        rad2deg(angle(fundamental[1] / fundamental[2])) : 0.0
    parity = [String(item.parity) for item in first(result.output_amplitudes)]
    right_amplitudes = reduce(hcat, [
        ComplexF64[item.right_amplitude for item in port]
        for port in result.output_amplitudes
    ])
    left_amplitudes = reduce(hcat, [
        ComplexF64[item.left_amplitude for item in port]
        for port in result.output_amplitudes
    ])
    JLD2.jldsave(
        result_path();
        format_version=1,
        variant=VARIANT,
        frequency_hz=FREQUENCY_HZ,
        incident_power_w_per_m=result.incident_power_w_per_m,
        reflected_power_w_per_m=result.reflected_power_w_per_m,
        transmitted_power_w_per_m=result.transmitted_power_w_per_m,
        incoming_output_power_w_per_m=result.incoming_output_power_w_per_m,
        fundamental_amplitude=fundamental,
        fundamental_power_w_per_m=useful_power,
        reflection_power_fraction=reflection,
        all_mode_power_transmission=all_mode_transmission,
        useful_fundamental_power_transmission=useful_transmission,
        power_balance_ratio=power_balance,
        small_to_large_fundamental_amplitude_ratio=amplitude_ratio,
        small_to_large_fundamental_phase_deg=relative_phase_deg,
        propagating_mode_parity=parity,
        output_right_amplitudes=right_amplitudes,
        output_left_amplitudes=left_amplitudes,
        port_basis_sizes=result.port_basis_sizes,
        port_gram_condition_numbers=result.port_gram_condition_numbers,
    )
    println("[+] $VARIANT R/T/Tfund=$reflection / $all_mode_transmission / $useful_transmission")
    println("[+] ratio/phase=$amplitude_ratio / $relative_phase_deg deg")
    println("[+] power balance=$power_balance")
    println("[+] $(result_path())")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
