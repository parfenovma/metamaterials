module SweepAluminiumHornTTDPulseCycles

using JLD2

include(joinpath(@__DIR__, "optimize_aluminium_horn_ttd_broadband_weights.jl"))
using .OptimizeAluminiumHornTTDBroadbandWeights
const O = OptimizeAluminiumHornTTDBroadbandWeights

const CYCLES = Float64[5, 6, 7, 8, 10, 12, 15, 20, 30, 50]

function main()
    matrices = JLD2.load.(O.matrix_path.(O.FREQUENCIES_HZ))
    references = JLD2.load.(O.reference_path.(O.FREQUENCIES_HZ))
    focus_response = [vec(item["focus_response_m"][1, :]) for item in matrices]
    power_matrix = [O.active_power_matrix(
        item["source_normal_displacement_integral_m3"], O.FREQUENCIES_HZ[index],
    ) for (index, item) in enumerate(matrices)]
    reference_focus = ComplexF64[item["focus_displacement_m"][1] for item in references]
    reference_power = abs.(Float64[item["active_input_power_w"] for item in references])
    carrier_weights = Float64.(JLD2.load(joinpath(
        O.LENS_ROOT, "15_power_optimal_weights.jld2",
    ))["pressure_weights"])
    five_cycle_weights = Float64.(JLD2.load(joinpath(
        O.LENS_ROOT, "15_broadband_optimal_weights.jld2",
    ))["pressure_weights"])
    single = [O.optimize_single_frequency(
        focus_response[index], power_matrix[index], reference_power[index],
    ) for index in eachindex(O.FREQUENCIES_HZ)]
    base_starts = vcat(
        [carrier_weights, five_cycle_weights, ones(8)],
        [item.weights for item in single],
    )
    rows = NamedTuple[]
    previous_weights = five_cycle_weights
    for cycles in CYCLES
        coarse = O.build_objective(
            focus_response, power_matrix, reference_focus, reference_power;
            cycles, samples_per_period=60, fft_length=16384,
        )
        starts = vcat([previous_weights], base_starts)
        searches = [O.pattern_search(coarse, start) for start in starts]
        optimum = searches[argmax(getproperty.(searches, :gain))]
        fine = O.build_objective(
            focus_response, power_matrix, reference_focus, reference_power;
            cycles, samples_per_period=80, fft_length=65536,
        )
        result = fine(optimum.weights; full=true)
        energy_fraction = O.spectral_energy_fraction(
            result.pulse;
            lower_hz=first(O.FREQUENCIES_HZ),
            upper_hz=last(O.FREQUENCIES_HZ),
        )
        passed = result.gain >= 2.0 && result.metrics.broadening_ratio <= 1.25 &&
                 result.metrics.pulse_correlation >= 0.90 &&
                 result.metrics.postcursor_ratio <= 0.12
        push!(rows, (;
            cycles,
            weights=optimum.weights,
            gain=result.gain,
            metrics=result.metrics,
            energy_fraction,
            equal_energy_scale=result.equal_energy_scale,
            evaluations=sum(getproperty.(searches, :evaluations)),
            passed,
        ))
        previous_weights = optimum.weights
        println("[+] cycles=$cycles: G=$(result.gain), " *
                "Bt=$(result.metrics.broadening_ratio), " *
                "rho=$(result.metrics.pulse_correlation), " *
                "post=$(result.metrics.postcursor_ratio), passed=$passed")
    end
    first_pass = findfirst(getproperty.(rows, :passed))
    selected = isnothing(first_pass) ? rows[argmax(getproperty.(rows, :gain))] :
               rows[first_pass]
    summary_path = joinpath(O.LENS_ROOT, "15_pulse_cycle_sweep.csv")
    open(summary_path, "w") do io
        println(io, "pulse_cycles,band_spectral_energy_fraction,optimized_anchor_impulse_gain,broadening_ratio,pulse_correlation,postcursor_ratio,equal_energy_scale,minimum_weight,maximum_weight,objective_evaluations,passed")
        for row in rows
            println(io, join((
                row.cycles, row.energy_fraction, row.gain,
                row.metrics.broadening_ratio, row.metrics.pulse_correlation,
                row.metrics.postcursor_ratio, row.equal_energy_scale,
                minimum(row.weights), maximum(row.weights), row.evaluations,
                row.passed,
            ), ','))
        end
    end
    design_path = joinpath(O.LENS_ROOT, "15_selected_pulse_drive.jld2")
    JLD2.jldsave(
        design_path;
        format_version=1,
        frequencies_hz=O.FREQUENCIES_HZ,
        pulse_cycles=selected.cycles,
        pressure_weights=selected.weights,
        predicted_anchor_impulse_gain=selected.gain,
        predicted_anchor_metrics=selected.metrics,
        band_spectral_energy_fraction=selected.energy_fraction,
        equal_energy_scale=selected.equal_energy_scale,
        passed=selected.passed,
    )
    weights_path = joinpath(O.LENS_ROOT, "15_selected_pulse_drive.csv")
    centers_mm = collect(0.0:8.2:57.4)
    open(weights_path, "w") do io
        println(io, "group,y_mm,pressure_weight")
        for index in eachindex(selected.weights)
            println(io, "$(index),$(centers_mm[index]),$(selected.weights[index])")
        end
    end
    println("[+] selected pulse=$(selected.cycles) cycles, gain=$(selected.gain), " *
            "passed=$(selected.passed)")
    println("[+] selected weights=$(selected.weights)")
    println("[+] $summary_path")
    println("[+] $design_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
