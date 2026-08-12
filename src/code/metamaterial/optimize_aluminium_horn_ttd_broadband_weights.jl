module OptimizeAluminiumHornTTDBroadbandWeights

using FFTW
using JLD2
using LinearAlgebra

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const LENS_ROOT = joinpath(
    PROJECT_ROOT, "tmp", "aluminium_horn_ttd_full_diffuser_aperture_242khz",
)
const REFERENCE_ROOT = joinpath(
    PROJECT_ROOT, "tmp", "aluminium_horn_ttd_full_aperture_242khz",
)
const FREQUENCIES_HZ = Float64[193.8e3, 242.0e3, 290.1e3]
const PRESSURE_PA = 1.0e6

include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
include(joinpath(@__DIR__, "spectral_analysis.jl"))
using .ImpulseRiskAnalysis
using .SpectralAnalysis

frequency_suffix(frequency_hz) = isapprox(frequency_hz, 242.0e3; atol=1e-6) ? "" :
    "_f$(round(Int, frequency_hz))hz"
matrix_path(frequency_hz) = joinpath(
    LENS_ROOT, "lens_half_response_matrix$(frequency_suffix(frequency_hz)).jld2",
)
reference_path(frequency_hz) = joinpath(
    REFERENCE_ROOT, "uniform_half_harmonic$(frequency_suffix(frequency_hz)).jld2",
)

function active_power_matrix(displacement_work, frequency_hz)
    raw = real.(0.5im * PRESSURE_PA * 2pi * frequency_hz .* conj.(displacement_work))
    Symmetric((raw + transpose(raw)) / 2)
end

active_power(matrix, weights) = real(dot(weights, matrix * weights))

function subset_solvers(power_matrix)
    count = size(power_matrix, 1)
    result = NamedTuple[]
    for mask in 1:(2^count - 1)
        indices = findall(index -> !iszero(mask & (1 << (index - 1))), 1:count)
        matrix = Matrix(power_matrix[indices, indices])
        minimum(eigvals(Symmetric(matrix))) > 0 || continue
        push!(result, (; indices, inverse=inv(matrix)))
    end
    result
end

function optimize_single_frequency(focus_response, power_matrix, power_budget;
                                   phase_samples=3600)
    solvers = subset_solvers(power_matrix)
    best = nothing
    for phase in range(0.0, 2pi; length=phase_samples + 1)[1:end-1]
        projection = real.(cis(-phase) .* focus_response)
        for solver in solvers
            local_weights = solver.inverse * projection[solver.indices]
            all(>=(0.0), local_weights) || continue
            weights = zeros(length(focus_response))
            weights[solver.indices] .= local_weights
            power = active_power(power_matrix, weights)
            power > 0 || continue
            weights .*= sqrt(power_budget / power)
            amplitude = abs(sum(focus_response .* weights))
            if isnothing(best) || amplitude > best.amplitude
                best = (; weights, amplitude)
            end
        end
    end
    isnothing(best) && error("single-frequency nonnegative optimum was not found")
    best
end

function unwrap_phase(phases)
    result = Float64[first(phases)]
    for phase in Iterators.drop(phases, 1)
        push!(result, last(result) + mod(phase - last(result) + pi, 2pi) - pi)
    end
    result
end

function linear_interpolate(x, y, query)
    right = searchsortedfirst(x, query)
    right == firstindex(x) && return y[right]
    right > lastindex(x) && return y[end]
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * y[left] + fraction * y[right]
end

function analytic_envelope(signal)
    count = length(signal)
    spectrum = fft(signal)
    multiplier = zeros(Float64, count)
    multiplier[1] = 1.0
    if iseven(count)
        multiplier[2:(count ÷ 2)] .= 2.0
        multiplier[count ÷ 2 + 1] = 1.0
    else
        multiplier[2:((count + 1) ÷ 2)] .= 2.0
    end
    abs.(ifft(spectrum .* multiplier))
end

function build_objective(
    focus_response,
    power_matrix,
    reference_focus,
    reference_power;
    cycles=5.0,
    samples_per_period=60,
    fft_length=16384,
)
    pulse = pulse_spectrum(PulseConfig(
        center_frequency_hz=242.0e3,
        cycles=Float64(cycles),
        samples_per_period=Int(samples_per_period),
        fft_length=Int(fft_length),
    ))
    active = (pulse.frequency_hz .>= first(FREQUENCIES_HZ)) .&
             (pulse.frequency_hz .<= last(FREQUENCIES_HZ))
    active_indices = findall(active)
    reference_transfer = zeros(ComplexF64, length(pulse.frequency_hz))
    reference_transfer[active] .= 1.0 + 0im
    reference_waveform = irfft(
        pulse.spectrum .* reference_transfer, length(pulse.signal),
    )
    reference_peak = maximum(analytic_envelope(reference_waveform))
    spectral_weight = abs2.(pulse.spectrum)

    function evaluate(weights; full=false)
        all(>=(0.0), weights) || return full ? nothing : -Inf
        any(>(0.0), weights) || return full ? nothing : -Inf
        lens_focus = ComplexF64[
            sum(focus_response[index] .* weights)
            for index in eachindex(FREQUENCIES_HZ)
        ]
        lens_power = Float64[
            active_power(power_matrix[index], weights)
            for index in eachindex(FREQUENCIES_HZ)
        ]
        all(>(0.0), lens_power) || return full ? nothing : -Inf
        relative = lens_focus ./ reference_focus
        phase = unwrap_phase(angle.(relative))
        log_amplitude = log.(abs.(relative))
        log_lens_power = log.(lens_power)
        log_reference_power = log.(reference_power)
        transfer = zeros(ComplexF64, length(pulse.frequency_hz))
        interpolated_lens_power = zeros(Float64, length(pulse.frequency_hz))
        interpolated_reference_power = zeros(Float64, length(pulse.frequency_hz))
        for index in active_indices
            frequency_hz = pulse.frequency_hz[index]
            transfer[index] = exp(linear_interpolate(
                FREQUENCIES_HZ, log_amplitude, frequency_hz,
            )) * cis(linear_interpolate(FREQUENCIES_HZ, phase, frequency_hz))
            interpolated_lens_power[index] = exp(linear_interpolate(
                FREQUENCIES_HZ, log_lens_power, frequency_hz,
            ))
            interpolated_reference_power[index] = exp(linear_interpolate(
                FREQUENCIES_HZ, log_reference_power, frequency_hz,
            ))
        end
        lens_energy = sum(spectral_weight[active] .* interpolated_lens_power[active])
        reference_energy = sum(
            spectral_weight[active] .* interpolated_reference_power[active],
        )
        equal_energy_scale = sqrt(reference_energy / lens_energy)
        transfer .*= equal_energy_scale
        waveform = irfft(pulse.spectrum .* transfer, length(pulse.signal))
        gain = maximum(analytic_envelope(waveform)) / reference_peak
        full || return gain
        metrics = pulse_metrics(
            waveform, reference_waveform, reference_waveform,
            pulse.time_s[2] - pulse.time_s[1],
        )
        (; gain, metrics, waveform, reference_waveform, transfer, pulse,
           lens_focus, lens_power, relative_phase_rad=phase, equal_energy_scale)
    end
    evaluate
end

function normalized_ratios(weights)
    weights = max.(Float64.(weights), 1e-8)
    weights ./ weights[1]
end

function pattern_search(evaluate, initial; initial_step=0.35, minimum_step=0.004)
    ratios = normalized_ratios(initial)
    variables = log.(ratios[2:end])
    best_weights = vcat(1.0, exp.(variables))
    best_gain = evaluate(best_weights)
    step = initial_step
    evaluations = 1
    while step >= minimum_step
        improved = false
        for index in eachindex(variables), direction in (-1.0, 1.0)
            candidate_variables = copy(variables)
            candidate_variables[index] += direction * step
            candidate_weights = vcat(1.0, exp.(candidate_variables))
            gain = evaluate(candidate_weights)
            evaluations += 1
            if gain > best_gain + 1e-8
                variables = candidate_variables
                best_weights = candidate_weights
                best_gain = gain
                improved = true
            end
        end
        improved || (step /= 2)
    end
    (; weights=best_weights, gain=best_gain, evaluations)
end

function main()
    matrices = JLD2.load.(matrix_path.(FREQUENCIES_HZ))
    references = JLD2.load.(reference_path.(FREQUENCIES_HZ))
    focus_response = [vec(item["focus_response_m"][1, :]) for item in matrices]
    power_matrix = [active_power_matrix(
        item["source_normal_displacement_integral_m3"], FREQUENCIES_HZ[index],
    ) for (index, item) in enumerate(matrices)]
    reference_focus = ComplexF64[item["focus_displacement_m"][1] for item in references]
    reference_power = abs.(Float64[item["active_input_power_w"] for item in references])
    carrier_design = JLD2.load(joinpath(LENS_ROOT, "15_power_optimal_weights.jld2"))
    carrier_weights = Float64.(carrier_design["pressure_weights"])
    single = [optimize_single_frequency(
        focus_response[index], power_matrix[index], reference_power[index],
    ) for index in eachindex(FREQUENCIES_HZ)]
    single_gain = [single[index].amplitude / abs(reference_focus[index])
                   for index in eachindex(FREQUENCIES_HZ)]
    evaluate = build_objective(
        focus_response, power_matrix, reference_focus, reference_power,
    )
    starts = vcat([carrier_weights, ones(8)], [item.weights for item in single])
    searches = [pattern_search(evaluate, start) for start in starts]
    optimum = searches[argmax(getproperty.(searches, :gain))]
    result = evaluate(optimum.weights; full=true)
    isnothing(result) && error("broadband optimum evaluation failed")
    passed = result.gain >= 2.0 && result.metrics.broadening_ratio <= 1.25 &&
             result.metrics.pulse_correlation >= 0.90 &&
             result.metrics.postcursor_ratio <= 0.12

    design_path = joinpath(LENS_ROOT, "15_broadband_optimal_weights.jld2")
    JLD2.jldsave(
        design_path;
        format_version=1,
        frequencies_hz=FREQUENCIES_HZ,
        pressure_weights=optimum.weights,
        predicted_anchor_impulse_gain=result.gain,
        predicted_anchor_metrics=result.metrics,
        single_frequency_upper_bound_gain=single_gain,
        carrier_weight_anchor_gain=evaluate(carrier_weights),
        passed,
    )
    weights_path = joinpath(LENS_ROOT, "15_broadband_optimal_weights.csv")
    centers_mm = collect(0.0:8.2:57.4)
    open(weights_path, "w") do io
        println(io, "group,y_mm,pressure_weight")
        for index in eachindex(optimum.weights)
            println(io, "$(index),$(centers_mm[index]),$(optimum.weights[index])")
        end
    end
    summary_path = joinpath(LENS_ROOT, "15_broadband_weight_optimization.csv")
    open(summary_path, "w") do io
        println(io, "carrier_weight_anchor_impulse_gain,optimized_anchor_impulse_gain,lower_single_frequency_upper_bound_gain,carrier_single_frequency_upper_bound_gain,upper_single_frequency_upper_bound_gain,broadening_ratio,pulse_correlation,postcursor_ratio,equal_energy_scale,minimum_weight,maximum_weight,total_objective_evaluations,passed")
        println(io, join((
            evaluate(carrier_weights), result.gain, single_gain...,
            result.metrics.broadening_ratio, result.metrics.pulse_correlation,
            result.metrics.postcursor_ratio, result.equal_energy_scale,
            minimum(optimum.weights), maximum(optimum.weights),
            sum(getproperty.(searches, :evaluations)), passed,
        ), ','))
    end
    println("[+] single-frequency equal-power upper bounds=$single_gain")
    println("[+] carrier-weight anchor impulse gain=$(evaluate(carrier_weights))")
    println("[+] optimized anchor impulse gain=$(result.gain), passed=$passed")
    println("[+] weights=$(optimum.weights)")
    println("[+] Bt=$(result.metrics.broadening_ratio), " *
            "rho=$(result.metrics.pulse_correlation), " *
            "post=$(result.metrics.postcursor_ratio)")
    println("[+] $summary_path")
    println("[+] $design_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
