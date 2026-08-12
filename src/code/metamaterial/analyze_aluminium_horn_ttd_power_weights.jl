module AnalyseAluminiumHornTTDPowerWeights

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_AL_HORN_TTD_FULL_OUTPUT",
    argument("root", joinpath(
        PROJECT_ROOT, "tmp", "aluminium_horn_ttd_full_aperture_242khz",
    )),
)
const BASELINE_ROOT = argument("baseline-root", OUTPUT_ROOT)
const FREQUENCY_HZ = 242.0e3
const PRESSURE_PA = 1.0e6

using JLD2
using LinearAlgebra

function active_power_matrix(displacement_work)
    raw = real.(0.5im * PRESSURE_PA * 2pi * FREQUENCY_HZ .* conj.(displacement_work))
    Symmetric((raw + transpose(raw)) / 2)
end

function active_power(matrix, weights)
    real(dot(weights, matrix * weights))
end

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

function optimize_nonnegative_weights(focus_response, power_matrix, power_budget;
                                      phase_samples=7200)
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
                best = (; weights, amplitude, phase)
            end
        end
    end
    isnothing(best) && error("no nonnegative equal-power aperture weights were found")
    best
end

function main()
    matrix = JLD2.load(joinpath(OUTPUT_ROOT, "lens_half_response_matrix.jld2"))
    matched_uniform = JLD2.load(joinpath(OUTPUT_ROOT, "uniform_half_harmonic.jld2"))
    uniform = JLD2.load(joinpath(BASELINE_ROOT, "uniform_half_harmonic.jld2"))
    lens = JLD2.load(joinpath(OUTPUT_ROOT, "lens_half_harmonic.jld2"))
    focus_response = vec(matrix["focus_response_m"][1, :])
    work = matrix["source_normal_displacement_integral_m3"]
    power_matrix = active_power_matrix(work)
    uniform_power_w = abs(Float64(uniform["active_input_power_w"]))
    uniform_focus_m = abs(uniform["focus_displacement_m"][1])
    matched_uniform_power_w = abs(Float64(matched_uniform["active_input_power_w"]))
    matched_uniform_focus_m = abs(matched_uniform["focus_displacement_m"][1])
    reconstructed_power_w = active_power(power_matrix, ones(length(focus_response)))
    stored_lens_power_w = abs(Float64(lens["active_input_power_w"]))
    power_reconstruction_error = reconstructed_power_w / stored_lens_power_w - 1
    minimum_power_eigenvalue = minimum(eigvals(power_matrix))
    optimum = optimize_nonnegative_weights(focus_response, power_matrix, uniform_power_w)
    optimized_gain = optimum.amplitude / uniform_focus_m
    current_equal_power_gain = abs(sum(focus_response)) / uniform_focus_m *
                               sqrt(uniform_power_w / reconstructed_power_w)
    current_matched_equal_power_gain = abs(sum(focus_response)) /
        matched_uniform_focus_m * sqrt(matched_uniform_power_w / reconstructed_power_w)
    passed = optimized_gain >= 2.0

    summary_path = joinpath(OUTPUT_ROOT, "15_power_weight_upper_bound.csv")
    open(summary_path, "w") do io
        println(io, "baseline_root,current_equal_power_gain,current_matched_equal_power_gain,optimized_nonnegative_equal_power_gain,baseline_uniform_power_w,matched_uniform_power_w,reconstructed_lens_power_w,stored_lens_power_w,power_reconstruction_relative_error,minimum_power_matrix_eigenvalue,minimum_weight,maximum_weight,passed")
        println(io, join((BASELINE_ROOT, current_equal_power_gain,
                          current_matched_equal_power_gain, optimized_gain, uniform_power_w,
                          matched_uniform_power_w,
                          reconstructed_power_w, stored_lens_power_w,
                          power_reconstruction_error, minimum_power_eigenvalue,
                          minimum(optimum.weights), maximum(optimum.weights), passed), ','))
    end
    weights_path = joinpath(OUTPUT_ROOT, "15_power_optimal_weights.csv")
    centers_mm = collect(0.0:8.2:57.4)
    open(weights_path, "w") do io
        println(io, "group,y_mm,pressure_weight")
        for index in eachindex(optimum.weights)
            println(io, "$(index),$(centers_mm[index]),$(optimum.weights[index])")
        end
    end
    weights_field_path = joinpath(OUTPUT_ROOT, "15_power_optimal_weights.jld2")
    JLD2.jldsave(
        weights_field_path;
        format_version=1,
        baseline_root=BASELINE_ROOT,
        frequency_hz=FREQUENCY_HZ,
        pressure_weights=optimum.weights,
        power_budget_w=uniform_power_w,
        predicted_focus_amplitude_m=optimum.amplitude,
        predicted_equal_power_gain=optimized_gain,
    )
    verdict_path = joinpath(OUTPUT_ROOT, "15_power_weight_verdict.txt")
    open(verdict_path, "w") do io
        println(io, passed ?
            "PASS: amplitude/impedance weighting has a >2 equal-power upper bound." :
            "STOP: even optimal nonnegative equal-power weighting cannot reach 2.")
        println(io, "current_equal_power_gain=$current_equal_power_gain")
        println(io, "current_matched_equal_power_gain=$current_matched_equal_power_gain")
        println(io, "optimized_equal_power_gain=$optimized_gain")
        println(io, "power_reconstruction_error=$power_reconstruction_error")
        println(io, "weights=$(join(optimum.weights, ';'))")
    end
    println("[+] current equal-power gain=$current_equal_power_gain")
    println("[+] current matched equal-power gain=$current_matched_equal_power_gain")
    println("[+] optimal nonnegative equal-power gain=$optimized_gain, passed=$passed")
    println("[+] weights=$(optimum.weights)")
    println("[+] power reconstruction error=$power_reconstruction_error")
    println("[+] $summary_path")
    println("[+] $weights_field_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
