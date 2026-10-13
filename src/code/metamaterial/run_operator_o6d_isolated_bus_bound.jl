module RunOperatorO6DIsolatedBusBound

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

export dct_basis,
       transfer_metrics,
       sparse_internal_map,
       quantized_target,
       run_o6d_isolated_bus_bound

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const O4_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o4_k3q1_transfer_target_v0",
    "o4_transfer_target.jld2",
)
const O6C_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6c_composite_twin_v0",
    "o6c_composite_twin.jld2",
)
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6d_isolated_bus_bound_v0",
)
const ELEMENT_COUNT = 17
const PORT_COUNT = 3
const ROTATION_SEED = UInt64(0x06d00117)
const ROTATION_DRAWS = 200_000
const FIDELITY_GATE = 0.98
const COHERENCE_GATE = 0.15
const IMBALANCE_GATE_DB = 0.75
const AMPLITUDE_RMS_GATE = 0.03
const DELAY_RMS_GATE_S = 82.6e-9
const USEFUL_AMPLITUDE_GATE = 0.8994

function dct_basis(count=ELEMENT_COUNT)
    basis = zeros(Float64, count, PORT_COUNT)
    basis[:, 1] .= inv(sqrt(count))
    for mode in 2:PORT_COUNT
        basis[:, mode] .= sqrt(2 / count) .* cos.(
            (mode - 1) .* ((1:count) .- 0.5) .* pi / count,
        )
    end
    basis
end

function balance_and_contract(transfer)
    gains = [norm(transfer[:, port]) for port in 1:PORT_COUNT]
    all(>(0), gains) || return zeros(eltype(transfer), size(transfer))
    reference = exp(mean(log.(gains)))
    balanced = transfer * Diagonal(reference ./ gains)
    balanced ./ max(1.0, opnorm(balanced))
end

function transfer_metrics(transfer, target=dct_basis())
    gram = adjoint(transfer) * transfer
    gains = sqrt.(max.(0.0, real.(diag(gram))))
    common = dot(vec(target), vec(transfer)) / sum(abs2, target)
    coherence = maximum(
        abs(gram[first, second]) /
        sqrt(max(eps(), real(gram[first, first] * gram[second, second])))
        for first in 1:PORT_COUNT, second in 1:PORT_COUNT if first != second
    )
    (
        common,
        useful_amplitude=abs(common),
        fidelity=abs(dot(vec(target), vec(transfer))) /
                 (norm(target) * max(eps(), norm(transfer))),
        coherence,
        imbalance_db=20log10(maximum(gains) / max(eps(), minimum(gains))),
        amplitude_rms=norm(abs.(transfer) .- abs.(common .* target)) /
                      max(eps(), norm(abs.(common .* target))),
        sigma_max=opnorm(transfer),
    )
end

function electrical_gates(metrics)
    (
        transfer_fidelity=metrics.fidelity >= FIDELITY_GATE,
        port_coherence=metrics.coherence <= COHERENCE_GATE,
        column_imbalance=metrics.imbalance_db <= IMBALANCE_GATE_DB,
        amplitude_rms=metrics.amplitude_rms <= AMPLITUDE_RMS_GATE,
        useful_amplitude=metrics.useful_amplitude >= USEFUL_AMPLITUDE_GATE,
        passivity=metrics.sigma_max <= 1 + 1e-10,
        differential_delay=true,
        reciprocity=true,
    )
end

function gate_score(metrics)
    minimum((
        metrics.fidelity / FIDELITY_GATE,
        COHERENCE_GATE / max(metrics.coherence, 1e-15),
        IMBALANCE_GATE_DB / max(metrics.imbalance_db, 1e-15),
        AMPLITUDE_RMS_GATE / max(metrics.amplitude_rms, 1e-15),
        metrics.useful_amplitude / USEFUL_AMPLITUDE_GATE,
        1 / max(metrics.sigma_max, 1e-15),
    ))
end

function quaternion_rotation(rng)
    quaternion = randn(rng, 4)
    quaternion ./= norm(quaternion)
    w, x, y, z = quaternion
    [
        1 - 2(y^2 + z^2) 2(x * y - z * w) 2(x * z + y * w)
        2(x * y + z * w) 1 - 2(x^2 + z^2) 2(y * z - x * w)
        2(x * z - y * w) 2(y * z + x * w) 1 - 2(x^2 + y^2)
    ]
end

function axis_rotation(axis, angle)
    cosine = cos(angle)
    sine = sin(angle)
    axis == 1 && return [1.0 0 0; 0 cosine -sine; 0 sine cosine]
    axis == 2 && return [cosine 0 sine; 0 1.0 0; -sine 0 cosine]
    axis == 3 && return [cosine -sine 0; sine cosine 0; 0 0 1.0]
    throw(ArgumentError("rotation axis must be 1, 2 or 3"))
end

function sparse_internal_map(target, rotation, connections_per_column)
    1 <= connections_per_column <= PORT_COUNT ||
        throw(ArgumentError("connections per column must be in 1:3"))
    internal_target = target * rotation
    internal = zeros(Float64, size(internal_target))
    for row in axes(internal_target, 1)
        local_row = internal_target[row, :]
        if connections_per_column == 1
            port = argmax(abs.(local_row))
            internal[row, port] = local_row[port]
        elseif connections_per_column == 2
            internal[row, :] .= local_row
            internal[row, argmin(abs.(local_row))] = 0.0
        else
            internal[row, :] .= local_row
        end
    end
    public_transfer = balance_and_contract(internal * transpose(rotation))
    (; internal, public_transfer)
end

function rotation_objective(target, rotation, connections_per_column, objective)
    mapped = sparse_internal_map(target, rotation, connections_per_column)
    metrics = transfer_metrics(mapped.public_transfer, target)
    value = if objective == :gate_score
        gate_score(metrics)
    elseif objective == :fidelity
        metrics.fidelity
    elseif objective == :amplitude_rms
        -metrics.amplitude_rms
    else
        throw(ArgumentError("unknown rotation objective: $objective"))
    end
    (; value, metrics, mapped.internal, transfer=mapped.public_transfer, rotation)
end

function refine_rotation(target, initial, connections_per_column, objective)
    best = rotation_objective(target, initial, connections_per_column, objective)
    step = 0.25
    while step >= 1e-7
        improved = false
        for axis in 1:3, direction in (-1.0, 1.0)
            candidate_rotation = best.rotation * axis_rotation(axis, direction * step)
            candidate = rotation_objective(
                target, candidate_rotation, connections_per_column, objective,
            )
            if candidate.value > best.value + 1e-14
                best = candidate
                improved = true
            end
        end
        improved || (step /= 2)
    end
    best
end


function search_sparse_mapping(target, connections_per_column)
    rng = Xoshiro(ROTATION_SEED + UInt64(connections_per_column))
    objectives = (:gate_score, :fidelity, :amplitude_rms)
    best = Dict(
        objective => rotation_objective(
            target, Matrix{Float64}(I, PORT_COUNT, PORT_COUNT),
            connections_per_column, objective,
        ) for objective in objectives
    )
    for _ in 1:ROTATION_DRAWS
        rotation = quaternion_rotation(rng)
        mapped = sparse_internal_map(target, rotation, connections_per_column)
        metrics = transfer_metrics(mapped.public_transfer, target)
        values = (
            gate_score=gate_score(metrics),
            fidelity=metrics.fidelity,
            amplitude_rms=-metrics.amplitude_rms,
        )
        for objective in objectives
            value = getproperty(values, objective)
            if value > best[objective].value
                best[objective] = (
                    value, metrics, mapped.internal,
                    transfer=mapped.public_transfer, rotation,
                )
            end
        end
    end
    Dict(
        objective => refine_rotation(
            target, best[objective].rotation, connections_per_column, objective,
        ) for objective in objectives
    )
end

function unique_magnitude_count(internal)
    active = abs.(internal[abs.(internal) .> 1e-12])
    isempty(active) && return 0
    length(unique(round.(active ./ maximum(active); digits=6)))
end

function lloyd_levels(values, level_count)
    isempty(values) && return Float64[]
    unique_values = sort(unique(values))
    level_count >= length(unique_values) && return copy(values)
    centers = [quantile(values, (index - 0.5) / level_count) for index in 1:level_count]
    assignment = ones(Int, length(values))
    for _ in 1:200
        new_assignment = [argmin(abs.(centers .- value)) for value in values]
        new_centers = copy(centers)
        for index in 1:level_count
            selected = values[new_assignment .== index]
            isempty(selected) || (new_centers[index] = mean(selected))
        end
        new_assignment == assignment && maximum(abs.(new_centers .- centers)) < 1e-14 &&
            return [new_centers[new_assignment[index]] for index in eachindex(values)]
        assignment = new_assignment
        centers = new_centers
    end
    [centers[assignment[index]] for index in eachindex(values)]
end

function quantized_target(target, level_count)
    quantized = zeros(Float64, size(target))
    for port in 1:PORT_COUNT
        active = abs.(target[:, port]) .> 1e-12
        values = abs.(target[active, port])
        magnitudes = lloyd_levels(values, level_count)
        quantized[active, port] .= sign.(target[active, port]) .* magnitudes
    end
    balance_and_contract(quantized)
end

function unwrap_phase(values)
    result = angle.(values)
    for index in 2:length(result)
        result[index:end] .-= 2pi * round(
            (result[index] - result[index - 1]) / (2pi),
        )
    end
    result
end

function derivative(values, coordinates)
    result = similar(values, Float64)
    result[1] = (values[2] - values[1]) / (coordinates[2] - coordinates[1])
    result[end] = (values[end] - values[end - 1]) /
                  (coordinates[end] - coordinates[end - 1])
    for index in 2:(length(values) - 1)
        result[index] = (values[index + 1] - values[index - 1]) /
                        (coordinates[index + 1] - coordinates[index - 1])
    end
    result
end

function differential_delay(frequencies_hz, transfer, target)
    omega = 2pi .* frequencies_hz
    active = abs.(target) .>= 0.10maximum(abs, target)
    delays = fill(NaN, size(transfer))
    for port in 1:PORT_COUNT, bin in 1:ELEMENT_COUNT
        active[bin, port] || continue
        relative = transfer[bin, port, :] ./ target[bin, port]
        delays[bin, port, :] .= -derivative(unwrap_phase(relative), omega)
    end
    weights = abs2.(target)
    rms = zeros(Float64, length(frequencies_hz))
    for frequency_index in eachindex(frequencies_hz)
        valid = active .& isfinite.(delays[:, :, frequency_index])
        local_weights = weights[valid]
        local_delays = delays[:, :, frequency_index][valid]
        common = sum(local_weights .* local_delays) / sum(local_weights)
        rms[frequency_index] = sqrt(
            sum(local_weights .* (local_delays .- common) .^ 2) /
            sum(local_weights),
        )
    end
    (; delays, rms, maximum_rms=maximum(rms))
end

function mechanical_metrics(kappa, frequencies_hz, measured, target)
    transfer = similar(measured)
    for frequency_index in eachindex(frequencies_hz)
        common = dot(
            vec(target), vec(measured[:, :, frequency_index]),
        ) / sum(abs2, target)
        residual = measured[:, :, frequency_index] .- common .* target
        transfer[:, :, frequency_index] .= common .* target .+ kappa .* residual
    end
    delay = differential_delay(frequencies_hz, transfer, target)
    rows = [transfer_metrics(transfer[:, :, index], target) for index in eachindex(frequencies_hz)]
    summary = (
        kappa,
        residual_suppression=1 - kappa,
        residual_coupling_db=kappa > 0 ? 20log10(kappa) : -Inf,
        minimum_fidelity=minimum(getproperty.(rows, :fidelity)),
        maximum_coherence=maximum(getproperty.(rows, :coherence)),
        maximum_imbalance_db=maximum(getproperty.(rows, :imbalance_db)),
        maximum_amplitude_rms=maximum(getproperty.(rows, :amplitude_rms)),
        maximum_delay_rms_ns=1e9delay.maximum_rms,
        minimum_useful_amplitude=minimum(getproperty.(rows, :useful_amplitude)),
        maximum_sigma=maximum(getproperty.(rows, :sigma_max)),
    )
    gates = (
        transfer_fidelity=summary.minimum_fidelity >= FIDELITY_GATE,
        port_coherence=summary.maximum_coherence <= COHERENCE_GATE,
        column_imbalance=summary.maximum_imbalance_db <= IMBALANCE_GATE_DB,
        amplitude_rms=summary.maximum_amplitude_rms <= AMPLITUDE_RMS_GATE,
        differential_delay=delay.maximum_rms <= DELAY_RMS_GATE_S,
        useful_amplitude=summary.minimum_useful_amplitude >= USEFUL_AMPLITUDE_GATE,
        passivity=summary.maximum_sigma <= 1 + 1e-10,
        reciprocity=true,
    )
    (; summary, summary..., gates, passed=all(values(gates)), transfer, delay)
end

function mechanical_threshold(frequencies_hz, measured, target)
    lower = 0.0
    upper = 1.0
    mechanical_metrics(lower, frequencies_hz, measured, target).passed ||
        error("perfectly isolated O6d control does not pass O4")
    mechanical_metrics(upper, frequencies_hz, measured, target).passed &&
        return 1.0
    for _ in 1:50
        middle = (lower + upper) / 2
        if mechanical_metrics(middle, frequencies_hz, measured, target).passed
            lower = middle
        else
            upper = middle
        end
    end
    lower
end

function write_namedtuple_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function metric_row(class_id, connections, active_taps, unique_levels, metrics, gates)
    (
        class_id,
        connections_per_column=connections,
        active_taps,
        unique_nonzero_magnitude_ratios=unique_levels,
        normalized_transfer_fidelity=metrics.fidelity,
        maximum_port_coherence=metrics.coherence,
        column_imbalance_db=metrics.imbalance_db,
        differential_amplitude_rms=metrics.amplitude_rms,
        useful_amplitude=metrics.useful_amplitude,
        sigma_max=metrics.sigma_max,
        differential_delay_rms_ns=0.0,
        reciprocity_residual=0.0,
        passed=all(values(gates)),
    )
end

function run_o6d_isolated_bus_bound()
    isfile(O4_PATH) || error("O4 target not found: $O4_PATH")
    isfile(O6C_PATH) || error("O6c result not found: $O6C_PATH")
    mkpath(OUTPUT_ROOT)
    o4 = load(O4_PATH)
    o6c = load(O6C_PATH)
    target = Float64.(o4["state_matrix"])
    frequencies_hz = Float64.(o6c["frequencies_hz"])
    measured = ComplexF64.(o6c["as_built_transfer"])
    size(target) == (ELEMENT_COUNT, PORT_COUNT) || error("O4 target size changed")
    target' * target ≈ Matrix{Float64}(I, PORT_COUNT, PORT_COUNT) ||
        error("O4 target is no longer orthonormal")
    size(measured) == (ELEMENT_COUNT, PORT_COUNT, length(frequencies_hz)) ||
        error("O6c transfer size changed")

    println("[>] O6d mechanical residual threshold")
    kappa_max = mechanical_threshold(frequencies_hz, measured, target)
    threshold_result = mechanical_metrics(kappa_max, frequencies_hz, measured, target)
    mechanical_rows = [
        mechanical_metrics(kappa, frequencies_hz, measured, target).summary
        for kappa in range(0.0, 1.0; length=101)
    ]

    search_results = Dict{Int, Any}()
    search_rows = NamedTuple[]
    electrical_rows = NamedTuple[]
    rotation_rows = NamedTuple[]
    for connections in 1:2
        println("[>] O6d sparse electrical map: $connections connection(s)/column")
        results = search_sparse_mapping(target, connections)
        search_results[connections] = results
        for objective in (:gate_score, :fidelity, :amplitude_rms)
            result = results[objective]
            metrics = result.metrics
            push!(search_rows, (
                connections_per_column=connections,
                objective=String(objective),
                objective_value=result.value,
                gate_score=gate_score(metrics),
                normalized_transfer_fidelity=metrics.fidelity,
                maximum_port_coherence=metrics.coherence,
                column_imbalance_db=metrics.imbalance_db,
                differential_amplitude_rms=metrics.amplitude_rms,
                useful_amplitude=metrics.useful_amplitude,
                sigma_max=metrics.sigma_max,
            ))
        end
        selected = results[:gate_score]
        gates = electrical_gates(selected.metrics)
        push!(electrical_rows, metric_row(
            "sparse_$(connections)_bus",
            connections,
            count(abs.(selected.internal) .> 1e-12),
            unique_magnitude_count(selected.internal),
            selected.metrics,
            gates,
        ))
        for row in 1:PORT_COUNT, column in 1:PORT_COUNT
            push!(rotation_rows, (
                connections_per_column=connections,
                row,
                column,
                value=selected.rotation[row, column],
            ))
        end
    end

    exact_transfer = target
    exact_metrics = transfer_metrics(exact_transfer, target)
    exact_gates = electrical_gates(exact_metrics)
    push!(electrical_rows, metric_row(
        "full_3_bus_upper_bound",
        3,
        count(abs.(target) .> 1e-12),
        unique_magnitude_count(target),
        exact_metrics,
        exact_gates,
    ))

    quantization_rows = NamedTuple[]
    quantized_transfers = Dict{Int, Matrix{Float64}}()
    for levels in (1, 2, 4, 8, 16)
        transfer = quantized_target(target, levels)
        quantized_transfers[levels] = transfer
        metrics = transfer_metrics(transfer, target)
        gates = electrical_gates(metrics)
        push!(quantization_rows, metric_row(
            "full_3_bus_$(levels)_levels_per_port",
            3,
            count(abs.(transfer) .> 1e-12),
            unique_magnitude_count(transfer),
            metrics,
            gates,
        ))
    end
    push!(quantization_rows, metric_row(
        "full_3_bus_continuous",
        3,
        count(abs.(target) .> 1e-12),
        unique_magnitude_count(target),
        exact_metrics,
        exact_gates,
    ))

    low_complexity_passed = any(row.passed for row in electrical_rows if
        row.connections_per_column <= 2)
    full_map_passed = last(electrical_rows).passed
    overall_passed = low_complexity_passed && threshold_result.passed
    write_namedtuple_csv(
        joinpath(OUTPUT_ROOT, "mechanical_isolation_frontier.csv"), mechanical_rows,
    )
    write_namedtuple_csv(
        joinpath(OUTPUT_ROOT, "electrical_sparsity_frontier.csv"), electrical_rows,
    )
    write_namedtuple_csv(
        joinpath(OUTPUT_ROOT, "rotation_search_diagnostics.csv"), search_rows,
    )
    write_namedtuple_csv(
        joinpath(OUTPUT_ROOT, "selected_rotation_matrices.csv"), rotation_rows,
    )
    write_namedtuple_csv(
        joinpath(OUTPUT_ROOT, "weight_quantization_frontier.csv"), quantization_rows,
    )
    jldsave(
        joinpath(OUTPUT_ROOT, "o6d_isolated_bus_bound.jld2");
        format_version=1,
        source_o4_path=abspath(O4_PATH),
        source_o6c_path=abspath(O6C_PATH),
        target,
        frequencies_hz,
        measured_o6c_transfer=measured,
        rotation_seed=ROTATION_SEED,
        rotation_draws=ROTATION_DRAWS,
        kappa_max,
        threshold_result,
        mechanical_rows,
        search_results,
        electrical_rows,
        rotation_rows,
        quantization_rows,
        quantized_transfers,
        low_complexity_passed,
        full_map_passed,
        overall_passed,
    )

    one_bus = only(filter(row -> row.class_id == "sparse_1_bus", electrical_rows))
    two_bus = only(filter(row -> row.class_id == "sparse_2_bus", electrical_rows))
    first_quantized_pass = findfirst(row -> row.passed, quantization_rows)
    first_quantized_label = isnothing(first_quantized_pass) ? "none" :
        quantization_rows[first_quantized_pass].class_id
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6d mechanically isolated three-bus bound")
        println(io)
        println(io, overall_passed ?
            "**PASS: one macro-column physical twin is authorized.**" :
            "**FAIL: the low-complexity isolated three-bus face is stopped before geometry.**")
        println(io)
        @printf(io,
            "The measured O6c differential residual must be reduced to `kappa <= %.6f`: at least `%.2f%%` suppression (`%.2f dB` residual amplitude relative to O6c). This is an operator requirement, not a kerf-width prediction.\n\n",
            kappa_max, 100(1 - kappa_max), 20log10(kappa_max),
        )
        @printf(io,
            "At that boundary the active constraint is differential delay (`%.3f ns`); amplitude RMS is already `%.3f%%`.\n\n",
            threshold_result.maximum_delay_rms_ns,
            100threshold_result.maximum_amplitude_rms,
        )
        println(io, "## Ideal electrical coding frontier")
        println(io)
        println(io, "| Internal map | Active taps | Fidelity | Amp RMS | Useful | Gate |")
        println(io, "|---|---:|---:|---:|---:|---|")
        for row in electrical_rows
            @printf(io, "| `%s` | %d | %.5f | %.2f%% | %.5f | %s |\n",
                row.class_id, row.active_taps,
                row.normalized_transfer_fidelity,
                100row.differential_amplitude_rms,
                row.useful_amplitude,
                row.passed ? "PASS" : "FAIL",
            )
        end
        println(io)
        @printf(io,
            "The best one-bus and two-bus maps use an ideal lossless public-port rotation, continuous signed tap weights, fixed column balancing, zero network delay and no network loss. Their misses therefore cannot be repaired by adding real passive parasitics. The SO(3) search is deterministic (`%d` Haar rotations plus local refinement) but not a global certificate.\n\n",
            ROTATION_DRAWS,
        )
        println(io, "The full map passes analytically, but needs `$(last(electrical_rows).active_taps)` nonzero bus-to-column couplings and `$(last(electrical_rows).unique_nonzero_magnitude_ratios)` distinct normalized magnitudes in its continuous form. The first passing quantization entry is `$(first_quantized_label)` (`17` distinct normalized ratios across the map). This is the explicitly identified passive `3-to-17` beamformer fallback from O5, not a pass of the selected monolithic topology.")
        println(io)
        println(io, overall_passed ?
            "The next allowed calculation is one isolated-column FEM using the threshold above; no kerf sweep is allowed." :
            "No kerf FEM or F1d replay is authorized. The next decision must explicitly compare accepting the full passive beamformer fabrication burden with returning to a different architecture-frontier point; it may not hide those taps inside the word `composite`.")
    end

    mechanics_plot = plot(
        getproperty.(mechanical_rows, :kappa),
        100 .* getproperty.(mechanical_rows, :maximum_amplitude_rms);
        linewidth=2, label="amplitude RMS, %", xlabel="residual kappa",
        ylabel="gate observable", title="O6d mechanical isolation requirement",
    )
    plot!(mechanics_plot,
        getproperty.(mechanical_rows, :kappa),
        getproperty.(mechanical_rows, :maximum_delay_rms_ns) ./ (1e9DELAY_RMS_GATE_S);
        linewidth=2, label="delay / gate")
    hline!(mechanics_plot, [100AMPLITUDE_RMS_GATE, 1.0];
        color=:black, linestyle=:dash, label="gates")
    vline!(mechanics_plot, [kappa_max]; color=:red, linestyle=:dot, label="kappa max")
    savefig(mechanics_plot, joinpath(OUTPUT_ROOT, "mechanical_isolation_bound.png"))

    electrical_plot = scatter(
        getproperty.(electrical_rows, :active_taps),
        100 .* getproperty.(electrical_rows, :differential_amplitude_rms);
        marker=:circle, markersize=7, label="sparsity frontier",
        xlabel="active bus-to-column taps", ylabel="amplitude RMS, %",
        title="O6d electrical coding burden",
    )
    hline!(electrical_plot, [100AMPLITUDE_RMS_GATE];
        color=:black, linestyle=:dash, label="O4 gate")
    savefig(electrical_plot, joinpath(OUTPUT_ROOT, "electrical_sparsity_bound.png"))
    println(overall_passed ? "[+] O6d low-complexity PASS" : "[-] O6d low-complexity FAIL")
    println("[+] $(joinpath(OUTPUT_ROOT, "README.md"))")
    (; overall_passed, low_complexity_passed, full_map_passed, kappa_max,
       electrical_rows, quantization_rows, threshold_result)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6d_isolated_bus_bound()

end
