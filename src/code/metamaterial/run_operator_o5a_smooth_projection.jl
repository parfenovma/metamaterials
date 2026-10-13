module RunOperatorO5ASmoothProjection

using JLD2
using LinearAlgebra
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o3_f1c_stochastic.jl"))

using .RunOperatorO3F1CStochastic

export run_o5a_validation

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o5a_smooth_projection_v0",
)
const STATIC_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "static_frontier_v0",
    "static_frontier.jld2",
)
const F1C_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "fullwave_f1c_stochastic_v0",
    "f1c_stochastic.jld2",
)
const SOURCE_ID = :trained_k4_q2_unconstrained
const MODE_COUNTS = collect(4:12)
const REFERENCE_ALPHA = 3.0

function dct_basis(count)
    basis = zeros(Float64, count, count)
    basis[:, 1] .= inv(sqrt(count))
    for mode in 2:count
        basis[:, mode] .= sqrt(2 / count) .* cos.(
            (mode - 1) .* ((1:count) .- 0.5) .* pi / count,
        )
    end
    basis
end

function smooth_states(source_states, mode_count)
    basis = dct_basis(size(first(source_states), 1))[:, 1:mode_count]
    [
        begin
            projected = basis * (transpose(basis) * state)
            Matrix(qr(projected).Q)[:, 1:size(state, 2)]
        end
        for state in source_states
    ]
end

function definition_for(states, mode_count)
    (
        id=Symbol("o5a_smooth_m$(mode_count)"),
        kind=:static,
        electronic_channels=4,
        physical_states=2,
        transmit_events=8,
        indices=Int[],
        states,
    )
end

function state_fidelity(source_states, candidate_states)
    values = Float64[]
    minimum_principal_cosine = 1.0
    for (source, candidate) in zip(source_states, candidate_states)
        push!(values, abs(dot(vec(source), vec(candidate))) / (norm(source) * norm(candidate)))
        minimum_principal_cosine = min(
            minimum_principal_cosine,
            minimum(svdvals(transpose(source) * candidate)),
        )
    end
    mean(values), minimum_principal_cosine
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

function run_o5a_validation()
    mkpath(OUTPUT_ROOT)
    source_states = JLD2.load(STATIC_PATH, "architecture_states")[SOURCE_ID]
    f1c = JLD2.load(F1C_PATH)
    electronics_sigma = Float64(f1c["electronics_sigma"])
    validation_h0_specs = RunOperatorO3F1CStochastic.nuisance_specs(
        RunOperatorO3F1CStochastic.SEEDS.validation_h0,
        [:validation_1],
        RunOperatorO3F1CStochastic.VALIDATION_H0_REPLICATES;
        with_targets=false,
    )
    validation_signal_specs = RunOperatorO3F1CStochastic.nuisance_specs(
        RunOperatorO3F1CStochastic.SEEDS.validation_signal,
        [:validation_1],
        RunOperatorO3F1CStochastic.VALIDATION_SIGNAL_REPLICATES;
        with_targets=true,
    )
    rows = NamedTuple[]
    state_map = Dict{Int,Vector{Matrix{Float64}}}()
    for mode_count in MODE_COUNTS
        println("[+] O5a validation M=$mode_count")
        states = smooth_states(source_states, mode_count)
        state_map[mode_count] = states
        definition = definition_for(states, mode_count)
        model = RunOperatorO3F1CStochastic.build_architecture_model(definition)
        base = RunOperatorO3F1CStochastic.load_architecture_scenes(
            definition, [:validation_1],
        )
        factor = RunOperatorO3F1CStochastic.noise_factor(model)
        h0 = RunOperatorO3F1CStochastic.score_components(
            model, base, validation_h0_specs, electronics_sigma,
        )
        signal = RunOperatorO3F1CStochastic.score_components(
            model, base, validation_signal_specs, electronics_sigma,
        )
        threshold = quantile(
            RunOperatorO3F1CStochastic.h0_statistics(
                Xoshiro(
                    RunOperatorO3F1CStochastic.SEEDS.validation_noise +
                    RunOperatorO3F1CStochastic.architecture_offset(definition.id),
                ),
                h0, factor, RunOperatorO3F1CStochastic.H0_VALIDATION_TRIALS,
            ),
            0.99,
        )
        metrics = [
            RunOperatorO3F1CStochastic.signal_metrics(
                Xoshiro(
                    RunOperatorO3F1CStochastic.SEEDS.validation_noise +
                    RunOperatorO3F1CStochastic.architecture_offset(definition.id) + UInt64(index),
                ),
                signal, factor, threshold, amplitude,
            )
            for (index, amplitude) in enumerate(RunOperatorO3F1CStochastic.AMPLITUDES)
        ]
        reference_index = only(findall(==(REFERENCE_ALPHA), RunOperatorO3F1CStochastic.AMPLITUDES))
        reference = metrics[reference_index]
        fidelity, minimum_principal_cosine = state_fidelity(source_states, states)
        alpha90 = RunOperatorO3F1CStochastic.crossing_amplitude(
            RunOperatorO3F1CStochastic.AMPLITUDES, getproperty.(metrics, :pd),
        )
        validation_passed = reference.pd >= 0.90 && reference.pdl >= 0.80 &&
            reference.minimum_case_pd >= 0.75 && alpha90 <= REFERENCE_ALPHA
        push!(rows, (
            mode_count,
            state_fidelity=fidelity,
            minimum_principal_cosine,
            nuisance_rank=model.nuisance_rank,
            template_stable_rank=model.stable_rank,
            maximum_template_coherence=model.maximum_template_coherence,
            threshold,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            minimum_case_pd_reference=reference.minimum_case_pd,
            minimum_case_pdl_reference=reference.minimum_case_pdl,
            alpha90_validation=alpha90,
            validation_passed,
        ))
    end
    passing = filter(row -> row.validation_passed, rows)
    selected_mode_count = isempty(passing) ? 0 : minimum(getproperty.(passing, :mode_count))
    selected_states = selected_mode_count == 0 ? Matrix{Float64}[] : state_map[selected_mode_count]
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "validation_summary.csv"), rows)
    jldsave(
        joinpath(OUTPUT_ROOT, "o5a_validation.jld2");
        format_version=1,
        source_id=SOURCE_ID,
        mode_counts=MODE_COUNTS,
        reference_alpha=REFERENCE_ALPHA,
        electronics_sigma,
        rows,
        selected_mode_count,
        selected_states,
        test_opened=false,
    )
    readme_path = joinpath(OUTPUT_ROOT, "README.md")
    open(readme_path, "w") do io
        println(io, "# O5a smooth operator projection — validation freeze")
        println(io)
        println(io, "> Train/validation only. No test response is loaded by this script.")
        println(io)
        println(io, "| DCT modes M | State fidelity | Min principal cosine | Rank | Coherence | Pd@3 | Pdl@3 | Min Pd@3 | alpha90 val | Gate |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
        for row in rows
            println(io, @sprintf(
                "| %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %s |",
                row.mode_count, row.state_fidelity, row.minimum_principal_cosine,
                row.template_stable_rank, row.maximum_template_coherence,
                row.pd_reference, row.pdl_reference, row.minimum_case_pd_reference,
                row.alpha90_validation, row.validation_passed ? "pass" : "fail",
            ))
        end
        println(io)
        if selected_mode_count == 0
            println(io, "No `M<=12` projection passes validation; geometry synthesis remains blocked.")
        else
            println(io, "Frozen smallest passing projection: **M=$(selected_mode_count)**. A fresh physical test context is required before geometry synthesis.")
        end
    end
    println("[+] O5a validation report: $readme_path")
    rows
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o5a_validation()

end
