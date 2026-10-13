module TopologyBendingMetalGRINMultimodeCheck

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV, "METAMATERIALS_BENDING_METAL_GRIN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "bending_metal_grin"),
)

include(joinpath(@__DIR__, "multimode_port_harmonic_solver.jl"))
using .MultimodePortHarmonicSolver
using .SinusoidalMaterialLens: aluminium_6061

function fundamental_amplitude(amplitudes, direction)
    symmetric = filter(a -> a.parity == :symmetric, amplitudes)
    length(symmetric) == 1 || error("expected one propagating symmetric mode")
    direction == :right ? only(symmetric).right_amplitude : only(symmetric).left_amplitude
end

function solve_case(variant, frequency_hz; evanescent_modes=8)
    mesh = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
    state = solve_multimode_port_harmonic(
        mesh, ["Microphone"], aluminium_6061();
        config=MultimodePortConfig(
            frequency_hz=frequency_hz,
            element_order=2,
            quadrature_degree=4,
            port_element_count=120,
            evanescent_modes_per_port=evanescent_modes,
        ),
    )
    incident = fundamental_amplitude(state.source_amplitudes, :right)
    reflected = fundamental_amplitude(state.source_amplitudes, :left)
    transmitted = fundamental_amplitude(only(state.output_amplitudes), :right)
    (
        variant,
        frequency_hz,
        evanescent_modes,
        incident_power=state.incident_power_w_per_m,
        reflected_power=state.reflected_power_w_per_m,
        transmitted_power=only(state.transmitted_power_w_per_m),
        R_all=state.reflected_power_w_per_m / state.incident_power_w_per_m,
        T_all=only(state.transmitted_power_w_per_m) / state.incident_power_w_per_m,
        power_sum=(state.reflected_power_w_per_m + only(state.transmitted_power_w_per_m)) /
                  state.incident_power_w_per_m,
        S11_complex=reflected / incident,
        S21_complex=transmitted / incident,
        source_basis_size=first(state.port_basis_sizes),
        output_basis_size=last(state.port_basis_sizes),
        max_gram_condition=maximum(state.port_gram_condition_numbers),
    )
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, c)) for c in columns), ','))
        end
    end
    println("[+] $path")
end

function main()
    frequencies = (162e3, 202e3, 242e3, 262e3, 282e3, 302e3, 322e3)
    rows = NamedTuple[]
    for variant in (:solid, :dense_c), frequency_hz in frequencies
        row = solve_case(variant, frequency_hz)
        push!(rows, row)
        println(
            "[+] multimode $variant $(frequency_hz / 1e3) kHz: T=",
            round(row.T_all; digits=4), ", R=", round(row.R_all; digits=4),
            ", sum=", round(row.power_sum; digits=5),
        )
    end
    # One deterministic evanescent-basis convergence anchor at the carrier.
    push!(rows, solve_case(:dense_c, 242e3; evanescent_modes=12))
    write_csv(joinpath(OUTPUT_ROOT, "g1_multimode_anchor_check.csv"), rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
