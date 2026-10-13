module RunOperatorO8MaterialFreeze

using JLD2
using LinearAlgebra

include(joinpath(@__DIR__, "operator_o8_plane_strain_solver.jl"))

using .OperatorO8PlaneStrainSolver

export run_o8_material_freeze

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const O8_ROOT = joinpath(PROJECT_ROOT, "results", "operator_space", "o8_plane_strain_v0")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "results", "operator_space", "o8_material_freeze_v0")
const MESH_PATH = joinpath(O8_ROOT, "meshes", "pzt_base.msh")
const CARRIER_HZ = 242.0e3
const PRESSURE_SPEED_M_S = 6122.102437409232
const APERTURE_HEIGHT_M = 8PRESSURE_SPEED_M_S / CARRIER_HZ
const BLOCK_LENGTH_M = 2PRESSURE_SPEED_M_S / CARRIER_HZ
const PATCH_LENGTH_M = 0.5 * (20.0819e-3 + 20.4241e-3)
const COUPON_THICKNESS_M = 8.0e-3
const PZT_THICKNESS_M = 1.0e-3
const EPSILON_0_F_M = 8.8541878128e-12
const TARGET_A_HZ = 258.989e3
const TARGET_B_HZ = 253.411e3
const TARGET_Q_A = 6.382
const TARGET_Q_B = 4.232
const FREQUENCY_RATIOS = collect(range(0.8, 1.2; length=7))
const BITS = ((0, 0), (1, 0), (0, 1), (1, 1))

const MATERIALS = [
    (
        key="cts_3195hd_bulk",
        description="CTS 3195HD bulk PZT",
        relative_permittivity=1900.0,
        thickness_coupling=0.49,
        mechanical_quality_factor=80.0,
        c33e_pa=117.0e9,
        procurement_class="custom bulk plate",
    ),
    (
        key="cts_3203hd_bulk",
        description="CTS 3203HD bulk PZT",
        relative_permittivity=3800.0,
        thickness_coupling=0.55,
        mechanical_quality_factor=50.0,
        c33e_pa=117.0e9,
        procurement_class="custom bulk plate",
    ),
    (
        key="pi_pic255_1_3_floor",
        description="PI PIC255 1-3 composite, vendor lower-bound coupling",
        relative_permittivity=0.69 * 1800.0,
        thickness_coupling=0.60,
        mechanical_quality_factor=20.0,
        c33e_pa=65.0e9,
        procurement_class="custom 1-3 composite",
    ),
]

function write_csv(path, rows)
    isempty(rows) && return
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function tuned_load(capacitance_per_patch_f, target_hz, target_q)
    pair_capacitance_f = 2capacitance_per_patch_f
    omega = 2pi * target_hz
    inductance_h = inv(omega^2 * pair_capacitance_f)
    resistance_ohm = omega * inductance_h / target_q
    O8Load(; inductance_h, resistance_ohm)
end

function material_definition(material)
    capacitance_per_patch_f = EPSILON_0_F_M * material.relative_permittivity *
        PATCH_LENGTH_M * COUPON_THICKNESS_M / PZT_THICKNESS_M
    pzt = O8PZTBoundary(
        thickness_m=PZT_THICKNESS_M,
        density_kg_m3=material.key == "pi_pic255_1_3_floor" ? 6300.0 : 7800.0,
        c33e_pa=material.c33e_pa,
        mechanical_quality_factor=material.mechanical_quality_factor,
        coupling_squared=material.thickness_coupling^2,
        capacitance_per_patch_f=capacitance_per_patch_f,
    )
    load_a = tuned_load(capacitance_per_patch_f, TARGET_A_HZ, TARGET_Q_A)
    load_b = tuned_load(capacitance_per_patch_f, TARGET_B_HZ, TARGET_Q_B)
    (; material, pzt, load_a, load_b, capacitance_per_patch_f)
end

function solve_material(definition)
    results = NamedTuple[]
    for (index, ratio) in enumerate(FREQUENCY_RATIOS)
        frequency_hz = CARRIER_HZ * ratio
        println("[>] material freeze $(definition.material.key) $(index)/$(length(FREQUENCY_RATIOS)): $(frequency_hz / 1e3) kHz")
        config = O8PlaneStrainConfig(
            ; frequency_hz, aperture_height_m=APERTURE_HEIGHT_M,
            block_length_m=BLOCK_LENGTH_M,
        )
        result = solve_o8_frequency(
            MESH_PATH; mechanism=:pzt, config, pzt=definition.pzt,
            load_a=definition.load_a, load_b=definition.load_b,
        )
        push!(results, merge((; frequency_index=index, frequency_ratio=ratio), result))
    end
    results
end

function run_o8_material_freeze()
    isfile(MESH_PATH) || error("missing frozen O8 mesh: $MESH_PATH")
    mkpath(OUTPUT_ROOT)
    summary_rows = NamedTuple[]
    state_rows = NamedTuple[]
    distance_rows = NamedTuple[]
    bundles = Dict{String,Any}()
    for material in MATERIALS
        definition = material_definition(material)
        results = solve_material(definition)
        for result in results, state in 1:4
            push!(state_rows, (
                material=material.key,
                frequency_hz=result.frequency_hz,
                frequency_ratio=result.frequency_ratio,
                state="q$(BITS[state][1])$(BITS[state][2])",
                amplitude_reference_ratio=result.amplitude_ratios[state],
                raw_mirror_error=result.raw_mirror_errors[state],
                reciprocity_error=result.reciprocity_errors[state],
            ))
        end
        pairs = [(first, second) for first in 1:4 for second in (first + 1):4]
        for result in results, (pair_index, pair) in enumerate(pairs)
            push!(distance_rows, (
                material=material.key,
                frequency_hz=result.frequency_hz,
                frequency_ratio=result.frequency_ratio,
                first_state="q$(BITS[pair[1]][1])$(BITS[pair[1]][2])",
                second_state="q$(BITS[pair[2]][1])$(BITS[pair[2]][2])",
                projector_distance=result.state_distances[pair_index],
            ))
        end
        carrier = results[4]
        minimum_amplitude = minimum(
            result.amplitude_ratios[state] for result in results for state in 1:4
        )
        maximum_amplitude = maximum(
            result.amplitude_ratios[state] for result in results for state in 1:4
        )
        carrier_minimum_state_distance = minimum(carrier.state_distances)
        build_pass = minimum_amplitude >= 0.70 && maximum_amplitude <= 1.50 &&
            carrier_minimum_state_distance >= 0.05 &&
            maximum(result.raw_mirror_errors[state] for result in results for state in 1:4) <= 0.02 &&
            maximum(result.reciprocity_errors[state] for result in results for state in 1:4) <= 1.0e-6
        push!(summary_rows, (
            material=material.key,
            description=replace(material.description, ',' => ';'),
            procurement_class=replace(material.procurement_class, ',' => ';'),
            relative_permittivity=material.relative_permittivity,
            thickness_coupling=material.thickness_coupling,
            mechanical_quality_factor=material.mechanical_quality_factor,
            capacitance_per_patch_nf=1e9 * definition.capacitance_per_patch_f,
            inductance_a_uh=1e6 * definition.load_a.inductance_h,
            resistance_a_ohm=definition.load_a.resistance_ohm,
            inductance_b_uh=1e6 * definition.load_b.inductance_h,
            resistance_b_ohm=definition.load_b.resistance_ohm,
            minimum_amplitude,
            maximum_amplitude,
            carrier_minimum_state_distance,
            build_pass,
        ))
        bundles[material.key] = (; definition, results, build_pass)
    end
    write_csv(joinpath(OUTPUT_ROOT, "material_summary.csv"), summary_rows)
    write_csv(joinpath(OUTPUT_ROOT, "state_band_metrics.csv"), state_rows)
    write_csv(joinpath(OUTPUT_ROOT, "state_distances.csv"), distance_rows)
    jldsave(joinpath(OUTPUT_ROOT, "material_freeze.jld2"); summary_rows, bundles)
    summary_rows
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    for row in run_o8_material_freeze()
        println(row)
    end
end

end
