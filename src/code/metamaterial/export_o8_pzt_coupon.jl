module ExportO8PZTCoupon

using Gmsh: gmsh
using JLD2
using LinearAlgebra

include(joinpath(@__DIR__, "operator_o8_boundary_mesher.jl"))
include(joinpath(@__DIR__, "operator_o8_plane_strain_solver.jl"))

using .OperatorO8BoundaryMesher
using .OperatorO8PlaneStrainSolver

export export_o8_pzt_coupon

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "manufacturing", "o8_pzt_coupon_v1")
const CARRIER_HZ = 242.0e3
const PRESSURE_SPEED_M_S = 6122.102437409232
const SHEAR_SPEED_M_S = 3083.810277185563
const BLOCK_HEIGHT_M = 8PRESSURE_SPEED_M_S / CARRIER_HZ
const BLOCK_LENGTH_M = 2PRESSURE_SPEED_M_S / CARRIER_HZ
const BLOCK_THICKNESS_M = 8.0e-3
const PZT_LENGTH_M = 20.0e-3
const PZT_WIDTH_M = 8.0e-3
const PZT_THICKNESS_M = 1.0e-3
const INNER_CENTER_M = 0.19835729692786852 * BLOCK_HEIGHT_M / 2
const OUTER_CENTER_M = BLOCK_HEIGHT_M / 2 - PZT_LENGTH_M / 2
const PZT_CENTERS_M = (-OUTER_CENTER_M, -INNER_CENTER_M, INNER_CENTER_M, OUTER_CENTER_M)
const PZT_NAMES = ("B-", "A-", "A+", "B+")
const EPSILON_0_F_M = 8.8541878128e-12
const RELATIVE_PERMITTIVITY = 1900.0
const PZT_CAPACITANCE_F = EPSILON_0_F_M * RELATIVE_PERMITTIVITY *
    PZT_LENGTH_M * PZT_WIDTH_M / PZT_THICKNESS_M
const TARGET_A_HZ = 258.989e3
const TARGET_B_HZ = 253.411e3
const TARGET_Q_A = 6.382
const TARGET_Q_B = 4.232
const FREQUENCY_RATIOS = collect(range(0.8, 1.2; length=7))
const BITS = ((0, 0), (1, 0), (0, 1), (1, 1))

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

function tuned_load(target_hz, target_q)
    omega = 2pi * target_hz
    inductance_h = inv(omega^2 * 2PZT_CAPACITANCE_F)
    resistance_ohm = omega * inductance_h / target_q
    O8Load(; inductance_h, resistance_ohm)
end

function write_step(path; assembly)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add(assembly ? "o8_pzt_coupon_assembly" : "o8_pzt_coupon_block")
        block = gmsh.model.occ.addBox(
            0.0, -BLOCK_HEIGHT_M / 2, -BLOCK_THICKNESS_M / 2,
            BLOCK_LENGTH_M, BLOCK_HEIGHT_M, BLOCK_THICKNESS_M,
        )
        gmsh.model.setEntityName(3, block, "Al_6061_T6_block")
        if assembly
            for (name, center) in zip(PZT_NAMES, PZT_CENTERS_M)
                patch = gmsh.model.occ.addBox(
                    -PZT_THICKNESS_M, center - PZT_LENGTH_M / 2,
                    -PZT_WIDTH_M / 2, PZT_THICKNESS_M, PZT_LENGTH_M, PZT_WIDTH_M,
                )
                gmsh.model.setEntityName(3, patch, "PZT_$name")
            end
        end
        gmsh.model.occ.synchronize()
        gmsh.write(path)
    finally
        gmsh.finalize()
    end
end

function write_dimensioned_svg(path)
    height_mm = 1e3BLOCK_HEIGHT_M
    length_mm = 1e3BLOCK_LENGTH_M
    inner_mm = 1e3INNER_CENTER_M
    outer_mm = 1e3OUTER_CENTER_M
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="760" viewBox="0 0 1200 760">
<defs><marker id="a" markerWidth="8" markerHeight="8" refX="4" refY="4" orient="auto"><path d="M8,4 L0,0 L0,8 Z" fill="#222"/></marker></defs>
<style>text{font-family:Arial,sans-serif;fill:#111}.t{font-size:25px;font-weight:700}.n{font-size:18px}.d{stroke:#222;stroke-width:1.5;marker-start:url(#a);marker-end:url(#a)}.e{stroke:#777;stroke-width:1}.c{stroke:#111;stroke-dasharray:8 6}.pA{fill:#2693ff;stroke:#111}.pB{fill:#f39c34;stroke:#111}.al{fill:#e1e4e8;stroke:#111;stroke-width:2}</style>
<text x="50" y="45" class="t">O8 PZT coupon v1 — all dimensions in mm</text>
<text x="50" y="85" class="n">Plane-strain x-y section (not to scale), origin at geometric center</text>
<rect x="150" y="110" width="32" height="500" class="al"/>
<line x1="115" y1="360" x2="220" y2="360" class="c"/>
<rect x="142" y="110" width="8" height="49.41" class="pB"/><text x="75" y="140" class="n">B−</text>
<rect x="142" y="285.63" width="8" height="49.41" class="pA"/><text x="75" y="317" class="n">A−</text>
<rect x="142" y="384.96" width="8" height="49.41" class="pA"/><text x="75" y="416" class="n">A+</text>
<rect x="142" y="560.59" width="8" height="49.41" class="pB"/><text x="75" y="592" class="n">B+</text>
<line x1="250" y1="110" x2="250" y2="610" class="d"/><line x1="182" y1="110" x2="265" y2="110" class="e"/><line x1="182" y1="610" x2="265" y2="610" class="e"/><text x="270" y="365" class="n">H = $(round(height_mm; digits=2))</text>
<text x="300" y="160" class="n">4 × CTS 3195HD, 20.00 × 8.00 × 1.00</text>
<text x="300" y="195" class="n">centers y: ±$(round(inner_mm; digits=2)), ±$(round(outer_mm; digits=2))</text>
<text x="300" y="230" class="n">B plates flush with the top/bottom edges</text>
<text x="650" y="85" class="n">Side view</text>
<rect x="650" y="250" width="405" height="64" class="al"/><rect x="642" y="250" width="8" height="64" class="pA"/>
<line x1="650" y1="360" x2="1055" y2="360" class="d"/><line x1="650" y1="314" x2="650" y2="375" class="e"/><line x1="1055" y1="314" x2="1055" y2="375" class="e"/><text x="790" y="393" class="n">L = $(round(length_mm; digits=2))</text>
<line x1="1095" y1="250" x2="1095" y2="314" class="d"/><text x="1110" y="288" class="n">8.00</text>
<text x="642" y="230" class="n">1.00 PZT</text>
<text x="650" y="455" class="n">Block: Al 6061-T6, $(round(length_mm; digits=2)) × $(round(height_mm; digits=2)) × 8.00</text>
<text x="650" y="490" class="n">Machining: ±0.05 overall; source face flatness 0.03; Ra ≤ 1.6 µm</text>
<text x="650" y="525" class="n">Bond line: EPO-TEK 301-2, target 25–50 µm</text>
<text x="50" y="690" class="n">Blue pair A: equal electrical polarity. Orange pair B: reverse one plate electrically.</text>
<text x="50" y="725" class="n">Do not infer electrode polarity from this drawing; verify every plate with the vendor certificate and mark the positive electrode before bonding.</text>
</svg>"""
    open(path, "w") do io
        write(io, svg)
    end
end

function write_wiring_svg(path, load_a, load_b)
    svg = """<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="650" viewBox="0 0 1200 650">
<style>text{font-family:Arial,sans-serif;fill:#111}.t{font-size:25px;font-weight:700}.n{font-size:18px}.wire{stroke:#111;stroke-width:3;fill:none}.box{fill:#f5f7f9;stroke:#111;stroke-width:2}.pA{fill:#d8ebff;stroke:#1777c8;stroke-width:2}.pB{fill:#fff0dc;stroke:#d87500;stroke-width:2}</style>
<text x="40" y="45" class="t">Two-bit passive shunt — bench wiring</text>
<rect x="55" y="105" width="165" height="120" class="pA"/><text x="83" y="145" class="n">A− ∥ A+</text><text x="78" y="178" class="n">same polarity</text>
<line x1="220" y1="135" x2="310" y2="135" stroke="#111" stroke-width="3"/><line x1="220" y1="195" x2="310" y2="195" stroke="#111" stroke-width="3"/>
<rect x="310" y="100" width="150" height="130" class="box"/><text x="333" y="140" class="n">K_A NO</text><text x="326" y="174" class="n">Omron</text><text x="326" y="202" class="n">G5V-1 5 V</text>
<line x1="460" y1="135" x2="525" y2="135" stroke="#111" stroke-width="3"/><rect x="525" y="105" width="125" height="60" class="box"/><text x="548" y="142" class="n">L_A</text><line x1="650" y1="135" x2="710" y2="135" stroke="#111" stroke-width="3"/><rect x="710" y="105" width="125" height="60" class="box"/><text x="733" y="142" class="n">R_A</text><path d="M835 135 H900 V195 H460" stroke="#111" stroke-width="3" fill="none"/>
<text x="505" y="255" class="n">target total: L=$(round(1e6load_a.inductance_h; digits=2)) µH, R=$(round(load_a.resistance_ohm; digits=2)) Ω</text>
<rect x="55" y="335" width="165" height="120" class="pB"/><text x="78" y="375" class="n">B− ∥ B+</text><text x="73" y="408" class="n">one reversed</text>
<line x1="220" y1="365" x2="310" y2="365" stroke="#111" stroke-width="3"/><line x1="220" y1="425" x2="310" y2="425" stroke="#111" stroke-width="3"/>
<rect x="310" y="330" width="150" height="130" class="box"/><text x="333" y="370" class="n">K_B NO</text><text x="326" y="404" class="n">Omron</text><text x="326" y="432" class="n">G5V-1 5 V</text>
<line x1="460" y1="365" x2="525" y2="365" stroke="#111" stroke-width="3"/><rect x="525" y="335" width="125" height="60" class="box"/><text x="548" y="372" class="n">L_B</text><line x1="650" y1="365" x2="710" y2="365" stroke="#111" stroke-width="3"/><rect x="710" y="335" width="125" height="60" class="box"/><text x="733" y="372" class="n">R_B</text><path d="M835 365 H900 V425 H460" stroke="#111" stroke-width="3" fill="none"/>
<text x="505" y="485" class="n">target total: L=$(round(1e6load_b.inductance_h; digits=2)) µH, R=$(round(load_b.resistance_ohm; digits=2)) Ω</text>
<rect x="955" y="105" width="195" height="350" class="box"/><text x="980" y="145" class="n">5 V coil side</text><text x="980" y="185" class="n">SW_A → K_A</text><text x="980" y="220" class="n">SW_B → K_B</text><text x="980" y="270" class="n">1N4148 across</text><text x="980" y="305" class="n">each coil</text><text x="980" y="355" class="n">galvanically</text><text x="980" y="390" class="n">isolated from</text><text x="980" y="425" class="n">PZT ports</text>
<text x="55" y="555" class="n">q00: both relays open; q10: A closed; q01: B closed; q11: both closed.</text>
<text x="55" y="590" class="n">L and R are total measured series values at the PZT terminals, including inductor ESR and relay contact resistance.</text>
</svg>"""
    open(path, "w") do io
        write(io, svg)
    end
end

function build_validation_mesh(path)
    half = BLOCK_HEIGHT_M / 2
    config = O8BoundaryGeometry(
        aperture_height_m=BLOCK_HEIGHT_M,
        block_length_m=BLOCK_LENGTH_M,
        position_a=INNER_CENTER_M / half,
        width_a=(PZT_LENGTH_M / 2) / half,
        position_b=OUTER_CENTER_M / half,
        width_b=(PZT_LENGTH_M / 2) / half,
        mesh_size_m=(SHEAR_SPEED_M_S / CARRIER_HZ) / 6,
    )
    build_o8_boundary_mesh(path; config)
end

function validate_frozen_coupon(mesh_path, load_a, load_b)
    pzt = O8PZTBoundary(
        thickness_m=PZT_THICKNESS_M,
        density_kg_m3=7800.0,
        c33e_pa=117.0e9,
        mechanical_quality_factor=80.0,
        coupling_squared=0.49^2,
        capacitance_per_patch_f=PZT_CAPACITANCE_F,
    )
    results = NamedTuple[]
    for (index, ratio) in enumerate(FREQUENCY_RATIOS)
        frequency_hz = CARRIER_HZ * ratio
        println("[>] coupon validation $index/$(length(FREQUENCY_RATIOS)): $(frequency_hz / 1e3) kHz")
        config = O8PlaneStrainConfig(
            ; frequency_hz, aperture_height_m=BLOCK_HEIGHT_M,
            block_length_m=BLOCK_LENGTH_M,
        )
        result = solve_o8_frequency(
            mesh_path; mechanism=:pzt, config, pzt, load_a, load_b,
        )
        push!(results, merge((; frequency_index=index, frequency_ratio=ratio), result))
    end
    carrier = results[4]
    minimum_amplitude = minimum(
        result.amplitude_ratios[state] for result in results for state in 1:4
    )
    maximum_amplitude = maximum(
        result.amplitude_ratios[state] for result in results for state in 1:4
    )
    carrier_minimum_state_distance = minimum(carrier.state_distances)
    maximum_mirror_error = maximum(
        result.raw_mirror_errors[state] for result in results for state in 1:4
    )
    maximum_reciprocity_error = maximum(
        result.reciprocity_errors[state] for result in results for state in 1:4
    )
    passed = minimum_amplitude >= 0.70 && maximum_amplitude <= 1.50 &&
        carrier_minimum_state_distance >= 0.05 && maximum_mirror_error <= 0.02 &&
        maximum_reciprocity_error <= 1.0e-6
    summary = (;
        material="CTS 3195HD", minimum_amplitude, maximum_amplitude,
        carrier_minimum_state_distance, maximum_mirror_error,
        maximum_reciprocity_error, passed,
    )
    rows = NamedTuple[]
    for result in results, state in 1:4
        push!(rows, (
            frequency_hz=result.frequency_hz,
            frequency_ratio=result.frequency_ratio,
            state="q$(BITS[state][1])$(BITS[state][2])",
            amplitude_reference_ratio=result.amplitude_ratios[state],
            raw_mirror_error=result.raw_mirror_errors[state],
            reciprocity_error=result.reciprocity_errors[state],
        ))
    end
    write_csv(joinpath(OUTPUT_ROOT, "model_validation_band.csv"), rows)
    jldsave(joinpath(OUTPUT_ROOT, "model_validation.jld2"); summary, results)
    summary
end

function write_placement_files(load_a, load_b)
    placement = [
        (
            patch=name, pair=first(name), center_y_mm=1e3center,
            lower_y_mm=1e3(center - PZT_LENGTH_M / 2),
            upper_y_mm=1e3(center + PZT_LENGTH_M / 2),
            length_mm=1e3PZT_LENGTH_M, width_z_mm=1e3PZT_WIDTH_M,
            thickness_x_mm=1e3PZT_THICKNESS_M,
            electrical_phase=startswith(name, "A") ? "same" :
                (name == "B-" ? "reference" : "reversed"),
        ) for (name, center) in zip(PZT_NAMES, PZT_CENTERS_M)
    ]
    write_csv(joinpath(OUTPUT_ROOT, "pzt_placement.csv"), placement)
    receiver_pitch = BLOCK_HEIGHT_M / 17
    receiver_rows = [
        (index, y_mm=1e3(-BLOCK_HEIGHT_M / 2 + (index - 0.5) * receiver_pitch))
        for index in 1:17
    ]
    write_csv(joinpath(OUTPUT_ROOT, "receiver_scan_points.csv"), receiver_rows)
    source_rows = [
        (port=1, y_mm=-50.0), (port=2, y_mm=0.0), (port=3, y_mm=50.0),
    ]
    write_csv(joinpath(OUTPUT_ROOT, "source_points.csv"), source_rows)
    shunt_rows = [
        (
            pair="A", target_resonance_hz=TARGET_A_HZ, target_q=TARGET_Q_A,
            nominal_patch_capacitance_nf=1e9PZT_CAPACITANCE_F,
            target_pair_capacitance_nf=2e9PZT_CAPACITANCE_F,
            target_total_inductance_uh=1e6load_a.inductance_h,
            target_total_series_resistance_ohm=load_a.resistance_ohm,
        ),
        (
            pair="B", target_resonance_hz=TARGET_B_HZ, target_q=TARGET_Q_B,
            nominal_patch_capacitance_nf=1e9PZT_CAPACITANCE_F,
            target_pair_capacitance_nf=2e9PZT_CAPACITANCE_F,
            target_total_inductance_uh=1e6load_b.inductance_h,
            target_total_series_resistance_ohm=load_b.resistance_ohm,
        ),
    ]
    write_csv(joinpath(OUTPUT_ROOT, "shunt_targets.csv"), shunt_rows)
end

function export_o8_pzt_coupon()
    mkpath(OUTPUT_ROOT)
    load_a = tuned_load(TARGET_A_HZ, TARGET_Q_A)
    load_b = tuned_load(TARGET_B_HZ, TARGET_Q_B)
    write_step(joinpath(OUTPUT_ROOT, "coupon_block.step"); assembly=false)
    write_step(joinpath(OUTPUT_ROOT, "coupon_assembly.step"); assembly=true)
    write_dimensioned_svg(joinpath(OUTPUT_ROOT, "dimensioned_drawing.svg"))
    write_wiring_svg(joinpath(OUTPUT_ROOT, "wiring.svg"), load_a, load_b)
    write_placement_files(load_a, load_b)
    mesh_path = joinpath(OUTPUT_ROOT, "coupon_validation.msh")
    build_validation_mesh(mesh_path)
    summary = validate_frozen_coupon(mesh_path, load_a, load_b)
    println("[+] coupon manufacturing gate: $(summary)")
    summary
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    export_o8_pzt_coupon()
end

end
