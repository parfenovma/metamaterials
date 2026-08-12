module DesignPassiveFeedMMI

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_mmi")

using LinearAlgebra
using JLD2
include(joinpath(@__DIR__, "port_mode_solver.jl"))
using .ElasticPortModes

const FREQUENCY_HZ = 242.0e3
const DENSITY_KG_M3 = 2700.0
const PRESSURE_SPEED_M_S = 6122.102437409232
const SHEAR_SPEED_M_S = 3083.810277185563
const INPUT_WIDTH_M = 7.0e-3
const MMI_WIDTH_M = 15.2e-3
const OUTPUT_GAP_M = 1.2e-3

function linear_interpolate(x, values, query)
    query <= first(x) && return values[1]
    query >= last(x) && return values[end]
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function trapezoidal_inner(y, first_x, first_y, second_x, second_y)
    values = conj.(first_x) .* second_x .+ conj.(first_y) .* second_y
    sum(
        (y[index + 1] - y[index]) * (values[index] + values[index + 1]) / 2
        for index in 1:(length(y) - 1)
    )
end

function embedded_mode(mode, y, centers)
    ux = zeros(ComplexF64, length(y))
    uy = zeros(ComplexF64, length(y))
    half_width = INPUT_WIDTH_M / 2
    for center in centers, index in eachindex(y)
        local_y = y[index] - center
        if -half_width <= local_y <= half_width
            ux[index] += linear_interpolate(mode.y_m, mode.displacement_x, local_y)
            uy[index] += linear_interpolate(mode.y_m, mode.displacement_y, local_y)
        end
    end
    ux, uy
end

function run()
    local_modes = solve_port_modes(PortModeConfig(
        height_m=INPUT_WIDTH_M,
        density=DENSITY_KG_M3,
        pressure_wave_speed=PRESSURE_SPEED_M_S,
        shear_wave_speed=SHEAR_SPEED_M_S,
        frequency_hz=FREQUENCY_HZ,
        element_count=140,
    ))
    local_fundamental = fundamental_quasi_longitudinal(local_modes)
    wide_modes = solve_port_modes(PortModeConfig(
        height_m=MMI_WIDTH_M,
        density=DENSITY_KG_M3,
        pressure_wave_speed=PRESSURE_SPEED_M_S,
        shear_wave_speed=SHEAR_SPEED_M_S,
        frequency_hz=FREQUENCY_HZ,
        element_count=240,
    ))
    symmetric = filter(
        mode -> mode.kind == :propagating && mode.direction == :right &&
                mode.parity == :symmetric,
        wide_modes,
    )
    sort!(symmetric; by=mode -> real(mode.wavenumber_per_m))
    y = collect(range(-MMI_WIDTH_M / 2, MMI_WIDTH_M / 2; length=1201))
    initial_x, initial_y = embedded_mode(local_fundamental, y, [0.0])
    output_center = (INPUT_WIDTH_M + OUTPUT_GAP_M) / 2
    target_x, target_y = embedded_mode(
        local_fundamental, y, [-output_center, output_center],
    )
    wide_x = [ComplexF64[
        linear_interpolate(mode.y_m, mode.displacement_x, value) for value in y
    ] for mode in symmetric]
    wide_y = [ComplexF64[
        linear_interpolate(mode.y_m, mode.displacement_y, value) for value in y
    ] for mode in symmetric]
    count = length(symmetric)
    gram = ComplexF64[
        trapezoidal_inner(y, wide_x[row], wide_y[row], wide_x[column], wide_y[column])
        for row in 1:count, column in 1:count
    ]
    initial_overlap = ComplexF64[
        trapezoidal_inner(y, wide_x[index], wide_y[index], initial_x, initial_y)
        for index in 1:count
    ]
    coefficients = gram \ initial_overlap
    target_norm = real(trapezoidal_inner(y, target_x, target_y, target_x, target_y))
    rows = NamedTuple[]
    for length_mm in range(0.0, 500.0; step=0.05)
        propagated = coefficients .* exp.(getproperty.(symmetric, :gamma_per_m) .* (length_mm * 1e-3))
        field_x = sum(propagated[index] .* wide_x[index] for index in 1:count)
        field_y = sum(propagated[index] .* wide_y[index] for index in 1:count)
        field_norm = real(trapezoidal_inner(y, field_x, field_y, field_x, field_y))
        overlap = trapezoidal_inner(y, target_x, target_y, field_x, field_y)
        coherence = abs(overlap) / sqrt(target_norm * field_norm)
        push!(rows, (; length_mm, coherence, overlap_phase_deg=rad2deg(angle(overlap))))
    end
    best = rows[argmax(getproperty.(rows, :coherence))]
    mkpath(OUTPUT_ROOT)
    open(joinpath(OUTPUT_ROOT, "mmi_modal_length_sweep.csv"), "w") do io
        println(io, "length_mm,target_coherence,overlap_phase_deg")
        for row in rows
            println(io, "$(row.length_mm),$(row.coherence),$(row.overlap_phase_deg)")
        end
    end
    JLD2.jldsave(
        joinpath(OUTPUT_ROOT, "mmi_modal_design.jld2");
        format_version=1,
        frequency_hz=FREQUENCY_HZ,
        input_width_mm=INPUT_WIDTH_M * 1e3,
        mmi_width_mm=MMI_WIDTH_M * 1e3,
        output_gap_mm=OUTPUT_GAP_M * 1e3,
        symmetric_wavenumber_per_m=getproperty.(symmetric, :wavenumber_per_m),
        excitation_coefficients=coefficients,
        best_length_mm=best.length_mm,
        best_target_coherence=best.coherence,
        best_overlap_phase_deg=best.overlap_phase_deg,
    )
    println("[+] MMI symmetric k=$(real.(getproperty.(symmetric, :wavenumber_per_m))) 1/m")
    println("[+] MMI coefficients=$coefficients")
    println("[+] best MMI length=$(best.length_mm) mm, coherence=$(best.coherence)")
    best
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
