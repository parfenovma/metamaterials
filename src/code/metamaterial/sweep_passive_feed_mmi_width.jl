module SweepPassiveFeedMMIWidth

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_mmi")

include(joinpath(@__DIR__, "design_passive_feed_mmi.jl"))
using .DesignPassiveFeedMMI
using .DesignPassiveFeedMMI.ElasticPortModes
using LinearAlgebra

function evaluate_width(local_fundamental, width_mm)
    width_m = width_mm * 1e-3
    modes = solve_port_modes(PortModeConfig(
        height_m=width_m,
        density=DesignPassiveFeedMMI.DENSITY_KG_M3,
        pressure_wave_speed=DesignPassiveFeedMMI.PRESSURE_SPEED_M_S,
        shear_wave_speed=DesignPassiveFeedMMI.SHEAR_SPEED_M_S,
        frequency_hz=DesignPassiveFeedMMI.FREQUENCY_HZ,
        element_count=max(140, round(Int, 12width_mm)),
    ))
    symmetric = filter(
        mode -> mode.kind == :propagating && mode.direction == :right &&
                mode.parity == :symmetric,
        modes,
    )
    sort!(symmetric; by=mode -> real(mode.wavenumber_per_m))
    y = collect(range(-width_m / 2, width_m / 2; length=1401))
    initial_x, initial_y = DesignPassiveFeedMMI.embedded_mode(
        local_fundamental, y, [0.0],
    )
    output_center = (
        DesignPassiveFeedMMI.INPUT_WIDTH_M + DesignPassiveFeedMMI.OUTPUT_GAP_M
    ) / 2
    target_x, target_y = DesignPassiveFeedMMI.embedded_mode(
        local_fundamental, y, [-output_center, output_center],
    )
    wide_x = [ComplexF64[
        DesignPassiveFeedMMI.linear_interpolate(
            mode.y_m, mode.displacement_x, value,
        ) for value in y
    ] for mode in symmetric]
    wide_y = [ComplexF64[
        DesignPassiveFeedMMI.linear_interpolate(
            mode.y_m, mode.displacement_y, value,
        ) for value in y
    ] for mode in symmetric]
    count = length(symmetric)
    gram = ComplexF64[
        DesignPassiveFeedMMI.trapezoidal_inner(
            y, wide_x[row], wide_y[row], wide_x[column], wide_y[column],
        ) for row in 1:count, column in 1:count
    ]
    initial_overlap = ComplexF64[
        DesignPassiveFeedMMI.trapezoidal_inner(
            y, wide_x[index], wide_y[index], initial_x, initial_y,
        ) for index in 1:count
    ]
    coefficients = gram \ initial_overlap
    target_norm = real(DesignPassiveFeedMMI.trapezoidal_inner(
        y, target_x, target_y, target_x, target_y,
    ))
    best_coherence = -Inf
    best_length_mm = NaN
    for length_mm in range(0.0, 600.0; step=0.2)
        propagated = coefficients .* exp.(
            getproperty.(symmetric, :gamma_per_m) .* (length_mm * 1e-3),
        )
        field_x = sum(propagated[index] .* wide_x[index] for index in 1:count)
        field_y = sum(propagated[index] .* wide_y[index] for index in 1:count)
        field_norm = real(DesignPassiveFeedMMI.trapezoidal_inner(
            y, field_x, field_y, field_x, field_y,
        ))
        overlap = DesignPassiveFeedMMI.trapezoidal_inner(
            y, target_x, target_y, field_x, field_y,
        )
        coherence = abs(overlap) / sqrt(target_norm * field_norm)
        if coherence > best_coherence
            best_coherence = coherence
            best_length_mm = length_mm
        end
    end
    (
        width_mm=Float64(width_mm),
        symmetric_mode_count=count,
        best_length_mm,
        best_coherence,
    )
end

function run()
    local_modes = solve_port_modes(PortModeConfig(
        height_m=DesignPassiveFeedMMI.INPUT_WIDTH_M,
        density=DesignPassiveFeedMMI.DENSITY_KG_M3,
        pressure_wave_speed=DesignPassiveFeedMMI.PRESSURE_SPEED_M_S,
        shear_wave_speed=DesignPassiveFeedMMI.SHEAR_SPEED_M_S,
        frequency_hz=DesignPassiveFeedMMI.FREQUENCY_HZ,
        element_count=140,
    ))
    local_fundamental = fundamental_quasi_longitudinal(local_modes)
    rows = [
        evaluate_width(local_fundamental, width_mm)
        for width_mm in range(15.2, 40.0; step=0.8)
    ]
    best = rows[argmax(getproperty.(rows, :best_coherence))]
    mkpath(OUTPUT_ROOT)
    open(joinpath(OUTPUT_ROOT, "mmi_width_sweep.csv"), "w") do io
        println(io, "mmi_width_mm,symmetric_mode_count,best_length_mm,best_target_coherence")
        for row in rows
            println(io, join((
                row.width_mm, row.symmetric_mode_count,
                row.best_length_mm, row.best_coherence,
            ), ','))
        end
    end
    println("[+] best MMI width=$best")
    best
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
