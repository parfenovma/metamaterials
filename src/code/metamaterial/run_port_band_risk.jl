module PortBandRisk

ENV["GKSwstype"] = "100"

using Plots
using Printf

include(joinpath(@__DIR__, "port_mode_solver.jl"))
using .ElasticPortModes
include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .ImpulseRiskAnalysis

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_PORT_BAND_RISK_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "port_band_risk"),
)
const ELEMENT_COUNT = parse(Int, get(ENV, "METAMATERIALS_PORT_BAND_ELEMENTS", "50"))
const HEIGHTS_MM = (2.8, 3.0, 3.2, 3.4, 3.6, 4.2)

function right_going_modes(frequency_hz, height_m)
    solve_port_modes(PortModeConfig(
        height_m=height_m,
        frequency_hz=frequency_hz,
        element_count=ELEMENT_COUNT,
    )) |> modes -> propagating_modes(modes; direction=:right)
end

symmetric_count(frequency_hz, height_m) = count(
    mode -> mode.parity == :symmetric,
    right_going_modes(frequency_hz, height_m),
)

function cuton_frequency_hz(height_m; lower_hz=100.0e3, upper_hz=600.0e3)
    symmetric_count(lower_hz, height_m) == 1 ||
        error("lower cut-on bracket is not single-symmetric-mode")
    symmetric_count(upper_hz, height_m) > 1 || return NaN
    lower = lower_hz
    upper = upper_hz
    while upper - lower > 50.0
        middle = (lower + upper) / 2
        if symmetric_count(middle, height_m) > 1
            upper = middle
        else
            lower = middle
        end
    end
    upper
end

function write_rows(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function reference_band_rows(spectrum, lower_hz, upper_hz)
    frequencies = collect(floor(lower_hz / 2e3) * 2e3:2e3:ceil(upper_hz / 2e3) * 2e3)
    rows = NamedTuple[]
    for frequency_hz in frequencies
        modes = right_going_modes(frequency_hz, 4.2e-3)
        for (mode_index, mode) in enumerate(modes)
            pulse_index = argmin(abs.(spectrum.frequency_hz .- frequency_hz))
            push!(rows, (
                frequency_hz,
                mode_index,
                parity=string(mode.parity),
                wavenumber_per_m=real(mode.wavenumber_per_m),
                phase_velocity_m_per_s=2pi * frequency_hz / abs(real(mode.wavenumber_per_m)),
                p_fraction=mode.p_fraction,
                axial_fraction=mode.axial_displacement_fraction,
                pulse_relative_amplitude=spectrum.relative_amplitude[pulse_index],
            ))
        end
    end
    rows
end

function track_reference_branch(lower_hz, upper_hz)
    frequencies = collect(10.0e3:2.0e3:ceil(upper_hz / 2e3) * 2e3)
    branch = track_fundamental_branch(
        PortModeConfig(frequency_hz=first(frequencies), element_count=ELEMENT_COUNT),
        frequencies,
    )
    wavenumbers = real.(getproperty.(branch.modes, :wavenumber_per_m))
    group_velocity = fill(NaN, length(frequencies))
    for index in 2:(length(frequencies) - 1)
        dk_df = (wavenumbers[index + 1] - wavenumbers[index - 1]) /
                (frequencies[index + 1] - frequencies[index - 1])
        group_velocity[index] = 2pi / dk_df
    end
    [(
        frequency_hz=frequencies[index],
        wavenumber_per_m=wavenumbers[index],
        phase_velocity_m_per_s=2pi * frequencies[index] / wavenumbers[index],
        group_velocity_m_per_s=group_velocity[index],
        p_fraction=branch.modes[index].p_fraction,
        axial_fraction=branch.modes[index].axial_displacement_fraction,
        overlap=branch.overlaps[index],
    ) for index in eachindex(frequencies) if frequencies[index] >= lower_hz]
end

function save_mode_plot(mode_rows, branch_rows, spectrum, band20, cuton_hz)
    panel = plot(
        xlabel="frequency, kHz",
        ylabel="wavenumber k, 1/m",
        title="Propagating modes of the 4.2 mm strip across B-20",
        legend=:topleft,
    )
    for (parity, color) in (("symmetric", :royalblue), ("antisymmetric", :darkorange))
        selected = filter(row -> row.parity == parity, mode_rows)
        scatter!(
            panel,
            getproperty.(selected, :frequency_hz) ./ 1e3,
            abs.(getproperty.(selected, :wavenumber_per_m));
            marker_z=getproperty.(selected, :axial_fraction),
            color=color,
            markersize=3,
            label=parity,
        )
    end
    vline!(panel, [cuton_hz / 1e3]; color=:black, linestyle=:dash, label="extra symmetric cut-on")
    vspan!(panel, [band20.lower_hz / 1e3, band20.upper_hz / 1e3];
           color=:gray, alpha=0.08, label="pulse B-20")

    group = plot(
        getproperty.(branch_rows, :frequency_hz) ./ 1e3,
        getproperty.(branch_rows, :group_velocity_m_per_s);
        xlabel="frequency, kHz",
        ylabel="group velocity, m/s",
        title="Tracked fundamental symmetric branch",
        linewidth=2,
        label="v_g",
    )
    vline!(group, [cuton_hz / 1e3]; color=:black, linestyle=:dash, label="cut-on")
    savefig(plot(panel, group; layout=(2, 1), size=(1050, 850)),
            joinpath(OUTPUT_ROOT, "port_band_map.png"))
end

function save_height_plot(rows, band20)
    panel = plot(
        getproperty.(rows, :height_mm),
        getproperty.(rows, :extra_symmetric_cuton_hz) ./ 1e3;
        marker=:circle,
        linewidth=2,
        xlabel="strip height, mm",
        ylabel="first extra symmetric cut-on, kHz",
        title="Single-symmetric-mode margin for the provisional pulse",
        label="computed cut-on",
    )
    hline!(panel, [band20.upper_hz / 1e3]; color=:firebrick, linestyle=:dash,
           label="B-20 upper edge")
    savefig(panel, joinpath(OUTPUT_ROOT, "port_height_cuton.png"))
end

function run()
    mkpath(OUTPUT_ROOT)
    spectrum = pulse_spectrum()
    band20 = spectral_band(spectrum, -20.0)
    height_rows = [begin
        cuton = cuton_frequency_hz(height_mm * 1e-3)
        (
            height_mm,
            extra_symmetric_cuton_hz=cuton,
            margin_above_band20_hz=cuton - band20.upper_hz,
            single_symmetric_mode_in_band20=cuton > band20.upper_hz,
        )
    end for height_mm in HEIGHTS_MM]
    write_rows(joinpath(OUTPUT_ROOT, "port_height_cuton.csv"), height_rows)

    reference_cuton = only(filter(row -> row.height_mm == 4.2, height_rows)).extra_symmetric_cuton_hz
    mode_rows = reference_band_rows(spectrum, band20.lower_hz, band20.upper_hz)
    branch_rows = track_reference_branch(band20.lower_hz, band20.upper_hz)
    write_rows(joinpath(OUTPUT_ROOT, "port_modes_across_pulse_band.csv"), mode_rows)
    write_rows(joinpath(OUTPUT_ROOT, "fundamental_branch_across_pulse_band.csv"), branch_rows)
    save_mode_plot(mode_rows, branch_rows, spectrum, band20, reference_cuton)
    save_height_plot(height_rows, band20)

    energy_multimode = spectral_energy_fraction(spectrum; lower_hz=reference_cuton)
    carrier_row = branch_rows[argmin(abs.(getproperty.(branch_rows, :frequency_hz) .- 242.0e3))]
    @printf("[+] B-20 upper edge: %.2f kHz\n", band20.upper_hz / 1e3)
    @printf("[+] 4.2 mm first extra symmetric cut-on: %.2f kHz\n", reference_cuton / 1e3)
    @printf("[+] pulse energy at/above cut-on: %.2f%%\n", 100energy_multimode)
    @printf("[+] tracked branch vg at carrier: %.1f m/s\n", carrier_row.group_velocity_m_per_s)
    for row in height_rows
        @printf(
            "[+] h=%.1f mm: cut-on %.2f kHz, margin %.2f kHz, pass=%s\n",
            row.height_mm,
            row.extra_symmetric_cuton_hz / 1e3,
            row.margin_above_band20_hz / 1e3,
            row.single_symmetric_mode_in_band20,
        )
    end
    println("[+] $OUTPUT_ROOT")
    (; spectrum, band20, height_rows, mode_rows, branch_rows, energy_multimode)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
