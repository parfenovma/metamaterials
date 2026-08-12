module PortModeCalibration

ENV["GKSwstype"] = "100"

using JLD2
using Plots
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_PORT_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "port_mode_calibration"),
)
const FREQUENCY_HZ = parse(Float64, get(ENV, "METAMATERIALS_PORT_FREQUENCY_HZ", "242000"))
const ELEMENT_COUNT = parse(Int, get(ENV, "METAMATERIALS_PORT_ELEMENTS", "80"))
const CONTINUATION_STEP_HZ = parse(
    Float64,
    get(ENV, "METAMATERIALS_PORT_CONTINUATION_STEP_HZ", "4000"),
)

include(joinpath(@__DIR__, "port_mode_solver.jl"))
using .ElasticPortModes

function magnitude_profiles(mode)
    values = vcat(mode.displacement_x, mode.displacement_y)
    scale = maximum(abs, values)
    abs.(mode.displacement_x) ./ scale,
    abs.(mode.displacement_y) ./ scale
end

function write_summary(path, modes)
    open(path, "w") do io
        println(io, "index,kind,direction,k_real_per_m,k_imag_per_m,gamma_real_per_m,gamma_imag_per_m,power_w_per_m,parity,parity_score,p_fraction,s_fraction,axial_fraction,residual")
        for (index, mode) in enumerate(modes)
            @printf(
                io,
                "%d,%s,%s,%.12g,%.12g,%.12g,%.12g,%.12g,%s,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                index,
                mode.kind,
                mode.direction,
                real(mode.wavenumber_per_m),
                imag(mode.wavenumber_per_m),
                real(mode.gamma_per_m),
                imag(mode.gamma_per_m),
                mode.power_w_per_m,
                mode.parity,
                mode.parity_score,
                mode.p_fraction,
                mode.s_fraction,
                mode.axial_displacement_fraction,
                mode.relative_residual,
            )
        end
    end
end

function plot_propagating_modes(path, modes, input_mode)
    right_going = propagating_modes(modes; direction=:right)
    panels = Any[]
    for mode in right_going
        ux, uy = magnitude_profiles(mode)
        selected = mode === input_mode ? " [input m=0]" : ""
        title = @sprintf(
            "k=%.1f  %s%s\nP=%.3f, axial=%.3f",
            real(mode.wavenumber_per_m),
            mode.parity,
            selected,
            mode.p_fraction,
            mode.axial_displacement_fraction,
        )
        push!(panels, plot(
            mode.y_m .* 1e3,
            ux;
            label="|ux|",
            xlabel="y, mm",
            ylabel="normalized displacement",
            title,
            linewidth=2,
        ))
        plot!(panels[end], mode.y_m .* 1e3, uy; label="|uy|", linewidth=2)
    end
    isempty(panels) && return nothing
    figure = plot(
        panels...;
        layout=(length(panels), 1),
        size=(900, 300length(panels)),
        plot_title=@sprintf("Right-going strip modes at %.1f kHz", FREQUENCY_HZ / 1e3),
    )
    savefig(figure, path)
    path
end

function write_branch_summary(path, frequencies_hz, branch)
    open(path, "w") do io
        println(io, "frequency_hz,k_per_m,phase_velocity_m_per_s,p_fraction,s_fraction,axial_fraction,overlap")
        for (frequency_hz, mode, overlap) in zip(frequencies_hz, branch.modes, branch.overlaps)
            @printf(
                io,
                "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                frequency_hz,
                real(mode.wavenumber_per_m),
                2pi * frequency_hz / real(mode.wavenumber_per_m),
                mode.p_fraction,
                mode.s_fraction,
                mode.axial_displacement_fraction,
                overlap,
            )
        end
    end
end

function plot_branch(path, frequencies_hz, branch)
    frequency_khz = frequencies_hz ./ 1e3
    wavenumber = real.(getproperty.(branch.modes, :wavenumber_per_m))
    top = plot(
        frequency_khz,
        wavenumber;
        label="tracked m=0",
        xlabel="frequency, kHz",
        ylabel="k, 1/m",
        linewidth=2,
        marker=:circle,
        markersize=2,
        title="Fundamental symmetric port branch",
    )
    bottom = plot(
        frequency_khz,
        getproperty.(branch.modes, :p_fraction);
        label="P diagnostic",
        xlabel="frequency, kHz",
        ylabel="fraction",
        linewidth=2,
        ylim=(0, 1),
    )
    plot!(
        bottom,
        frequency_khz,
        getproperty.(branch.modes, :axial_displacement_fraction);
        label="axial displacement",
        linewidth=2,
    )
    figure = plot(top, bottom; layout=(2, 1), size=(900, 700))
    savefig(figure, path)
    path
end

function run()
    mkpath(OUTPUT_ROOT)
    config = PortModeConfig(frequency_hz=FREQUENCY_HZ, element_count=ELEMENT_COUNT)
    modes = solve_port_modes(config)
    right_going = propagating_modes(modes; direction=:right)
    input_mode = fundamental_quasi_longitudinal(modes)
    input_index = findfirst(mode -> mode === input_mode, modes)
    result_path = joinpath(OUTPUT_ROOT, "port_modes.jld2")
    summary_path = joinpath(OUTPUT_ROOT, "port_modes.csv")
    figure_path = joinpath(OUTPUT_ROOT, "right_going_modes.png")
    continuation_frequencies = collect(10.0e3:CONTINUATION_STEP_HZ:FREQUENCY_HZ)
    last(continuation_frequencies) == FREQUENCY_HZ || push!(continuation_frequencies, FREQUENCY_HZ)
    continuation_config = PortModeConfig(
        frequency_hz=first(continuation_frequencies),
        element_count=min(ELEMENT_COUNT, 50),
    )
    branch = track_fundamental_branch(continuation_config, continuation_frequencies)
    branch_summary_path = joinpath(OUTPUT_ROOT, "fundamental_branch.csv")
    branch_figure_path = joinpath(OUTPUT_ROOT, "fundamental_branch.png")
    @save result_path config modes input_index
    write_summary(summary_path, modes)
    plot_propagating_modes(figure_path, modes, input_mode)
    write_branch_summary(branch_summary_path, continuation_frequencies, branch)
    plot_branch(branch_figure_path, continuation_frequencies, branch)

    println("[+] propagating modes: $(length(right_going)) right-going / $(length(propagating_modes(modes))) total")
    @printf(
        "[+] input m=0: k=%.6g 1/m, parity=%s, P-content=%.4f, axial=%.4f\n",
        real(input_mode.wavenumber_per_m),
        input_mode.parity,
        input_mode.p_fraction,
        input_mode.axial_displacement_fraction,
    )
    println("[+] $summary_path")
    println("[+] $figure_path")
    println("[+] $branch_figure_path")
    (; config, modes, input_index, result_path, summary_path, figure_path,
       branch_summary_path, branch_figure_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
