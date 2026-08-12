module AnalyseAluminiumHornTTDTemporalCheck

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_AL_HORN_TTD_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_pilot_242khz"),
)

ENV["GKSwstype"] = "100"
using JLD2
using Plots
if !isdefined(Main, :AluminiumHornTTDPilot)
    Base.include(Main, joinpath(@__DIR__, "run_aluminium_horn_ttd_pilot.jl"))
end
include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
include(joinpath(@__DIR__, "spectral_analysis.jl"))
using .ImpulseRiskAnalysis
using .SpectralAnalysis

function probe_signal(data, name, component=:x)
    index = findfirst(==(name), data["probe_names"])
    isnothing(index) && error("probe $name is absent")
    key = component == :x ? "probe_velocity_x_m_per_s" : "probe_velocity_y_m_per_s"
    vec(data[key][index, :])
end

function check_row(label, maximum_data, equal_data)
    time_s = maximum_data["time_s"]
    time_s == equal_data["time_s"] || error("maximum/equal time grids differ")
    maximum_signal = probe_signal(maximum_data, "target_axis")
    equal_signal = probe_signal(equal_data, "target_axis")
    transverse_signal = probe_signal(maximum_data, "target_axis", :y)
    dt_s = time_s[2] - time_s[1]
    pulse = pulse_metrics(maximum_signal, equal_signal, equal_signal, dt_s)
    maximum_envelope = analytic_envelope(maximum_signal)
    equal_envelope = analytic_envelope(equal_signal)
    tail_mask = time_s .>= max(time_s[end] - 10e-6, time_s[1])
    (
        resolution=label,
        final_time_us=time_s[end] * 1e6,
        samples_per_period=round(Int, 1 / (242e3 * dt_s)),
        peak_over_equal_path=pulse.gain_peak,
        envelope_peak_over_equal_path=maximum(maximum_envelope) / maximum(equal_envelope),
        energy_over_equal_path=sum(abs2, maximum_signal) / sum(abs2, equal_signal),
        broadening_ratio=pulse.broadening_ratio,
        pulse_correlation=pulse.pulse_correlation,
        postcursor_ratio=pulse.postcursor_ratio,
        transverse_energy_ratio=sum(abs2, transverse_signal) / sum(abs2, maximum_signal),
        maximum_minus_equal_envelope_delay_us=envelope_delay(
            time_s, equal_signal, maximum_signal,
        ) * 1e6,
        maximum_tail_envelope_over_peak=maximum(maximum_envelope[tail_mask]) /
                                        maximum(maximum_envelope),
        equal_tail_envelope_over_peak=maximum(equal_envelope[tail_mask]) /
                                      maximum(equal_envelope),
    )
end

function write_rows(path, rows)
    fields = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(fields), ','))
        for row in rows
            println(io, join((getproperty(row, field) for field in fields), ','))
        end
    end
end

function main()
    signal_root = joinpath(OUTPUT_ROOT, "signals")
    coarse_maximum = JLD2.load(joinpath(signal_root, "maximum.jld2"))
    coarse_equal = JLD2.load(joinpath(signal_root, "equal_path.jld2"))
    refined_maximum = JLD2.load(joinpath(signal_root, "maximum_time_refined.jld2"))
    refined_equal = JLD2.load(joinpath(signal_root, "equal_path_time_refined.jld2"))
    rows = [
        check_row("100us_30spp", coarse_maximum, coarse_equal),
        check_row("130us_40spp", refined_maximum, refined_equal),
    ]
    summary_path = joinpath(OUTPUT_ROOT, "aluminium_horn_ttd_temporal_check.csv")
    write_rows(summary_path, rows)

    refined = rows[2]
    passed = refined.peak_over_equal_path >= 0.8 &&
             refined.broadening_ratio <= 1.25 &&
             refined.pulse_correlation >= 0.90 &&
             refined.postcursor_ratio <= 0.10 &&
             refined.transverse_energy_ratio <= 0.10
    relative_peak_change = refined.peak_over_equal_path / rows[1].peak_over_equal_path - 1

    panels = Any[]
    for (title, maximum_data, equal_data) in (
        ("100 μs, 30 samples/period", coarse_maximum, coarse_equal),
        ("130 μs, 40 samples/period", refined_maximum, refined_equal),
    )
        time_us = maximum_data["time_s"] .* 1e6
        maximum_envelope = analytic_envelope(probe_signal(maximum_data, "target_axis"))
        equal_envelope = analytic_envelope(probe_signal(equal_data, "target_axis"))
        panel = plot(
            time_us, equal_envelope;
            label="straight equal path", color=:seagreen, linewidth=2,
            xlabel="time, μs", ylabel="|analytic vₓ|, m/s", title,
            gridalpha=0.25,
        )
        plot!(panel, time_us, maximum_envelope;
              label="maximum sin⁴", color=:royalblue, linewidth=2.4)
        push!(panels, panel)
    end
    figure_path = joinpath(OUTPUT_ROOT, "aluminium_horn_ttd_temporal_check.png")
    savefig(plot(panels...; layout=(2, 1), size=(900, 760)), figure_path)

    verdict_path = joinpath(OUTPUT_ROOT, "temporal_check_verdict.txt")
    open(verdict_path, "w") do io
        println(io, passed ? "PASS" : "STOP")
        println(io, "relative_peak_gate_change=$(relative_peak_change)")
        for field in propertynames(refined)
            println(io, "$(field)=$(getproperty(refined, field))")
        end
    end
    println("[+] refined peak/equal=$(refined.peak_over_equal_path), " *
            "change=$(100relative_peak_change)%")
    println("[+] refined Bt=$(refined.broadening_ratio), " *
            "rho=$(refined.pulse_correlation), post=$(refined.postcursor_ratio), " *
            "transverse=$(refined.transverse_energy_ratio), passed=$passed")
    println("[+] $summary_path")
    println("[+] $figure_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
