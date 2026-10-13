module TopologyHybridMonostaticAnalysis

using FFTW
using JLD2
using LinearAlgebra

include(joinpath(@__DIR__, "ndt_benchmark.jl"))
using .NDTBenchmark
include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
using .TopologyHybridTTDGRINDesign

export AnalysisConfig,
       matched_filter_trace,
       analyze_signal_pair,
       analyze_benchmark

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RAW_ROOT = get(
    ENV,
    "METAMATERIALS_HYBRID_MONOSTATIC_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_monostatic_g3"),
)
const RESULT_ROOT = joinpath(
    PROJECT_ROOT, "results", "topology_benchmark", "hybrid_monostatic_g3",
)
const TOPOLOGIES = (:bare, :straight_footprint, :pure_ttd, :hybrid)

Base.@kwdef struct AnalysisConfig
    frequency_hz::Float64 = 242.0e3
    pulse_cycles::Float64 = 5.0
    bare_echo_window_us::Tuple{Float64, Float64} = (24.0, 40.0)
    device_echo_window_us::Tuple{Float64, Float64} = (45.0, 68.0)
    dead_zone_threshold_db::Float64 = -40.0
    dead_zone_hold_cycles::Float64 = 1.0
    pressure_wave_speed_m_s::Float64 = 6122.102437409232
    calibrated_guide_speed_m_s::Float64 = 4793.0
end

signal_path(level, topology, case_id) = joinpath(
    RAW_ROOT, "signals", "signal_$(level)_$(topology)_$(case_id).jld2",
)

function load_signal(level, topology, case_id)
    path = signal_path(level, topology, case_id)
    isfile(path) || error("signal missing: $path")
    names = (
        "time_s",
        "source_normal_velocity_m_s",
        "source_tangent_velocity_m_s",
        "source_drive_pa",
        "incident_energy_j_per_m",
        "dt_s",
        "elapsed_s",
        "dof_count",
        "source_length_m",
    )
    data = JLD2.load(path, names...)
    NamedTuple{Symbol.(names)}(Tuple(data))
end

function pulse_template(dt_s, config::AnalysisConfig)
    duration = config.pulse_cycles / config.frequency_hz
    time = collect(0.0:dt_s:duration)
    values = [
        t < duration ?
        0.5 * (1 - cospi(2t / duration)) * sinpi(2config.frequency_hz * t) : 0.0
        for t in time
    ]
    (; time, values, duration)
end

"""Causal real matched filter; returned times refer to the template centre."""
function matched_filter_trace(time, signal, template)
    length(time) == length(signal) || throw(DimensionMismatch("signal time mismatch"))
    dt = time[2] - time[1]
    count = length(signal) - length(template.values) + 1
    count > 0 || throw(DimensionMismatch("template longer than signal"))
    normalization = sqrt(sum(abs2, template.values) * dt)
    values = zeros(Float64, count)
    for start in 1:count
        values[start] = dot(
            @view(signal[start:(start + length(template.values) - 1)]),
            template.values,
        ) * dt / normalization
    end
    centre_time = time[1:count] .+ template.duration / 2
    (; time_s=centre_time, values)
end

function analytic_envelope(signal)
    count = length(signal)
    spectrum = fft(signal)
    multiplier = zeros(Float64, count)
    multiplier[1] = 1.0
    if iseven(count)
        multiplier[2:(count ÷ 2)] .= 2.0
        multiplier[count ÷ 2 + 1] = 1.0
    else
        multiplier[2:((count + 1) ÷ 2)] .= 2.0
    end
    abs.(ifft(spectrum .* multiplier))
end

function echo_window(topology, config)
    topology == :bare ? config.bare_echo_window_us : config.device_echo_window_us
end

function device_two_way_offset_s(topology, config)
    topology == :bare && return 0.0
    straight = 55e-3 / config.calibrated_guide_speed_m_s
    topology == :straight_footprint && return 2straight
    common = 1.5e-3 / config.calibrated_guide_speed_m_s
    design = HybridTTDGRINConfig()
    relative = topology == :pure_ttd ? maximum_delay_s(design) : 1.56771e-6
    2(straight + common + relative)
end

function dead_zone_depth_mm(time, baseline, topology, config)
    envelope = analytic_envelope(baseline)
    drive_end = config.pulse_cycles / config.frequency_hz
    drive_mask = time .<= drive_end
    threshold = maximum(envelope[drive_mask]) * 10.0^(config.dead_zone_threshold_db / 20)
    hold_count = max(1, round(Int, config.dead_zone_hold_cycles /
                                  (config.frequency_hz * (time[2] - time[1]))))
    first_search = searchsortedfirst(time, drive_end)
    clearing_index = length(time)
    for index in first_search:(length(time) - hold_count + 1)
        if maximum(@view envelope[index:(index + hold_count - 1)]) <= threshold
            clearing_index = index
            break
        end
    end
    effective_time = max(0.0, time[clearing_index] - device_two_way_offset_s(topology, config))
    0.5 * config.pressure_wave_speed_m_s * effective_time * 1e3
end

function window_peak(trace, window_us)
    selected = (trace.time_s .>= window_us[1] * 1e-6) .&
               (trace.time_s .<= window_us[2] * 1e-6)
    any(selected) || error("matched-filter window contains no samples")
    maximum(abs, trace.values[selected])
end

function window_energy(time, signal, window_us)
    selected = (time .>= window_us[1] * 1e-6) .&
               (time .<= window_us[2] * 1e-6)
    any(selected) || error("energy window contains no samples")
    sum(abs2, signal[selected]) * (time[2] - time[1])
end

function analyze_signal_pair(level, topology; config=AnalysisConfig())
    baseline = load_signal(level, topology, :none)
    defect = load_signal(level, topology, :void_r2_focus)
    baseline.time_s == defect.time_s || error("paired time grids differ")
    delta = defect.source_normal_velocity_m_s - baseline.source_normal_velocity_m_s
    template = pulse_template(baseline.dt_s, config)
    echo_trace = matched_filter_trace(baseline.time_s, delta, template)
    ringing_trace = matched_filter_trace(
        baseline.time_s, baseline.source_normal_velocity_m_s, template,
    )
    window = echo_window(topology, config)
    echo_mf = window_peak(echo_trace, window)
    ringing_mf = window_peak(ringing_trace, window)
    echo_energy = window_energy(baseline.time_s, delta, window)
    ringing_energy = window_energy(
        baseline.time_s, baseline.source_normal_velocity_m_s, window,
    )
    normal_energy = sum(abs2, baseline.source_normal_velocity_m_s) * baseline.dt_s
    tangent_energy = sum(abs2, baseline.source_tangent_velocity_m_s) * baseline.dt_s
    peak_index = argmax(abs.(echo_trace.values) .* (
        (echo_trace.time_s .>= window[1] * 1e-6) .&
        (echo_trace.time_s .<= window[2] * 1e-6)
    ))
    (
        level,
        topology,
        echo_window_start_us=window[1],
        echo_window_stop_us=window[2],
        echo_mf,
        ringing_mf,
        echo_to_ringing_db=20log10(echo_mf / ringing_mf),
        echo_energy,
        ringing_energy,
        energy_echo_to_ringing_db=10log10(echo_energy / ringing_energy),
        matched_filter_peak_time_us=echo_trace.time_s[peak_index] * 1e6,
        dead_zone_depth_mm=dead_zone_depth_mm(
            baseline.time_s, baseline.source_normal_velocity_m_s, topology, config,
        ),
        receiver_tangent_fraction=tangent_energy / (normal_energy + tangent_energy),
        input_energy_j_per_m=baseline.incident_energy_j_per_m,
        dt_us=baseline.dt_s * 1e6,
        dof_none=baseline.dof_count,
        dof_void=defect.dof_count,
        solve_elapsed_none_s=baseline.elapsed_s,
        solve_elapsed_void_s=defect.elapsed_s,
        source_length_mm=baseline.source_length_m * 1e3,
    )
end

function csv_value(value)
    value isa Symbol && return String(value)
    value isa AbstractString && return value
    value isa Bool && return lowercase(string(value))
    string(value)
end

function write_rows(path, rows)
    mkpath(dirname(path))
    isempty(rows) && error("cannot write empty table")
    names = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(names), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, name)) for name in names), ','))
        end
    end
end

function old_isolated_metrics()
    path = joinpath(
        PROJECT_ROOT, "results", "topology_benchmark", "hybrid_ttd_grin",
        "isolated_summary.csv",
    )
    lines = readlines(path)
    header = Symbol.(split(first(lines), ','))
    Dict(begin
        values = split(line, ',')
        row = Dict(header[index] => values[index] for index in eachindex(header))
        Symbol(row[:variant]) => row
    end for line in Iterators.drop(lines, 1))
end

parse_metric(row, name) = parse(Float64, row[name])

function verdict_from_results(coarse, refined)
    coarse_by = Dict(row.topology => row for row in coarse)
    coarse_advantage = 20log10(coarse_by[:hybrid].echo_mf / coarse_by[:pure_ttd].echo_mf)
    if isnothing(refined)
        return (
            verdict=:conditional,
            reason="refinement_missing",
            coarse_advantage_db=coarse_advantage,
            refined_advantage_db=NaN,
            sign_stable=false,
        )
    end
    refined_by = Dict(row.topology => row for row in refined)
    refined_advantage = 20log10(refined_by[:hybrid].echo_mf / refined_by[:pure_ttd].echo_mf)
    sign_stable = coarse_advantage * refined_advantage > 0
    hybrid = refined_by[:hybrid]
    pure = refined_by[:pure_ttd]
    no_dead_zone_worsening = hybrid.dead_zone_depth_mm <= pure.dead_zone_depth_mm + 1e-9
    no_substantial_ringing_worsening = hybrid.echo_to_ringing_db >= pure.echo_to_ringing_db - 3.0
    if refined_advantage >= 2.0 && hybrid.echo_to_ringing_db >= 0 &&
       no_dead_zone_worsening && sign_stable
        verdict, reason = :continue, "hybrid_exceeds_pure_by_2db"
    elseif 0 <= refined_advantage < 2.0 && hybrid.echo_to_ringing_db >= 0 &&
           hybrid.dead_zone_depth_mm <= pure.dead_zone_depth_mm + 5.0 &&
           no_substantial_ringing_worsening && sign_stable
        verdict, reason = :conditional, "sub_2db_gain_with_40pct_path_reduction"
    else
        verdict, reason = :stop, "exact_g3_decision_gate_failed"
    end
    (
        verdict,
        reason,
        coarse_advantage_db=coarse_advantage,
        refined_advantage_db=refined_advantage,
        sign_stable,
    )
end

function write_common_metrics(coarse, verdict)
    bare = only(filter(row -> row.topology == :bare, coarse))
    isolated = old_isolated_metrics()
    rows = NamedTuple[]
    for result in coarse
        variant = result.topology == :hybrid ? :hybrid :
                  result.topology == :pure_ttd ? :pure_same_depth : :straight
        g1 = isolated[variant]
        echo_gain = result.echo_mf / bare.echo_mf
        energy_gain = result.echo_energy / bare.echo_energy
        push!(rows, (
            topology=result.topology,
            case_id=:void_r2_focus,
            reference=:bare,
            frequency_hz=242000.0,
            input_energy_j_per_m=result.input_energy_j_per_m,
            echo_mf_gain_bare=echo_gain,
            echo_energy_gain_bare=energy_gain,
            echo_to_ringing_db=result.echo_to_ringing_db,
            equal_noise_snr_proxy_db=20log10(echo_gain),
            dead_zone_depth_mm=result.dead_zone_depth_mm,
            focus_peak_gain_bare=NaN,
            focus_energy_gain_bare=NaN,
            fwhm_mm=NaN,
            dof_minus3db_mm=NaN,
            sidelobe_ratio=NaN,
            pulse_broadening_ratio=parse_metric(g1, :pulse_broadening_ratio),
            pulse_correlation=parse_metric(g1, :pulse_correlation),
            postcursor_ratio=parse_metric(g1, :postcursor_ratio),
            weighted_transmission=result.topology == :bare ? 1.0 :
                parse_metric(g1, :weighted_transmission),
            weighted_reflection=result.topology == :bare ? 0.0 :
                parse_metric(g1, :weighted_reflection),
            weighted_mode_conversion=result.topology == :bare ? 0.0 :
                parse_metric(g1, :weighted_mode_conversion),
            minimum_feature_mm=result.topology == :hybrid ?
                design_summary().minimum_insert_width_mm :
                result.topology == :pure_ttd ? 3.2 : NaN,
            verdict=result.topology == :hybrid ? verdict.verdict : :reference,
        ))
    end
    write_rows(joinpath(RESULT_ROOT, "common_metrics.csv"), rows)
end

function analyze_benchmark()
    coarse = [analyze_signal_pair(:coarse, topology) for topology in TOPOLOGIES]
    refined_available = all(
        isfile(signal_path(:refined, topology, case_id))
        for topology in (:pure_ttd, :hybrid), case_id in (:none, :void_r2_focus)
    )
    refined = refined_available ?
        [analyze_signal_pair(:refined, topology) for topology in (:pure_ttd, :hybrid)] : nothing
    verdict = verdict_from_results(coarse, refined)
    write_rows(joinpath(RESULT_ROOT, "g3_metrics.csv"), coarse)
    if !isnothing(refined)
        write_rows(joinpath(RESULT_ROOT, "convergence_metrics.csv"), vcat(
            filter(row -> row.topology in (:pure_ttd, :hybrid), coarse), refined,
        ))
    end
    write_rows(joinpath(RESULT_ROOT, "decision.csv"), [verdict])
    write_common_metrics(coarse, verdict)
    println("[+] verdict=$(verdict.verdict), coarse hybrid-pure=$(round(verdict.coarse_advantage_db;digits=3)) dB, refined=$(round(verdict.refined_advantage_db;digits=3)) dB")
    (; coarse, refined, verdict)
end

if abspath(PROGRAM_FILE) == @__FILE__
    analyze_benchmark()
end

end
