module TopologyHybridRobustnessSweep

using FFTW
using JLD2
using LinearAlgebra
using Statistics
ENV["GKSwstype"] = "100"
using Plots

include(joinpath(@__DIR__, "topology_hybrid_monostatic_analysis.jl"))
using .TopologyHybridMonostaticAnalysis
include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
using .TopologyHybridTTDGRINDesign

export RobustnessPoint,
       aperture_error_transfer,
       apply_robustness_transfer,
       evaluate_point,
       run_robustness_sweep

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RAW_ROOT = get(
    ENV,
    "METAMATERIALS_HYBRID_MONOSTATIC_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_monostatic_g3"),
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HYBRID_ROBUSTNESS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_robustness_sweep"),
)
const RESULT_ROOT = joinpath(
    PROJECT_ROOT, "results", "topology_benchmark", "hybrid_robustness_sweep",
)

const NOMINAL_MODE_CONVERSION = 0.036018748758661887
const FREQUENCY_LOW_HZ = 162.2e3
const FREQUENCY_HIGH_HZ = 321.7e3

Base.@kwdef struct RobustnessPoint
    one_way_excess_loss_db::Float64 = 0.0
    maximum_profile_delay_error_ns::Float64 = 0.0
    differential_phase_ripple_deg::Float64 = 0.0
    phase_ripple_phase_deg::Float64 = 0.0
    absolute_mode_conversion_fraction::Float64 = NOMINAL_MODE_CONVERSION
    extra_ringing_penalty_db::Float64 = 0.0
end

function validate(point::RobustnessPoint)
    point.one_way_excess_loss_db >= 0 || throw(ArgumentError("loss must be non-negative"))
    point.differential_phase_ripple_deg >= 0 ||
        throw(ArgumentError("phase ripple must be non-negative"))
    0 <= point.absolute_mode_conversion_fraction < 1 ||
        throw(ArgumentError("mode conversion must lie in [0,1)"))
    point.extra_ringing_penalty_db >= 0 ||
        throw(ArgumentError("ringing penalty must be non-negative"))
    point
end

function signed_fft_frequencies(count::Integer, dt_s::Real)
    scale = 1 / (count * Float64(dt_s))
    [
        (index <= fld(count, 2) ? index : index - count) * scale
        for index in 0:(count - 1)
    ]
end

function radial_delay_fractions(design::HybridTTDGRINConfig=HybridTTDGRINConfig())
    maximum_delay = maximum_delay_s(design)
    [geometric_delay_s(y, design) / maximum_delay for y in channel_centers_mm(design)]
end

function phase_ripple_rad(frequency_hz, peak_deg, phase_offset_deg)
    frequency_hz == 0 && return 0.0
    span = FREQUENCY_HIGH_HZ - FREQUENCY_LOW_HZ
    phase = 2pi * (frequency_hz - FREQUENCY_LOW_HZ) / span
    deg2rad(peak_deg) * sin(phase + deg2rad(phase_offset_deg))
end

"""
Additional two-way useful-mode transfer applied to the exact refined hybrid
defect trace. Loss is specified per one-way traversal. Delay error and phase
ripple are the edge-to-centre differential errors; equal channel weights are
used for the deliberately conservative reduced aperture model. Conversion is
absolute useful-to-unwanted power, normalised to the already-modelled nominal
3.602% value so that RobustnessPoint() is exactly the identity.
"""
function aperture_error_transfer(
    frequencies_hz::AbstractVector{<:Real},
    point::RobustnessPoint;
    fractions=radial_delay_fractions(),
)
    validate(point)
    useful_amplitude_ratio = sqrt(
        (1 - point.absolute_mode_conversion_fraction) /
        (1 - NOMINAL_MODE_CONVERSION),
    )
    one_way_amplitude = 10.0^(-point.one_way_excess_loss_db / 20) *
                        useful_amplitude_ratio
    delay_s = point.maximum_profile_delay_error_ns * 1e-9
    ComplexF64[map(frequencies_hz) do frequency
        absolute_frequency = abs(Float64(frequency))
            ripple = phase_ripple_rad(
                absolute_frequency,
                point.differential_phase_ripple_deg,
                point.phase_ripple_phase_deg,
        )
        coherence = mean(cis.(-(
            2pi * absolute_frequency * delay_s .* fractions .+
            ripple .* fractions
        )))
        one_way = one_way_amplitude * coherence
        value = one_way^2
        frequency < 0 ? conj(value) : value
    end...]
end

function apply_robustness_transfer(
    signal::AbstractVector{<:Real},
    dt_s::Real,
    point::RobustnessPoint,
)
    count = length(signal)
    padded_count = nextpow(2, 2count)
    padded = zeros(Float64, padded_count)
    padded[1:count] .= signal
    frequencies = signed_fft_frequencies(padded_count, dt_s)
    transfer = aperture_error_transfer(frequencies, point)
    filtered = real.(ifft(fft(padded) .* transfer))
    filtered[1:count]
end

function signal_data(topology, case_id)
    path = joinpath(
        RAW_ROOT, "signals", "signal_refined_$(topology)_$(case_id).jld2",
    )
    isfile(path) || error("missing refined signal: $path")
    names = ("time_s", "source_normal_velocity_m_s", "dt_s")
    Dict(name => value for (name, value) in zip(names, JLD2.load(path, names...)))
end

function duration_90(time, signal, window_us)
    selected = findall(
        (time .>= window_us[1] * 1e-6) .& (time .<= window_us[2] * 1e-6),
    )
    isempty(selected) && error("empty duration window")
    envelope_energy = TopologyHybridMonostaticAnalysis.analytic_envelope(signal[selected]).^2
    total = sum(envelope_energy)
    total > 0 || return 0.0
    target = 0.9total
    best = length(selected)
    right = 1
    accumulated = 0.0
    for left in eachindex(selected)
        while right <= length(selected) && accumulated < target
            accumulated += envelope_energy[right]
            right += 1
        end
        accumulated >= target && (best = min(best, right - left))
        accumulated -= envelope_energy[left]
    end
    best * (time[2] - time[1])
end

function aligned_correlation(time, candidate, reference, window_us, frequency_hz=242e3)
    indices = findall(
        (time .>= window_us[1] * 1e-6) .& (time .<= window_us[2] * 1e-6),
    )
    isempty(indices) && error("empty correlation window")
    maximum_lag = round(Int, 2 / (frequency_hz * (time[2] - time[1])))
    best = -1.0
    for lag in (-maximum_lag):maximum_lag
        shifted = indices .+ lag
        valid = (shifted .>= 1) .& (shifted .<= length(time))
        count(valid) > 8 || continue
        a = candidate[shifted[valid]]
        b = reference[indices[valid]]
        denominator = norm(a) * norm(b)
        denominator > 0 || continue
        best = max(best, abs(dot(a, b)) / denominator)
    end
    best
end

function refined_context(; config=AnalysisConfig())
    hybrid_none = signal_data(:hybrid, :none)
    hybrid_void = signal_data(:hybrid, :void_r2_focus)
    pure = analyze_signal_pair(:refined, :pure_ttd; config)
    hybrid = analyze_signal_pair(:refined, :hybrid; config)
    time = hybrid_none["time_s"]
    dt_s = hybrid_none["dt_s"]
    delta = hybrid_void["source_normal_velocity_m_s"] .-
            hybrid_none["source_normal_velocity_m_s"]
    template = TopologyHybridMonostaticAnalysis.pulse_template(dt_s, config)
    window = TopologyHybridMonostaticAnalysis.echo_window(:hybrid, config)
    nominal_duration = duration_90(time, delta, window)
    (; time, dt_s, delta, template, window, pure, hybrid, nominal_duration)
end

function evaluate_point(point::RobustnessPoint, context=refined_context())
    validate(point)
    filtered = apply_robustness_transfer(context.delta, context.dt_s, point)
    trace = matched_filter_trace(context.time, filtered, context.template)
    echo_mf = TopologyHybridMonostaticAnalysis.window_peak(trace, context.window)
    echo_energy = TopologyHybridMonostaticAnalysis.window_energy(
        context.time, filtered, context.window,
    )
    ringing_mf = context.hybrid.ringing_mf *
                 10.0^(point.extra_ringing_penalty_db / 20)
    broadening = duration_90(context.time, filtered, context.window) /
                 context.nominal_duration
    correlation = aligned_correlation(
        context.time, filtered, context.delta, context.window,
    )
    advantage_db = 20log10(echo_mf / context.pure.echo_mf)
    echo_to_ringing_db = 20log10(echo_mf / ringing_mf)
    (
        one_way_excess_loss_db=point.one_way_excess_loss_db,
        maximum_profile_delay_error_ns=point.maximum_profile_delay_error_ns,
        differential_phase_ripple_deg=point.differential_phase_ripple_deg,
        phase_ripple_phase_deg=point.phase_ripple_phase_deg,
        absolute_mode_conversion_fraction=point.absolute_mode_conversion_fraction,
        extra_ringing_penalty_db=point.extra_ringing_penalty_db,
        echo_mf,
        echo_energy,
        hybrid_to_pure_advantage_db=advantage_db,
        echo_to_ringing_db,
        pulse_broadening_ratio=broadening,
        pulse_correlation=correlation,
        conditional_gate_pass=advantage_db >= 0 && echo_to_ringing_db >= 0 &&
                              broadening <= 1.25 && correlation >= 0.90,
        continue_2db_gate_pass=advantage_db >= 2 && echo_to_ringing_db >= 0 &&
                               broadening <= 1.25 && correlation >= 0.90,
    )
end

csv_value(value::Bool) = lowercase(string(value))
csv_value(value) = string(value)

function write_rows(path, rows)
    isempty(rows) && error("cannot write empty rows")
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, column)) for column in columns), ','))
        end
    end
    path
end

function last_passing_limit(rows, coordinate; transform=identity)
    ordered = sort(rows; by=row -> transform(getproperty(row, coordinate)))
    passed = filter(row -> row.conditional_gate_pass, ordered)
    isempty(passed) && return (limit=NaN, censored=false)
    value = maximum(transform(getproperty(row, coordinate)) for row in passed)
    maximum_tested = maximum(transform(getproperty(row, coordinate)) for row in ordered)
    (limit=value, censored=value == maximum_tested)
end

function symmetric_last_passing_limit(rows, coordinate)
    magnitudes = sort(unique(abs(getproperty(row, coordinate)) for row in rows))
    passing = Float64[]
    for magnitude in magnitudes
        matching = filter(row -> abs(getproperty(row, coordinate)) == magnitude, rows)
        length(matching) == (magnitude == 0 ? 1 : 2) || continue
        all(row -> row.conditional_gate_pass, matching) && push!(passing, magnitude)
    end
    isempty(passing) && return (limit=NaN, censored=false)
    limit = maximum(passing)
    (; limit, censored=limit == maximum(magnitudes))
end

worst_row(rows) = rows[argmin(getproperty.(rows, :hybrid_to_pure_advantage_db))]

function write_specification(axis_rows)
    loss = last_passing_limit(axis_rows[:loss], :one_way_excess_loss_db)
    delay = symmetric_last_passing_limit(axis_rows[:delay], :maximum_profile_delay_error_ns)
    ripple = last_passing_limit(axis_rows[:ripple], :differential_phase_ripple_deg)
    conversion = last_passing_limit(axis_rows[:conversion], :absolute_mode_conversion_fraction)
    ringing = last_passing_limit(axis_rows[:ringing], :extra_ringing_penalty_db)
    rows = [
        # The proposed gates are one combined qualification point, not five
        # independent 80%-of-break-even limits. Applying the independent
        # limits simultaneously would spend the same 0.97 dB margin repeatedly.
        (parameter="one_way_excess_loss_db", break_even_limit=loss.limit, censored=loss.censored, proposed_coupon_gate=0.10),
        (parameter="maximum_profile_delay_error_ns", break_even_limit=delay.limit, censored=delay.censored, proposed_coupon_gate=100.0),
        (parameter="differential_phase_ripple_deg", break_even_limit=ripple.limit, censored=ripple.censored, proposed_coupon_gate=5.0),
        (parameter="absolute_mode_conversion_fraction", break_even_limit=conversion.limit, censored=conversion.censored, proposed_coupon_gate=0.05),
        (parameter="extra_ringing_penalty_db", break_even_limit=ringing.limit, censored=ringing.censored, proposed_coupon_gate=6.0),
    ]
    write_rows(joinpath(RESULT_ROOT, "coupon_specification.csv"), rows)
    rows
end


function write_combined_qualification(context)
    cases = (
        (name="qualification_target", loss=0.10, delay=100.0, ripple=5.0, conversion=0.05, ringing=6.0),
        (name="relaxed_screen", loss=0.15, delay=100.0, ripple=10.0, conversion=0.06, ringing=6.0),
    )
    rows = NamedTuple[]
    for case in cases
        candidates = [
            evaluate_point(RobustnessPoint(
                one_way_excess_loss_db=case.loss,
                maximum_profile_delay_error_ns=sign * case.delay,
                differential_phase_ripple_deg=case.ripple,
                phase_ripple_phase_deg=phase,
                absolute_mode_conversion_fraction=case.conversion,
                extra_ringing_penalty_db=case.ringing,
            ), context)
            for sign in (-1.0, 1.0), phase in (0.0, 90.0, 180.0, 270.0)
        ]
        worst = worst_row(vec(candidates))
        push!(rows, merge((qualification=case.name,), worst))
    end
    write_rows(joinpath(RESULT_ROOT, "combined_qualification.csv"), rows)
    rows
end

function save_plots(loss_delay, loss_ripple, loss_conversion, loss_ringing)
    function matrix(rows, xname, yname, value)
        xs = sort(unique(getproperty(row, xname) for row in rows))
        ys = sort(unique(getproperty(row, yname) for row in rows))
        z = fill(NaN, length(ys), length(xs))
        xi, yi = Dict(x => i for (i, x) in enumerate(xs)), Dict(y => i for (i, y) in enumerate(ys))
        for row in rows
            z[yi[getproperty(row, yname)], xi[getproperty(row, xname)]] = getproperty(row, value)
        end
        xs, ys, z
    end
    panels = Any[]
    for (rows, xname, yname, ylabel) in (
        (loss_delay, :one_way_excess_loss_db, :maximum_profile_delay_error_ns, "profile delay error, ns"),
        (loss_ripple, :one_way_excess_loss_db, :differential_phase_ripple_deg, "phase ripple, deg"),
        (loss_conversion, :one_way_excess_loss_db, :absolute_mode_conversion_fraction, "absolute mode conversion"),
    )
        x, y, z = matrix(rows, xname, yname, :hybrid_to_pure_advantage_db)
        push!(panels, heatmap(
            x, y, z; xlabel="one-way excess loss, dB", ylabel,
            color=:balance, clims=(-2, 2), colorbar_title="hybrid-pure, dB",
            contour=false,
        ))
        contour!(panels[end], x, y, z; levels=[0.0], color=:black, linewidth=2, label=false)
    end
    x, y, z = matrix(
        loss_ringing, :one_way_excess_loss_db, :extra_ringing_penalty_db,
        :echo_to_ringing_db,
    )
    push!(panels, heatmap(
        x, y, z; xlabel="one-way excess loss, dB", ylabel="extra ringing penalty, dB",
        color=:viridis, colorbar_title="echo/ringing, dB",
    ))
    contour!(panels[end], x, y, z; levels=[0.0], color=:white, linewidth=2, label=false)
    figure = plot(
        panels...; layout=(2, 2), size=(1500, 1050), margin=4Plots.mm,
        plot_title="Hybrid monostatic robustness sweep (exact refined G3 trace + reduced error transfer)",
    )
    mkpath(OUTPUT_ROOT)
    path = joinpath(OUTPUT_ROOT, "hybrid_robustness_maps.png")
    savefig(figure, path)
    path
end

function run_robustness_sweep()
    mkpath(OUTPUT_ROOT)
    mkpath(RESULT_ROOT)
    context = refined_context()
    nominal = evaluate_point(RobustnessPoint(), context)

    axis_rows = Dict{Symbol, Vector{NamedTuple}}()
    axis_rows[:loss] = [evaluate_point(RobustnessPoint(one_way_excess_loss_db=value), context) for value in 0.0:0.005:0.70]
    axis_rows[:delay] = [evaluate_point(RobustnessPoint(maximum_profile_delay_error_ns=value), context) for value in -700.0:5.0:700.0]
    axis_rows[:ripple] = [worst_row([
        evaluate_point(RobustnessPoint(
            differential_phase_ripple_deg=value,
            phase_ripple_phase_deg=phase,
        ), context)
        for phase in (0.0, 90.0, 180.0, 270.0)
    ]) for value in 0.0:0.5:60.0]
    axis_rows[:conversion] = [evaluate_point(RobustnessPoint(absolute_mode_conversion_fraction=value), context) for value in NOMINAL_MODE_CONVERSION:0.001:0.16]
    axis_rows[:ringing] = [evaluate_point(RobustnessPoint(extra_ringing_penalty_db=value), context) for value in 0.0:0.05:18.0]
    for (name, rows) in axis_rows
        write_rows(joinpath(RESULT_ROOT, "axis_$(name).csv"), rows)
    end

    loss_values = 0.0:0.05:1.0
    loss_delay = [worst_row([
        evaluate_point(RobustnessPoint(
            one_way_excess_loss_db=loss,
            maximum_profile_delay_error_ns=sign * delay,
        ), context)
        for sign in (delay == 0 ? (1.0,) : (-1.0, 1.0))
    ]) for loss in loss_values for delay in 0.0:50.0:1000.0]
    loss_ripple = [worst_row([
        evaluate_point(RobustnessPoint(
            one_way_excess_loss_db=loss,
            differential_phase_ripple_deg=ripple,
            phase_ripple_phase_deg=phase,
        ), context)
        for phase in (0.0, 90.0, 180.0, 270.0)
    ]) for loss in loss_values for ripple in 0.0:5.0:90.0]
    loss_conversion = [evaluate_point(RobustnessPoint(one_way_excess_loss_db=loss, absolute_mode_conversion_fraction=conversion), context) for loss in loss_values for conversion in range(NOMINAL_MODE_CONVERSION, 0.25; length=22)]
    loss_ringing = [evaluate_point(RobustnessPoint(one_way_excess_loss_db=loss, extra_ringing_penalty_db=penalty), context) for loss in loss_values for penalty in 0.0:1.0:20.0]
    write_rows(joinpath(RESULT_ROOT, "surface_loss_delay.csv"), loss_delay)
    write_rows(joinpath(RESULT_ROOT, "surface_loss_ripple.csv"), loss_ripple)
    write_rows(joinpath(RESULT_ROOT, "surface_loss_conversion.csv"), loss_conversion)
    write_rows(joinpath(RESULT_ROOT, "surface_loss_ringing.csv"), loss_ringing)
    specification = write_specification(axis_rows)
    combined = write_combined_qualification(context)
    plot_path = save_plots(loss_delay, loss_ripple, loss_conversion, loss_ringing)
    println("[+] nominal hybrid-pure=$(round(nominal.hybrid_to_pure_advantage_db; digits=4)) dB")
    specification_path = joinpath(RESULT_ROOT, "coupon_specification.csv")
    println("[+] specification=$specification_path")
    println("[+] maps=$plot_path")
    (; nominal, axis_rows, specification, combined, plot_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_robustness_sweep()
end

end
