module AnalyzeSinusoidalWaveguideTemporalReversal

ENV["GKSwstype"] = "100"

using FFTW
using JLD2
using LinearAlgebra
using Plots
using Statistics

export analyze_sinusoidal_waveguide_temporal_reversal,
       analytic_envelope,
       temporal_refocus,
       temporal_focus_metrics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PACKAGE_ROOT = joinpath(
    PROJECT_ROOT, "results", "sinusoidal_material_transfer_p6_242khz",
)
const INPUT_PATH = joinpath(
    PACKAGE_ROOT, "broadband", "aluminium_broadband_impulse.jld2",
)
const OUTPUT_ROOT = joinpath(PACKAGE_ROOT, "time_reversal_broadband_scalar")
const FREQUENCY_HZ = 242.0e3
const PULSE_CYCLES = 5.0

function analytic_envelope(signal::AbstractVector{<:Real})
    count = length(signal)
    count >= 3 || throw(ArgumentError("signal is too short"))
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

function temporal_refocus(
    transfer::AbstractVector{<:Complex},
    pulse::AbstractVector{<:Real},
)
    pulse_spectrum = rfft(pulse)
    length(transfer) == length(pulse_spectrum) || throw(DimensionMismatch(
        "transfer must be defined on the pulse rFFT grid",
    ))
    received_spectrum = transfer .* pulse_spectrum
    reversed_record = irfft(conj.(received_spectrum), length(pulse))
    replay_norm = norm(reversed_record)
    replay_norm > 0 || error("received record vanished")
    refocused = irfft(
        transfer .* conj.(received_spectrum) ./ replay_norm,
        length(pulse),
    )
    envelope = analytic_envelope(refocused)
    centre = length(pulse) ÷ 2 + 1
    shift = centre - argmax(envelope)
    shifted = circshift(refocused, shift)
    (
        waveform=shifted,
        envelope=analytic_envelope(shifted),
        replay_norm_before_normalization=replay_norm,
        circular_shift_samples=shift,
    )
end

function temporal_focus_metrics(result, dt_s; focus_window_s=PULSE_CYCLES / FREQUENCY_HZ)
    dt_s > 0 || throw(ArgumentError("sample interval must be positive"))
    focus_window_s > 0 || throw(ArgumentError("focus window must be positive"))
    envelope = result.envelope
    waveform = result.waveform
    centre = length(envelope) ÷ 2 + 1
    peak = envelope[centre]
    peak > 0 || error("time-reversal focus vanished")
    threshold = peak / sqrt(2)
    left = centre
    right = centre
    while left > firstindex(envelope) && envelope[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(envelope) && envelope[right + 1] >= threshold
        right += 1
    end
    time_s = ((1:length(envelope)) .- centre) .* dt_s
    far = abs.(time_s) .> focus_window_s / 2
    total_energy = sum(abs2, waveform)
    (
        peak_envelope=peak,
        half_power_width_s=(right - left + 1) * dt_s,
        half_power_energy_fraction=sum(abs2, waveform[left:right]) / total_energy,
        ringing_energy_fraction=sum(abs2, waveform[far]) / total_energy,
        peak_far_envelope_db=20log10(maximum(envelope[far]) / peak),
        replay_norm_before_normalization=result.replay_norm_before_normalization,
        focus_window_s,
        time_s,
    )
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function make_plot(results, metrics)
    time_us = metrics[:sinusoidal].time_s .* 1.0e6
    visible = abs.(time_us) .<= 30.0
    envelope_plot = plot(
        xlabel="time from focus / microseconds", ylabel="normalized envelope",
        title="Equal-energy scalar time reversal", size=(960, 580),
        legend=:topright,
    )
    for (geometry, label, color) in (
        (:sinusoidal, "sinusoidal", :darkorange),
        (:straight, "straight", :steelblue),
    )
        normalized = results[geometry].envelope ./ maximum(results[geometry].envelope)
        plot!(envelope_plot, time_us[visible], normalized[visible];
            linewidth=2.5, label, color)
    end
    half_window_us = metrics[:sinusoidal].focus_window_s * 1.0e6 / 2
    vline!(envelope_plot, [-half_window_us, half_window_us];
        color=:black, linestyle=:dash, label="one input-burst window")
    savefig(envelope_plot, joinpath(OUTPUT_ROOT, "temporal_refocus_envelope.png"))

    comparison = bar(
        ["peak ratio", "FWHM ratio", "ringing ratio"],
        [
            metrics[:sinusoidal].peak_envelope / metrics[:straight].peak_envelope,
            metrics[:sinusoidal].half_power_width_s / metrics[:straight].half_power_width_s,
            metrics[:sinusoidal].ringing_energy_fraction /
                metrics[:straight].ringing_energy_fraction,
        ];
        ylabel="sinusoidal / straight", title="Broadband time-reversal comparison",
        color=[:darkorange, :steelblue, :mediumpurple], label=false,
        size=(820, 520),
    )
    hline!(comparison, [1.0]; color=:black, linestyle=:dash, label="equal")
    savefig(comparison, joinpath(OUTPUT_ROOT, "temporal_refocus_comparison.png"))
end

function write_readme(rows, comparison)
    values = Dict(row.geometry => row for row in rows)
    lens = values[:sinusoidal]
    straight = values[:straight]
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# Sinusoidal aluminium waveguides: scalar broadband time reversal")
        println(io)
        println(io, "## Result")
        println(io)
        println(io, "The saved 17-anchor full-FEM channel was time reversed with equal replay energy. The sinusoidal channel gives a `$(round(comparison.peak_ratio; digits=3))x` larger refocused envelope peak than the matched straight strips. Its half-power width is `$(round(1.0e6lens.half_power_width_s; digits=3)) microseconds` versus `$(round(1.0e6straight.half_power_width_s; digits=3)) microseconds`.")
        println(io)
        println(io, "| geometry | equal-energy peak / m | half-power width / us | energy in half-power lobe | ringing energy outside 5-cycle window | far peak / dB |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for row in rows
            println(io, "| $(row.geometry) | $(row.peak_envelope_m) | $(round(1.0e6row.half_power_width_s; digits=3)) | $(round(100row.half_power_energy_fraction; digits=2))% | $(round(100row.ringing_energy_fraction; digits=2))% | $(round(row.peak_far_envelope_db; digits=2)) |")
        end
        println(io)
        println(io, "The sinusoidal path leaves `$(round(100lens.ringing_energy_fraction; digits=2))%` of refocused energy outside the centred five-cycle window, compared with `$(round(100straight.ringing_energy_fraction; digits=2))%` for the straight control. Thus its advantage in this SISO experiment is temporal matched-filter concentration, even though the separate carrier `8 x 8` audit found worse spatial-channel conditioning.")
        println(io)
        println(io, "## Boundary")
        println(io)
        println(io, "This is a reciprocal scalar record-and-replay calculation from the common uniform input to the design focus. It neither demonstrates added spatial modes nor estimates modal Q. The 17-point response was already validated for the five-cycle waveform, but a broadband spatial time-reversal claim still requires the full response matrix at every frequency and odd as well as even source groups.")
    end
end

function analyze_sinusoidal_waveguide_temporal_reversal()
    isfile(INPUT_PATH) || error("missing frozen broadband response: $INPUT_PATH")
    mkpath(OUTPUT_ROOT)
    data = JLD2.load(INPUT_PATH)
    pulse = Float64.(data["pulse_input"])
    time_s = Float64.(data["pulse_time_s"])
    dt_s = time_s[2] - time_s[1]
    transfers = Dict(
        :sinusoidal => ComplexF64.(data["selected_transfer_m"]),
        :straight => ComplexF64.(data["uniform_transfer_m"]),
    )
    results = Dict(
        geometry => temporal_refocus(transfer, pulse)
        for (geometry, transfer) in transfers
    )
    metrics = Dict(
        geometry => temporal_focus_metrics(result, dt_s)
        for (geometry, result) in results
    )
    rows = [
        let value=metrics[geometry]
            (
                geometry,
                peak_envelope_m=value.peak_envelope,
                half_power_width_s=value.half_power_width_s,
                half_power_energy_fraction=value.half_power_energy_fraction,
                ringing_energy_fraction=value.ringing_energy_fraction,
                peak_far_envelope_db=value.peak_far_envelope_db,
                replay_norm_before_normalization=
                    value.replay_norm_before_normalization,
            )
        end
        for geometry in (:sinusoidal, :straight)
    ]
    comparison = (
        peak_ratio=metrics[:sinusoidal].peak_envelope /
            metrics[:straight].peak_envelope,
        half_power_width_ratio=metrics[:sinusoidal].half_power_width_s /
            metrics[:straight].half_power_width_s,
        ringing_energy_ratio=metrics[:sinusoidal].ringing_energy_fraction /
            metrics[:straight].ringing_energy_fraction,
        far_peak_change_db=metrics[:sinusoidal].peak_far_envelope_db -
            metrics[:straight].peak_far_envelope_db,
    )
    write_csv(joinpath(OUTPUT_ROOT, "temporal_refocus_summary.csv"), rows)
    write_csv(joinpath(OUTPUT_ROOT, "temporal_refocus_comparison.csv"), [comparison])
    visible = abs.(metrics[:sinusoidal].time_s) .<= 50.0e-6
    waveform_rows = [
        (
            time_s=metrics[:sinusoidal].time_s[index],
            sinusoidal_waveform_m=results[:sinusoidal].waveform[index],
            sinusoidal_envelope_m=results[:sinusoidal].envelope[index],
            straight_waveform_m=results[:straight].waveform[index],
            straight_envelope_m=results[:straight].envelope[index],
        )
        for index in findall(visible)
    ]
    write_csv(joinpath(OUTPUT_ROOT, "temporal_refocus_waveforms.csv"), waveform_rows)
    make_plot(results, metrics)
    write_readme(rows, comparison)
    jldsave(
        joinpath(OUTPUT_ROOT, "sinusoidal_temporal_reversal.jld2");
        frequency_hz=FREQUENCY_HZ, pulse_cycles=PULSE_CYCLES, dt_s,
        results, metrics, rows, comparison,
    )
    println("[+] Equal-energy peak ratio: $(comparison.peak_ratio)")
    println("[+] Half-power width ratio: $(comparison.half_power_width_ratio)")
    println("[+] Ringing-energy ratio: $(comparison.ringing_energy_ratio)")
    rows
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    analyze_sinusoidal_waveguide_temporal_reversal()
end

end
