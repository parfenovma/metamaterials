module TopologyHybridTTDGRINAperture

ENV["GKSwstype"] = "100"

using FFTW
using JLD2
using LinearAlgebra
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HYBRID_TTD_GRIN_BENCHMARK_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_ttd_grin", "aperture"),
)
const ISOLATED_ROOT = get(
    ENV,
    "METAMATERIALS_HYBRID_TTD_GRIN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_ttd_grin", "isolated"),
)

include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
using .TopologyHybridTTDGRINDesign
include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .ImpulseRiskAnalysis
include(joinpath(@__DIR__, "ndt_benchmark.jl"))
using .NDTBenchmark

const DESIGN = HybridTTDGRINConfig()
const STAND = NDTStandConfig()
const TOPOLOGIES = (:bare, :straight_footprint, :pure_ttd, :hybrid)

function unwrap_phase(values)
    result = Float64[first(values)]
    for value in Iterators.drop(values, 1)
        push!(result, last(result) + mod(value - last(result) + pi, 2pi) - pi)
    end
    result
end

function linear_interpolate(x, values, query)
    query <= first(x) && return first(values)
    query >= last(x) && return last(values)
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function load_isolated_transfer()
    path = joinpath(ISOLATED_ROOT, "isolated_analysis.jld2")
    isfile(path) || error("isolated G1 result is missing: $path")
    data = JLD2.load(path)
    rows = data["rows"]
    selected(variant) = sort(
        filter(row -> row.variant == variant, rows);
        by=row -> row.frequency_hz,
    )
    straight = selected(:straight)
    hybrid = selected(:hybrid)
    pure = selected(:pure_same_depth)
    frequencies = getproperty.(straight, :frequency_hz)
    hybrid_relative = getproperty.(hybrid, :S21) ./ getproperty.(straight, :S21)
    pure_relative = getproperty.(pure, :S21) ./ getproperty.(straight, :S21)
    (
        frequencies_hz=frequencies,
        hybrid_relative_amplitude=abs.(hybrid_relative),
        hybrid_relative_phase=unwrap_phase(angle.(hybrid_relative)),
        pure_relative_amplitude=abs.(pure_relative),
        pure_relative_phase=unwrap_phase(angle.(pure_relative)),
        hybrid_s11=getproperty.(hybrid, :S11),
        pure_s11=getproperty.(pure, :S11),
        straight_s11=getproperty.(straight, :S11),
    )
end

function complex_extrapolated(frequency, x, amplitude, phase)
    if frequency < first(x)
        delay = -(phase[2] - phase[1]) / (2pi * (x[2] - x[1]))
        return first(amplitude) * cis(first(phase) - 2pi * (frequency - first(x)) * delay)
    elseif frequency > last(x)
        delay = -(phase[end] - phase[end - 1]) / (2pi * (x[end] - x[end - 1]))
        return last(amplitude) * cis(last(phase) - 2pi * (frequency - last(x)) * delay)
    end
    linear_interpolate(x, amplitude, frequency) *
    cis(linear_interpolate(x, phase, frequency))
end

function amplitude_phase_extrapolated(frequency, x, amplitude, phase)
    if frequency < first(x)
        delay = -(phase[2] - phase[1]) / (2pi * (x[2] - x[1]))
        return first(amplitude), first(phase) - 2pi * (frequency - first(x)) * delay
    elseif frequency > last(x)
        delay = -(phase[end] - phase[end - 1]) / (2pi * (x[end] - x[end - 1]))
        return last(amplitude), last(phase) - 2pi * (frequency - last(x)) * delay
    end
    linear_interpolate(x, amplitude, frequency), linear_interpolate(x, phase, frequency)
end

function complex_interpolated(frequency, x, values; outside=0.0 + 0.0im)
    first(x) <= frequency <= last(x) || return ComplexF64(outside)
    linear_interpolate(x, real.(values), frequency) +
    im * linear_interpolate(x, imag.(values), frequency)
end

function source_elements(topology)
    if topology == :bare
        count = 301
        width = DESIGN.aperture_height_mm / count
        positions = collect(range(
            -DESIGN.aperture_height_mm / 2 + width / 2,
            DESIGN.aperture_height_mm / 2 - width / 2;
            length=count,
        ))
        return [(
            y_mm=position,
            width_mm=width,
            source_amplitude=inv(sqrt(DESIGN.aperture_height_mm)),
            delay_fraction=0.0,
        ) for position in positions]
    end
    subdivisions = 5
    active_width = DESIGN.channel_count * DESIGN.guide_width_mm
    elements = NamedTuple[]
    maximum_delay = maximum_delay_s(DESIGN)
    for center in channel_centers_mm(DESIGN)
        delay_fraction = geometric_delay_s(center, DESIGN) / maximum_delay
        for offset in range(
            -DESIGN.guide_width_mm / 2 + DESIGN.guide_width_mm / (2subdivisions),
            DESIGN.guide_width_mm / 2 - DESIGN.guide_width_mm / (2subdivisions);
            length=subdivisions,
        )
            push!(elements, (
                y_mm=center + offset,
                width_mm=DESIGN.guide_width_mm / subdivisions,
                source_amplitude=inv(sqrt(active_width)),
                delay_fraction,
            ))
        end
    end
    elements
end

function device_transfer(frequency_hz, topology, fraction, isolated)
    topology == :bare && return 1.0 + 0.0im
    base = cis(-2pi * frequency_hz * DESIGN.device_depth_mm * 1e-3 /
               DESIGN.calibrated_aluminium_guide_speed_m_s)
    topology == :straight_footprint && return base
    amplitude_values = topology == :hybrid ? isolated.hybrid_relative_amplitude :
                       isolated.pure_relative_amplitude
    phase_values = topology == :hybrid ? isolated.hybrid_relative_phase :
                   isolated.pure_relative_phase
    maximum_amplitude, maximum_phase = amplitude_phase_extrapolated(
        frequency_hz,
        isolated.frequencies_hz,
        amplitude_values,
        phase_values,
    )
    base * maximum_amplitude^fraction * cis(fraction * maximum_phase)
end

function field_transfer(
    frequency_hz,
    topology,
    observation_x_mm,
    observation_y_mm,
    isolated,
)
    frequency_hz <= 0 && return 0.0 + 0.0im
    wavenumber = 2pi * frequency_hz / STAND.pressure_wave_speed_m_s
    sum(source_elements(topology)) do element
        dy_m = element.width_mm * 1e-3
        dx_m = observation_x_mm * 1e-3
        transverse_m = (observation_y_mm - element.y_mm) * 1e-3
        distance_m = hypot(dx_m, transverse_m)
        obliquity = dx_m / distance_m
        device = device_transfer(
            frequency_hz,
            topology,
            element.delay_fraction,
            isolated,
        )
        element.source_amplitude * dy_m * obliquity /
        sqrt(distance_m) * cis(-wavenumber * distance_m) * device
    end
end

function spectrum_transfer(pulse, topology, isolated; two_way=false)
    threshold = maximum(abs, pulse.spectrum) * 1e-5
    result = zeros(ComplexF64, length(pulse.spectrum))
    for index in eachindex(result)
        abs(pulse.spectrum[index]) >= threshold || continue
        frequency = pulse.frequency_hz[index]
        field = field_transfer(
            frequency,
            topology,
            STAND.nominal_focus_depth_mm,
            0.0,
            isolated,
        )
        result[index] = two_way ? field^2 * (frequency / DESIGN.frequency_hz)^2 : field
    end
    result
end

function ringing_transfer(pulse, topology, isolated)
    topology == :bare && return zeros(ComplexF64, length(pulse.spectrum))
    source = topology == :hybrid ? isolated.hybrid_s11 :
             topology == :pure_ttd ? isolated.pure_s11 : isolated.straight_s11
    ComplexF64[
        complex_interpolated(
            frequency,
            isolated.frequencies_hz,
            source;
            outside=0.0 + 0.0im,
        )
        for frequency in pulse.frequency_hz
    ]
end

analytic_envelope(signal) = abs.(ImpulseRiskAnalysis.analytic_signal(signal))

function contiguous_fwhm(coordinate, amplitude)
    peak = argmax(amplitude)
    threshold = amplitude[peak] / sqrt(2)
    left, right = peak, peak
    while left > 1 && amplitude[left - 1] >= threshold
        left -= 1
    end
    while right < length(amplitude) && amplitude[right + 1] >= threshold
        right += 1
    end
    (width=coordinate[right] - coordinate[left], left, right, peak)
end

function topology_result(topology, pulse, isolated, bare_echo_norm)
    focus_transfer = spectrum_transfer(pulse, topology, isolated)
    echo_transfer = spectrum_transfer(pulse, topology, isolated; two_way=true)
    ring_transfer = ringing_transfer(pulse, topology, isolated)
    focus = irfft(pulse.spectrum .* focus_transfer, length(pulse.signal))
    echo = irfft(pulse.spectrum .* echo_transfer, length(pulse.signal))
    ringing = irfft(pulse.spectrum .* ring_transfer, length(pulse.signal))
    echo_norm = norm(echo)
    echo_mf_gain = echo_norm / bare_echo_norm
    echo_energy_gain = echo_mf_gain^2
    echo_peak = argmax(analytic_envelope(echo))
    half_window = round(Int, 5e-6 / (pulse.time_s[2] - pulse.time_s[1]))
    window = max(1, echo_peak - half_window):min(length(echo), echo_peak + half_window)
    ringing_peak = maximum(analytic_envelope(ringing)[window])
    echo_peak_value = maximum(analytic_envelope(echo)[window])
    echo_to_ringing_db = ringing_peak > 0 ?
        20log10(echo_peak_value / ringing_peak) : Inf

    transverse_mm = collect(-35.0:0.25:35.0)
    carrier = [
        abs(field_transfer(
            DESIGN.frequency_hz,
            topology,
            STAND.nominal_focus_depth_mm,
            y,
            isolated,
        ))
        for y in transverse_mm
    ]
    fwhm = contiguous_fwhm(transverse_mm, carrier)
    exclusion = max(2, fwhm.right - fwhm.left + 1)
    sidelobe_indices = [
        index for index in eachindex(carrier)
        if abs(index - fwhm.peak) > exclusion
    ]
    sidelobe = isempty(sidelobe_indices) ? NaN :
               maximum(carrier[sidelobe_indices]) / maximum(carrier)
    (
        topology,
        focus,
        echo,
        ringing,
        focus_transfer,
        echo_transfer,
        echo_mf_gain_bare=echo_mf_gain,
        echo_energy_gain_bare=echo_energy_gain,
        echo_to_ringing_db,
        equal_noise_snr_proxy_db=20log10(echo_mf_gain),
        fwhm_mm=fwhm.width,
        sidelobe_ratio=sidelobe,
        carrier_focus_amplitude=maximum(carrier),
        transverse_mm,
        carrier_profile=carrier,
    )
end

function write_common_metrics(path, results, isolated_summary)
    isolated_by_variant = Dict(row.variant => row for row in isolated_summary)
    open(path, "w") do io
        println(io, join(string.(metric_columns()), ','))
        for result in results
            isolated_variant = result.topology == :hybrid ? :hybrid :
                               result.topology == :pure_ttd ? :pure_same_depth : :straight
            local_metrics = isolated_by_variant[isolated_variant]
            values = (
                topology=result.topology,
                case_id=:void_r2_focus,
                reference=:bare,
                frequency_hz=DESIGN.frequency_hz,
                input_energy_j_per_m=1.0,
                echo_mf_gain_bare=NaN,
                echo_energy_gain_bare=NaN,
                echo_to_ringing_db=NaN,
                equal_noise_snr_proxy_db=NaN,
                dead_zone_depth_mm=NaN,
                focus_peak_gain_bare=NaN,
                focus_energy_gain_bare=NaN,
                fwhm_mm=NaN,
                dof_minus3db_mm=NaN,
                sidelobe_ratio=NaN,
                pulse_broadening_ratio=local_metrics.pulse_broadening_ratio,
                pulse_correlation=local_metrics.pulse_correlation,
                postcursor_ratio=local_metrics.postcursor_ratio,
                weighted_transmission=result.topology == :bare ? 1.0 :
                    local_metrics.weighted_transmission,
                weighted_reflection=result.topology == :bare ? 0.0 :
                    local_metrics.weighted_reflection,
                weighted_mode_conversion=result.topology == :bare ? 0.0 :
                    local_metrics.weighted_mode_conversion,
                minimum_feature_mm=result.topology == :hybrid ?
                    design_summary(DESIGN).minimum_insert_width_mm :
                    result.topology == :pure_ttd ? DESIGN.guide_width_mm : NaN,
                verdict=result.topology == :hybrid ? :conditional : :control,
            )
            println(io, join((getproperty(values, column) for column in metric_columns()), ','))
        end
    end
end

function write_proxy_metrics(path, results)
    open(path, "w") do io
        println(io, "model_class,defect_proxy,topology,echo_mf_gain_bare,echo_energy_gain_bare,echo_to_device_ringing_db,equal_noise_snr_proxy_db,carrier_focus_gain_bare,fwhm_mm,sidelobe_ratio")
        bare_focus = first(results).carrier_focus_amplitude
        for result in results
            println(io, join((
                "reciprocal_2d_scalar_strict_proxy",
                "point_scatterer_ka_approximately_0p5_not_explicit_void",
                result.topology,
                result.echo_mf_gain_bare,
                result.echo_energy_gain_bare,
                result.echo_to_ringing_db,
                result.equal_noise_snr_proxy_db,
                result.carrier_focus_amplitude / bare_focus,
                result.fwhm_mm,
                result.sidelobe_ratio,
            ), ','))
        end
    end
end

function main()
    isolated = load_isolated_transfer()
    isolated_data = JLD2.load(joinpath(ISOLATED_ROOT, "isolated_analysis.jld2"))
    isolated_summary = isolated_data["summary_rows"]
    pulse = pulse_spectrum(PulseConfig(
        center_frequency_hz=DESIGN.frequency_hz,
        cycles=5.0,
        samples_per_period=80,
        fft_length=16384,
    ))
    bare_transfer = spectrum_transfer(pulse, :bare, isolated; two_way=true)
    bare_echo = irfft(pulse.spectrum .* bare_transfer, length(pulse.signal))
    bare_norm = norm(bare_echo)
    results = [topology_result(topology, pulse, isolated, bare_norm) for topology in TOPOLOGIES]
    mkpath(OUTPUT_ROOT)
    common_path = joinpath(OUTPUT_ROOT, "common_metrics.csv")
    write_common_metrics(common_path, results, isolated_summary)
    proxy_path = joinpath(OUTPUT_ROOT, "strict_proxy_metrics.csv")
    write_proxy_metrics(proxy_path, results)
    jldsave(joinpath(OUTPUT_ROOT, "aperture_proxy.jld2");
            results, pulse, isolated_summary, model_class="reciprocal_2D_scalar_strict_proxy")

    echo_panel = plot(
        xlabel="time, μs", ylabel="normalized defect echo",
        title="Reciprocal point-scatterer pulse-echo proxy",
    )
    profile_panel = plot(
        xlabel="transverse y, mm", ylabel="carrier amplitude / bare peak",
        title="Carrier transverse field at 60 mm",
    )
    metrics_panel = bar(
        String.(getproperty.(results, :topology)),
        getproperty.(results, :echo_mf_gain_bare);
        ylabel="white-noise matched-filter gain / bare",
        title="Primary G3 proxy", label=false,
        color=[:gray50, :darkorange, :firebrick, :royalblue],
    )
    ringing_panel = bar(
        String.(getproperty.(results, :topology)),
        [isfinite(row.echo_to_ringing_db) ? row.echo_to_ringing_db : 60.0 for row in results];
        ylabel="echo / device ringing, dB", title="G3 ringing proxy",
        label=false, color=[:gray50, :darkorange, :firebrick, :royalblue],
    )
    colors = Dict(:bare => :black, :straight_footprint => :darkorange,
                  :pure_ttd => :firebrick, :hybrid => :royalblue)
    for result in results
        envelope = analytic_envelope(result.echo)
        plot!(echo_panel, pulse.time_s .* 1e6, envelope ./ maximum(analytic_envelope(first(results).echo));
              color=colors[result.topology], linewidth=2, label=String(result.topology),
              xlims=(10, 70))
        plot!(profile_panel, result.transverse_mm,
              result.carrier_profile ./ maximum(first(results).carrier_profile);
              color=colors[result.topology], linewidth=2, label=String(result.topology))
    end
    savefig(plot(
        echo_panel, profile_panel, metrics_panel, ringing_panel;
        layout=(2, 2), size=(1500, 950), margin=5Plots.mm,
    ), joinpath(OUTPUT_ROOT, "aperture_defect_proxy.png"))
    for row in results
        println(
            "[+] $(row.topology): echo MF gain=",
            round(row.echo_mf_gain_bare; digits=3),
            ", echo/ringing=", round(row.echo_to_ringing_db; digits=2),
            " dB, FWHM=", row.fwhm_mm,
            " mm, sidelobe=", round(row.sidelobe_ratio; digits=3),
        )
    end
    println("[+] $common_path")
    println("[+] $proxy_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
