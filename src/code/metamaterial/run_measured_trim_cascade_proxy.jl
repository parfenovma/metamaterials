module MeasuredTrimCascadeProxy

ENV["GKSwstype"] = "100"

using FFTW
using JLD2
using Plots

include(joinpath(@__DIR__, "run_dispersive_delay_proxy.jl"))
using .DispersiveDelayProxy

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SOURCE_ROOT = get(
    ENV,
    "METAMATERIALS_MEASURED_TRIM_SOURCE",
    joinpath(PROJECT_ROOT, "tmp", "tapered_shunt_trim_pilot"),
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MEASURED_TRIM_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "measured_trim_cascade_proxy"),
)
const MAX_CASCADE_COUNT = 4

point_path(stage, variant, frequency_hz) = joinpath(
    SOURCE_ROOT,
    String(stage),
    "$(variant)_$(round(Int, frequency_hz))hz.jld2",
)

const FEM_FREQUENCIES_HZ = collect(162.0e3:4.0e3:322.0e3)

function load_calibrated_transfer(stage)
    transfer = ComplexF64[]
    for frequency_hz in FEM_FREQUENCIES_HZ
        solid = JLD2.load(point_path(stage, :solid, frequency_hz))["row"]
        trim = JLD2.load(point_path(stage, :max, frequency_hz))["row"]
        push!(transfer, trim.transmission_amplitude / solid.transmission_amplitude)
    end
    raw_phase = angle.(transfer)
    phase = Float64[first(raw_phase)]
    for index in 2:length(raw_phase)
        push!(phase, last(phase) + mod(raw_phase[index] - raw_phase[index - 1] + pi, 2pi) - pi)
    end
    (; amplitude=abs.(transfer), phase)
end

function interpolate_transfer(data, frequency_hz)
    amplitude = DispersiveDelayProxy.interpolate_linear(
        FEM_FREQUENCIES_HZ,
        data.amplitude,
        frequency_hz,
    )
    phase = DispersiveDelayProxy.interpolate_linear(
        FEM_FREQUENCIES_HZ,
        data.phase,
        frequency_hz,
    )
    isfinite(amplitude) && isfinite(phase) ? amplitude * cis(phase) : 0.0 + 0.0im
end

circular_error(first, second) = mod(first - second + pi, 2pi) - pi

function select_cascade_counts(required_trim_rad, unit_phase_rad)
    best = nothing
    for common_offset_rad in range(0.0, 2pi; length=1441)
        counts = Int[]
        errors = Float64[]
        for required in required_trim_rad
            candidates = collect(0:MAX_CASCADE_COUNT)
            candidate_errors = [
                circular_error(count * unit_phase_rad, required + common_offset_rad)
                for count in candidates
            ]
            index = argmin(abs2.(candidate_errors))
            push!(counts, candidates[index])
            push!(errors, candidate_errors[index])
        end
        rms = sqrt(sum(abs2, errors) / length(errors))
        candidate = (; common_offset_rad, counts, errors, rms)
        (best === nothing || candidate.rms < best.rms) && (best = candidate)
    end
    best
end

function measured_focus_transfer(spectrum, branch, delays_s, trim_data; path_loss_scale)
    config = DispersiveDelayProxy.APERTURE
    centers = DispersiveDelayProxy.aperture_centers_m(config)
    distances = hypot.(config.focal_distance_m, centers)
    frequency0 = spectrum.config.center_frequency_hz
    omega0 = 2pi * frequency0
    k0 = DispersiveDelayProxy.interpolate_linear(
        branch.frequencies_hz,
        branch.wavenumber,
        frequency0,
    )
    vg0 = DispersiveDelayProxy.interpolate_linear(
        branch.frequencies_hz,
        branch.group_velocity,
        frequency0,
    )
    path_lengths_m = vg0 .* delays_s
    required_trim_rad = k0 .* path_lengths_m .- omega0 .* delays_s
    trim_at_f0 = interpolate_transfer(trim_data, frequency0)
    trim_phase_rad = angle(trim_at_f0)
    trim_phase_rad > 0 || error("measured trim must provide positive carrier phase")
    selection = select_cascade_counts(required_trim_rad, trim_phase_rad)
    cascade_count = selection.counts
    implemented_trim_rad = cascade_count .* trim_phase_rad
    phase_error_rad = selection.errors

    transfer = zeros(ComplexF64, length(spectrum.frequency_hz))
    for (frequency_index, frequency_hz) in enumerate(spectrum.frequency_hz)
        iszero(spectrum.spectrum[frequency_index]) && continue
        omega = 2pi * frequency_hz
        k = DispersiveDelayProxy.interpolate_linear(
            branch.frequencies_hz,
            branch.wavenumber,
            frequency_hz,
        )
        vg = DispersiveDelayProxy.interpolate_linear(
            branch.frequencies_hz,
            branch.group_velocity,
            frequency_hz,
        )
        isfinite(k) && isfinite(vg) || continue
        trim = interpolate_transfer(trim_data, frequency_hz)
        wavelength_m = config.propagation_speed_m_per_s / frequency_hz
        decay_rate = path_loss_scale * (
            config.rayleigh_alpha_per_s + config.rayleigh_beta_s * omega^2
        ) / 2
        transfer[frequency_index] = sum(eachindex(centers)) do index
            propagation = config.element_width_m /
                          sqrt(wavelength_m * distances[index]) *
                          cis(-omega * distances[index] / config.propagation_speed_m_per_s)
            modal_delay = path_lengths_m[index] / vg
            path = exp(-decay_rate * modal_delay) * cis(-k * path_lengths_m[index])
            propagation * path * trim^cascade_count[index]
        end
    end
    (
        transfer=transfer,
        path_lengths_m=path_lengths_m,
        required_trim_rad=required_trim_rad,
        implemented_trim_rad=implemented_trim_rad,
        phase_error_rad=phase_error_rad,
        cascade_count=cascade_count,
        trim_phase_rad=trim_phase_rad,
        common_phase_offset_rad=selection.common_offset_rad,
    )
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

function run()
    mkpath(OUTPUT_ROOT)
    raw_spectrum = DispersiveDelayProxy.pulse_spectrum()
    band20 = DispersiveDelayProxy.spectral_band(raw_spectrum, -20.0)
    spectrum = DispersiveDelayProxy.truncated_spectrum(raw_spectrum, band20)
    branch = DispersiveDelayProxy.branch_data(band20.upper_hz)
    delays_s = DispersiveDelayProxy.quantize_delays(
        DispersiveDelayProxy.ideal_delays_s(DispersiveDelayProxy.APERTURE),
        6,
    )
    uniform_transfer = DispersiveDelayProxy.uniform_transfer(spectrum)
    uniform = irfft(spectrum.spectrum .* uniform_transfer, length(spectrum.signal))
    ideal = DispersiveDelayProxy.proxy_transfer(
        spectrum,
        branch,
        DispersiveDelayProxy.ideal_delays_s(DispersiveDelayProxy.APERTURE);
        model=:ideal_ttd,
        loss_scale=0.0,
    )
    ideal_waveform = irfft(spectrum.spectrum .* ideal.transfer, length(spectrum.signal))
    dt_s = spectrum.time_s[2] - spectrum.time_s[1]

    rows = NamedTuple[]
    waveforms = Dict{Symbol, Vector{Float64}}(:uniform => uniform, :ideal => ideal_waveform)
    for (label, stage, path_loss_scale) in (
        (:lossless, :fine_lossless, 0.0),
        (:rayleigh_x1, :fine_lossy, 1.0),
    )
        trim_data = load_calibrated_transfer(stage)
        proxy = measured_focus_transfer(
            spectrum,
            branch,
            delays_s,
            trim_data;
            path_loss_scale,
        )
        waveform = irfft(spectrum.spectrum .* proxy.transfer, length(spectrum.signal))
        waveforms[label] = waveform
        metrics = DispersiveDelayProxy.pulse_metrics(
            waveform,
            uniform,
            ideal_waveform,
            dt_s,
        )
        push!(rows, (
            case=label,
            delay_states=6,
            maximum_cascade_count=maximum(proxy.cascade_count),
            distinct_cascade_counts=join(sort(unique(proxy.cascade_count)), ';'),
            maximum_path_mm=maximum(proxy.path_lengths_m) * 1e3,
            maximum_required_trim_rad=maximum(proxy.required_trim_rad),
            maximum_implemented_trim_rad=maximum(proxy.implemented_trim_rad),
            trim_unit_phase_rad=proxy.trim_phase_rad,
            common_phase_offset_rad=proxy.common_phase_offset_rad,
            rms_carrier_phase_error_rad=sqrt(sum(abs2, proxy.phase_error_rad) / length(proxy.phase_error_rad)),
            gain_peak=metrics.gain_peak,
            broadening_ratio=metrics.broadening_ratio,
            pulse_correlation=metrics.pulse_correlation,
            postcursor_ratio=metrics.postcursor_ratio,
            pass_all=metrics.gain_peak >= 2.0 &&
                     metrics.broadening_ratio <= 1.25 &&
                     metrics.pulse_correlation >= 0.90 &&
                     metrics.postcursor_ratio <= 0.10,
        ))
    end

    csv_path = joinpath(OUTPUT_ROOT, "measured_trim_cascade_proxy.csv")
    write_rows(csv_path, rows)
    ideal_peak = argmax(abs.(ideal_waveform))
    time_us = spectrum.time_s .* 1e6
    panel = plot(
        time_us,
        uniform;
        xlim=(max(0.0, time_us[ideal_peak] - 22), time_us[ideal_peak] + 42),
        linewidth=1.5,
        color=:black,
        label="uniform",
        xlabel="time, μs",
        ylabel="focus observable, a.u.",
        title="Measured shunt trim cascaded onto group-matched path proxy",
        gridalpha=0.25,
    )
    plot!(panel, time_us, ideal_waveform; linewidth=2, color=:seagreen, label="ideal TTD")
    plot!(panel, time_us, waveforms[:lossless]; linewidth=2, color=:royalblue,
          label="measured trim, lossless")
    plot!(panel, time_us, waveforms[:rayleigh_x1]; linewidth=2, color=:darkorange,
          label="measured trim, Rayleigh x1")
    png_path = joinpath(OUTPUT_ROOT, "measured_trim_cascade_proxy.png")
    savefig(panel, png_path)
    println("[+] $csv_path")
    println("[+] $png_path")
    for row in rows
        println(
            "[+] $(row.case): counts=", row.distinct_cascade_counts,
            ", trim=", round(row.maximum_implemented_trim_rad; digits=3), " rad",
            ", G=", round(row.gain_peak; digits=3),
            ", Bt=", round(row.broadening_ratio; digits=3),
            ", rho=", round(row.pulse_correlation; digits=4),
            ", post=", round(row.postcursor_ratio; digits=4),
            ", pass=", row.pass_all,
        )
    end
    (; rows, waveforms)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
