include(joinpath(@__DIR__, "step4_wavelet_analyzer.jl"))
include(joinpath(@__DIR__, "spectral_analysis.jl"))
using .SpectralAnalysis

const CONTROL_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", "tmp", "pilot_220khz"))
const CONTROL_SIGNAL_DIR = joinpath(CONTROL_ROOT, "signals")
const CONTROL_OUTPUT_DIR = joinpath(CONTROL_ROOT, "wavelets")

const WAVELET_CONTROL_CASES = [
    (slug="sin_A_2.5_N_2_STG", carrier_khz=220.0),
    (slug="pow_A_2.5_N_2_P_4.0_STG", carrier_khz=220.0),
    (slug="exp_A_2.5_N_2_K_3.0_STG", carrier_khz=240.0),
    (slug="A_2.5", carrier_khz=240.0),
]

function signal_path(slug, carrier_khz)
    joinpath(CONTROL_SIGNAL_DIR, "data_$(slug)_F_$(carrier_khz).jld2")
end

function write_summary(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function fft_control_point(sample_path, reference_path, carrier_frequency_hz)
    sample = load_wavelet_dataset(sample_path)
    reference = load_wavelet_dataset(reference_path)
    sample_transfer = analyze_transfer(
        sample.time_s,
        sample.drive_mpa,
        sample.right_pressure_mpa,
    )
    reference_transfer = analyze_transfer(
        reference.time_s,
        reference.drive_mpa,
        reference.right_pressure_mpa,
    )
    calibrated = relative_transfer(sample_transfer, reference_transfer)
    SpectralAnalysis.value_at_frequency(calibrated, carrier_frequency_hz)
end

function run_wavelet_controls()
    mkpath(CONTROL_OUTPUT_DIR)
    rows = NamedTuple[]
    for control in WAVELET_CONTROL_CASES
        sample_path = signal_path(control.slug, control.carrier_khz)
        reference_path = signal_path("sin_A_0.0_N_2", control.carrier_khz)
        output_path = save_wavelet_analysis(
            sample_path,
            reference_path;
            output_dir=CONTROL_OUTPUT_DIR,
        )
        saved = JLD2.load(output_path)
        point = wavelet_value_at_frequency(
            saved["comparison"],
            control.carrier_khz * 1000.0,
        )
        fft_point = fft_control_point(
            sample_path,
            reference_path,
            control.carrier_khz * 1000.0,
        )
        row = (
            slug=control.slug,
            carrier_frequency_hz=control.carrier_khz * 1000.0,
            cwt_frequency_hz=point.frequency_hz,
            cwt_energy_amplitude_ratio=point.amplitude_ratio,
            centroid_delay_s=point.centroid_delay_s,
            peak_delay_s=point.peak_delay_s,
            correlation_delay_s=point.correlation_delay_s,
            envelope_correlation=point.envelope_correlation,
            reference_energy_relative=point.reference_energy_relative,
            cwt_valid=point.valid,
            fft_amplitude=fft_point.amplitude,
            fft_group_delay_s=fft_point.group_delay_s,
            fft_valid=fft_point.valid,
            combined_confident=point.valid && fft_point.valid && fft_point.amplitude >= 0.05,
            result_path=abspath(output_path),
        )
        push!(rows, row)
        println(
            "[+] $(control.slug) @ $(control.carrier_khz) kHz: ",
            "amplitude_ratio=$(round(point.amplitude_ratio; digits=4)), ",
            "centroid_delay_us=$(round(1e6 * point.centroid_delay_s; digits=3)), ",
            "correlation_delay_us=$(round(1e6 * point.correlation_delay_s; digits=3)), ",
            "fft_delay_us=$(round(1e6 * fft_point.group_delay_s; digits=3)), ",
            "confident=$(row.combined_confident)",
        )
    end
    summary_path = joinpath(CONTROL_OUTPUT_DIR, "wavelet_control_summary.csv")
    write_summary(summary_path, rows)
    println("[+] Wavelet control summary: $summary_path")
    summary_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_wavelet_controls()
end
