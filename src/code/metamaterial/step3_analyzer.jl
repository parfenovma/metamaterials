using JLD2

include(joinpath(@__DIR__, "spectral_analysis.jl"))
using .SpectralAnalysis

const ANALYSIS_FORMAT_VERSION = 1

function require_v2(data, path)
    get(data, "data_format_version", 0) >= 2 ||
        error("$path does not contain v2 port observables; rerun step2_solver.jl")
end

function analyze_dataset(path::AbstractString; config::SpectrumConfig=SpectrumConfig())
    data = jldopen(path, "r") do file
        Dict(
            "data_format_version" => file["data_format_version"],
            "time_s" => file["time_s"],
            "source_drive_mpa" => file["source_drive_mpa"],
            "left_normal_traction_mpa" => file["left_normal_traction_mpa"],
            "right_normal_traction_mpa" => file["right_normal_traction_mpa"],
            "profile_type" => file["profile_type"],
            "profile_parameters" => file["profile_parameters"],
            "frequency_hz" => file["frequency_hz"],
        )
    end
    require_v2(data, path)

    time_s = data["time_s"]
    drive = data["source_drive_mpa"]
    # Positive pressure is defined as minus outward normal traction.
    pressure_left = .-data["left_normal_traction_mpa"]
    pressure_right = .-data["right_normal_traction_mpa"]

    (
        metadata=(
            source_path=abspath(path),
            profile_type=data["profile_type"],
            profile_parameters=data["profile_parameters"],
            frequency_hz=data["frequency_hz"],
        ),
        drive_to_left=analyze_transfer(time_s, drive, pressure_left; config),
        drive_to_right=analyze_transfer(time_s, drive, pressure_right; config),
        envelope_delay_right_s=envelope_delay(time_s, drive, pressure_right),
    )
end

function save_analysis(
    sample_path::AbstractString;
    reference_path::Union{Nothing, AbstractString}=nothing,
    output_dir::AbstractString="4_characteristics",
    config::SpectrumConfig=SpectrumConfig(),
)
    sample = analyze_dataset(sample_path; config)
    reference = isnothing(reference_path) ? nothing : analyze_dataset(reference_path; config)
    calibrated_transmission = isnothing(reference) ?
                              nothing :
                              relative_transfer(sample.drive_to_right, reference.drive_to_right)
    excess_envelope_delay_s = isnothing(reference) ?
                              nothing :
                              sample.envelope_delay_right_s - reference.envelope_delay_right_s

    mkpath(output_dir)
    stem = splitext(basename(sample_path))[1]
    output_path = joinpath(output_dir, "characteristics_$(stem).jld2")
    jldsave(
        output_path;
        analysis_format_version=ANALYSIS_FORMAT_VERSION,
        spectrum_config=config,
        metadata=sample.metadata,
        drive_to_left=sample.drive_to_left,
        drive_to_right=sample.drive_to_right,
        envelope_delay_right_s=sample.envelope_delay_right_s,
        reference_path=isnothing(reference_path) ? nothing : abspath(reference_path),
        calibrated_transmission,
        excess_envelope_delay_s,
    )
    output_path
end

function main(args=ARGS)
    isempty(args) && error(
        "usage: julia step3_analyzer.jl SAMPLE.jld2 [REFERENCE_A0.jld2] [OUTPUT_DIR]",
    )
    sample_path = args[1]
    reference_path = length(args) >= 2 ? args[2] : nothing
    output_dir = length(args) >= 3 ? args[3] : "4_characteristics"
    output_path = save_analysis(sample_path; reference_path, output_dir)
    println("[+] Characteristics saved: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
