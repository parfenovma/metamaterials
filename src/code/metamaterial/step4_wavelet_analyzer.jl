using JLD2

include(joinpath(@__DIR__, "wavelet_analysis.jl"))
using .WaveletAnalysis

const WAVELET_ANALYSIS_FORMAT_VERSION = 1

function load_wavelet_dataset(path::AbstractString)
    data = jldopen(path, "r") do file
        Dict(
            "data_format_version" => file["data_format_version"],
            "time_s" => file["time_s"],
            "source_drive_mpa" => file["source_drive_mpa"],
            "right_normal_traction_mpa" => file["right_normal_traction_mpa"],
            "profile_type" => file["profile_type"],
            "profile_parameters" => file["profile_parameters"],
            "frequency_hz" => file["frequency_hz"],
        )
    end
    get(data, "data_format_version", 0) >= 2 ||
        error("$path does not contain v2 port observables; rerun step2_solver.jl")
    (
        time_s=Float64.(data["time_s"]),
        drive_mpa=Float64.(data["source_drive_mpa"]),
        # Positive pressure is minus outward normal traction.
        right_pressure_mpa=.-Float64.(data["right_normal_traction_mpa"]),
        metadata=(
            source_path=abspath(path),
            profile_type=data["profile_type"],
            profile_parameters=data["profile_parameters"],
            carrier_frequency_hz=Float64(data["frequency_hz"]),
        ),
    )
end

function config_description(config::MorletConfig)
    Dict(string(name) => getfield(config, name) for name in fieldnames(typeof(config)))
end

function analyze_wavelet_pair(
    sample_path::AbstractString,
    reference_path::AbstractString;
    config::Union{Nothing, MorletConfig}=nothing,
)
    sample = load_wavelet_dataset(sample_path)
    reference = load_wavelet_dataset(reference_path)
    sample.time_s == reference.time_s ||
        error("sample and reference time grids differ; use the matching carrier reference")
    sample.metadata.carrier_frequency_hz == reference.metadata.carrier_frequency_hz ||
        error("sample and reference carrier frequencies differ")

    effective_config = isnothing(config) ?
                       recommended_morlet_config(reference.time_s) : config

    drive_cwt = continuous_wavelet_transform(
        reference.time_s,
        reference.drive_mpa;
        config=effective_config,
    )
    reference_cwt = continuous_wavelet_transform(
        reference.time_s,
        reference.right_pressure_mpa;
        config=effective_config,
    )
    sample_cwt = continuous_wavelet_transform(
        sample.time_s,
        sample.right_pressure_mpa;
        config=effective_config,
    )
    comparison = compare_wavelets(
        reference_cwt,
        sample_cwt;
        config=effective_config,
    )
    (
        sample_metadata=sample.metadata,
        reference_metadata=reference.metadata,
        effective_config,
        drive_cwt,
        reference_cwt,
        sample_cwt,
        comparison,
    )
end

function save_wavelet_analysis(
    sample_path::AbstractString,
    reference_path::AbstractString;
    output_dir::AbstractString="5_wavelets",
    config::Union{Nothing, MorletConfig}=nothing,
)
    result = analyze_wavelet_pair(sample_path, reference_path; config)
    mkpath(output_dir)
    stem = splitext(basename(sample_path))[1]
    output_path = joinpath(output_dir, "wavelet_$(stem).jld2")
    jldsave(
        output_path;
        wavelet_analysis_format_version=WAVELET_ANALYSIS_FORMAT_VERSION,
        morlet_config=config_description(result.effective_config),
        sample_metadata=result.sample_metadata,
        reference_metadata=result.reference_metadata,
        drive_cwt=result.drive_cwt,
        reference_cwt=result.reference_cwt,
        sample_cwt=result.sample_cwt,
        comparison=result.comparison,
    )
    output_path
end

function main(args=ARGS)
    length(args) >= 2 || error(
        "usage: julia step4_wavelet_analyzer.jl SAMPLE.jld2 REFERENCE.jld2 [OUTPUT_DIR]",
    )
    output_dir = length(args) >= 3 ? args[3] : "5_wavelets"
    output_path = save_wavelet_analysis(args[1], args[2]; output_dir)
    println("[+] Wavelet analysis saved: $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
