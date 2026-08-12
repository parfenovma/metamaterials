using Test
using JLD2

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "step3_analyzer.jl"))

@testset "step3 JLD2 integration" begin
    mktempdir() do directory
        sample_count = 512
        dt = 1.0e-6
        time = collect(0:(sample_count - 1)) .* dt
        drive = zeros(sample_count)
        drive[50] = 1.0

        function write_fixture(path, delay_samples, gain, amplitude)
            right_pressure = zeros(sample_count)
            right_pressure[50 + delay_samples] = gain
            left_pressure = zeros(sample_count)
            left_pressure[50] = amplitude
            jldsave(
                path;
                data_format_version=2,
                time_s=time,
                source_drive_mpa=drive,
                left_normal_traction_mpa=-left_pressure,
                right_normal_traction_mpa=-right_pressure,
                profile_type="SyntheticProfile",
                profile_parameters=Dict("gain" => gain),
                frequency_hz=100e3,
            )
        end

        sample_path = joinpath(directory, "sample.jld2")
        reference_path = joinpath(directory, "reference.jld2")
        write_fixture(sample_path, 9, 0.4, 0.1)
        write_fixture(reference_path, 4, 0.8, 0.0)

        output_path = save_analysis(
            sample_path;
            reference_path,
            output_dir=joinpath(directory, "characteristics"),
            config=SpectrumConfig(
                zero_padding_factor=1,
                input_floor_relative=0.0,
                regularization_relative=0.0,
            ),
        )

        @test isfile(output_path)
        saved = load(output_path)
        @test saved["analysis_format_version"] == 1
        calibrated = saved["calibrated_transmission"]
        interior = 3:(length(calibrated.frequency_hz) - 2)
        @test maximum(abs.(calibrated.amplitude[interior] .- 0.5)) < 1.0e-12
        @test maximum(abs.(calibrated.group_delay_s[interior] .- 5.0e-6)) < 1.0e-12
    end
end
