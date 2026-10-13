using DelimitedFiles
using LinearAlgebra
using Test

const ROOT = @__DIR__
const PROFILE_PATH = joinpath(ROOT, "data", "receiver_scattered_profiles.csv")
const SUMMARY_PATH = joinpath(ROOT, "data", "material_lens_defect_summary.csv")

profile_data, profile_header = readdlm(PROFILE_PATH, ',', header=true)
summary_data, summary_header = readdlm(SUMMARY_PATH, ',', header=true)

profile_columns = vec(String.(profile_header))
summary_columns = vec(String.(summary_header))
column(name, data, columns) = Float64.(data[:, only(findall(==(name), columns))])

y_mm = column("y_mm", profile_data, profile_columns)
sinusoidal = column("selected_scattered_ux_m", profile_data, profile_columns)
straight = column("uniform_scattered_ux_m", profile_data, profile_columns)

summary(name) = Float64(summary_data[1, only(findall(==(name), summary_columns))])
rms(values) = norm(values) / sqrt(length(values))

@testset "sinusoidal defectoscope data provenance" begin
    @test length(y_mm) == 31
    @test y_mm == collect(-30.0:2.0:30.0)
    @test all(>(0), sinusoidal)
    @test all(>(0), straight)

    @test rms(sinusoidal) ≈ summary("selected_receiver_scattered_ux_rms_m") rtol=1e-12
    @test rms(straight) ≈ summary("uniform_receiver_scattered_ux_rms_m") rtol=1e-12

    scattering_gain = rms(sinusoidal) / rms(straight)
    snr_improvement_db = 20 * log10(scattering_gain)
    @test scattering_gain ≈ summary("scattering_gain_equal_noise_snr_gain") rtol=1e-12
    @test snr_improvement_db ≈ summary("snr_improvement_db") rtol=1e-12
    @test summary("local_energy_gain") > 1.9

    javascript = read(joinpath(ROOT, "demo_data.js"), String)
    function javascript_array(name)
        found = match(Regex(name * raw": Object\.freeze\(\[([^\]]+)\]\)"), javascript)
        @assert found !== nothing
        return parse.(Float64, strip.(split(String(found.captures[1]), ',')))
    end

    @test occursin("source: \"full-FEM paired defect-minus-reference profile, 242 kHz\"", javascript)
    @test occursin(string(summary("scattering_gain_equal_noise_snr_gain")), javascript)
    @test occursin(string(summary("snr_improvement_db")), javascript)
    @test javascript_array("yMm") == y_mm
    @test javascript_array("sinusoidalUxNm") ≈ 1e9 .* sinusoidal rtol=1e-15
    @test javascript_array("straightUxNm") ≈ 1e9 .* straight rtol=1e-15

    html = read(joinpath(ROOT, "index.html"), String)
    @test occursin("FEM + синтетический шум", html)
    @test occursin("не измеренная POD", html)
    @test occursin("demo_data.js", html)
end

println("[+] demo verified: scattering gain = $(round(rms(sinusoidal) / rms(straight); digits=4)), " *
        "equal-noise SNR proxy = $(round(20 * log10(rms(sinusoidal) / rms(straight)); digits=3)) dB")
