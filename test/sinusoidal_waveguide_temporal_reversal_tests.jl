module SinusoidalWaveguideTemporalReversalTests

using FFTW
using LinearAlgebra
using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "analyze_sinusoidal_waveguide_temporal_reversal.jl",
))
using .AnalyzeSinusoidalWaveguideTemporalReversal

@testset "scalar temporal time reversal" begin
    count = 256
    pulse = zeros(Float64, count)
    pulse[20:24] .= [0.0, 0.5, 1.0, 0.5, 0.0]
    transfer = ones(ComplexF64, count ÷ 2 + 1)
    result = temporal_refocus(transfer, pulse)
    metrics = temporal_focus_metrics(result, 1.0; focus_window_s=20.0)
    @test argmax(result.envelope) == count ÷ 2 + 1
    @test norm(irfft(conj.(rfft(pulse)), count)) ≈
        metrics.replay_norm_before_normalization
    @test metrics.peak_envelope > 0
    @test metrics.half_power_width_s > 0
    @test 0 <= metrics.half_power_energy_fraction <= 1
    @test 0 <= metrics.ringing_energy_fraction <= 1

    delayed = transfer .* exp.(-2pi * im .* collect(0:(count ÷ 2)) .* 17 / count)
    delayed_result = temporal_refocus(delayed, pulse)
    @test delayed_result.envelope ≈ result.envelope atol=1.0e-10

    @test_throws DimensionMismatch temporal_refocus(transfer[1:end-1], pulse)
    @test_throws ArgumentError temporal_focus_metrics(result, 0.0)
end

end
