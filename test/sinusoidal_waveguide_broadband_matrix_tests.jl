module SinusoidalWaveguideBroadbandMatrixTests

using LinearAlgebra
using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "analyze_sinusoidal_waveguide_broadband_matrix.jl",
))
using .AnalyzeSinusoidalWaveguideBroadbandMatrix

@testset "broadband response-matrix analysis" begin
    matrix = reshape(ComplexF64.(1:9), 3, 3)
    projected = reflection_project(matrix)
    @test reverse(reverse(projected; dims=1); dims=2) ≈ projected
    @test reflection_project(projected) ≈ projected

    even, odd = parity_bases(5)
    @test size(even) == (5, 3)
    @test size(odd) == (5, 2)
    @test even' * even ≈ Matrix{Float64}(I, 3, 3)
    @test odd' * odd ≈ Matrix{Float64}(I, 2, 2)
    @test even' * odd ≈ zeros(3, 2) atol=1.0e-14
    @test_throws ArgumentError parity_bases(4)

    frequencies = [0.8, 0.9, 1.0, 1.2]
    weights, quadrature = trapezoid_pulse_weights(frequencies; carrier_hz=1.0)
    @test norm(weights) ≈ 1.0
    @test quadrature ≈ [0.05, 0.10, 0.15, 0.10]

    power = Diagonal([2.0, 3.0])
    omega = 2pi
    work = -im .* (2 / omega) .* power
    whitened = power_whiten(Matrix{ComplexF64}(work), 1.0)
    @test whitened.power ≈ power
    @test whitened.inverse_root' * power * whitened.inverse_root ≈
        Matrix{ComplexF64}(I, 2, 2)
    @test whitened.regularized_count == 0
end

end
