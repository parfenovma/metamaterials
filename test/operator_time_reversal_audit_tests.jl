module OperatorTimeReversalAuditTests

using LinearAlgebra
using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "operator_time_reversal_audit.jl",
))
using .OperatorTimeReversalAudit

@testset "time-reversal modal audit" begin
    weights = normalized_pulse_weights(collect(range(0.8, 1.2; length=7)))
    @test norm(weights) ≈ 1.0
    @test all(>(0), weights)
    @test weights ≈ reverse(weights) rtol=5.0e-3

    repeated = zeros(ComplexF64, 2, 1, 4, 1)
    repeated[:, 1, :, 1] .= reshape([1.0, 0.0], 2, 1, 1)
    q00 = weighted_operator(repeated, [1], [1.0])
    q4 = weighted_operator(repeated, 1:4, [1.0])
    @test time_reversal_kernel(q4) ≈ time_reversal_kernel(q00)

    complementary = zeros(ComplexF64, 2, 1, 2, 1)
    complementary[:, 1, 1, 1] .= [1.0, 0.0]
    complementary[:, 1, 2, 1] .= [0.0, 1.0]
    one_state = weighted_operator(complementary, [1], [1.0])
    two_states = weighted_operator(complementary, 1:2, [1.0])
    @test singular_summary(one_state).stable_rank ≈ 1.0
    @test singular_summary(two_states).stable_rank ≈ 2.0
    @test singular_summary(two_states).observable_modes[20.0] == 2
    @test subspace_novelty(one_state, complementary[:, :, 2, 1], 20.0) ≈ 1.0

    kernel = Matrix{ComplexF64}(I, 3, 3)
    metrics = focus_metrics(kernel, 2; bin_spacing=0.5)
    @test metrics.peak_bin == 2
    @test metrics.target_is_global_maximum
    @test metrics.fwhm_bins == 1
    @test metrics.fwhm == 0.5
    @test maximum(coherence_kernel(kernel) - Matrix{Float64}(I, 3, 3)) <= 1e-14

    @test_throws DimensionMismatch weighted_operator(
        zeros(ComplexF64, 2, 1, 1, 2), [1], [1.0],
    )
    @test_throws ArgumentError observable_mode_count([1.0], 0.0)
end

end
