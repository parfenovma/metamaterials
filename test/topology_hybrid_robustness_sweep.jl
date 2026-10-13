using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "topology_hybrid_robustness_sweep.jl",
))
using .TopologyHybridRobustnessSweep

@testset "hybrid robustness transfer" begin
    frequency = [-242e3, 0.0, 242e3]
    identity_transfer = aperture_error_transfer(frequency, RobustnessPoint())
    @test identity_transfer ≈ ones(ComplexF64, 3) atol=1e-12
    @test identity_transfer[1] ≈ conj(identity_transfer[3])

    loss = aperture_error_transfer(
        frequency,
        RobustnessPoint(one_way_excess_loss_db=0.5),
    )
    @test abs(loss[3]) ≈ 10.0^(-1.0 / 20) rtol=1e-12

    conversion = aperture_error_transfer(
        frequency,
        RobustnessPoint(absolute_mode_conversion_fraction=0.15),
    )
    @test abs(conversion[3]) < 1

    delayed = aperture_error_transfer(
        frequency,
        RobustnessPoint(maximum_profile_delay_error_ns=500.0),
    )
    @test abs(delayed[3]) < 1
    @test delayed[1] ≈ conj(delayed[3]) atol=1e-12
end

@testset "hybrid exact trace robustness metrics" begin
    nominal = evaluate_point(RobustnessPoint())
    lossy = evaluate_point(RobustnessPoint(one_way_excess_loss_db=0.6))
    @test nominal.hybrid_to_pure_advantage_db > 0
    @test nominal.echo_to_ringing_db > 0
    @test lossy.hybrid_to_pure_advantage_db < nominal.hybrid_to_pure_advantage_db
    @test lossy.echo_mf < nominal.echo_mf
    @test nominal.pulse_correlation > 0.999
end
