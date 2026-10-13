using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "topology_fresnel_rings_model.jl",
))
using .TopologyFresnelRingsModel

@testset "Fresnel rings analytic design" begin
    design = design_fresnel_rings()
    @test design.config.level_count == 3
    @test design.maximum_delay_s * 1e6 ≈ 1.611 atol=0.002
    @test design.carrier_phase_span_rad < pi
    @test design.full_wrap_count == 0
    @test isempty(halfwave_zone_radii_mm(design.config))
    @test design.first_halfwave_radius_mm > design.config.stand.aperture_height_mm / 2
    @test design.total_device_depth_mm < design.config.stand.maximum_device_depth_mm
    @test all(level -> level.peek_depth_mm >= 0, design.levels)
    @test minimum(
        level.outer_radius_mm - level.inner_radius_mm for level in design.levels
    ) > design.config.minimum_feature_mm
    @test all(
        level -> level.target_delay_s ≈ level.wrapped_delay_s,
        design.levels,
    )
end

@testset "lossless transfer conserves power" begin
    design = design_fresnel_rings()
    level = design.levels[2]
    result = level_transfer(level, design, 242.0e3; loss_scale=0.0)
    @test result.power_transmission + result.power_reflection ≈ 1 atol=1e-10
    @test 0 <= result.power_transmission <= 1
end

@testset "scalar controls are finite" begin
    design = design_fresnel_rings()
    for topology in (:bare, :straight_footprint, :unwrapped_continuous, :fresnel_quantized)
        response = scalar_aperture_transfer(design, 242.0e3, topology; radial_samples=24)
        @test isfinite(real(response))
        @test isfinite(imag(response))
    end
end
