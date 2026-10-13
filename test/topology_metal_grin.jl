using Test

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "topology_metal_grin_mesher.jl",
))
using .TopologyMetalGRINMesher

@testset "metal GRIN waisted-ring parameterization" begin
    config = MetalGRINCellConfig()
    @test validate_config(config) === nothing
    @test total_length_mm(config) == 55.0
    @test config.active_length_mm == config.period_count * config.period_mm
    @test neck_height_mm(config, :solid) == config.port_height_mm
    @test neck_height_mm(config, :neck_0p8) == 0.8
    @test minimum_feature_mm(config, :neck_0p8) == 0.8
    @test local_height_mm(config, :neck_1p6, 0.0) == config.port_height_mm
    @test local_height_mm(config, :neck_1p6, config.lead_length_mm) == config.port_height_mm
    @test local_height_mm(
        config,
        :neck_1p6,
        config.lead_length_mm + config.period_mm / 2,
    ) ≈ 1.6 atol=1.0e-12
    @test local_height_mm(
        config,
        :neck_1p6,
        config.lead_length_mm + config.period_mm,
    ) ≈ config.port_height_mm atol=1.0e-12
    @test minimum_curvature_radius_mm(config, :solid) == Inf
    @test minimum_curvature_radius_mm(config, :neck_0p8) > 0.16
    @test all(
        minimum_feature_mm(config, variant) >= config.nominal_minimum_feature_mm
        for variant in METAL_GRIN_VARIANTS
    )
    @test_throws ArgumentError neck_height_mm(config, :unknown)
    @test_throws ArgumentError validate_config(MetalGRINCellConfig(period_count=9))
    @test_throws ArgumentError validate_config(MetalGRINCellConfig(
        nominal_minimum_feature_mm=0.9,
    ))
    @test_throws ArgumentError validate_config(MetalGRINCellConfig(
        transition_fraction=0.5,
    ))
end
