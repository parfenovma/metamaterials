using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "topology_bending_metal_grin_mesher.jl",
))
using .TopologyBendingMetalGRINMesher

@testset "bending metal GRIN parameterization" begin
    config = BendingMetalGRINConfig()
    @test validate_config(config) === nothing
    @test total_length_mm(config) == 55.0
    @test length(BENDING_METAL_GRIN_VARIANTS) == 11
    @test variant_parameters(config, :dense_c).ligament_mm == 2.70
    @test variant_parameters(config, :bd_a).ligament_mm == 1.70
    @test variant_parameters(config, :bd_g).ligament_mm == 0.80
    @test variant_parameters(config, :bd_g).mass_height_mm >
          variant_parameters(config, :bd_a).mass_height_mm
    @test minimum_feature_mm(config, :bd_g) == 0.8
    @test all(
        minimum_feature_mm(config, variant) >= config.nominal_minimum_feature_mm
        for variant in BENDING_METAL_GRIN_VARIANTS
    )
    @test_throws ArgumentError variant_parameters(config, :unknown)
    @test validate_config(BendingMetalGRINConfig(period_count=1)) === nothing
    @test_throws ArgumentError validate_config(BendingMetalGRINConfig(period_count=0))
    @test_throws ArgumentError validate_config(BendingMetalGRINConfig(
        boundary_node_width_mm=0.7,
    ))
end
