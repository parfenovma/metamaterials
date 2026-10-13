using Test

include(joinpath(@__DIR__, "topology_composite_rings_model.jl"))
using .TopologyCompositeRingsModel

@testset "composite-rings measured-property model" begin
    anchors = measured_dp190_w_epoxy()
    @test length(anchors) == 4
    @test interpolate_composite(0.0992).pressure_wave_speed_m_s == 1982.0
    midpoint = interpolate_composite((0.0992 + 0.2022) / 2)
    @test 1698.0 < midpoint.pressure_wave_speed_m_s < 1982.0
    @test 5.9 < midpoint.impedance_mrayl < 8.5
    @test_throws ArgumentError interpolate_composite(0.5)
end

@testset "seven manufacturable concentric rings" begin
    config = CompositeRingsConfig()
    zones = design_ring_zones(config)
    @test length(zones) == 7
    @test all(diff(getproperty.(zones, :target_delay_s)) .< 0)
    @test all(0.0992 .<= getproperty.(getproperty.(zones, :property), :tungsten_volume_fraction) .<= 0.2022)
    @test minimum(getproperty.(zones, :outer_radius_mm) .- getproperty.(zones, :inner_radius_mm)) > 5.0
    @test config.core_thickness_mm + 2config.face_thickness_mm <= 55.0
end

@testset "passive layer transfer and bond sensitivity" begin
    config = CompositeRingsConfig(aperture_samples=101)
    zones = design_ring_zones(config)
    lossless = LossCase(label="lossless", core_db_per_cm_at_f0=0.0)
    perfect = BondCase(label="perfect", thickness_mm=0.0)
    response = zone_transfer(zones[1].property, config.center_frequency_hz, config, lossless, perfect)
    @test response.power_transmission <= 1 + 1.0e-10
    @test abs2(response.reflection) + response.power_transmission ≈ 1 atol=1.0e-9

    lossy = LossCase(label="lossy", core_db_per_cm_at_f0=2.0)
    lossy_response = zone_transfer(zones[1].property, config.center_frequency_hz, config, lossy, perfect)
    @test lossy_response.power_transmission < response.power_transmission
end

println("composite-rings topology tests passed")
