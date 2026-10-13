using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "topology_hybrid_ttd_grin_design.jl",
))
using .TopologyHybridTTDGRINDesign

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "topology_hybrid_ttd_grin_transfer.jl",
))
using .TopologyHybridTTDGRINTransfer

@testset "hybrid TTD--GRIN design" begin
    config = HybridTTDGRINConfig()
    summary = design_summary(config)
    @test validate(config) === config
    @test length(channel_centers_mm(config)) == 15
    @test first(channel_centers_mm(config)) ≈ -34.2 atol=1e-12
    @test last(channel_centers_mm(config)) ≈ 34.2 atol=1e-12
    @test maximum_delay_s(config) * 1e6 ≈ 1.48030524 rtol=1e-7
    @test summary.extra_path_reduction_fraction ≈ 0.40
    @test summary.hybrid_extra_path_mm / config.calibrated_path_delay_speed_m_s * 1e-3 ≈
          config.path_delay_fraction * maximum_delay_s(config) rtol=1e-9
    @test summary.hybrid_inner_radius_mm >= config.minimum_inner_radius_mm
    @test summary.minimum_insert_width_mm >= config.minimum_feature_mm
    @test summary.hybrid_manufacturable
    @test summary.hybrid_amplitude_mm < summary.pure_curvature_compliant_amplitude_mm

    zero_insert = slow_insert_profile(0.0, config)
    @test isempty(zero_insert.segments)
    @test zero_insert.total_length_mm == 0.0

    insert = slow_insert_profile(maximum_delay_s(config), config)
    @test insert.maximum_copper_fraction == 1.0
    @test insert.plateau_length_mm > 0
    al = guided_material(0.0, config)
    realized_model_delay = sum(insert.segments) do segment
        material = guided_material(segment.copper_fraction, config)
        segment.length_mm * 1e-3 *
        (1 / material.guided_speed_m_s - 1 / al.guided_speed_m_s)
    end
    expected_physical_delay = (1 - config.path_delay_fraction) * maximum_delay_s(config)
    @test realized_model_delay * config.slow_delay_calibration_factor ≈
          expected_physical_delay rtol=0.12

    reference_line_impedance = al.line_impedance_per_area * config.guide_width_mm
    for fraction in 0.0:0.1:1.0
        material = guided_material(fraction, config)
        @test material.line_impedance_per_area * guide_width_mm(fraction, config) ≈
              reference_line_impedance rtol=1e-12
    end
end

@testset "hybrid transfer-matrix conservation" begin
    config = HybridTTDGRINConfig()
    for frequency in (162.2e3, 242e3, 321.7e3)
        row = channel_scattering(frequency, maximum_delay_s(config), config)
        @test row.transmission + row.reflection ≈ 1.0 atol=1e-10
        @test row.transmission > 0.999999
    end
    gate = run_transfer_gate(config)
    @test gate.passed
    @test gate.pulse.metrics.pulse_correlation > 0.99
end
