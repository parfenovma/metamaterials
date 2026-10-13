using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "topology_hybrid_monostatic_mesher.jl",
))
using .TopologyHybridMonostaticMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "topology_hybrid_monostatic_transient.jl",
))
using .TopologyHybridMonostaticTransient

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "topology_hybrid_monostatic_analysis.jl",
))
using .TopologyHybridMonostaticAnalysis

@testset "hybrid monostatic G3 geometry" begin
    straight = HybridMonostaticMeshConfig(topology=:straight_footprint)
    pure = HybridMonostaticMeshConfig(topology=:pure_ttd)
    hybrid = HybridMonostaticMeshConfig(topology=:hybrid)

    straight_center = topology_channel_geometry(0.0, straight)
    pure_center = topology_channel_geometry(0.0, pure)
    hybrid_center = topology_channel_geometry(0.0, hybrid)
    @test straight_center.total_extra_path_mm == 0
    @test pure_center.relative_extra_path_mm > hybrid_center.relative_extra_path_mm > 0
    @test pure_center.total_extra_path_mm > hybrid_center.total_extra_path_mm
    @test pure_center.amplitude_mm > hybrid_center.amplitude_mm
    @test hybrid_center.insert.total_length_mm > 0

    edge = first(TopologyHybridMonostaticMesher.TopologyHybridTTDGRINDesign.channel_centers_mm())
    @test topology_channel_geometry(edge, hybrid).delay_fraction ≈ 0 atol=1e-12
    @test topology_channel_geometry(edge, pure).relative_extra_path_mm ≈ 0 atol=1e-12

    # Equal-energy pressure normalization: p^2*aperture is invariant.
    bare_length = straight.stand.aperture_height_mm
    device_length = straight.design.channel_count * straight.design.guide_width_mm
    @test source_pressure_scale(bare_length, straight.stand)^2 * bare_length ≈
          source_pressure_scale(device_length, straight.stand)^2 * device_length
end

@testset "hybrid monostatic waveform and receiver metrics" begin
    transient = HybridMonostaticTransientConfig()
    duration = transient.pulse_cycles / transient.frequency_hz
    @test hann_burst(-1e-6, transient) == 0
    @test hann_burst(duration + 1e-6, transient) == 0
    @test abs(hann_burst(duration / 2, transient)) <= 1
    @test nominal_incident_energy_j_per_m(0.0716, 1.0, transient) > 0

    dt = 1 / (18 * 242e3)
    analysis = AnalysisConfig()
    template = TopologyHybridMonostaticAnalysis.pulse_template(dt, analysis)
    time = collect(0.0:dt:80e-6)
    signal = zeros(length(time))
    start = 80
    signal[start:(start + length(template.values) - 1)] .= 2 .* template.values
    trace = matched_filter_trace(time, signal, template)
    expected = time[start] + template.duration / 2
    @test trace.time_s[argmax(trace.values)] ≈ expected atol=dt
    @test maximum(trace.values) > 0
end

@testset "hybrid monostatic exact result contract" begin
    result_root = joinpath(
        @__DIR__, "..", "results", "topology_benchmark", "hybrid_monostatic_g3",
    )
    for name in ("g3_metrics.csv", "convergence_metrics.csv", "common_metrics.csv", "decision.csv")
        @test isfile(joinpath(result_root, name))
    end
    decision = readlines(joinpath(result_root, "decision.csv"))
    @test length(decision) == 2
    @test startswith(decision[1], "verdict,reason")
    @test startswith(decision[2], "conditional,")
end
