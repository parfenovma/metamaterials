using Test
using Gridap: VectorValue
using LinearAlgebra: I, ⋅, Diagonal, Hermitian, diag, isposdef

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "ndt_benchmark.jl"))
using .NDTBenchmark

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "ndt_fullwave_harmonic_solver.jl",
))
using .NDTFullWaveHarmonicSolver

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "ndt_fullwave_transient_solver.jl",
))
using .NDTFullWaveTransientSolver: NDTFullWaveTransientConfig, hann_burst

@testset "shared NDT topology benchmark" begin
    config = NDTStandConfig()
    @test validate(config) === config
    @test pressure_wavelength_mm(config) ≈ 25.298 rtol=1.0e-4
    @test shear_wavelength_mm(config) ≈ 12.743 rtol=1.0e-4
    @test burst_duration_s(config) ≈ 5 / 242.0e3
    @test round_trip_time_s(60.0, config) ≈ 2 * 0.060 / config.pressure_wave_speed_m_s
    @test nominal_dead_zone_depth_mm(config) > config.nominal_focus_depth_mm
    @test Symbol.(getproperty.(required_defects(config), :id)) == [:none, :void_r2_focus]
    @test length(all_defects(config)) == 10
    @test :echo_mf_gain_bare in metric_columns()
    @test topology_result_directory("metal_grin") ==
          joinpath("results", "topology_benchmark", "metal_grin")
    @test_throws ArgumentError topology_result_directory("Metal GRIN")
end

@testset "NDT full-wave direct transient" begin
    config = NDTFullWaveTransientConfig()
    @test NDTFullWaveTransientSolver.validate(config) === config
    @test hann_burst(0.0, config.frequency_hz, config.pulse_cycles) == 0.0
    @test hann_burst(config.pulse_cycles / config.frequency_hz, config.frequency_hz, config.pulse_cycles) ≈ 0.0 atol=1.0e-14
    @test hann_burst(6 / config.frequency_hz, config.frequency_hz, config.pulse_cycles) == 0.0
    @test hann_burst(1.25 / config.frequency_hz, config.frequency_hz, config.pulse_cycles) ≈ 0.5 atol=1.0e-12
    @test_throws ArgumentError NDTFullWaveTransientSolver.validate(
        NDTFullWaveTransientConfig(samples_per_cycle=5),
    )
end

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "run_operator_o3_f1c_scenes.jl",
))
using .RunOperatorO3F1CScenes

@testset "O3 F1c structural-scene contract" begin
    contexts = RunOperatorO3F1CScenes.context_rows()
    @test getproperty.(contexts, :id) == [:train_1, :train_2, :validation_1, :test_1]
    @test getproperty.(contexts, :split) == [:train, :train, :validation, :test]
    @test length(unique(getproperty.(contexts, :seed))) == length(contexts)
    for context in contexts
        first_draw = RunOperatorO3F1CScenes.generate_backgrounds(context.seed)
        second_draw = RunOperatorO3F1CScenes.generate_backgrounds(context.seed)
        @test first_draw == second_draw
        @test length(first_draw) == 5
        @test all(defect -> defect.kind == :circular_void, first_draw)
        @test all(defect -> 0.6 <= defect.radius_mm <= 1.1, first_draw)
        @test all(
            defect -> RunOperatorO3F1CScenes.valid_background(
                defect, filter(other -> other.id != defect.id, first_draw),
            ),
            first_draw,
        )
    end
end

@testset "NDT full-wave carrier contract" begin
    config = NDTFullWaveHarmonicConfig(element_count=17)
    @test NDTFullWaveHarmonicSolver.validate(config) === config
    @test_throws ArgumentError NDTFullWaveHarmonicSolver.validate(
        NDTFullWaveHarmonicConfig(pressure_wave_speed_m_s=2000.0, shear_wave_speed_m_s=3000.0),
    )

    include(joinpath(
        @__DIR__, "..", "src", "code", "metamaterial", "run_operator_o3_fullwave_carrier.jl",
    ))
    identity_response = Matrix{ComplexF64}(I, 5, 5)
    a_scan = (
        id=:a, kind=:a_scan, electronic_channels=1, physical_states=1,
        transmit_events=1, indices=Int[], states=(),
    )
    motion = (
        id=:m, kind=:motion, electronic_channels=1, physical_states=1,
        transmit_events=3, indices=[1, 3, 5], states=(),
    )
    @test RunOperatorO3FullWaveCarrier.project_operator(identity_response, a_scan) ≈ [1.0 + 0im]
    @test RunOperatorO3FullWaveCarrier.project_operator(identity_response, motion) ≈
          fill(inv(sqrt(3)) + 0im, 3)
    @test RunOperatorO3FullWaveCarrier.relative_residual(identity_response, identity_response) == 0
end

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "measurement_operator_space.jl",
))
using .MeasurementOperatorSpace

@testset "measurement-operator design space" begin
    architectures = architecture_space()
    @test MeasurementOperatorSpace.validate.(architectures) == architectures
    @test first(architectures).id == :a_scan
    @test last(architectures).id == :paut_oracle
    @test count(architecture -> architecture.role == :candidate, architectures) == 8

    identity_dictionary = Matrix{Float64}(I, 3, 3)
    metrics = operator_metrics(identity_dictionary)
    @test metrics.minimum_detectability ≈ 1.0
    @test metrics.median_detectability ≈ 1.0
    @test metrics.maximum_pair_coherence ≈ 0.0 atol=1.0e-14
    @test metrics.stable_rank ≈ 3.0
    @test metrics.numerical_rank == 3
    @test metrics.information_logdet_nats ≈ 3 * log(2)

    ambiguous_dictionary = [1.0 1.0; 0.0 0.0]
    ambiguous = operator_metrics(ambiguous_dictionary)
    @test ambiguous.maximum_pair_coherence ≈ 1.0
    @test ambiguous.stable_rank ≈ 1.0
    @test ambiguous.numerical_rank == 1

    complex_dictionary = ComplexF64[1 im; im 1]
    noise_covariance = [2.0 0.0; 0.0 2.0]
    whitened = whiten_dictionary(complex_dictionary; noise_covariance)
    @test whitened ≈ complex_dictionary / sqrt(2)
    @test normalized_gram(whitened)[1, 1] ≈ 1.0

    @test_throws ArgumentError normalized_gram([1.0 0.0; 0.0 0.0])
    @test_throws ArgumentError operator_metrics(
        identity_dictionary; noise_covariance=zeros(3, 3),
    )
end

module OperatorO6CCompositeTests
using Test
using LinearAlgebra

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6c_composite_twin.jl",
))

@testset "O6c frozen homogenized composite" begin
    metadata = RunOperatorO6CCompositeTwin.composite_metadata()
    @test metadata.ceramic_volume_fraction == 0.65
    @test metadata.density_kg_m3 == 5461.0
    @test metadata.axial_impedance_rayl == 17.0e6
    @test metadata.positive_definite
    @test metadata.axial_stiffness_pa * metadata.lateral_stiffness_pa >
          metadata.cross_stiffness_pa^2
    @test metadata.nominal_half_wave_frequency_hz ≈ 259415.2475 rtol=1.0e-9
    material = RunOperatorO6CCompositeTwin.composite_material()
    @test material.quality_factor == 35.0
    basis = RunOperatorO6CCompositeTwin.dct_basis()
    @test basis' * basis ≈ Matrix{Float64}(I, 3, 3) atol=1.0e-12
end

end

module OperatorO6DIsolatedBusTests
using Test
using LinearAlgebra

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6d_isolated_bus_bound.jl",
))

@testset "O6d isolated three-bus bound" begin
    target = RunOperatorO6DIsolatedBusBound.dct_basis()
    @test transpose(target) * target ≈ Matrix{Float64}(I, 3, 3) atol=1.0e-12
    identity_rotation = Matrix{Float64}(I, 3, 3)
    one = RunOperatorO6DIsolatedBusBound.sparse_internal_map(
        target, identity_rotation, 1,
    )
    two = RunOperatorO6DIsolatedBusBound.sparse_internal_map(
        target, identity_rotation, 2,
    )
    full = RunOperatorO6DIsolatedBusBound.sparse_internal_map(
        target, identity_rotation, 3,
    )
    @test count(abs.(one.internal) .> 1.0e-12) == 17
    @test count(abs.(two.internal) .> 1.0e-12) <= 34
    @test full.public_transfer ≈ target atol=1.0e-12
    four_level = RunOperatorO6DIsolatedBusBound.quantized_target(target, 4)
    eight_level = RunOperatorO6DIsolatedBusBound.quantized_target(target, 8)
    @test RunOperatorO6DIsolatedBusBound.transfer_metrics(
        four_level, target,
    ).amplitude_rms > 0.03
    @test RunOperatorO6DIsolatedBusBound.transfer_metrics(
        eight_level, target,
    ).amplitude_rms < 0.03
end

end

module OperatorO6EResourceForkTests
using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6e_resource_fork.jl",
))

@testset "O6e resource fork statistics" begin
    module_under_test = RunOperatorO6EResourceFork
    @test module_under_test.wilson_lower_bound(50_000, 50_000) > 0.999
    @test module_under_test.wilson_lower_bound(49_538, 50_000) > 0.99
    @test module_under_test.wilson_lower_bound(49_500, 50_000) < 0.99
    target = module_under_test.RunOperatorO6DIsolatedBusBound.dct_basis()
    common_delay = fill(20.0e-9, size(target))
    @test module_under_test.delay_rms(common_delay, target) ≈ 0.0 atol=1.0e-20
    differential = copy(common_delay)
    differential[1, 1] += 100.0e-9
    @test module_under_test.delay_rms(differential, target) > 0
end

end

module OperatorO6FBeamformerTests
using Test
using LinearAlgebra

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6f_beamformer_network.jl",
))

@testset "O6f reciprocal Givens lattice" begin
    module_under_test = RunOperatorO6FBeamformerNetwork
    seed = Matrix(qr(reshape(sin.(1.0:51.0), 17, 3)).Q[:, 1:3])
    lattice = module_under_test.synthesize_isometry_lattice(seed)
    @test length(lattice.stages) == 45
    @test lattice.isometry_embedding * lattice.core_rotation ≈ seed atol=1e-11
    transfer = module_under_test.realized_transfer(
        lattice.stages, lattice.core_rotation, 242.0e3,
    )
    @test real.(transfer) ≈ seed atol=1e-11
    @test maximum(abs, imag.(transfer)) <= 1e-12
end

end

module OperatorO6FMotionRiskTests
using Test
using LinearAlgebra

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6f_motion_risk.jl",
))

@testset "O6f arbitrary-position motion compressor" begin
    module_under_test = RunOperatorO6FMotionRisk
    config = module_under_test.SurrogateConfig(
        frequency_ratios=[0.9, 1.0, 1.1], element_count=5,
    )
    hypotheses = module_under_test.localization_hypotheses(
        depths=[3.0], lateral_positions=[-0.5, 0.5],
    )
    positions = module_under_test.element_positions(config)
    features = module_under_test.scan_features_at_positions(
        hypotheses, config, positions,
    )
    architecture = module_under_test.motion_architecture(
        :test_motion, collect(1:5), :test,
    )
    expected = module_under_test.architecture_dictionary(
        module_under_test.response_library(hypotheses, config),
        architecture,
        config,
    )
    @test features ./ sqrt(5) ≈ expected rtol=1e-12
    delayed = module_under_test.scan_features_at_positions(
        hypotheses,
        config,
        positions;
        delay_errors_cycles=fill(0.05, 5),
    )
    @test norm(delayed - features) > 0
end

end

module OperatorO6GNetworkBudgetTests
using Test
using LinearAlgebra

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6g_network_budget.jl",
))

@testset "O6g complete signed q8 factorization" begin
    module_under_test = RunOperatorO6GNetworkBudget
    left = Matrix(qr(reshape(cos.(1.0:51.0), 17, 3)).Q[:, 1:3])
    right = Matrix(qr(reshape(sin.(1.0:9.0), 3, 3)).Q)
    singular_values = [1.0, 0.99, 0.97]
    output = module_under_test.F.synthesize_isometry_lattice(left)
    first_input = module_under_test.factor_orthogonal_transform(right')
    second_input = module_under_test.factor_orthogonal_transform(
        output.core_rotation,
    )
    transfer = module_under_test.complete_network_transfer(
        output,
        first_input,
        second_input,
        singular_values,
        242.0e3,
    )
    @test transfer ≈ left * Diagonal(singular_values) * right' atol=1e-11
    @test length(output.stages) + length(first_input.stages) +
          length(second_input.stages) == 51
    lossy = module_under_test.complete_network_transfer(
        output,
        first_input,
        second_input,
        singular_values,
        242.0e3;
        coupler_insertion_loss_db=0.01,
    )
    @test opnorm(lossy) <= 1 + 1e-12
end

end

module OperatorO6HLoadedPhysicsTests
using Test
using LinearAlgebra

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6h_bvd_loading.jl",
))

@testset "O6h reciprocal scattering completion" begin
    module_under_test = RunOperatorO6HBVDLoading
    isometry = Matrix(qr(reshape(cos.(1.0:51.0), 17, 3)).Q[:, 1:3])
    transfer = isometry * Diagonal([1.0, 0.99, 0.97])
    completion = module_under_test.scattering_completion(transfer)
    @test completion.forward' * completion.forward ≈
          Matrix{Float64}(I, 3, 3) atol=1e-11
    @test completion.full_scattering' * completion.full_scattering ≈
          Matrix{Float64}(I, 23, 23) atol=1e-10
    loaded = module_under_test.loaded_transfer(
        completion, 1.0, zeros(17), ones(17), ones(17),
    )
    @test loaded.mechanical ≈ transfer atol=1e-12
    @test norm(loaded.public_reflection) <= 1e-12
    @test module_under_test.bvd_impedance(1.0, 35.0, 1.0, 1.0) == 1 + 0im
end

end

module OperatorO6HContactPhysicsTests
using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6h_angular_contact.jl",
))

@testset "O6h exact layered-contact positive controls" begin
    module_under_test = RunOperatorO6HAngularContact
    contact = module_under_test.C
    frequency = contact.CARRIER_FREQUENCY_HZ
    one_dimensional = contact.layer_power_transmission(
        frequency,
        contact.NOMINAL_THICKNESS_M,
        contact.COUPLANT_IMPEDANCE_RAYL,
        contact.COUPLANT_SPEED_M_S,
    )
    angular = module_under_test.elastic_layer_scattering(
        frequency,
        0.0,
        contact.NOMINAL_THICKNESS_M,
        contact.COUPLANT_IMPEDANCE_RAYL,
        contact.COUPLANT_SPEED_M_S,
    )
    @test angular.energy_balance ≈ 1.0 atol=1e-12
    @test angular.transmitted_power ≈ abs2(one_dimensional) atol=1e-12
    direct = module_under_test.direct_interface_scattering(frequency, 30.0)
    @test direct.energy_balance ≈ 1.0 atol=1e-12
end

end

module OperatorO6ISolidContactTests
using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o6i_mitigation_map.jl",
))

@testset "O6i shear-transmitting solid contact" begin
    module_under_test = RunOperatorO6IMitigationMap
    contact = module_under_test.C
    material = module_under_test.A.solid_material(
        2053.6068837506846, 6122.1024, 1530.5256,
    )
    result = module_under_test.solid_layer_scattering(
        contact.CARRIER_FREQUENCY_HZ,
        45.0,
        505.9588760330579e-6,
        material,
    )
    @test result.energy_balance ≈ 1.0 atol=1e-12
    @test result.transmitted_power > 0.9
    @test 0 <= result.transmitted_s_fraction <= 1
end

end

module OperatorO7GlobalStateTests
using Test
using LinearAlgebra
using Random

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "run_operator_o7_global_states.jl",
))

@testset "O7 stationary two-bit modal boundary" begin
    config = SurrogateConfig(frequency_ratios=[0.8, 1.0, 1.2])
    parameters = o7_random_parameters(Xoshiro(707), :pzt)
    first_states, first_scale, first_raw_sigma = o7_transfer_states(
        :pzt, parameters, config,
    )
    second_states, second_scale, second_raw_sigma = o7_transfer_states(
        :pzt, parameters, config,
    )
    @test first_states == second_states
    @test first_scale == second_scale
    @test first_raw_sigma == second_raw_sigma
    architecture = o7_architecture(
        :test_pzt, :pzt_2bit, first_states, :test; mechanism=:pzt,
    )
    checks = o7_physical_checks(architecture, parameters, config)
    @test checks.passed
    @test checks.maximum_sigma <= 0.98 + 1.0e-12
    @test checks.reflection_error <= 1.0e-10
    @test checks.reciprocity_error <= 1.0e-10
    @test checks.minimum_state_distance > 1.0e-4
    @test checks.zero_control_collapse == 0.0
    @test all(
        imag(o7_feedback(
            :pzt, 1, frequency, parameters.strength_a,
            parameters.resonance_a, parameters.quality_a, parameters.loss_a,
        )) >= 0
        for frequency in config.frequency_ratios
    )

    mechanical_parameters = o7_random_parameters(Xoshiro(708), :mechanical)
    mechanical_states, _, _ = o7_transfer_states(
        :mechanical, mechanical_parameters, config,
    )
    mechanical = o7_architecture(
        :test_mechanical, :mechanical_2bit, mechanical_states, :test;
        mechanism=:mechanical,
    )
    @test o7_physical_checks(mechanical, mechanical_parameters, config).passed

    hypotheses = localization_hypotheses(
        depths=[3.0], lateral_positions=[-0.5, 0.5],
    )
    library = response_library(hypotheses, config)
    dictionary = o7_compress_library(library, architecture, config)
    @test size(dictionary) == (4 * 3^2 * 3, 2)
end

end

module OperatorO8PlaneStrainBoundaryTests
using Test
using LinearAlgebra
using Random

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "operator_o8_boundary_mesher.jl",
))
include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "operator_o8_plane_strain_solver.jl",
))

using .OperatorO8BoundaryMesher
using .OperatorO8PlaneStrainSolver

@testset "O8 frozen plane-strain boundary physics" begin
    geometry = O8BoundaryGeometry(
        aperture_height_m=0.20,
        block_length_m=0.05,
        position_a=0.20,
        width_a=0.10,
        position_b=0.90,
        width_b=0.11,
        mesh_size_m=0.002,
    )
    @test OperatorO8BoundaryMesher.validate(geometry) === geometry
    intervals = control_intervals(geometry)
    @test intervals.a_minus[1] ≈ -intervals.a_plus[2]
    @test intervals.a_minus[2] ≈ -intervals.a_plus[1]
    @test intervals.b_minus[1] ≈ -intervals.b_plus[2]
    @test intervals.b_plus[2] == geometry.aperture_height_m / 2

    pzt = O8PZTBoundary()
    load = O8Load(inductance_h=95e-6, resistance_ohm=25.0)
    @test OperatorO8PlaneStrainSolver.validate(pzt) === pzt
    @test OperatorO8PlaneStrainSolver.validate(load) === load
    omega = 2pi * 242e3
    capacitance = 4e-9
    open_coefficient = OperatorO8PlaneStrainSolver.electrical_coefficient(
        0, omega, capacitance, load,
    )
    loaded_coefficient = OperatorO8PlaneStrainSolver.electrical_coefficient(
        1, omega, capacitance, load,
    )
    @test open_coefficient ≈ inv(capacitance)
    @test imag(loaded_coefficient) > 0

    mechanical = O8MechanicalBoundary(
        stiffness_a_n_m3=1e12,
        stiffness_b_n_m3=2e12,
        loss_a=0.03,
        loss_b=0.10,
    )
    @test OperatorO8PlaneStrainSolver.validate(mechanical) === mechanical
    rng = Xoshiro(808)
    transfer = randn(rng, ComplexF64, 17, 3)
    @test state_projector_distance(transfer, transfer) <= 1e-12
end

end

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "elastic_operator_surrogate.jl",
))
using .ElasticOperatorSurrogate

@testset "dimensionless elastic operator surrogate" begin
    config = SurrogateConfig(frequency_ratios=[0.9, 1.1], element_count=5)
    @test ElasticOperatorSurrogate.validate(config) === config
    @test element_positions(config) == [-4.0, -2.0, 0.0, 2.0, 4.0]

    hypotheses = localization_hypotheses(
        depths=[3.0, 4.0], lateral_positions=[-1.0, 0.0, 1.0],
    )
    @test length(hypotheses) == 6
    @test mirror_hypothesis_index(hypotheses) == [3, 2, 1, 6, 5, 4]

    library = response_library(hypotheses, config)
    @test size(library) == (5, 5, 2, 6)
    @test maximum_reciprocity_error(library) <= 1.0e-12
    @test maximum_mirror_error(library, hypotheses) <= 1.0e-12

    @test sum(abs2, uniform_port(config)) ≈ 1.0
    @test sum(abs2, fixed_focus_port(1.0, 3.5, 0.0, config)) ≈ 1.0
    @test sum(abs2, dispersive_port(1.0, config; symmetry=:even)) ≈ 1.0
    ports = cosine_ports(2, config)
    @test ports' * ports ≈ Matrix{Float64}(I, 2, 2)

    one_port = one_port_dictionary(library, _ -> ComplexF64.(uniform_port(config)), config)
    two_port = multiport_dictionary(library, ports, config)
    full = fmc_dictionary(library, config)
    @test size(one_port) == (2, 6)
    @test size(two_port) == (8, 6)
    @test size(full) == (50, 6)

    compressed_1, energy_1, _ = ideal_compression_dictionary(full, 1)
    compressed_2, energy_2, _ = ideal_compression_dictionary(full, 2)
    @test size(compressed_1) == (1, 6)
    @test size(compressed_2) == (2, 6)
    @test 0 < energy_1 <= energy_2 <= 1

    @test_throws ArgumentError ElasticOperatorSurrogate.validate(
        SurrogateConfig(element_count=4),
    )
end

module CausalOperatorScriptTests
using Test
using LinearAlgebra
using Random

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "optimize_operator_multistate_codebooks.jl",
))

@testset "causal and multistate operator constraints" begin
    frequencies = [0.8, 1.0, 1.2]
    spectrum = hann_burst_spectrum(frequencies; sample_count=2001)
    @test length(spectrum) == 3
    @test spectrum[2] > spectrum[1]
    @test spectrum[2] > spectrum[3]

    covariance = frequency_covariance(frequencies, 2; model=:correlated)
    @test size(covariance) == (6, 6)
    @test isposdef(Hermitian(covariance))

    rng = Xoshiro(17)
    raw_ports = randn(rng, 7, 4)
    ports = causal_orthonormal_columns(raw_ports, 4)
    static_ports, _, static_sigma = frequency_ports(
        ports, 0.8; model=:constant_contractive,
    )
    @test opnorm(static_ports) <= 1 + 1.0e-12
    @test static_sigma <= 1 + 1.0e-12

    complex_ports = ports .* exp.(2pi * im .* rand(rng, 7, 4))
    delay_ports, sigma_before, sigma_after = frequency_ports(
        complex_ports, 0.8; model=:delay_contractive,
    )
    @test sigma_before >= sigma_after
    @test opnorm(delay_ports) <= 1 + 1.0e-12
    @test sigma_after <= 1 + 1.0e-12

    hypotheses = [
        ElasticOperatorSurrogate.DefectHypothesis(:left, 3.0, -1.0, Matrix{Float64}(I, 2, 2)),
        ElasticOperatorSurrogate.DefectHypothesis(:right, 3.0, 1.0, Matrix{Float64}(I, 2, 2)),
    ]
    first_state = ComplexF64[1 0; 0 1; 1 1]
    second_state = ComplexF64[1 1; 1 -1; 0 1]
    base_covariance = Matrix{Float64}(I, 3, 3)
    _, combined, combined_covariance = combined_state_operator(
        [first_state, second_state], base_covariance, hypotheses,
    )
    @test combined == vcat(first_state, second_state) / sqrt(2)
    @test combined_covariance == Matrix{Float64}(I, 6, 6)
    @test vec(sum(abs2, combined; dims=1)) ≈
          (vec(sum(abs2, first_state; dims=1)) +
           vec(sum(abs2, second_state; dims=1))) / 2
end

end

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "operator_o3_statistics.jl",
))
using .OperatorO3Statistics

@testset "O3 paired-clutter decision statistics" begin
    dictionary = ComplexF64[1 0; 0 1; 1 1]
    empty_clutter = zeros(ComplexF64, 3, 0)
    @test precision_dictionary(dictionary, empty_clutter) == dictionary

    clutter = ComplexF64[1; 0; 0;;]
    precision_applied = precision_dictionary(dictionary, clutter)
    @test precision_applied ≈ Diagonal([0.5, 1.0, 1.0]) * dictionary

    output_covariance = matched_output_covariance(precision_applied, clutter)
    @test size(output_covariance) == (2, 2)
    @test isposdef(Hermitian(
        output_covariance + 1.0e-9 * Matrix{Float64}(I, 2, 2),
    ))
    factor = OperatorO3Statistics.covariance_factor(output_covariance)
    @test size(factor) == (2, 2)

    means = matched_mean(precision_applied, dictionary)
    scores, indices = normalized_scores(means, max.(real.(diag(means)), eps()))
    @test length(scores) == 2
    @test indices == [1, 2]

    @test 0.01 < wilson_upper_bound(100, 10_000) < 0.0125
    @test crossing_amplitude([1.0, 2.0, 4.0], [0.1, 0.8, 0.95], 0.9) ≈
          exp(log(2.0) + (0.1 / 0.15) * (log(4.0) - log(2.0)))
    @test crossing_amplitude([1.0, 2.0], [0.1, 0.2], 0.9) == Inf
end

include(joinpath(@__DIR__, "topology_fresnel_rings.jl"))

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "profiles.jl"))
using .MetamaterialProfiles

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "step2_solver.jl"))

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "side_mass_mesher.jl"))
using .SideMassMesher

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "huygens_pair_mesher.jl"))
using .HuygensPairMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "run_huygens_pair_pilot.jl",
))

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "balanced_huygens_mesher.jl",
))

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "run_balanced_huygens_pilot.jl",
))

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "modal_solver.jl"))
using .ConservativeElasticModes

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "port_mode_solver.jl"))
using .ElasticPortModes

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "modal_projection.jl"))
using .ElasticModalProjection

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "modal_harmonic_solver.jl"))

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "uniform_reference_mesher.jl",
))

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "baseline_r0_mesher.jl"))

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "inline_mass_component_mesher.jl",
))
using .InlineMassComponentMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "stiffness_component_mesher.jl",
))
using .StiffnessComponentMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "combined_kmk_mesher.jl",
))
using .CombinedKMKMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "distributed_slow_wave_mesher.jl",
))
using .DistributedSlowWaveMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "tapered_shunt_trim_mesher.jl",
))
using .TaperedShuntTrimMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "perforated_channel_mesher.jl",
))
using .PerforatedChannelMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "folded_path_mesher.jl",
))
using .FoldedPathMesher

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "spectral_analysis.jl"))
using .SpectralAnalysis

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "wavelet_analysis.jl"))
using .WaveletAnalysis

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "impulse_risk_analysis.jl",
))
using .ImpulseRiskAnalysis

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "passive_feed_splitter_mesher.jl",
))
import .PassiveFeedSplitterMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "passive_feed_manifold_mesher.jl",
))
import .PassiveFeedManifoldMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "design_aluminium_horn_ttd_passive_feed.jl",
))
import .DesignAluminiumHornTTDPassiveFeed

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "run_pilot_sweep.jl"))
include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "run_working_band_boundary_study.jl",
))
include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "run_period_length_pilot.jl",
))

include(joinpath(@__DIR__, "..", "src", "code", "metamaterial", "run_lens_prototype.jl"))
include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "run_gap_pilot.jl",
))
include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "run_ideal_huygens_aperture.jl",
))
using .LensDesign

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "build_target_frequency_library.jl",
))
include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "analyze_phase_efficiency.jl",
))

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "sinusoidal_material_lens.jl",
))
using .SinusoidalMaterialLens

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "material_lens_mesher.jl",
))
using .MaterialLensMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "horn_point_radiator_mesher.jl",
))
using .HornPointRadiatorMesher

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "monotonic_horn_lens.jl",
))
using .MonotonicHornLens

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "horn_monotonic_array_3d_mesher.jl",
))

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "horn_neighbor_3d_mesher.jl",
))

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "run_horn_monotonic_array_3d_convergence.jl",
))

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "run_horn_monotonic_array_3d_half_symmetry.jl",
))

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "run_horn_monotonic_response_matrix.jl",
))

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "run_horn_monotonic_group_delay.jl",
))

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "analyze_horn_monotonic_v1_carrier.jl",
))

include(joinpath(
    @__DIR__,
    "..",
    "src",
    "code",
    "metamaterial",
    "analyze_horn_monotonic_broadband_coupling.jl",
))

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "point_radiator_aperture_mesher.jl",
))
using .PointRadiatorApertureMesher

@testset "sinusoidal material lens reduced model" begin
    polymer = photopolymer()
    aluminium = aluminium_6061()
    matching = geometric_matching_material(polymer, aluminium)
    frequency = 242.0e3
    straight = SinusoidalCell(length_mm=17.0)

    @test cell_gap_mm(straight) == 7.0
    @test 3.0 < quarter_wave_thickness_mm(matching, frequency) < 5.0

    homogeneous = cell_transfer(straight, frequency, aluminium, aluminium)
    @test abs(homogeneous.stress) ≈ 1.0 atol=1.0e-10
    @test abs(homogeneous.velocity) ≈ 1.0 atol=1.0e-10

    direct = cell_transfer(straight, frequency, polymer, aluminium)
    matched = cell_transfer(
        straight,
        frequency,
        polymer,
        aluminium;
        matching_material=matching,
        matching_thickness_mm=quarter_wave_thickness_mm(matching, frequency),
    )
    @test abs(matched.velocity) > abs(direct.velocity)
end

@testset "material lens defect mesh configuration" begin
    cells = fill(SinusoidalCell(length_mm=28.0), 3)
    config = MaterialLensMeshConfig(
        defect_radius_mm=2.0,
        defect_center_x_mm=40.0,
        defect_center_y_mm=0.0,
    )
    @test MaterialLensMesher.validate(cells, 0.0, config) === nothing
    @test_throws ArgumentError MaterialLensMesher.validate(
        cells,
        0.0,
        MaterialLensMeshConfig(defect_radius_mm=-1.0),
    )
end

@testset "rounded-notch lens geometry" begin
    cell = RoundedNotchCell(
        length_mm=20.0,
        minimum_gap_mm=3.0,
        notch_count=4,
        notch_width_mm=1.6,
    )
    @test notch_centers_mm(cell) ≈ [3.4, 7.8, 12.2, 16.6]
    @test SinusoidalMaterialLens.gap_at(cell, 3.4) ≈ 3.0
    @test SinusoidalMaterialLens.gap_at(cell, 5.0) ≈ 7.0
    @test cell_gap_mm(cell) == 3.0
    @test_throws ArgumentError SinusoidalMaterialLens.validate(RoundedNotchCell(
        length_mm=8.0,
        notch_count=5,
        notch_width_mm=2.0,
    ))
end

@testset "horn-fed point-radiator geometry" begin
    base = HornPointRadiatorConfig()
    gentle = HornPointRadiatorConfig(rounded_axial_length_mm=26.0)
    long_smooth = HornPointRadiatorConfig(
        straight_guide_length_mm=42.578058,
        rounded_axial_length_mm=38.158058,
    )
    @test horn_height_mm(base, 0.0) ≈ 3.2
    @test horn_height_mm(base, base.horn_length_mm) ≈ 1.6
    @test adiabatic_parameter(base) <= 0.05
    @test rounded_guide_length_mm(base, rounded_guide_amplitude_mm(base)) ≈
          base.straight_guide_length_mm atol=1e-8
    @test HornPointRadiatorMesher.minimum_rounded_inner_radius_mm(gentle) > 3.9
    @test outlet_x_mm(gentle, :collector_gentle) == 51.0
    @test :collector_smooth in HORN_POINT_VARIANTS
    @test smooth_guide_length_mm(long_smooth, smooth_guide_amplitude_mm(long_smooth)) ≈
          long_smooth.straight_guide_length_mm atol=1e-8
    @test HornPointRadiatorMesher.smooth_centerline_slope(
        long_smooth, 0.0, smooth_guide_amplitude_mm(long_smooth),
    ) ≈ 0.0 atol=1e-12
    @test HornPointRadiatorMesher.smooth_centerline_slope(
        long_smooth, long_smooth.rounded_axial_length_mm,
        smooth_guide_amplitude_mm(long_smooth),
    ) ≈ 0.0 atol=1e-12
    aluminium_maximum = HornPointRadiatorConfig(
        input_height_mm=7.0,
        throat_height_mm=3.5,
        horn_length_mm=65.41,
        lower_band_frequency_hz=162.2e3,
        pressure_wave_speed_m_s=6122.102437409232,
        straight_guide_length_mm=155.53,
        rounded_axial_length_mm=132.5,
        receiver_length_mm=85.0,
        receiver_half_height_mm=45.0,
        target_distance_mm=60.0,
        minimum_inner_radius_mm=1.0,
    )
    @test adiabatic_parameter(aluminium_maximum) <= 0.05
    @test smooth_guide_length_mm(
        aluminium_maximum, smooth_guide_amplitude_mm(aluminium_maximum),
    ) ≈ aluminium_maximum.straight_guide_length_mm atol=1e-8
    @test minimum_smooth_inner_radius_mm(aluminium_maximum) >= 1.0
    @test isnothing(HornPointRadiatorMesher.validate_config(
        aluminium_maximum; variant=:collector_smooth,
    ))
    aluminium_diffuser = HornPointRadiatorConfig(
        input_height_mm=7.0,
        throat_height_mm=3.5,
        horn_length_mm=65.41,
        lower_band_frequency_hz=162.2e3,
        pressure_wave_speed_m_s=6122.102437409232,
        straight_guide_length_mm=155.53,
        rounded_axial_length_mm=132.5,
        receiver_length_mm=85.0,
        receiver_half_height_mm=45.0,
        target_distance_mm=60.0,
        diffuser_output_height_mm=7.0,
        diffuser_length_mm=65.41,
    )
    @test diffuser_height_mm(aluminium_diffuser, 0.0) ≈ 3.5
    @test diffuser_height_mm(aluminium_diffuser, 65.41) ≈ 7.0
    @test diffuser_adiabatic_parameter(aluminium_diffuser) <= 0.05
    @test isnothing(HornPointRadiatorMesher.validate_config(
        aluminium_diffuser; variant=:collector_smooth_diffuser,
    ))
    aluminium_input_expander = HornPointRadiatorConfig(
        input_height_mm=1.4,
        throat_height_mm=3.5,
        horn_length_mm=112.2,
        lower_band_frequency_hz=162.2e3,
        pressure_wave_speed_m_s=6122.102437409232,
        straight_guide_length_mm=155.53,
        rounded_axial_length_mm=132.5,
        receiver_length_mm=85.0,
        receiver_half_height_mm=45.0,
        target_distance_mm=60.0,
        diffuser_output_height_mm=7.0,
        diffuser_length_mm=65.41,
    )
    @test horn_height_mm(aluminium_input_expander, 0.0) ≈ 1.4
    @test horn_height_mm(aluminium_input_expander, 112.2) ≈ 3.5
    @test adiabatic_parameter(aluminium_input_expander) <= 0.05
    @test isnothing(HornPointRadiatorMesher.validate_config(
        aluminium_input_expander; variant=:collector_straight_diffuser,
    ))
    aluminium_output_contractor = HornPointRadiatorConfig(
        input_height_mm=7.0,
        throat_height_mm=3.5,
        horn_length_mm=65.41,
        lower_band_frequency_hz=162.2e3,
        pressure_wave_speed_m_s=6122.102437409232,
        straight_guide_length_mm=155.53,
        rounded_axial_length_mm=132.5,
        receiver_length_mm=85.0,
        receiver_half_height_mm=45.0,
        target_distance_mm=60.0,
        diffuser_output_height_mm=1.4,
        diffuser_length_mm=112.2,
    )
    @test diffuser_height_mm(aluminium_output_contractor, 0.0) ≈ 3.5
    @test diffuser_height_mm(aluminium_output_contractor, 112.2) ≈ 1.4
    @test diffuser_adiabatic_parameter(aluminium_output_contractor) <= 0.05
    @test isnothing(HornPointRadiatorMesher.validate_config(
        aluminium_output_contractor; variant=:collector_straight_diffuser,
    ))
    @test outlet_x_mm(aluminium_diffuser, :collector_smooth_diffuser) ≈ 263.32
    probes = probe_points_mm(gentle, :collector_gentle)
    @test "target_axis" in probes.names
    @test length(probes.names) == length(probes.x_mm) == length(probes.y_mm)

    aperture = PointRadiatorApertureConfig()
    @test radiator_centers_mm(aperture) == collect(-7:7) .* 4.8
    @test maximum(abs, radiator_centers_mm(aperture)) == 33.6

    neighbor = HornNeighbor3DMesher.HornNeighbor3DConfig(
        device_guide_length_mm=38.158058,
        gentle_guide_axial_length_mm=38.158058,
        gentle_guide_path_length_mm=42.578058,
    )
    @test HornNeighbor3DMesher.channel_centers_mm(neighbor, 3) == [-4.8, 0.0, 4.8]
    @test HornNeighbor3DMesher.smooth_guide_amplitude_mm(neighbor) ≈ 7.797432875 atol=1e-8
    @test HornNeighbor3DMesher.outlet_x_mm(neighbor) ≈ 63.158058
    @test isnothing(HornNeighbor3DMesher.validate_config(neighbor, [:device, :smooth, :device]))
    neighbor_probes = HornNeighbor3DMesher.probe_points_mm(
        neighbor, [:device, :smooth, :device],
    )
    @test length(neighbor_probes.names) == 7
    @test neighbor_probes.y_mm[1:2:end-1] == [-4.8, 0.0, 4.8]
end

@testset "passive aluminium feed geometry" begin
    splitter = PassiveFeedSplitterMesher.PassiveFeedSplitterConfig()
    @test isnothing(PassiveFeedSplitterMesher.validate_config(splitter))
    @test PassiveFeedSplitterMesher.small_junction_width_mm(splitter) ≈ 1.011222055969761
    @test PassiveFeedSplitterMesher.large_junction_width_mm(splitter) ≈ 5.9887779440302396
    @test PassiveFeedSplitterMesher.target_output_amplitude_ratio(splitter) ≈ 0.4109170499805894
    @test PassiveFeedSplitterMesher.adiabatic_parameter(
        splitter, PassiveFeedSplitterMesher.small_junction_width_mm(splitter),
    ) <= 0.051
    @test PassiveFeedSplitterMesher.output_x_mm(splitter) ≈ 171.00271326349065

    calibrated = PassiveFeedSplitterMesher.PassiveFeedSplitterConfig(
        small_power_fraction=0.1955777738982686,
        transformer_length_mm=144.4780677722518,
    )
    @test isnothing(PassiveFeedSplitterMesher.validate_config(calibrated))
    @test PassiveFeedSplitterMesher.small_junction_width_mm(calibrated) ≈
          1.3690444172878802

    rounded = PassiveFeedSplitterMesher.PassiveFeedSplitterConfig(
        small_power_fraction=0.499999,
        transformer_length_mm=131.0,
        output_center_offset_mm=4.5,
        splitter_tip_radius_mm=0.75,
        profile_samples=600,
    )
    @test isnothing(PassiveFeedSplitterMesher.validate_config(rounded))

    manifold = PassiveFeedManifoldMesher.PassiveFeedManifoldConfig()
    @test isnothing(PassiveFeedManifoldMesher.validate_config(manifold))
    @test PassiveFeedManifoldMesher.aperture_width_mm(manifold) ≈ 121.8
    @test PassiveFeedManifoldMesher.required_transformer_length_mm(manifold) ≈
          193.4127058472977
    @test PassiveFeedManifoldMesher.output_center_y_mm(manifold, 1) ≈ -57.4
    @test PassiveFeedManifoldMesher.output_center_y_mm(manifold, 8) ≈ 0.0 atol=1e-12
    @test PassiveFeedManifoldMesher.output_center_y_mm(manifold, 15) ≈ 57.4
end

@testset "passive feed topology ceilings" begin
    weights = [
        0.11593941601987037, 0.4500315894067127, 0.3990175013586375,
        0.7163447381726632, 0.7008448451782485, 0.7281942059110046,
        0.8212932885576412, 1.0, 0.8212932885576412,
        0.7281942059110046, 0.7008448451782485, 0.7163447381726632,
        0.3990175013586375, 0.4500315894067127, 0.11593941601987037,
    ]
    eta_equal = 0.9420469907840849
    eta_other = 0.9686039905811457
    _, unrestricted_nodes, unrestricted =
        DesignAluminiumHornTTDPassiveFeed.efficiency_optimal_tree(
            abs2.(weights), eta_equal, eta_other,
        )
    _, contiguous_nodes, contiguous =
        DesignAluminiumHornTTDPassiveFeed.efficiency_optimal_tree(
            abs2.(weights), eta_equal, eta_other; contiguous=true,
        )

    @test length(unrestricted_nodes) == length(weights) - 1
    @test length(contiguous_nodes) == length(weights) - 1
    @test unrestricted ≈ 0.8915994814751932 atol=1e-12
    @test contiguous ≈ 0.8905605108729353 atol=1e-12
    @test unrestricted >= contiguous
    @test_throws ArgumentError begin
        DesignAluminiumHornTTDPassiveFeed.efficiency_optimal_tree(
            [1.0, 0.0], eta_equal, eta_other,
        )
    end
end

@testset "monotonic horn lens geometry" begin
    config = MonotonicHornLensConfig()
    centers = MonotonicHornLens.aperture_centers_mm(config)
    delays = geometric_delays_s(config)
    fractions = normalized_delays(config)
    @test centers == collect(-7:7) .* 4.8
    @test delays == reverse(delays)
    @test fractions[1] == fractions[end] == 0.0
    @test fractions[8] == 1.0
    @test issorted(fractions[1:8])
    @test maximum(delays) * 1e6 ≈ 5.778 atol=0.01

    geometry = synthesize_monotonic_geometry(10.0; config)
    @test geometry.extra_path_mm == reverse(geometry.extra_path_mm)
    @test geometry.amplitude_mm == reverse(geometry.amplitude_mm)
    @test geometry.path_length_mm .- geometry.axial_length_mm ≈ geometry.extra_path_mm
    @test minimum(geometry.inner_radius_mm) >= config.minimum_inner_radius_mm - 1e-8
    @test smooth_arc_length_mm(
        geometry.axial_length_mm,
        maximum(geometry.amplitude_mm);
        samples=config.arc_samples,
    ) ≈ geometry.axial_length_mm + 10.0 atol=1e-8
    @test smooth_slope(0.0, geometry.axial_length_mm, maximum(geometry.amplitude_mm)) ≈ 0.0
    @test smooth_slope(
        geometry.axial_length_mm,
        geometry.axial_length_mm,
        maximum(geometry.amplitude_mm),
    ) ≈ 0.0 atol=1e-12
end

@testset "full monotonic 3D array geometry" begin
    geometry = synthesize_monotonic_geometry(9.724)
    config = HornMonotonicArray3DMesher.HornMonotonicArray3DConfig(
        guide_axial_length_mm=geometry.axial_length_mm,
        bend_amplitude_mm=geometry.amplitude_mm,
    )
    @test isnothing(HornMonotonicArray3DMesher.validate_config(config))
    @test HornMonotonicArray3DMesher.channel_centers_mm(config) == collect(-7:7) .* 4.8
    @test HornMonotonicArray3DMesher.outlet_x_mm(config) ≈ 75.6639442088 atol=1e-8
    @test HornMonotonicArray3DMesher.focus_x_mm(config) ≈ 110.6639442088 atol=1e-8
    @test HornMonotonicArray3DMesher.channel_center_z_mm(
        config, 8, config.guide_axial_length_mm / 2,
    ) ≈ maximum(config.bend_amplitude_mm)
    @test HornMonotonicArray3DMesher.channel_center_z_mm(
        config, 1, config.guide_axial_length_mm / 2,
    ) ≈ 0.0
    probes = HornMonotonicArray3DMesher.probe_points_mm(config)
    @test length(probes.names) == 31
    @test probes.names[end] == "focus"
    uniform = HornMonotonicArray3DMesher.HornMonotonicArray3DConfig(
        guide_axial_length_mm=geometry.axial_length_mm,
        bend_amplitude_mm=zeros(15),
    )
    @test isnothing(HornMonotonicArray3DMesher.validate_config(uniform))

    aluminium_small = HornMonotonicArray3DMesher.HornMonotonicArray3DConfig(
        input_height_mm=7.0,
        throat_height_mm=3.5,
        channel_depth_mm=7.0,
        channel_pitch_mm=28.7,
        horn_length_mm=65.41,
        guide_axial_length_mm=132.5,
        bend_amplitude_mm=[0.0, 24.66569557952898, 29.60200249816748,
                           24.66569557952898, 0.0],
        receiver_length_mm=85.0,
        receiver_half_width_mm=64.0,
        receiver_half_height_mm=6.0,
        focal_distance_mm=60.0,
    )
    @test isnothing(HornMonotonicArray3DMesher.validate_config(aluminium_small))
    @test HornMonotonicArray3DMesher.channel_centers_mm(aluminium_small) ≈
          [-57.4, -28.7, 0.0, 28.7, 57.4]
    aluminium_small_diffuser = HornMonotonicArray3DMesher.HornMonotonicArray3DConfig(
        input_height_mm=7.0,
        throat_height_mm=3.5,
        channel_depth_mm=7.0,
        channel_pitch_mm=28.7,
        horn_length_mm=65.41,
        guide_axial_length_mm=132.5,
        bend_amplitude_mm=aluminium_small.bend_amplitude_mm,
        diffuser_output_height_mm=7.0,
        diffuser_length_mm=65.41,
        receiver_length_mm=85.0,
        receiver_half_width_mm=64.0,
        receiver_half_height_mm=6.0,
        focal_distance_mm=60.0,
    )
    @test isnothing(HornMonotonicArray3DMesher.validate_config(aluminium_small_diffuser))
    @test HornMonotonicArray3DMesher.outlet_x_mm(aluminium_small_diffuser) ≈ 263.32
end

@testset "full monotonic 3D array convergence analysis" begin
    using .HornMonotonicArray3DConvergence
    using .DimensionlessWaveScaling

    @test getproperty.(LEVELS, :id) == ["coarse", "reference", "refined"]
    @test CARRIER_LAMBDA_P_MM ≈ 2340.0e3 / 242.0e3
    @test CARRIER_LAMBDA_S_MM ≈ 1170.0e3 / 242.0e3
    @test LEVELS[2].size_path_mm ≈ 0.90
    @test LEVELS[2].path_h_over_lambda_s ≈ 0.90 / CARRIER_LAMBDA_S_MM
    resolution = mesh_resolution(0.70, CARRIER_WAVE_SCALE)
    @test resolution.h_over_lambda_s ≈ 0.70 / CARRIER_LAMBDA_S_MM
    @test resolution.k_s_h ≈ 2pi * 0.70 / CARRIER_LAMBDA_S_MM
    @test resolution.shear_elements_per_wavelength ≈ CARRIER_LAMBDA_S_MM / 0.70
    @test physical_shear_length_mm(
        length_over_shear_wavelength(3.2, CARRIER_WAVE_SCALE),
        CARRIER_WAVE_SCALE,
    ) ≈ 3.2
    @test delay_cycles(5.0e-6, 242.0e3) ≈ 1.21
    options = parse_options(["--stage=solve", "--level=coarse", "--variant=lens"])
    @test options == (stage="solve", level="coarse", variant="lens")
    @test_throws ArgumentError parse_options(["--stage=solve", "--level=fine"])

    geometry = synthesize_monotonic_geometry(9.724)
    config = HornMonotonicArray3DMesher.HornMonotonicArray3DConfig(
        guide_axial_length_mm=geometry.axial_length_mm,
        bend_amplitude_mm=geometry.amplitude_mm,
    )
    probes = convergence_probe_points_mm(config)
    @test length(probes.channel_groups) == 15
    @test all(length.(probes.channel_groups) .== 9)
    @test length(probes.names) == 15 * 9 + 1
    @test probes.names[probes.focus_index] == "focus"

    displacement = zeros(ComplexF64, length(probes.names), 3)
    for (channel_index, group) in enumerate(probes.channel_groups)
        displacement[group, 1] .= channel_index + 2im
    end
    @test average_channel_outputs(displacement, probes.channel_groups) ==
          ComplexF64[channel + 2im for channel in 1:15]

    symmetric = ComplexF64[complex(index <= 8 ? index : 16 - index) for index in 1:15]
    mirror = mirror_pair_metrics(symmetric)
    @test mirror.maximum_amplitude_mismatch == 0.0
    @test mirror.maximum_phase_mismatch_deg == 0.0
    symmetric[end] *= 1.02cis(deg2rad(2.0))
    perturbed = mirror_pair_metrics(symmetric)
    @test perturbed.maximum_amplitude_mismatch ≈ 2 * 0.02 / 2.02
    @test perturbed.maximum_phase_mismatch_deg ≈ 2.0

    gate = convergence_gate(1.45, 1.472, 5.0, 5.0, 0.02, 2.0)
    @test gate.passed
    @test !convergence_gate(1.2, 1.472, 5.0, 7.0, 0.04, 4.0).passed
end

@testset "full monotonic 3D half-domain reconstruction" begin
    using .HornMonotonicArray3DHalfSymmetry

    phase1 = only(filter(
        level -> level.id == "phase1",
        HornMonotonicArray3DHalfSymmetry.HALF_LEVELS,
    ))
    @test phase1.path_h_over_lambda_s == 0.12
    @test phase1.focus_h_over_lambda_s == 0.20
    @test phase1.receiver_h_over_lambda_s == 0.30
    @test phase1.size_path_mm ≈ 0.5801652892561983
    @test phase1.size_focus_mm ≈ 0.9669421487603306
    @test phase1.size_receiver_mm ≈ 1.4504132231404957

    phase2 = only(filter(
        level -> level.id == "phase2",
        HornMonotonicArray3DHalfSymmetry.HALF_LEVELS,
    ))
    @test phase2.path_h_over_lambda_s == 0.10
    @test phase2.focus_h_over_lambda_s == 0.16
    @test phase2.receiver_h_over_lambda_s == 0.25
    @test phase2.size_path_mm ≈ 0.4834710743801653
    @test phase2.size_focus_mm ≈ 0.7735537190082644
    @test phase2.size_receiver_mm ≈ 1.2086776859504131

    phase3 = only(filter(
        level -> level.id == "phase3",
        HornMonotonicArray3DHalfSymmetry.HALF_LEVELS,
    ))
    @test phase3.path_h_over_lambda_s == 0.08
    @test phase3.focus_h_over_lambda_s == 0.16
    @test phase3.receiver_h_over_lambda_s == 0.25
    @test phase3.size_path_mm ≈ 0.3867768595041322
    @test phase3.size_focus_mm ≈ 0.7735537190082644
    @test phase3.size_receiver_mm ≈ 1.2086776859504131

    phase4 = only(filter(
        level -> level.id == "phase4",
        HornMonotonicArray3DHalfSymmetry.HALF_LEVELS,
    ))
    @test phase4.path_h_over_lambda_s == 0.07
    @test phase4.focus_h_over_lambda_s == 0.16
    @test phase4.receiver_h_over_lambda_s == 0.25
    @test phase4.size_path_mm ≈ 0.3384297520661157
    @test phase4.size_focus_mm ≈ 0.7735537190082644
    @test phase4.size_receiver_mm ≈ 1.2086776859504131

    nonnegative = ComplexF64.(8:15)
    reconstructed = reconstruct_full_outputs(nonnegative)
    @test reconstructed == ComplexF64[15, 14, 13, 12, 11, 10, 9, 8, 9, 10, 11, 12, 13, 14, 15]
    @test reconstructed == reverse(reconstructed)

    profile = reconstruct_full_profile(collect(0.0:2.0), [3.0, 2.0, 1.0])
    @test profile.coordinate == collect(-2.0:1.0:2.0)
    @test profile.profile == [1.0, 2.0, 3.0, 2.0, 1.0]
end

@testset "monotonic response-matrix precompensation" begin
    using .HornMonotonicResponseMatrix

    multiplicity = [1.0, 2.0]
    normalized = energy_normalize([1.0, 1.0], multiplicity)
    @test normalized == [1.0, 1.0]
    @test sum(multiplicity .* abs2.(normalized)) ≈ 3.0

    optimal = focus_optimal_weights(ComplexF64[1.0, 1.0], multiplicity)
    @test optimal[1] / optimal[2] ≈ 2.0 rtol=1e-3
    @test sum(multiplicity .* abs2.(optimal)) ≈ 3.0

    profile = profile_metrics(-2.0:1.0:2.0, ComplexF64[0.1, 0.5, 1.0, 0.5, 0.1] .* 1e-9)
    @test profile.peak_y_mm == 0.0
    @test profile.fwhm_mm == 0.0
    @test profile.sidelobe_amplitude_ratio ≈ 0.1

    power = input_power_w(reshape(ComplexF64[im], 1, 1), [1.0], 1.0, 1.0)
    @test power.active_w ≈ pi

    scan = ComplexF64[1.0 1.0; 0.5 0.5; 0.1 0.1] .* 1e-9
    design = optimize_focus_weights(
        ComplexF64[1.0, 1.0] .* 1e-9,
        [0.0, 1.0, 2.0],
        scan,
        multiplicity;
        blend_samples=21,
    )
    @test design.focus_amplitude_m >= 2e-9
    @test design.profile.fwhm_mm <= 6.0
end

@testset "monotonic broadband group delay" begin
    using .HornMonotonicGroupDelay

    unwrapped = unwrap_phase(deg2rad.([170.0, -170.0, -150.0]))
    @test rad2deg.(unwrapped) ≈ [170.0, 190.0, 210.0]

    frequencies_hz = [193.8e3, 242.0e3, 290.1e3]
    known_delay_s = 5.5e-6
    transfer = cis.(-2pi .* frequencies_hz .* known_delay_s)
    fit = linear_group_delay_s(frequencies_hz, transfer)
    @test fit.delay_s ≈ known_delay_s atol=1e-15
    @test fit.maximum_residual_deg < 1e-10

    @test isotonic_nonincreasing([4.0, 5.0, 2.0, 0.0]) ≈ [4.5, 4.5, 2.0, 0.0]
    old_path_mm = [9.0, 4.0, 0.0]
    measured_delay_s = old_path_mm .* 1e-3 ./ 1800.0
    @test effective_group_speed_m_s(old_path_mm, measured_delay_s) ≈ 1800.0
    target_delay_s = [5.5e-6, 2.5e-6, 0.0]
    correction = corrected_monotonic_paths_mm(
        old_path_mm,
        measured_delay_s,
        target_delay_s,
        1800.0,
    )
    @test correction.corrected_mm ≈ [9.9, 4.5, 0.0]
    @test all(diff(correction.corrected_mm) .<= 0)

    gated = corrected_monotonic_paths_mm(
        old_path_mm,
        [-2.0e-6, measured_delay_s[2], measured_delay_s[3]],
        target_delay_s,
        1800.0;
        reliable=[false, true, true],
    )
    @test gated.corrected_mm ≈ [9.9, 4.5, 0.0]

    options = HornMonotonicResponseMatrix.parse_options([
        "--stage=solve",
        "--variant=lens",
        "--source=1",
        "--frequency-khz=193.8",
    ])
    @test options.frequency_hz == 193.8e3
    @test endswith(
        HornMonotonicResponseMatrix.basis_path("lens", 1, options.frequency_hz),
        "basis_lens_phase2_source1_193p8khz.jld2",
    )
end

@testset "monotonic v1 carrier gate observables" begin
    using .HornMonotonicV1CarrierAnalysis

    @test wrap_phase_rad(3pi) ≈ -pi
    frequency_hz = 242.0e3
    target_delay_s = [0.0, 1.0e-6, 2.0e-6]
    common_phase = 0.37
    ideal_output = cis.(common_phase .- 2pi .* frequency_hz .* target_delay_s)
    @test maximum(abs, output_phase_error_deg(
        ideal_output,
        frequency_hz,
        target_delay_s,
    )) < 1e-10
end

@testset "monotonic broadband coupling observables" begin
    using .HornMonotonicBroadbandCoupling

    matrix = ComplexF64[1.0 0.3; 0.4 2.0]
    @test column_crosstalk_ratio(matrix, 1) ≈ 0.4
    @test column_crosstalk_ratio(matrix, 2) ≈ 0.15
end

@testset "impulse lens risk analysis" begin
    pulse = pulse_spectrum(PulseConfig(fft_length=16384))
    band6 = spectral_band(pulse, -6.0)
    band20 = spectral_band(pulse, -20.0)
    @test band6.lower_hz < 200.0e3 < band6.upper_hz
    @test band6.lower_hz < 289.0e3 < band6.upper_hz
    @test band20.lower_hz < band6.lower_hz
    @test band20.upper_hz > band6.upper_hz

    aperture = ApertureConfig()
    delays = ideal_delays_s(aperture)
    @test minimum(delays) == 0.0
    @test maximum(delays) > 5.0e-6
    @test maximum(delays) < 6.5e-6
    @test length(unique(quantize_delays(delays, 4))) <= 4

    response = focus_waveform(pulse, aperture, delays)
    metrics = pulse_metrics(
        response.focused,
        response.uniform,
        response.focused,
        pulse.time_s[2] - pulse.time_s[1],
    )
    @test metrics.gain_peak > 1.0
    @test metrics.pulse_correlation ≈ 1.0 atol=1.0e-10
    coordinates = collect(30.0e-3:2.5e-3:40.0e-3)
    profile = peak_profile(pulse, aperture, delays, coordinates; direction=:axial)
    @test length(profile) == length(coordinates)
    @test contiguous_width(coordinates, profile).peak_coordinate in coordinates
end

@testset "wall profiles" begin
    amplitude = 2.5

    sinusoidal = SinusoidalProfile(amplitude; periods=2)
    @test indentation(sinusoidal, 0.0) == 0.0
    @test indentation(sinusoidal, 0.25) ≈ amplitude
    @test indentation(sinusoidal, 1.0) ≈ 0.0 atol=1.0e-14

    exponential_limit = ExponentialProfile(amplitude; periods=2, sharpness=1.0e-12)
    for ξ in range(0.0, 1.0; length=25)
        @test indentation(exponential_limit, ξ) ≈ indentation(sinusoidal, ξ) atol=1.0e-10
    end

    exponential = ExponentialProfile(amplitude; periods=2, sharpness=3.0)
    @test indentation(exponential, 0.125) < indentation(sinusoidal, 0.125)
    @test indentation(exponential, 0.25) ≈ amplitude

    power = PowerProfile(amplitude; periods=2, power=2.0)
    @test indentation(power, 0.125) < indentation(sinusoidal, 0.125)
    @test indentation(power, 0.25) ≈ amplitude
    @test_throws ArgumentError PowerProfile(amplitude; power=0.5)

    single = CoupledConstrictionProfile([2.0], [2.0], [8.5])
    single_geometry = GeometryConfig(length_mm=17.0)
    @test minimum_gap_mm(single, single_geometry) ≈ 2.0 atol=1.0e-10
    @test MetamaterialProfiles.depth_at_x(single, single_geometry, 8.5) ≈ 2.5
    @test MetamaterialProfiles.depth_at_x(single, single_geometry, 7.5) ≈ 0.0 atol=1.0e-14

    pair = CoupledConstrictionProfile([2.0, 2.4], [2.0, 3.0], [6.0, 12.0])
    pair_geometry = GeometryConfig(length_mm=18.0)
    bottom, top = generate_boundary_points(pair, pair_geometry)
    @test minimum_gap_mm(pair, pair_geometry) ≈ 2.0 atol=1.0e-10
    @test minimum(y for (x, y) in top if x == 6.0) -
          maximum(y for (x, y) in bottom if x == 6.0) ≈ 2.0
    @test_throws ArgumentError CoupledConstrictionProfile([2.0, 2.0], [3.0, 3.0], [8.0, 9.0]) |>
                               profile -> validate_geometry(profile, GeometryConfig(length_mm=17.0))
end

@testset "transient right boundary configuration" begin
    @test SimulationConfig().right_boundary_condition == :absorbing
    @test SimulationConfig(right_boundary_condition=:free_reflecting).right_boundary_condition ==
          :free_reflecting
    @test_throws ArgumentError run_acoustic_simulation(
        SinusoidalProfile(0.0),
        220.0e3;
        simulation=SimulationConfig(right_boundary_condition=:invalid),
    )
end

@testset "side-mass geometry" begin
    config = SideMassConfig()
    @test isnothing(SideMassMesher.validate_config(config))
    @test length(SideMassMesher.solid_rectangles(config)) == 11
    @test length(SideMassMesher.solid_rectangles(config; variant=:b)) == 7
    @test length(SideMassMesher.solid_rectangles(config; variant=:r0)) == 3
    @test length(SideMassMesher.solid_rectangles(config; resonators=false)) == 3
    @test_throws ArgumentError SideMassMesher.validate_config(
        SideMassConfig(bright_neck_width_mm=0.0),
    )
    @test_throws ArgumentError SideMassMesher.validate_config(
        SideMassConfig(dark_center_x_mm=12.0),
    )
    @test_throws ArgumentError SideMassMesher.solid_rectangles(config; variant=:unknown)
    detuned = SideMassConfig(dark_mass_length_mm=3.6, dark_mass_height_mm=0.9)
    @test isnothing(SideMassMesher.validate_config(detuned))
    @test SideMassMesher.solid_rectangles(detuned)[8][3:4] == (3.6, 0.9)
end

@testset "Huygens pair geometry" begin
    config = HuygensPairConfig()
    @test isnothing(HuygensPairMesher.validate_config(config))
    @test length(HuygensPairMesher.solid_rectangles(config; variant=:r0)) == 3
    @test length(HuygensPairMesher.solid_rectangles(config; variant=:front)) == 7
    @test length(HuygensPairMesher.solid_rectangles(config; variant=:rear)) == 7
    @test length(HuygensPairMesher.solid_rectangles(config; variant=:pair)) == 11
    @test_throws ArgumentError HuygensPairMesher.validate_config(
        HuygensPairConfig(spine_height_mm=3.0),
    )
    @test_throws ArgumentError HuygensPairMesher.validate_config(
        HuygensPairConfig(front_mass_length_mm=12.0),
    )
    @test length(HuygensPairPilot.pair_cases()) == 20
    @test length(unique(case.id for case in HuygensPairPilot.pair_cases())) == 20
end

@testset "balanced Huygens geometry" begin
    config = BalancedHuygensMesher.BalancedHuygensConfig()
    @test isnothing(BalancedHuygensMesher.validate_config(config))
    @test length(BalancedHuygensMesher.solid_rectangles(config; variant=:r0)) == 3
    @test length(BalancedHuygensMesher.solid_rectangles(config; variant=:series)) == 5
    @test length(BalancedHuygensMesher.solid_rectangles(config; variant=:shunt)) == 7
    @test length(BalancedHuygensMesher.solid_rectangles(config; variant=:balanced)) == 9
    @test_throws ArgumentError BalancedHuygensMesher.validate_config(
        BalancedHuygensMesher.BalancedHuygensConfig(series_height_mm=1.3),
    )
    @test_throws ArgumentError BalancedHuygensMesher.validate_config(
        BalancedHuygensMesher.BalancedHuygensConfig(shunt_mass_length_mm=13.0),
    )
    @test length(BalancedHuygensPilot.balanced_cases()) == 16
    @test length(BalancedHuygensPilot.model_specs()) == 25
end

@testset "shift-invert modal solver" begin
    stiffness = [4.0 0.0 0.0; 0.0 9.0 0.0; 0.0 0.0 16.0]
    mass = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    eigenvalues, eigenvectors = shift_invert_eigenpairs(
        stiffness,
        mass,
        5.0;
        count=1,
        tolerance=1.0e-12,
    )
    @test only(eigenvalues) ≈ 4.0 atol=1.0e-10
    @test abs(eigenvectors[1, 1]) ≈ 1.0 atol=1.0e-10
    @test maximum(abs, eigenvectors[2:3, 1]) < 1.0e-10
end

@testset "inline mass component geometry" begin
    config = InlineMassConfig()
    @test isnothing(InlineMassComponentMesher.validate_config(config))
    @test InlineMassComponentMesher.component_length_mm(config) ≈ 3.0
    rectangles = InlineMassComponentMesher.solid_rectangles(config)
    @test length(rectangles) == 3
    @test rectangles[1] == (0.0, -0.3, 0.5, 0.6)
    @test rectangles[2] == (0.5, -1.5, 2.0, 3.0)
    @test rectangles[3] == (2.5, -0.3, 0.5, 0.6)
    @test_throws ArgumentError InlineMassComponentMesher.validate_config(
        InlineMassConfig(neck_height_mm=0.3),
    )
    @test_throws ArgumentError InlineMassComponentMesher.validate_config(
        InlineMassConfig(mass_height_mm=4.3),
    )
end

@testset "stiffness component geometry" begin
    config = StiffnessConfig()
    @test isnothing(StiffnessComponentMesher.validate_config(config))
    @test StiffnessComponentMesher.slot_height_mm(config) ≈ 2.6 / 3
    @test StiffnessComponentMesher.ligament_area_fraction(config) ≈ 1.6 / 4.2
    boxes = StiffnessComponentMesher.slot_boxes(config)
    @test length(boxes) == 3
    @test boxes[2][2] ≈ -StiffnessComponentMesher.slot_height_mm(config) / 2
    @test boxes[1][2] ≈ -boxes[3][2] - boxes[3][4]
    @test_throws ArgumentError StiffnessComponentMesher.validate_config(
        StiffnessConfig(ligament_height_mm=0.30),
    )
    @test_throws ArgumentError StiffnessComponentMesher.validate_config(
        StiffnessConfig(end_wall_mm=2.0),
    )
end

@testset "combined KMK geometry" begin
    config = CombinedKMKConfig()
    @test isnothing(CombinedKMKMesher.validate_config(config))
    @test CombinedKMKMesher.cell_length_mm(config) ≈ 5.3
    @test CombinedKMKMesher.total_length_mm(config) ≈ 29.3
    @test length(CombinedKMKMesher.solid_rectangles(config; variant=:solid)) == 1
    @test length(CombinedKMKMesher.solid_rectangles(config; variant=:k_only)) == 1
    @test length(CombinedKMKMesher.solid_rectangles(config; variant=:m_only)) == 5
    @test length(CombinedKMKMesher.solid_rectangles(config; variant=:combined)) == 5
    @test_throws ArgumentError CombinedKMKMesher.solid_rectangles(
        config;
        variant=:unknown,
    )
end

@testset "distributed slow-wave geometry" begin
    config = DistributedSlowWaveConfig()
    @test isnothing(DistributedSlowWaveMesher.validate_config(config))
    @test DistributedSlowWaveMesher.active_length_mm(config) ≈ 9.5
    @test DistributedSlowWaveMesher.total_length_mm(config) ≈ 33.5
    @test DistributedSlowWaveMesher.section_starts_mm(config) ≈ [12.0, 14.5, 17.0, 19.5]
    @test DistributedSlowWaveMesher.active_section_indices(config, :solid) == Int[]
    @test DistributedSlowWaveMesher.active_section_indices(config, :mid) == [2, 3]
    @test DistributedSlowWaveMesher.active_section_indices(config, :max) == [1, 2, 3, 4]
    @test length(DistributedSlowWaveMesher.slot_boxes(config, :mid)) == 6
    @test length(DistributedSlowWaveMesher.slot_boxes(config, :max)) == 12
    @test_throws ArgumentError DistributedSlowWaveMesher.validate_config(
        DistributedSlowWaveConfig(section_count=3),
    )
    @test_throws ArgumentError DistributedSlowWaveMesher.active_section_indices(config, :unknown)
end

@testset "tapered shunt-trim geometry" begin
    config = TaperedShuntTrimConfig()
    @test isnothing(TaperedShuntTrimMesher.validate_config(config))
    @test TaperedShuntTrimMesher.total_length_mm(config) ≈ 36.0
    @test TaperedShuntTrimMesher.taper_half_height_mm(config, 0.0) ≈ 1.6
    @test TaperedShuntTrimMesher.taper_half_height_mm(config, 3.0) ≈ 0.8
    @test TaperedShuntTrimMesher.taper_half_height_mm(config, 6.0) ≈ 0.8
    @test TaperedShuntTrimMesher.taper_half_height_mm(config, 12.0) ≈ 1.6
    mid = TaperedShuntTrimMesher.shunt_rectangles(config, :mid)
    max_state = TaperedShuntTrimMesher.shunt_rectangles(config, :max)
    @test length(mid) == 4
    @test mid[1][3] ≈ 2.4
    @test max_state[1][3] ≈ 4.0
    @test_throws ArgumentError TaperedShuntTrimMesher.validate_config(
        TaperedShuntTrimConfig(neck_width_mm=0.3),
    )
end

@testset "perforated straight-channel geometry" begin
    config = PerforatedChannelConfig()
    @test isnothing(PerforatedChannelMesher.validate_config(config))
    @test PerforatedChannelMesher.total_length_mm(config) ≈ 36.0
    @test length(PerforatedChannelMesher.hole_centers_mm(config)) == 8
    @test all(
        isapprox(
            PerforatedChannelMesher.void_area_mm2(config, variant),
            PerforatedChannelMesher.void_area_mm2(config, :round),
        )
        for variant in (:ellipse_longitudinal, :ellipse_transverse)
    )
    @test all(
        minimum(values(PerforatedChannelMesher.minimum_ligaments_mm(config, variant))) >=
        config.minimum_feature_mm
        for variant in (:round, :ellipse_longitudinal, :ellipse_transverse)
    )
    @test_throws ArgumentError PerforatedChannelMesher.validate_config(
        PerforatedChannelConfig(row_center_mm=0.9),
    )
end

@testset "equal-length global-fold geometry" begin
    config = FoldedPathConfig()
    @test isnothing(FoldedPathMesher.validate_config(config))
    @test FoldedPathMesher.unfolded_active_length_mm(config) ≈ 30.42
    @test FoldedPathMesher.total_axial_length_mm(config, :straight_device) ≈ 44.0
    @test FoldedPathMesher.total_axial_length_mm(config, :straight_unfolded) ≈ 54.42
    @test FoldedPathMesher.sharp_fold_amplitude_mm(config) ≈ 11.4605453622 atol=1e-9
    rounded_amplitude = FoldedPathMesher.rounded_fold_amplitude_mm(config)
    @test rounded_amplitude ≈ 10.75518175 atol=1e-6
    @test FoldedPathMesher.rounded_centerline_length_mm(config, rounded_amplitude) ≈
          FoldedPathMesher.unfolded_active_length_mm(config) atol=1e-8
    centerline = FoldedPathMesher.rounded_centerline_mm(config)
    @test all(isapprox.(first(centerline)[2:3], (0.0, 0.0); atol=1e-12))
    @test all(isapprox.(last(centerline)[2:3], (0.0, 0.0); atol=1e-12))
    @test FoldedPathMesher.minimum_rounded_inner_radius_mm(config) >=
          config.minimum_inner_radius_mm
    @test_throws ArgumentError FoldedPathMesher.validate_config(
        FoldedPathConfig(extra_path_length_mm=14.0),
    )
end

@testset "elastic strip port modes" begin
    low_frequency_modes = solve_port_modes(PortModeConfig(
        frequency_hz=10.0e3,
        element_count=16,
    ))
    right_going = propagating_modes(low_frequency_modes; direction=:right)
    left_going = propagating_modes(low_frequency_modes; direction=:left)
    @test length(right_going) == length(left_going)
    @test length(right_going) >= 2
    @test all(mode -> isapprox(mode.power_w_per_m, 1.0; atol=1.0e-9), right_going)
    @test all(mode -> isapprox(mode.power_w_per_m, -1.0; atol=1.0e-9), left_going)
    @test maximum(mode.relative_residual for mode in low_frequency_modes) < 1.0e-7

    input_mode = fundamental_quasi_longitudinal(low_frequency_modes)
    @test input_mode.parity == :symmetric
    @test input_mode.p_fraction > 0.9
    @test input_mode.axial_displacement_fraction > 0.9

    branch = track_fundamental_branch(
        PortModeConfig(frequency_hz=10.0e3, element_count=16),
        [10.0e3, 40.0e3, 80.0e3],
    )
    @test length(branch.modes) == 3
    @test all(mode -> mode.parity == :symmetric, branch.modes)
    @test minimum(branch.overlaps) > 0.9

    target_modes = solve_port_modes(PortModeConfig(element_count=20))
    target_right = propagating_modes(target_modes; direction=:right)
    @test length(target_right) >= 3
    @test any(mode -> mode.parity == :symmetric, target_right)
    @test fundamental_quasi_longitudinal(target_modes).parity == :symmetric
end

@testset "uniform reference geometry" begin
    config = UniformReferenceMesher.UniformReferenceConfig()
    @test config.height_mm == 4.2
    @test config.length_mm == 24.0
    @test isnothing(UniformReferenceMesher.validate_config(config))
    @test_throws ArgumentError UniformReferenceMesher.validate_config(
        UniformReferenceMesher.UniformReferenceConfig(height_mm=0.0),
    )
end

@testset "baseline R0 geometry" begin
    config = BaselineR0Mesher.BaselineR0Config()
    @test BaselineR0Mesher.total_length_mm(config) == 42.0
    @test length(BaselineR0Mesher.solid_rectangles(config)) == 3
    @test isnothing(BaselineR0Mesher.validate_config(config))
    @test_throws ArgumentError BaselineR0Mesher.validate_config(
        BaselineR0Mesher.BaselineR0Config(spine_height_mm=4.2),
    )
end

@testset "biorthogonal modal projection" begin
    modes = solve_port_modes(PortModeConfig(frequency_hz=40.0e3, element_count=20))
    right_mode = fundamental_quasi_longitudinal(modes)
    left_mode = counterpropagating_mode(right_mode, modes)
    right_amplitude = 0.7 - 0.2im
    left_amplitude = -0.1 + 0.3im
    sample = ElasticModalProjection.BoundaryFieldSample(
        copy(right_mode.y_m),
        right_amplitude .* right_mode.displacement_x .+
            left_amplitude .* left_mode.displacement_x,
        right_amplitude .* right_mode.displacement_y .+
            left_amplitude .* left_mode.displacement_y,
        right_amplitude .* right_mode.traction_xx_pa .+
            left_amplitude .* left_mode.traction_xx_pa,
        right_amplitude .* right_mode.traction_xy_pa .+
            left_amplitude .* left_mode.traction_xy_pa,
    )
    decomposition = decompose_mode_pair(sample, right_mode, left_mode)
    @test decomposition.right_amplitude ≈ right_amplitude atol=1.0e-10
    @test decomposition.left_amplitude ≈ left_amplitude atol=1.0e-10
end

@testset "rank-one modal boundary" begin
    displacement = VectorValue(1.0 + 0.5im, -0.2 + 0.7im)
    traction = VectorValue(-3.0 + 2.0im, 0.4 - 1.0im)
    operator = ModalHarmonicElasticity.rank_one_traction_map(displacement, traction)
    @test operator ⋅ displacement ≈ traction
end

@testset "phase efficiency" begin
    @test PhaseEfficiencyAnalysis.circular_error_deg(359.0, 1.0) == 2.0
    rows = [
        (valid=true, phase_deg=0.0, amplitude=0.4),
        (valid=true, phase_deg=8.0, amplitude=0.6),
        (valid=true, phase_deg=40.0, amplitude=0.9),
    ]
    best = PhaseEfficiencyAnalysis.best_for_phase(rows, 0.0, 10.0)
    @test best.amplitude == 0.6
    pareto = PhaseEfficiencyAnalysis.pareto_candidates(rows, 0.0)
    @test length(pareto) == 3
end

@testset "lens depth of focus" begin
    x = collect(1.0:1.0:7.0)
    y = [-1.0, 0.0, 1.0]
    amplitude = repeat(reshape([0.1, 0.5, 0.8, 1.0, 0.8, 0.5, 0.1], 1, :), 3, 1)
    depth = axial_depth_of_focus(amplitude, x, y)
    @test depth.start_mm == 3.0
    @test depth.stop_mm == 5.0
    @test depth.depth_mm == 2.0
end

@testset "ideal Huygens aperture" begin
    @test IdealHuygensAperture.nearest_odd_count(74.0, 4.8, 0.6) == 15
    @test IdealHuygensAperture.nearest_odd_count(74.0, 8.2, 1.2) == 9
    configs = IdealHuygensAperture.candidate_configs()
    @test length(configs) == 10
    @test all(isodd(config.element_count) for config in configs)
    fine_short = only(filter(
        config -> config.element_count == 15 &&
                  config.focal_distance_mm == 35.0 &&
                  config.element_width_mm + config.slot_width_mm == 4.8,
        configs,
    ))
    selection = IdealHuygensAperture.ideal_selection(fine_short)
    @test length(selection.entries) == 15
    @test abs(selection.focus_field) > 2.9
    @test IdealHuygensAperture.aperture_width_mm(fine_short) ≈ 71.4
end

@testset "target-frequency library" begin
    candidates = TargetFrequencyLibrary.collect_candidates()
    @test TargetFrequencyLibrary.TARGET_FREQUENCIES_HZ == [242.0e3, 244.0e3]
    @test length(unique(TargetFrequencyLibrary.geometry_key.(candidates))) == length(candidates)
    @test all(occursin("_STG_F_220.0.jld2", item.characteristic_path) for item in candidates)
    coverage = TargetFrequencyLibrary.phase_coverage([
        (spectral_valid=true, amplitude=0.5, phase_rad=0.0),
        (spectral_valid=true, amplitude=0.5, phase_rad=pi),
    ], 0.4)
    @test coverage.count == 2
    @test coverage.covered_arc_rad ≈ pi
end

@testset "pilot sweep definition" begin
    profiles = pilot_profiles()
    @test length(profiles) == 10
    @test SWEEP_FREQUENCIES_HZ == collect(160.0e3:20.0e3:280.0e3)
    @test PILOT_FREQUENCY_HZ in SWEEP_FREQUENCIES_HZ
    @test length(unique(profile_slug.(profiles))) == length(profiles)
    @test amplitude_mm(first(profiles)) == 0.0
    @test all(
        minimum_gap_mm(profile) ≈ 4.5
        for profile in profiles[2:9]
    )

    rows = [
        (slug="sample", requested_frequency_hz=160e3, phase_rad=3.0),
        (slug="sample", requested_frequency_hz=180e3, phase_rad=-3.0),
        (slug="other", requested_frequency_hz=160e3, phase_rad=-1.0),
    ]
    unwrapped = unwrap_frequency_rows(rows)
    @test unwrapped[2].phase_rad ≈ 2.0 * pi - 3.0
    @test unwrapped[3].phase_rad == -1.0
end

@testset "working-band boundary study definition" begin
    study = WorkingBandBoundaryStudy
    @test study.BAND_SCAN_FREQUENCIES_HZ == collect(120.0e3:10.0e3:160.0e3)
    @test length(study.TOPOLOGIES) == 10
    @test length(unique(getproperty.(study.TOPOLOGIES, :id))) == 10
    @test study.contiguous_band_indices([0.2, 0.8, 1.0, 0.75, 0.1], 3, 0.7) == 2:4
    edges = study.interpolated_band_edges(100.0:10.0:140.0, [0.2, 0.8, 1.0, 0.75, 0.1], 3, 0.7)
    @test edges.estimated_low_hz ≈ 108.33333333333333
    @test edges.estimated_high_hz ≈ 130.76923076923077
end

@testset "period and length pilot definition" begin
    period_cases = PeriodLengthPilot.period_sweep_cases()
    length_cases_n1 = PeriodLengthPilot.length_sweep_cases(1)
    length_cases_n2 = PeriodLengthPilot.length_sweep_cases(2)
    all_cases = PeriodLengthPilot.all_cases()

    @test getproperty.(getproperty.(period_cases, :profile), :periods) == [1, 2, 3, 4]
    @test getproperty.(getproperty.(length_cases_n1, :geometry), :length_mm) ==
          PeriodLengthPilot.LENGTH_VALUES_MM
    @test PeriodLengthPilot.FREQUENCY_MAP_HZ == collect(180.0e3:2.0e3:280.0e3)
    @test getproperty.(getproperty.(length_cases_n2, :geometry), :length_mm) ==
          PeriodLengthPilot.ORIGINAL_LENGTH_VALUES_MM
    @test length(PeriodLengthPilot.sample_cases()) ==
          12 + length(PeriodLengthPilot.EXTENDED_LENGTH_VALUES_MM)
    @test length(all_cases) ==
          length(PeriodLengthPilot.sample_cases()) + length(PeriodLengthPilot.LENGTH_VALUES_MM)
    @test length(unique(getproperty.(all_cases, :id))) == length(all_cases)
    @test all(
        PeriodLengthPilot.MetamaterialProfiles.minimum_gap_mm(case.profile, case.geometry) ≈ 4.5
        for case in PeriodLengthPilot.sample_cases()
    )
    @test all(
        PeriodLengthPilot.matching_reference(case).geometry.length_mm == case.geometry.length_mm
        for case in PeriodLengthPilot.sample_cases()
    )

    phase_rows = [(phase_rad=3.0,), (phase_rad=-3.0,), (phase_rad=-2.5,)]
    unwrapped = PeriodLengthPilot.unwrap_parameter_phase(phase_rows)
    @test unwrapped[2] ≈ 2.0 * pi - 3.0
    @test unwrapped[3] ≈ 2.0 * pi - 2.5
end

@testset "minimum-gap pilot definition" begin
    cases = GapPilot.gap_cases()
    new_cases = GapPilot.new_cases()

    @test length(cases) ==
          4 * length(GapPilot.ORIGINAL_GAP_VALUES_MM) +
          3 * length(GapPilot.TARGET_GAP_VALUES_MM)
    @test length(new_cases) == length(cases)
    @test length(unique(getproperty.(cases, :id))) == length(cases)
    @test sort(unique(getproperty.(cases, :gap_mm))) == GapPilot.GAP_VALUES_MM
    @test GapPilot.TARGET_FREQUENCIES_HZ == collect(225.0e3:2.0e3:250.0e3)
    @test all(
        GapPilot.MetamaterialProfiles.minimum_gap_mm(case.profile, case.geometry) ≈ case.gap_mm
        for case in cases
    )
    @test all(isfile(GapPilot.reference_signal_path(case)) for case in cases)

    phase_rows = [(phase_rad=3.0,), (phase_rad=-3.0,)]
    @test GapPilot.unwrap_parameter_phase(phase_rows)[2] ≈ 2.0 * pi - 3.0
end

@testset "independent-element lens design" begin
    library_path = joinpath(
        @__DIR__,
        "..",
        "tmp",
        "gap_pilot_220khz",
        "gap_pilot.csv",
    )
    entries = load_gap_library(library_path)
    config = LensConfig()
    centers = lens_centers_mm(config)
    selection = select_lens(entries; config)
    phase_values = phase_design_values(selection; config)

    @test length(entries) == length(GapPilot.gap_cases())
    @test centers == -reverse(centers)
    @test length(selection.entries) == config.element_count
    @test getproperty.(selection.entries, :case_id) == reverse(getproperty.(selection.entries, :case_id))
    @test required_phase_span(config) < 3.56
    @test isfinite(abs(selection.focus_field))
    @test maximum(abs, phase_values.errors) <= pi

    x = [config.focal_distance_mm]
    y = [-1.0, 0.0, 1.0]
    field = field_grid(x, y, selection; config)
    reference = reference_field_grid(x, y; config)
    @test size(field) == (3, 1)
    @test size(reference) == size(field)
    @test field[1, 1] ≈ field[3, 1]

    compact_entries = entries[1:min(4, length(entries))]
    objective = LensObjectiveConfig(
        axial_offsets_mm=[-15.0, 15.0],
        transverse_offsets_mm=[-12.0, 12.0],
    )
    contrast_selection = select_lens_contrast(compact_entries; config, objective)
    metrics = selection_objective_metrics(contrast_selection; config, objective)
    @test length(contrast_selection.entries) == config.element_count
    @test isfinite(metrics.score)
end

@testset "spectral transfer analysis" begin
    dt = 1.0e-6
    sample_count = 1024
    time = collect(0:(sample_count - 1)) .* dt
    input = zeros(sample_count)
    input[100] = 1.0

    delay_samples = 7
    gain = 0.4
    output = zeros(sample_count)
    output[(100 + delay_samples)] = gain

    config = SpectrumConfig(
        window=:rectangular,
        zero_padding_factor=1,
        input_floor_relative=0.0,
        regularization_relative=0.0,
    )
    result = analyze_transfer(time, input, output; config)
    interior = 3:(length(result.frequency_hz) - 2)

    @test all(result.valid)
    @test maximum(abs.(result.amplitude[interior] .- gain)) < 1.0e-12
    @test maximum(abs.(result.group_delay_s[interior] .- delay_samples * dt)) < 1.0e-12

    point = value_at_frequency(result, 100e3)
    @test point.valid
    @test point.amplitude ≈ gain
    @test envelope_delay(time, input, output) ≈ delay_samples * dt
end

@testset "reference-calibrated transfer" begin
    dt = 1.0e-6
    sample_count = 1024
    time = collect(0:(sample_count - 1)) .* dt
    input = zeros(sample_count)
    input[100] = 1.0

    reference_output = zeros(sample_count)
    reference_output[104] = 0.8
    sample_output = zeros(sample_count)
    sample_output[110] = 0.5

    config = SpectrumConfig(
        window=:rectangular,
        zero_padding_factor=1,
        input_floor_relative=0.0,
        regularization_relative=0.0,
    )
    reference = analyze_transfer(time, input, reference_output; config)
    sample = analyze_transfer(time, input, sample_output; config)
    calibrated = relative_transfer(sample, reference)
    interior = 3:(length(calibrated.frequency_hz) - 2)

    @test maximum(abs.(calibrated.amplitude[interior] .- 0.5 / 0.8)) < 1.0e-12
    @test maximum(abs.(calibrated.group_delay_s[interior] .- 6.0e-6)) < 1.0e-12
end

@testset "Morlet CWT timing" begin
    pilot_time = collect(range(0.0, 120.0e-6; length=801))
    pilot_grid = recommended_frequency_grid_hz(pilot_time)
    @test first(pilot_grid) == 50.0e3
    @test last(pilot_grid) == 1.0e6
    @test all(diff(pilot_grid) .== 5.0e3)
    pilot_config = recommended_morlet_config(pilot_time)
    @test pilot_config.frequencies_hz == pilot_grid

    dt = 1.0e-6
    sample_count = 512
    time = collect(0:(sample_count - 1)) .* dt
    frequency_hz = 80.0e3
    width_s = 24.0e-6
    reference_center_s = 180.0e-6
    delay_s = 12.0e-6
    packet(center_s) = @. exp(-0.5 * ((time - center_s) / width_s)^2) *
                           cos(2.0 * pi * frequency_hz * (time - center_s))
    reference_signal = packet(reference_center_s)
    sample_signal = 0.4 .* packet(reference_center_s + delay_s)
    config = MorletConfig(
        frequencies_hz=[frequency_hz],
        reference_energy_floor_relative=0.0,
        amplitude_ratio_floor=0.0,
        correlation_floor=0.0,
    )

    reference = continuous_wavelet_transform(time, reference_signal; config)
    sample = continuous_wavelet_transform(time, sample_signal; config)
    comparison = compare_wavelets(reference, sample; config)
    point = wavelet_value_at_frequency(comparison, frequency_hz)

    @test size(reference.coefficients) == (1, sample_count)
    @test !reference.cone_of_influence_valid[1, 1]
    @test reference.cone_of_influence_valid[1, sample_count ÷ 2]
    @test point.valid
    @test point.amplitude_ratio ≈ 0.4 rtol=2.0e-3
    @test point.centroid_delay_s ≈ delay_s atol=dt
    @test point.peak_delay_s ≈ delay_s atol=dt
    @test point.correlation_delay_s ≈ delay_s atol=dt

    self_comparison = compare_wavelets(reference, reference; config)
    self = wavelet_value_at_frequency(self_comparison, frequency_hz)
    @test self.centroid_delay_s ≈ 0.0 atol=eps(Float64)
    @test self.correlation_delay_s ≈ 0.0 atol=dt
end

@testset "geometry validation and boundaries" begin
    config = GeometryConfig(samples=200)
    profile = SinusoidalProfile(2.5)
    bottom, top = generate_boundary_points(profile, config)

    @test first(bottom) == (0.0, 0.0)
    @test last(bottom) == (config.length_mm, 0.0)
    @test first(top) == (config.length_mm, config.height_mm)
    @test last(top) == (0.0, config.height_mm)
    @test length(bottom) == config.samples + 3
    @test_throws ArgumentError validate_geometry(SinusoidalProfile(7.0), config)
    @test minimum_gap_mm(profile, config) ≈ 4.5

    for current in (
        profile,
        ExponentialProfile(2.5; sharpness=3.0),
        PowerProfile(2.5; power=4.0),
    )
        current_bottom, current_top = generate_boundary_points(current, config)
        bottom_y = Dict(x => y for (x, y) in current_bottom)
        top_y = Dict(x => y for (x, y) in current_top)
        common_x = intersect(keys(bottom_y), keys(top_y))
        interior_gaps = [
            top_y[x] - bottom_y[x]
            for x in common_x
            if config.end_margin_mm + 3.3 <= x <= config.length_mm - config.end_margin_mm - 3.3
        ]
        @test !isempty(interior_gaps)
        @test all(isapprox.(interior_gaps, 4.5; atol=1.0e-10))
        @test endswith(profile_slug(current), "_STG")
    end
end

@testset "rounded U-notch wall profile" begin
    config = GeometryConfig(samples=200)
    profile = RoundedNotchProfile(
        2.5;
        notch_count=4,
        notch_width_mm=1.4,
        end_margin_mm=1.2,
    )
    centers = MetamaterialProfiles.notch_centers_mm(profile, config)
    @test centers ≈ [3.025, 6.675, 10.325, 13.975]
    @test minimum_gap_mm(profile, config) ≈ 4.5
    @test MetamaterialProfiles.depth_at_x(profile, config, centers[1]) ≈ 2.5
    @test MetamaterialProfiles.top_depth_at_x(profile, config, centers[1]) ≈ 0.0
    @test MetamaterialProfiles.depth_at_x(profile, config, centers[2]) ≈ 0.0
    @test MetamaterialProfiles.top_depth_at_x(profile, config, centers[2]) ≈ 2.5
    bottom, top = generate_boundary_points(profile, config)
    @test first(bottom) == (0.0, 0.0)
    @test last(bottom) == (config.length_mm, 0.0)
    @test first(top) == (config.length_mm, config.height_mm)
    @test last(top) == (0.0, config.height_mm)
    @test endswith(profile_slug(profile), "_STG")
    @test_throws ArgumentError RoundedNotchProfile(0.5; notch_width_mm=1.4)
    @test_throws ArgumentError validate_geometry(
        RoundedNotchProfile(2.5; notch_count=4, notch_width_mm=4.0),
        config,
    )
end

@testset "legacy geometry compatibility" begin
    config = GeometryConfig(samples=200)
    amplitude = 2.8
    correction = amplitude * 4.0 / 100.0
    bottom, _ = generate_boundary_points(LegacySinusoidalProfile(amplitude), config)

    for i in (0, 50, 100, 150, 200)
        angle_deg = 4.0 * 180.0 * i / (config.samples + 1.0)
        expected = if i <= config.samples / 4 || i >= 3 * config.samples / 4
            amplitude / 2.0 * (1.0 - cosd(angle_deg))
        else
            (amplitude / 2.0 - correction / 2.0) * (1.0 - cosd(angle_deg)) + correction
        end
        @test bottom[i + 2][2] ≈ expected
    end


    # The legacy top and bottom profiles are staggered, so A > height/2 can
    # remain geometrically open and must not be rejected by W - 2A.
    @test minimum_gap_mm(LegacySinusoidalProfile(3.6), config) > 0
    @test isnothing(validate_geometry(LegacySinusoidalProfile(3.6), config))
end

include("operator_time_reversal_audit_tests.jl")
include("sinusoidal_waveguide_temporal_reversal_tests.jl")
include("material_lens_response_matrix_tests.jl")
include("sinusoidal_waveguide_broadband_matrix_tests.jl")
