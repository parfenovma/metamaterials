module TopologyHybridTTDGRINHarmonic

using Gridap

if !isdefined(parentmodule(@__MODULE__), :TopologyHybridTTDGRINDesign)
    include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
    using .TopologyHybridTTDGRINDesign
else
    using ..TopologyHybridTTDGRINDesign
end
include(joinpath(@__DIR__, "modal_harmonic_solver.jl"))
using .ModalHarmonicElasticity
include(joinpath(@__DIR__, "modal_projection.jl"))
using .ElasticModalProjection

export HybridHarmonicConfig, solve_hybrid_channel_scattering

Base.@kwdef struct HybridHarmonicConfig
    element_order::Int = 1
    quadrature_degree::Int = 2
    port_element_count::Int = 80
end

function bulk_speeds(material)
    rho = material.density_kg_m3
    young = material.young_modulus_pa
    nu = material.poisson_ratio
    mu = young / (2(1 + nu))
    lambda = young * nu / ((1 + nu) * (1 - 2nu))
    (pressure=sqrt((lambda + 2mu) / rho), shear=sqrt(mu / rho))
end

function stress_function(material)
    young = material.young_modulus_pa
    nu = material.poisson_ratio
    mu = young / (2(1 + nu))
    lambda = young * nu / ((1 + nu) * (1 - 2nu))
    identity = one(TensorValue{2, 2, ComplexF64})
    strain -> lambda * tr(strain) * identity + 2mu * strain
end

function project_all(state, tag)
    right_modes = filter(
        mode -> mode.kind == :propagating && mode.direction == :right,
        state.modes,
    )
    [begin
        left_candidates = filter(
            mode -> mode.kind == :propagating &&
                    mode.direction == :left &&
                    mode.parity == right_mode.parity,
            state.modes,
        )
        left_mode = first(sort(left_candidates; by=mode -> abs(
            mode.wavenumber_per_m + right_mode.wavenumber_per_m,
        )))
        projection = project_boundary_mode_pair(
            state,
            tag,
            right_mode,
            left_mode;
            quadrature_degree=6,
        )
        (; right_mode, left_mode, projection...)
    end for right_mode in right_modes]
end

function solve_hybrid_channel_scattering(
    mesh_path::AbstractString,
    frequency_hz::Real,
    material_fractions::AbstractVector{<:Real};
    design::HybridTTDGRINConfig=HybridTTDGRINConfig(),
    config::HybridHarmonicConfig=HybridHarmonicConfig(),
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    model = DiscreteModelFromFile(mesh_path)
    omega = 2pi * Float64(frequency_hz)
    al_constituent = aluminium_6061()
    al_speeds = bulk_speeds(al_constituent)
    port_config = ModalHarmonicElasticity.HarmonicElasticity.HarmonicConfig(
        density=al_constituent.density_kg_m3,
        pressure_wave_speed=al_speeds.pressure,
        shear_wave_speed=al_speeds.shear,
        rayleigh_alpha=0.0,
        rayleigh_beta=0.0,
        element_order=config.element_order,
        quadrature_degree=config.quadrature_degree,
    )
    modes, right_mode, left_mode = ModalHarmonicElasticity.matched_port_modes(
        port_config,
        frequency_hz,
        design.guide_width_mm * 1e-3,
        config.port_element_count,
    )

    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model,
        reference;
        conformity=:H1,
        vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)
    left_port = BoundaryTriangulation(model; tags=["Source"])
    right_port = BoundaryTriangulation(model; tags=["Microphone"])
    left_measure = Measure(left_port, config.quadrature_degree)
    right_measure = Measure(right_port, config.quadrature_degree)
    left_operator = ModalHarmonicElasticity.modal_boundary_fields(left_port, left_mode, -1.0)
    right_operator = ModalHarmonicElasticity.modal_boundary_fields(right_port, right_mode, 1.0)
    source_traction = ModalHarmonicElasticity.incident_source_field(
        left_port,
        right_mode,
        left_operator,
    )

    contributions = NamedTuple[]
    for raw_fraction in sort(unique(Float64.(material_fractions)))
        key = round(Int, 1000raw_fraction)
        tag = "MaterialF$(lpad(string(key), 4, '0'))"
        domain = Triangulation(model; tags=[tag])
        measure = Measure(domain, config.quadrature_degree)
        guided = guided_material(raw_fraction, design)
        material = (
            density_kg_m3=guided.density_kg_m3,
            young_modulus_pa=(1 - raw_fraction) * aluminium_6061().young_modulus_pa +
                             raw_fraction * copper_c110().young_modulus_pa,
            poisson_ratio=(1 - raw_fraction) * aluminium_6061().poisson_ratio +
                          raw_fraction * copper_c110().poisson_ratio,
        )
        push!(contributions, (
            fraction=raw_fraction,
            measure,
            density=material.density_kg_m3,
            sigma=stress_function(material),
        ))
    end
    isempty(contributions) && error("no material domains were supplied")

    function bilinear(u, v)
        first_contribution = first(contributions)
        result = ∫(
            first_contribution.sigma(ε(u)) ⊙ ε(v) -
            first_contribution.density * omega^2 * (u ⋅ v),
        )first_contribution.measure
        for contribution in Iterators.drop(contributions, 1)
            result = result + ∫(
                contribution.sigma(ε(u)) ⊙ ε(v) -
                contribution.density * omega^2 * (u ⋅ v),
            )contribution.measure
        end
        result -
        ∫((left_operator ⋅ u) ⋅ v)left_measure -
        ∫((right_operator ⋅ u) ⋅ v)right_measure
    end
    linear(v) = ∫(source_traction ⋅ v)left_measure
    displacement = solve(AffineFEOperator(bilinear, linear, trial_space, test_space))
    sigma_al = stress_function((
        density_kg_m3=aluminium_6061().density_kg_m3,
        young_modulus_pa=aluminium_6061().young_modulus_pa,
        poisson_ratio=aluminium_6061().poisson_ratio,
    ))
    state = (
        model,
        displacement,
        stress=sigma_al ∘ ε(displacement),
        omega,
        modes,
        right_mode,
        left_mode,
    )
    left = project_all(state, "Source")
    right = project_all(state, "Microphone")
    input_index = findfirst(pair -> pair.right_mode.parity == :symmetric, left)
    input_index === nothing && error("symmetric input mode was not found")
    incident = abs2(left[input_index].right_amplitude)
    reflection = abs2(left[input_index].left_amplitude) / incident
    transmission = abs2(right[input_index].right_amplitude) / incident
    conversion = sum(
        abs2(left_pair.left_amplitude) + abs2(right_pair.right_amplitude)
        for (index, (left_pair, right_pair)) in enumerate(zip(left, right))
        if index != input_index
    ) / incident
    transfer = right[input_index].right_amplitude / left[input_index].right_amplitude
    (
        frequency_hz=Float64(frequency_hz),
        S21=ComplexF64(transfer),
        S11=ComplexF64(left[input_index].left_amplitude /
                       left[input_index].right_amplitude),
        T00=transmission,
        R00=reflection,
        C0=conversion,
        power_sum=transmission + reflection + conversion,
        propagating_mode_count=length(left),
    )
end

end
