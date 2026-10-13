module OperatorO8PlaneStrainSolver

using Gridap
using GridapGmsh
using LinearAlgebra
using SparseArrays
using Statistics

export O8PlaneStrainConfig,
       O8PZTBoundary,
       O8MechanicalBoundary,
       O8Load,
       validate,
       solve_o8_frequency,
       state_projector_distance

Base.@kwdef struct O8PlaneStrainConfig
    frequency_hz::Float64 = 242.0e3
    aperture_height_m::Float64
    block_length_m::Float64
    density_kg_m3::Float64 = 2700.0
    pressure_speed_m_s::Float64 = 6122.102437409232
    shear_speed_m_s::Float64 = 3083.810277185563
    material_quality_factor::Float64 = 200.0
    source_pressure_scale_pa_sqrt_m::Float64 = 1.0
    output_bins::Int = 17
    sample_points_per_bin::Int = 7
    element_order::Int = 1
    quadrature_degree::Int = 2
end

Base.@kwdef struct O8PZTBoundary
    thickness_m::Float64 = 1.0e-3
    density_kg_m3::Float64 = 7750.0
    c33e_pa::Float64 = 117.0e9
    mechanical_quality_factor::Float64 = 55.0
    coupling_squared::Float64 = 0.45
    capacitance_per_patch_f::Float64 = 2.0e-9
    pair_a_parity::Int = 1
    pair_b_parity::Int = -1
end

Base.@kwdef struct O8MechanicalBoundary
    actuator_surface_mass_kg_m2::Float64 = 3.9
    off_stiffness_fraction::Float64 = 0.02
    stiffness_a_n_m3::Float64
    stiffness_b_n_m3::Float64
    loss_a::Float64
    loss_b::Float64
end

Base.@kwdef struct O8Load
    inductance_h::Float64
    resistance_ohm::Float64
end

function validate(config::O8PlaneStrainConfig)
    config.frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    config.aperture_height_m > 0 || throw(ArgumentError("aperture must be positive"))
    config.block_length_m > 0 || throw(ArgumentError("block length must be positive"))
    config.pressure_speed_m_s > config.shear_speed_m_s > 0 ||
        throw(ArgumentError("wave speeds must satisfy cP > cS > 0"))
    config.density_kg_m3 > 0 || throw(ArgumentError("density must be positive"))
    config.material_quality_factor > 0 || throw(ArgumentError("Q must be positive"))
    config.output_bins >= 3 || throw(ArgumentError("at least three output bins"))
    config.sample_points_per_bin >= 1 || throw(ArgumentError("samples per bin must be positive"))
    config
end

function validate(boundary::O8PZTBoundary)
    boundary.thickness_m > 0 || throw(ArgumentError("PZT thickness must be positive"))
    boundary.density_kg_m3 > 0 || throw(ArgumentError("PZT density must be positive"))
    boundary.c33e_pa > 0 || throw(ArgumentError("PZT stiffness must be positive"))
    boundary.mechanical_quality_factor > 0 || throw(ArgumentError("PZT Q must be positive"))
    0 < boundary.coupling_squared < 1 || throw(ArgumentError("0 < kt^2 < 1 required"))
    boundary.capacitance_per_patch_f > 0 || throw(ArgumentError("capacitance must be positive"))
    abs(boundary.pair_a_parity) == 1 || throw(ArgumentError("A parity must be +/-1"))
    abs(boundary.pair_b_parity) == 1 || throw(ArgumentError("B parity must be +/-1"))
    boundary
end

function validate(boundary::O8MechanicalBoundary)
    boundary.actuator_surface_mass_kg_m2 >= 0 || throw(ArgumentError("mass must be nonnegative"))
    0 <= boundary.off_stiffness_fraction <= 1 || throw(ArgumentError("invalid off stiffness"))
    boundary.stiffness_a_n_m3 > 0 || throw(ArgumentError("A stiffness must be positive"))
    boundary.stiffness_b_n_m3 > 0 || throw(ArgumentError("B stiffness must be positive"))
    boundary.loss_a >= 0 || throw(ArgumentError("A loss must be nonnegative"))
    boundary.loss_b >= 0 || throw(ArgumentError("B loss must be nonnegative"))
    boundary
end

function validate(load::O8Load)
    load.inductance_h > 0 || throw(ArgumentError("inductance must be positive"))
    load.resistance_ohm > 0 || throw(ArgumentError("resistance must be positive"))
    load
end

function stress(config, strain)
    mu = config.density_kg_m3 * config.shear_speed_m_s^2
    lambda = config.density_kg_m3 * config.pressure_speed_m_s^2 - 2mu
    identity = one(TensorValue{2, 2, ComplexF64})
    lambda * tr(strain) * identity + 2mu * strain
end

function absorbing_displacement(displacement, normal, config)
    normal_part = (displacement ⋅ normal) * normal
    tangent_part = displacement - normal_part
    config.density_kg_m3 * config.pressure_speed_m_s * normal_part +
        config.density_kg_m3 * config.shear_speed_m_s * tangent_part
end

function source_profile(port, y, height)
    port == 1 && return inv(sqrt(height))
    port == 2 && return sqrt(2 / height) * cos(2pi * y / height)
    port == 3 && return sqrt(2 / height) * sin(pi * y / height)
    throw(ArgumentError("only three frozen ports are allowed"))
end

function solve_matrix(matrix, right_hand_sides, trial_space)
    solver = LUSolver()
    numerical = numerical_setup(symbolic_setup(solver, matrix), matrix)
    coefficients = Matrix{ComplexF64}(undef, size(matrix, 1), length(right_hand_sides))
    temporary = similar(first(right_hand_sides))
    for index in eachindex(right_hand_sides)
        fill!(temporary, 0)
        solve!(temporary, numerical, right_hand_sides[index])
        coefficients[:, index] .= temporary
    end
    coefficients
end

function output_transfer(coefficients, trial_space, config)
    output = zeros(ComplexF64, config.output_bins, size(coefficients, 2))
    bin_height = config.aperture_height_m / config.output_bins
    x = config.block_length_m - 1.0e-10
    for port in axes(coefficients, 2)
        displacement = FEFunction(trial_space, coefficients[:, port])
        for bin in 1:config.output_bins
            lower = -config.aperture_height_m / 2 + (bin - 1) * bin_height
            y = range(
                lower + bin_height / (2config.sample_points_per_bin),
                lower + bin_height - bin_height / (2config.sample_points_per_bin);
                length=config.sample_points_per_bin,
            )
            output[bin, port] = mean(
                ComplexF64(displacement(Point(x, coordinate))[1]) for coordinate in y
            )
        end
    end
    output
end

function boundary_vector(test_space, measure, normal, length_m)
    vector = assemble_vector(v -> ∫(normal ⋅ v)measure, test_space)
    ComplexF64.(vector) ./ length_m
end

function boundary_unit_matrix(trial_space, test_space, measure, normal)
    assemble_matrix(
        (u, v) -> ∫((u ⋅ normal) * (v ⋅ normal))measure,
        trial_space, test_space,
    )
end

function sparse_outer(vector)
    sparse_vector = sparsevec(vector)
    sparse_vector * transpose(sparse_vector)
end

function electrical_coefficient(bit, omega, capacitance, load)
    bit == 0 && return inv(capacitance) + 0im
    impedance = load.resistance_ohm + im * omega * load.inductance_h
    admittance = inv(impedance)
    im * omega / (admittance + im * omega * capacitance)
end

function pair_electrical_data(
    plus_vector, minus_vector, patch_length, parity, pzt,
)
    surface_stiffness = pzt.c33e_pa / pzt.thickness_m
    short_stiffness = surface_stiffness * patch_length
    coupling_ratio = pzt.coupling_squared / (1 - pzt.coupling_squared)
    theta = sqrt(
        coupling_ratio * short_stiffness * pzt.capacitance_per_patch_f,
    )
    coupling_vector = theta .* (plus_vector .+ parity .* minus_vector)
    (
        coupling_vector,
        outer=sparse_outer(coupling_vector),
        pair_capacitance=2pzt.capacitance_per_patch_f,
        theta,
        short_stiffness,
    )
end

function state_projector_distance(first, second)
    first_svd = svd(first; full=false)
    second_svd = svd(second; full=false)
    first_q = first_svd.U[:, 1:size(first, 2)]
    second_q = second_svd.U[:, 1:size(second, 2)]
    norm(first_q * first_q' - second_q * second_q') / sqrt(2size(first, 2))
end

function solve_o8_frequency(
    mesh_path::AbstractString;
    mechanism::Symbol,
    config::O8PlaneStrainConfig,
    pzt::O8PZTBoundary=O8PZTBoundary(),
    load_a::Union{Nothing,O8Load}=nothing,
    load_b::Union{Nothing,O8Load}=nothing,
    mechanical::Union{Nothing,O8MechanicalBoundary}=nothing,
    zero_coupling::Bool=false,
    visualization_port::Union{Nothing,Int}=nothing,
)
    validate(config)
    isnothing(visualization_port) || 1 <= visualization_port <= 3 ||
        throw(ArgumentError("visualization port must be 1, 2, or 3"))
    mechanism in (:pzt, :mechanical) || throw(ArgumentError("unknown mechanism"))
    mechanism == :pzt && begin
        validate(pzt)
        isnothing(load_a) && error("PZT load A is required")
        isnothing(load_b) && error("PZT load B is required")
        validate(load_a)
        validate(load_b)
    end
    mechanism == :mechanical && begin
        isnothing(mechanical) && error("mechanical boundary is required")
        validate(mechanical)
    end
    isfile(mesh_path) || error("mesh not found: $mesh_path")

    model = GmshDiscreteModel(mesh_path)
    omega = 2pi * config.frequency_hz
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, config.element_order)
    test_space = TestFESpace(
        model, reference; conformity=:H1, vector_type=Vector{ComplexF64},
    )
    trial_space = TrialFESpace(test_space)
    domain = Triangulation(model; tags=["Block"])
    d_domain = Measure(domain, config.quadrature_degree)
    source = BoundaryTriangulation(model; tags=["Source"])
    radiation = BoundaryTriangulation(model; tags=["RadiationBoundary"])
    d_source = Measure(source, config.quadrature_degree)
    d_radiation = Measure(radiation, config.quadrature_degree)
    n_source = get_normal_vector(source)
    n_radiation = get_normal_vector(radiation)
    zone_names = ("ZoneAMinus", "ZoneAPlus", "ZoneBMinus", "ZoneBPlus")
    zones = [BoundaryTriangulation(model; tags=[name]) for name in zone_names]
    zone_measures = [Measure(zone, config.quadrature_degree) for zone in zones]
    zone_normals = get_normal_vector.(zones)
    zone_lengths = Float64[sum(∫(1.0)measure) for measure in zone_measures]

    function host_form(u, v)
        sigma = stress(config, ε(u))
        ∫(
            (1 + im / config.material_quality_factor) * (sigma ⊙ ε(v)) -
            config.density_kg_m3 * omega^2 * (u ⋅ v)
        )d_domain +
        ∫(im * omega * absorbing_displacement(u, n_source, config) ⋅ v)d_source +
        ∫(im * omega * absorbing_displacement(u, n_radiation, config) ⋅ v)d_radiation
    end
    host_matrix = assemble_matrix(host_form, trial_space, test_space)
    coordinate = get_physical_coordinate(source)
    source_function(port) = (point -> source_profile(
        port, point[2], config.aperture_height_m,
    )) ∘ coordinate
    source_linear(port) = v -> ∫(
        -config.source_pressure_scale_pa_sqrt_m * source_function(port) * (n_source ⋅ v)
    )d_source
    right_hand_sides = [
        assemble_vector(source_linear(port), test_space) for port in 1:3
    ]

    reference_coefficients = solve_matrix(host_matrix, right_hand_sides, trial_space)
    raw_reference_output = output_transfer(reference_coefficients, trial_space, config)
    reference_work = transpose(hcat(right_hand_sides...)) * reference_coefficients

    unit_matrices = [
        boundary_unit_matrix(
            trial_space, test_space, zone_measures[index], zone_normals[index],
        ) for index in eachindex(zones)
    ]
    device_base = copy(host_matrix)
    electrical_rows = NamedTuple[]
    state_matrices = SparseMatrixCSC{ComplexF64,Int}[]

    if mechanism == :pzt
        surface_stiffness = pzt.c33e_pa / pzt.thickness_m
        surface_mass = pzt.density_kg_m3 * pzt.thickness_m
        patch_coefficient = surface_stiffness * (1 + im / pzt.mechanical_quality_factor) -
            omega^2 * surface_mass
        for matrix in unit_matrices
            device_base += patch_coefficient .* matrix
        end
        vectors = [
            boundary_vector(
                test_space, zone_measures[index], zone_normals[index], zone_lengths[index],
            ) for index in eachindex(zones)
        ]
        pair_a = pair_electrical_data(
            vectors[2], vectors[1], mean(zone_lengths[1:2]), pzt.pair_a_parity, pzt,
        )
        pair_b = pair_electrical_data(
            vectors[4], vectors[3], mean(zone_lengths[3:4]), pzt.pair_b_parity, pzt,
        )
        for bits in ((0, 0), (1, 0), (0, 1), (1, 1))
            coefficient_a = zero_coupling ? 0.0 + 0im : electrical_coefficient(
                bits[1], omega, pair_a.pair_capacitance, load_a,
            )
            coefficient_b = zero_coupling ? 0.0 + 0im : electrical_coefficient(
                bits[2], omega, pair_b.pair_capacitance, load_b,
            )
            push!(state_matrices, device_base +
                coefficient_a .* pair_a.outer + coefficient_b .* pair_b.outer)
            push!(electrical_rows, (
                state="q$(bits[1])$(bits[2])",
                coefficient_a_real=real(coefficient_a),
                coefficient_a_imag=imag(coefficient_a),
                coefficient_b_real=real(coefficient_b),
                coefficient_b_imag=imag(coefficient_b),
                theta_a=pair_a.theta, theta_b=pair_b.theta,
            ))
        end
    else
        surface_mass = mechanical.actuator_surface_mass_kg_m2
        for matrix in unit_matrices
            device_base += (-omega^2 * surface_mass) .* matrix
        end
        matrix_a = unit_matrices[1] + unit_matrices[2]
        matrix_b = unit_matrices[3] + unit_matrices[4]
        for bits in ((0, 0), (1, 0), (0, 1), (1, 1))
            fraction_a = zero_coupling ? 0.0 :
                (bits[1] == 0 ? mechanical.off_stiffness_fraction : 1.0)
            fraction_b = zero_coupling ? 0.0 :
                (bits[2] == 0 ? mechanical.off_stiffness_fraction : 1.0)
            coefficient_a = fraction_a * mechanical.stiffness_a_n_m3 *
                (1 + im * mechanical.loss_a)
            coefficient_b = fraction_b * mechanical.stiffness_b_n_m3 *
                (1 + im * mechanical.loss_b)
            push!(state_matrices, device_base +
                coefficient_a .* matrix_a + coefficient_b .* matrix_b)
            push!(electrical_rows, (
                state="q$(bits[1])$(bits[2])",
                coefficient_a_real=real(coefficient_a),
                coefficient_a_imag=imag(coefficient_a),
                coefficient_b_real=real(coefficient_b),
                coefficient_b_imag=imag(coefficient_b),
                theta_a=NaN, theta_b=NaN,
            ))
        end
    end

    raw_state_outputs = zeros(ComplexF64, config.output_bins, 3, 4)
    port_work = zeros(ComplexF64, 3, 3, 4)
    field_coordinates = nothing
    field_cells = nothing
    field_ux = nothing
    field_uy = nothing
    for state in 1:4
        coefficients = solve_matrix(state_matrices[state], right_hand_sides, trial_space)
        raw_state_outputs[:, :, state] .= output_transfer(coefficients, trial_space, config)
        port_work[:, :, state] .= transpose(hcat(right_hand_sides...)) * coefficients
        if !isnothing(visualization_port)
            displacement = FEFunction(trial_space, coefficients[:, visualization_port])
            data = only(Gridap.Visualization.visualization_data(
                domain, "o8_state_field";
                order=1,
                cellfields=Dict("displacement" => displacement),
            ))
            coordinates = collect(Gridap.Geometry.get_node_coordinates(data.grid))
            values = collect(data.nodaldata["displacement"])
            if isnothing(field_coordinates)
                field_coordinates = coordinates
                field_cells = [
                    Int.(collect(ids))
                    for ids in Gridap.Geometry.get_cell_node_ids(data.grid)
                ]
                field_ux = zeros(ComplexF64, length(coordinates), 4)
                field_uy = zeros(ComplexF64, length(coordinates), 4)
            end
            field_ux[:, state] .= ComplexF64[value[1] for value in values]
            field_uy[:, state] .= ComplexF64[value[2] for value in values]
        end
    end
    parity = Diagonal(ComplexF64[1, 1, -1])
    project(operator) = 0.5 .* (
        operator + reverse(operator; dims=1) * parity
    )
    reference_output = project(raw_reference_output)
    state_outputs = similar(raw_state_outputs)
    for state in 1:4
        state_outputs[:, :, state] .= project(raw_state_outputs[:, :, state])
    end
    raw_mirror_errors = [
        norm(reverse(raw_state_outputs[:, :, state]; dims=1) -
            raw_state_outputs[:, :, state] * parity) /
            max(norm(raw_state_outputs[:, :, state]), eps(Float64))
        for state in 1:4
    ]
    projected_mirror_errors = [
        norm(reverse(state_outputs[:, :, state]; dims=1) -
            state_outputs[:, :, state] * parity) /
            max(norm(state_outputs[:, :, state]), eps(Float64))
        for state in 1:4
    ]
    reciprocity_errors = [
        norm(port_work[:, :, state] - transpose(port_work[:, :, state])) /
            max(norm(port_work[:, :, state]), eps(Float64))
        for state in 1:4
    ]
    amplitude_ratios = [
        norm(state_outputs[:, :, state]) / norm(reference_output) for state in 1:4
    ]
    state_distances = [
        state_projector_distance(
            state_outputs[:, :, first], state_outputs[:, :, second],
        ) for first in 1:4 for second in (first + 1):4
    ]
    accepted_work = [
        -imag(port_work[port, port, state])
        for port in 1:3, state in 1:4
    ]
    field_data = isnothing(visualization_port) ? nothing : (
        port=visualization_port,
        x_m=Float64[point[1] for point in field_coordinates],
        y_m=Float64[point[2] for point in field_coordinates],
        cells=field_cells,
        ux=field_ux,
        uy=field_uy,
    )
    (
        frequency_hz=config.frequency_hz, mechanism, raw_reference_output,
        reference_output, reference_work, raw_state_outputs, state_outputs,
        port_work, raw_mirror_errors, projected_mirror_errors,
        reciprocity_errors, amplitude_ratios, state_distances, accepted_work,
        electrical_rows, zone_lengths,
        field_data,
        matrix_size=size(host_matrix), matrix_nonzeros=nnz(host_matrix),
    )
end

end
