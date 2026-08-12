module ConservativeElasticModes

using ArnoldiMethod
using Gridap
using LinearAlgebra
using SparseArrays

export ModeConfig,
       ElasticMode,
       solve_conservative_modes,
       side_mass_mode_metrics,
       inline_mass_mode_metrics,
       select_inline_mass_mode,
       shift_invert_eigenpairs

Base.@kwdef struct ModeConfig
    density::Float64 = 1210.0
    pressure_wave_speed::Float64 = 2340.0
    shear_wave_speed::Float64 = 1170.0
    element_order::Int = 1
    quadrature_degree::Int = 2
    target_frequency_hz::Float64 = 242.0e3
    mode_count::Int = 12
    tolerance::Float64 = 1.0e-9
    clamp_tags::Vector{String} = ["Source", "Microphone"]
end

struct ElasticMode
    frequency_hz::Float64
    eigenvalue_rad2_s2::Float64
    relative_residual::Float64
    parity_y::Float64
    longitudinal_fraction::Float64
    node_x_m::Vector{Float64}
    node_y_m::Vector{Float64}
    node_area_m2::Vector{Float64}
    cell_node_ids::Vector{Vector{Int}}
    displacement_x::Vector{Float64}
    displacement_y::Vector{Float64}
end

struct ShiftInvertMap{TF, TM, TV}
    factorization::TF
    mass::TM
    temporary::TV
end

Base.size(operator::ShiftInvertMap) = size(operator.mass)
Base.size(operator::ShiftInvertMap, dimension::Integer) = size(operator.mass, dimension)
Base.eltype(operator::ShiftInvertMap) = eltype(operator.mass)

function LinearAlgebra.mul!(output, operator::ShiftInvertMap, input)
    mul!(operator.temporary, operator.mass, input)
    ldiv!(output, operator.factorization, operator.temporary)
    output
end

"""
Return generalized eigenpairs of `K x = lambda M x` closest to `target_lambda`.

The implementation applies Arnoldi iteration to
`(K-target_lambda*M)^(-1) M`; it therefore keeps the assembled FEM matrices
sparse and never forms a dense inverse.
"""
function shift_invert_eigenpairs(
    stiffness::AbstractMatrix{<:Real},
    mass::AbstractMatrix{<:Real},
    target_lambda::Real;
    count::Integer=8,
    tolerance::Real=1.0e-9,
)
    size(stiffness) == size(mass) || throw(DimensionMismatch("K and M must have equal size"))
    size(stiffness, 1) > count || throw(ArgumentError("mode count must be smaller than matrix size"))
    shifted = stiffness - Float64(target_lambda) * mass
    operator = ShiftInvertMap(
        lu(shifted),
        mass,
        zeros(eltype(stiffness), size(stiffness, 1)),
    )
    krylov_min = min(size(stiffness, 1), max(2Int(count) + 4, 20))
    krylov_max = min(size(stiffness, 1), max(2krylov_min, 40))
    decomposition, history = partialschur(
        operator;
        nev=Int(count),
        which=:LM,
        tol=Float64(tolerance),
        mindim=krylov_min,
        maxdim=krylov_max,
        restarts=300,
    )
    history.converged || @warn "Arnoldi iteration did not fully converge" history
    inverse_eigenvalues, eigenvectors = partialeigen(decomposition)
    eigenvalues = Float64(target_lambda) .+ 1.0 ./ inverse_eigenvalues
    order = sortperm(real.(eigenvalues))
    real.(eigenvalues[order]), real.(eigenvectors[:, order])
end

function assemble_elastic_matrices(model, config::ModeConfig)
    rho = config.density
    mu = rho * config.shear_wave_speed^2
    lambda = rho * config.pressure_wave_speed^2 - 2mu
    identity_tensor = one(TensorValue{2, 2, Float64})
    sigma(strain) = lambda * tr(strain) * identity_tensor + 2mu * strain

    reference_element = ReferenceFE(
        lagrangian,
        VectorValue{2, Float64},
        config.element_order,
    )
    test_space = TestFESpace(
        model,
        reference_element;
        conformity=:H1,
        dirichlet_tags=config.clamp_tags,
    )
    trial_space = TrialFESpace(test_space)
    domain = Triangulation(model)
    domain_measure = Measure(domain, config.quadrature_degree)
    zero_vector = VectorValue(0.0, 0.0)

    stiffness_form(u, v) = ∫(sigma(ε(u)) ⊙ ε(v))domain_measure
    mass_form(u, v) = ∫(rho * (u ⋅ v))domain_measure
    zero_form(v) = ∫(zero_vector ⋅ v)domain_measure

    stiffness = get_matrix(AffineFEOperator(
        stiffness_form,
        zero_form,
        trial_space,
        test_space,
    ))
    mass = get_matrix(AffineFEOperator(
        mass_form,
        zero_form,
        trial_space,
        test_space,
    ))
    stiffness, mass, trial_space, domain
end

function polygon_area(coordinates, node_ids)
    area_twice = 0.0
    for local_index in eachindex(node_ids)
        first_node = coordinates[node_ids[local_index]]
        second_node = coordinates[node_ids[mod1(local_index + 1, length(node_ids))]]
        area_twice += first_node[1] * second_node[2] - second_node[1] * first_node[2]
    end
    abs(area_twice) / 2
end

function nodal_area_weights(coordinates, cells)
    weights = zeros(Float64, length(coordinates))
    for node_ids in cells
        area_share = polygon_area(coordinates, node_ids) / length(node_ids)
        for node_id in node_ids
            weights[node_id] += area_share
        end
    end
    weights
end

function parity_score(x, y, ux, uy, weights)
    tolerance = max(maximum(x) - minimum(x), maximum(y) - minimum(y)) * 1.0e-8
    center_y = (minimum(y) + maximum(y)) / 2
    coordinate_key(x_value, y_value) = (
        round(Int, x_value / tolerance),
        round(Int, y_value / tolerance),
    )
    lookup = Dict(coordinate_key(x[index], y[index]) => index for index in eachindex(x))
    numerator = 0.0
    denominator = 0.0
    for index in eachindex(x)
        mirror_index = get(lookup, coordinate_key(x[index], 2center_y - y[index]), 0)
        mirror_index == 0 && continue
        weight = weights[index]
        numerator += weight * (ux[index] * ux[mirror_index] - uy[index] * uy[mirror_index])
        denominator += weight * (ux[index]^2 + uy[index]^2)
    end
    denominator > 0 ? clamp(numerator / denominator, -1.0, 1.0) : NaN
end

function visualization_mode(trial_space, domain, free_values)
    field = FEFunction(trial_space, free_values)
    data = only(Gridap.Visualization.visualization_data(
        domain,
        "elastic_mode";
        order=1,
        cellfields=Dict("displacement" => field),
    ))
    coordinates = collect(Gridap.Geometry.get_node_coordinates(data.grid))
    cells = [Int.(collect(ids)) for ids in Gridap.Geometry.get_cell_node_ids(data.grid)]
    displacement = collect(data.nodaldata["displacement"])
    x = Float64[coordinate[1] for coordinate in coordinates]
    y = Float64[coordinate[2] for coordinate in coordinates]
    ux = Float64[value[1] for value in displacement]
    uy = Float64[value[2] for value in displacement]
    weights = nodal_area_weights(coordinates, cells)
    total = sum(weights .* (ux .^ 2 .+ uy .^ 2))
    longitudinal = total > 0 ? sum(weights .* ux .^ 2) / total : NaN
    (; x, y, ux, uy, weights, cells, longitudinal,
       parity=parity_score(x, y, ux, uy, weights))
end

function solve_conservative_modes(
    model_path::AbstractString;
    config::ModeConfig=ModeConfig(),
)
    isfile(model_path) || error("Gridap model not found: $model_path")
    model = DiscreteModelFromFile(model_path)
    stiffness, mass, trial_space, domain = assemble_elastic_matrices(model, config)
    target_lambda = (2pi * config.target_frequency_hz)^2
    eigenvalues, eigenvectors = shift_invert_eigenpairs(
        stiffness,
        mass,
        target_lambda;
        count=config.mode_count,
        tolerance=config.tolerance,
    )

    modes = ElasticMode[]
    for column in axes(eigenvectors, 2)
        eigenvalue = eigenvalues[column]
        eigenvalue > 0 || continue
        free_values = copy(eigenvectors[:, column])
        mass_norm = sqrt(abs(dot(free_values, mass * free_values)))
        mass_norm > 0 || continue
        free_values ./= mass_norm
        largest_index = argmax(abs.(free_values))
        free_values[largest_index] < 0 && (free_values .*= -1)
        residual_vector = stiffness * free_values - eigenvalue * (mass * free_values)
        denominator = norm(stiffness * free_values) + abs(eigenvalue) * norm(mass * free_values)
        residual = denominator > 0 ? norm(residual_vector) / denominator : NaN
        visual = visualization_mode(trial_space, domain, free_values)
        push!(modes, ElasticMode(
            sqrt(eigenvalue) / (2pi),
            eigenvalue,
            residual,
            visual.parity,
            visual.longitudinal,
            visual.x,
            visual.y,
            visual.weights,
            visual.cells,
            visual.ux,
            visual.uy,
        ))
    end
    sort!(modes; by=mode -> mode.frequency_hz)
end

function rectangle_mask(mode::ElasticMode, center_x_mm, length_mm, y_min_mm, y_max_mm)
    x_mm = mode.node_x_m .* 1e3
    y_mm = mode.node_y_m .* 1e3
    half_length = length_mm / 2
    (x_mm .>= center_x_mm - half_length) .&
        (x_mm .<= center_x_mm + half_length) .&
        (y_mm .>= y_min_mm) .&
        (y_mm .<= y_max_mm)
end

function weighted_fraction(mode::ElasticMode, mask)
    intensity = mode.node_area_m2 .* (mode.displacement_x .^ 2 .+ mode.displacement_y .^ 2)
    denominator = sum(intensity)
    denominator > 0 ? sum(intensity[mask]) / denominator : NaN
end

"""Diagnostics for a central inline mass joined to two axial necks."""
function inline_mass_mode_metrics(mode::ElasticMode, config)
    x_mm = mode.node_x_m .* 1e3
    y_mm = mode.node_y_m .* 1e3
    tolerance = 1.0e-6
    mass_mask =
        (x_mm .>= config.neck_length_mm - tolerance) .&
        (x_mm .<= config.neck_length_mm + config.mass_length_mm + tolerance) .&
        (abs.(y_mm) .<= config.mass_height_mm / 2 + tolerance)
    weights = mode.node_area_m2
    ux = mode.displacement_x
    uy = mode.displacement_y
    intensity = weights .* (ux .^ 2 .+ uy .^ 2)
    total = sum(intensity)
    mass_fraction = total > 0 ? sum(intensity[mass_mask]) / total : NaN
    mass_weight = sum(weights[mass_mask])
    ux_rms = mass_weight > 0 ?
             sqrt(sum(weights[mass_mask] .* ux[mass_mask] .^ 2) / mass_weight) : 0.0
    ux_mean = mass_weight > 0 ?
              sum(weights[mass_mask] .* ux[mass_mask]) / mass_weight : 0.0
    translation_coherence = ux_rms > 0 ? clamp(abs(ux_mean) / ux_rms, 0.0, 1.0) : 0.0
    even_factor = clamp((mode.parity_y + 1) / 2, 0.0, 1.0)
    score = mass_fraction * translation_coherence * mode.longitudinal_fraction * even_factor
    (; mass_fraction, translation_coherence, score)
end

"""Select the y-even, longitudinal, approximately rigid translation of an inline mass."""
function select_inline_mass_mode(modes, metrics)
    length(modes) == length(metrics) ||
        throw(DimensionMismatch("modes and metrics must have equal length"))
    candidates = findall(eachindex(modes)) do index
        mode = modes[index]
        metric = metrics[index]
        mode.parity_y >= 0.75 &&
            mode.longitudinal_fraction >= 0.55 &&
            metric.translation_coherence >= 0.65
    end
    isempty(candidates) && (candidates = collect(eachindex(modes)))
    candidates[argmax(metrics[index].score for index in candidates)]
end

"""Approximate kinetic-energy fractions in the side masses from nodal quadrature."""
function side_mass_mode_metrics(mode::ElasticMode, config; variant::Symbol=:bd)
    variant in (:b, :bd) || throw(ArgumentError("variant must be :b or :bd"))
    bright_top_min = config.height_mm - config.wall_margin_mm - config.bright_mass_height_mm
    bright_top_max = config.height_mm - config.wall_margin_mm
    bright_bottom_min = config.wall_margin_mm
    bright_bottom_max = config.wall_margin_mm + config.bright_mass_height_mm
    bright_mask = rectangle_mask(
        mode,
        config.bright_center_x_mm,
        config.bright_mass_length_mm,
        bright_top_min,
        bright_top_max,
    ) .| rectangle_mask(
        mode,
        config.bright_center_x_mm,
        config.bright_mass_length_mm,
        bright_bottom_min,
        bright_bottom_max,
    )
    dark_fraction = 0.0
    if variant == :bd
        dark_top_min = config.height_mm - config.wall_margin_mm - config.dark_mass_height_mm
        dark_top_max = config.height_mm - config.wall_margin_mm
        dark_bottom_min = config.wall_margin_mm
        dark_bottom_max = config.wall_margin_mm + config.dark_mass_height_mm
        dark_mask = rectangle_mask(
            mode,
            config.dark_center_x_mm,
            config.dark_mass_length_mm,
            dark_top_min,
            dark_top_max,
        ) .| rectangle_mask(
            mode,
            config.dark_center_x_mm,
            config.dark_mass_length_mm,
            dark_bottom_min,
            dark_bottom_max,
        )
        dark_fraction = weighted_fraction(mode, dark_mask)
    end
    (
        bright_fraction=weighted_fraction(mode, bright_mask),
        dark_fraction=dark_fraction,
        side_mass_fraction=weighted_fraction(mode, bright_mask) + dark_fraction,
    )
end

end
