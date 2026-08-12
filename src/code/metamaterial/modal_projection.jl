module ElasticModalProjection

using Gridap

export BoundaryFieldSample,
       sample_boundary_field,
       counterpropagating_mode,
       reciprocity_product,
       decompose_mode_pair,
       project_boundary_mode_pair,
       decompose_propagating_modes

struct BoundaryFieldSample
    y_m::Vector{Float64}
    displacement_x::Vector{ComplexF64}
    displacement_y::Vector{ComplexF64}
    traction_xx_pa::Vector{ComplexF64}
    traction_xy_pa::Vector{ComplexF64}
end

"""
Sample distributed displacement and traction `sigma*e_x` on an external port.

The returned transverse coordinate is centred at zero so it uses the same
coordinate convention as `ElasticPortModes.PortMode`.
"""
function sample_boundary_field(state, tag::AbstractString; visualization_order::Integer=1)
    boundary = BoundaryTriangulation(state.model; tags=[String(tag)])
    traction_x = state.stress ⋅ VectorValue(1.0, 0.0)
    data = only(Gridap.Visualization.visualization_data(
        boundary,
        "modal_port";
        order=visualization_order,
        cellfields=Dict("displacement" => state.displacement, "traction_x" => traction_x),
    ))
    coordinates = collect(Gridap.Geometry.get_node_coordinates(data.grid))
    displacement = collect(data.nodaldata["displacement"])
    traction = collect(data.nodaldata["traction_x"])
    order = sortperm(coordinates; by=coordinate -> coordinate[2])
    y = Float64[coordinates[index][2] for index in order]
    center = (first(y) + last(y)) / 2
    BoundaryFieldSample(
        y .- center,
        ComplexF64[displacement[index][1] for index in order],
        ComplexF64[displacement[index][2] for index in order],
        ComplexF64[traction[index][1] for index in order],
        ComplexF64[traction[index][2] for index in order],
    )
end

function linear_interpolate(x, values, query)
    first(x) - 10eps(Float64) <= query <= last(x) + 10eps(Float64) ||
        throw(ArgumentError("projection point lies outside the modal cross-section"))
    query <= first(x) && return values[1]
    query >= last(x) && return values[end]
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function sample_mode(mode, y)
    (
        displacement_x=ComplexF64[linear_interpolate(mode.y_m, mode.displacement_x, value) for value in y],
        displacement_y=ComplexF64[linear_interpolate(mode.y_m, mode.displacement_y, value) for value in y],
        traction_xx=ComplexF64[linear_interpolate(mode.y_m, mode.traction_xx_pa, value) for value in y],
        traction_xy=ComplexF64[linear_interpolate(mode.y_m, mode.traction_xy_pa, value) for value in y],
    )
end

function trapezoidal_integral(y, values)
    result = zero(eltype(values))
    for index in 1:(length(y) - 1)
        result += (y[index + 1] - y[index]) * (values[index] + values[index + 1]) / 2
    end
    result
end

"""Unconjugated elastic reciprocity product B(first, second)."""
function reciprocity_product(y, first, second)
    integrand = (
        first.displacement_x .* second.traction_xx .+
        first.displacement_y .* second.traction_xy .-
        first.traction_xx .* second.displacement_x .-
        first.traction_xy .* second.displacement_y
    )
    trapezoidal_integral(y, integrand)
end

function mode_fields(mode, y)
    sampled = sample_mode(mode, y)
    (
        displacement_x=sampled.displacement_x,
        displacement_y=sampled.displacement_y,
        traction_xx=sampled.traction_xx,
        traction_xy=sampled.traction_xy,
    )
end

function boundary_fields(sample::BoundaryFieldSample)
    (
        displacement_x=sample.displacement_x,
        displacement_y=sample.displacement_y,
        traction_xx=sample.traction_xx_pa,
        traction_xy=sample.traction_xy_pa,
    )
end

"""Find the opposite-going propagating partner with matching parity and |k|."""
function counterpropagating_mode(mode, modes)
    target_direction = mode.direction == :right ? :left : :right
    candidates = filter(
        candidate -> candidate.kind == :propagating &&
                     candidate.direction == target_direction &&
                     candidate.parity == mode.parity,
        modes,
    )
    isempty(candidates) && error("no counterpropagating partner for k=$(mode.wavenumber_per_m)")
    first(sort(candidates; by=candidate -> abs(
        candidate.wavenumber_per_m + mode.wavenumber_per_m,
    )))
end

"""
Decompose a cross-sectional field into one right/left mode pair.

The formula is the biorthogonal Lorentz-reciprocity projection. No complex
conjugation is used in `B`; this is distinct from the Hermitian power product.
"""
function decompose_mode_pair(sample::BoundaryFieldSample, right_mode, left_mode)
    y = sample.y_m
    total = boundary_fields(sample)
    plus = mode_fields(right_mode, y)
    minus = mode_fields(left_mode, y)
    normalization = reciprocity_product(y, plus, minus)
    abs(normalization) > 100eps(Float64) || error("degenerate reciprocity normalization")
    right_amplitude = reciprocity_product(y, total, minus) / normalization
    left_amplitude = reciprocity_product(y, plus, total) / normalization
    (
        right_amplitude=ComplexF64(right_amplitude),
        left_amplitude=ComplexF64(left_amplitude),
        normalization=ComplexF64(normalization),
    )
end

function boundary_center_y(state, tag)
    boundary = BoundaryTriangulation(state.model; tags=[String(tag)])
    data = only(Gridap.Visualization.visualization_data(boundary, "port_coordinates"; order=1))
    coordinates = collect(Gridap.Geometry.get_node_coordinates(data.grid))
    y = getindex.(coordinates, 2)
    (minimum(y) + maximum(y)) / 2
end

function modal_cell_fields(boundary, mode, center_y_m)
    coordinate = get_physical_coordinate(boundary)
    displacement_function = point -> VectorValue(
        linear_interpolate(mode.y_m, mode.displacement_x, point[2] - center_y_m),
        linear_interpolate(mode.y_m, mode.displacement_y, point[2] - center_y_m),
    )
    traction_function = point -> VectorValue(
        linear_interpolate(mode.y_m, mode.traction_xx_pa, point[2] - center_y_m),
        linear_interpolate(mode.y_m, mode.traction_xy_pa, point[2] - center_y_m),
    )
    displacement_function ∘ coordinate, traction_function ∘ coordinate
end

"""
Project Gridap boundary fields directly by quadrature.

Unlike `sample_boundary_field`, this method does not nodally average the FEM
stress and is therefore the preferred path for quantitative S-parameters.
"""
function project_boundary_mode_pair(
    state,
    tag::AbstractString,
    right_mode,
    left_mode;
    quadrature_degree::Integer=4,
)
    boundary = BoundaryTriangulation(state.model; tags=[String(tag)])
    measure = Measure(boundary, quadrature_degree)
    center_y_m = boundary_center_y(state, tag)
    plus_u, plus_t = modal_cell_fields(boundary, right_mode, center_y_m)
    minus_u, minus_t = modal_cell_fields(boundary, left_mode, center_y_m)
    traction_x = state.stress ⋅ VectorValue(1.0, 0.0)
    normalization = sum(∫(plus_u ⋅ minus_t - plus_t ⋅ minus_u)measure)
    right_numerator = sum(∫(state.displacement ⋅ minus_t - traction_x ⋅ minus_u)measure)
    left_numerator = sum(∫(plus_u ⋅ traction_x - plus_t ⋅ state.displacement)measure)
    abs(normalization) > 100eps(Float64) || error("degenerate reciprocity normalization")
    (
        right_amplitude=ComplexF64(right_numerator / normalization),
        left_amplitude=ComplexF64(left_numerator / normalization),
        normalization=ComplexF64(normalization),
    )
end

function decompose_propagating_modes(sample::BoundaryFieldSample, modes)
    right_modes = filter(
        mode -> mode.kind == :propagating && mode.direction == :right,
        modes,
    )
    [
        (
            right_mode=right_mode,
            left_mode=counterpropagating_mode(right_mode, modes),
            decompose_mode_pair(
                sample,
                right_mode,
                counterpropagating_mode(right_mode, modes),
            )...,
        )
        for right_mode in right_modes
    ]
end

end
