module LensDesign

export LibraryEntry,
       LensConfig,
       LensObjectiveConfig,
       LensSelection,
       load_gap_library,
       lens_centers_mm,
       required_phase_span,
       select_lens,
       select_lens_contrast,
       selection_objective_metrics,
       field_at,
       field_grid,
       reference_field_grid,
       wrapped_phase_difference,
       phase_design_values

struct LibraryEntry
    case_id::String
    periods::Int
    length_mm::Float64
    gap_mm::Float64
    amplitude_mm::Float64
    transfer::ComplexF64
end

Base.@kwdef struct LensConfig
    frequency_hz::Float64 = 220.0e3
    longitudinal_speed_m_per_s::Float64 = 2340.0
    element_count::Int = 7
    element_width_mm::Float64 = 7.0
    slot_width_mm::Float64 = 1.2
    focal_distance_mm::Float64 = 60.0
    aperture_quadrature_points::Int = 9
end

Base.@kwdef struct LensObjectiveConfig
    axial_offsets_mm::Vector{Float64} = [-30.0, -24.0, -18.0, -14.0, 14.0, 18.0, 24.0, 30.0]
    transverse_offsets_mm::Vector{Float64} = [-24.0, -16.0, -12.0, 12.0, 16.0, 24.0]
    axial_weight::Float64 = 0.35
    sidelobe_weight::Float64 = 0.15
    minimum_focus_amplitude::Float64 = 0.0
end

struct LensSelection
    centers_mm::Vector{Float64}
    entries::Vector{LibraryEntry}
    focus_field::ComplexF64
end

wavelength_mm(config::LensConfig) =
    1000.0 * config.longitudinal_speed_m_per_s / config.frequency_hz

wavenumber_per_mm(config::LensConfig) = 2.0 * pi / wavelength_mm(config)

function validate_config(config::LensConfig)
    isodd(config.element_count) || throw(ArgumentError("element_count must be odd"))
    config.element_count >= 3 || throw(ArgumentError("at least three elements are required"))
    config.element_width_mm > 0 || throw(ArgumentError("element width must be positive"))
    config.slot_width_mm >= 0 || throw(ArgumentError("slot width must be non-negative"))
    config.focal_distance_mm > 0 || throw(ArgumentError("focal distance must be positive"))
    config.aperture_quadrature_points > 0 ||
        throw(ArgumentError("aperture quadrature count must be positive"))
    nothing
end

function load_gap_library(path::AbstractString)
    isfile(path) || error("gap library not found: $path")
    lines = readlines(path)
    length(lines) >= 2 || error("gap library is empty: $path")
    header = split(first(lines), ',')
    column = Dict(name => index for (index, name) in enumerate(header))
    required = (
        "case_id",
        "periods",
        "length_mm",
        "gap_mm",
        "amplitude_mm",
        "H_real",
        "H_imag",
        "valid",
    )
    all(haskey(column, name) for name in required) ||
        error("gap library has an incompatible schema")

    entries = LibraryEntry[]
    for line in Iterators.drop(lines, 1)
        isempty(strip(line)) && continue
        values = split(line, ',')
        lowercase(values[column["valid"]]) == "true" || continue
        push!(entries, LibraryEntry(
            values[column["case_id"]],
            parse(Int, values[column["periods"]]),
            parse(Float64, values[column["length_mm"]]),
            parse(Float64, values[column["gap_mm"]]),
            parse(Float64, values[column["amplitude_mm"]]),
            complex(
                parse(Float64, values[column["H_real"]]),
                parse(Float64, values[column["H_imag"]]),
            ),
        ))
    end
    isempty(entries) && error("gap library contains no valid entries")
    entries
end

function lens_centers_mm(config::LensConfig=LensConfig())
    validate_config(config)
    pitch_mm = config.element_width_mm + config.slot_width_mm
    half_count = (config.element_count - 1) ÷ 2
    collect((-half_count):half_count) .* pitch_mm
end

function required_phase_span(config::LensConfig=LensConfig())
    centers = lens_centers_mm(config)
    path_difference_mm = maximum(
        hypot(config.focal_distance_mm, center) - config.focal_distance_mm
        for center in centers
    )
    wavenumber_per_mm(config) * path_difference_mm
end

wrapped_phase_difference(a::Real, b::Real) = mod(a - b + pi, 2.0 * pi) - pi

function propagation_kernel(
    x_mm::Real,
    source_y_mm::Real,
    observation_y_mm::Real,
    config::LensConfig,
)
    x_mm > 0 || throw(ArgumentError("observation x must be positive"))
    distance_mm = hypot(x_mm, observation_y_mm - source_y_mm)
    cis(wavenumber_per_mm(config) * distance_mm - pi / 4.0) /
        sqrt(wavelength_mm(config) * distance_mm)
end

function element_kernel(
    center_y_mm::Real,
    x_mm::Real,
    observation_y_mm::Real,
    config::LensConfig,
)
    count = config.aperture_quadrature_points
    dy_mm = config.element_width_mm / count
    first_y_mm = center_y_mm - config.element_width_mm / 2.0 + dy_mm / 2.0
    dy_mm * sum(
        propagation_kernel(
            x_mm,
            first_y_mm + (index - 1) * dy_mm,
            observation_y_mm,
            config,
        )
        for index in 1:count
    )
end

function select_lens(entries::AbstractVector{LibraryEntry}; config::LensConfig=LensConfig())
    validate_config(config)
    centers = lens_centers_mm(config)
    absolute_centers = sort(unique(abs.(centers)))
    multiplicities = [iszero(center) ? 1 : 2 for center in absolute_centers]
    kernels = [
        element_kernel(center, config.focal_distance_mm, 0.0, config)
        for center in absolute_centers
    ]
    choices = ntuple(_ -> eachindex(entries), length(absolute_centers))

    best_amplitude = -Inf
    best_indices = nothing
    best_field = 0.0 + 0.0im
    for indices in Iterators.product(choices...)
        focus_field = sum(
            multiplicities[index] * entries[indices[index]].transfer * kernels[index]
            for index in eachindex(absolute_centers)
        )
        amplitude = abs(focus_field)
        if amplitude > best_amplitude
            best_amplitude = amplitude
            best_indices = Tuple(indices)
            best_field = focus_field
        end
    end

    selected_entries = [
        entries[best_indices[findfirst(==(abs(center)), absolute_centers)]]
        for center in centers
    ]
    LensSelection(centers, selected_entries, best_field)
end

function objective_observations(config::LensConfig, objective::LensObjectiveConfig)
    axial = [
        (config.focal_distance_mm + offset, 0.0)
        for offset in objective.axial_offsets_mm
        if config.focal_distance_mm + offset > 0
    ]
    transverse = [(config.focal_distance_mm, offset) for offset in objective.transverse_offsets_mm]
    vcat([(config.focal_distance_mm, 0.0)], axial, transverse), length(axial)
end

function symmetric_group_kernel(center, x, y, config)
    iszero(center) ?
    element_kernel(0.0, x, y, config) :
    element_kernel(center, x, y, config) + element_kernel(-center, x, y, config)
end

function select_lens_contrast(
    entries::AbstractVector{LibraryEntry};
    config::LensConfig=LensConfig(),
    objective::LensObjectiveConfig=LensObjectiveConfig(),
)
    validate_config(config)
    isempty(entries) && throw(ArgumentError("element library must not be empty"))
    objective.axial_weight >= 0 || throw(ArgumentError("axial weight must be non-negative"))
    objective.sidelobe_weight >= 0 || throw(ArgumentError("sidelobe weight must be non-negative"))
    objective.minimum_focus_amplitude >= 0 ||
        throw(ArgumentError("minimum focus amplitude must be non-negative"))

    centers = lens_centers_mm(config)
    absolute_centers = sort(unique(abs.(centers)))
    observations, axial_count = objective_observations(config, objective)
    contribution = Array{ComplexF64}(undef, length(observations), length(absolute_centers), length(entries))
    for (observation_index, (x, y)) in enumerate(observations)
        for (center_index, center) in enumerate(absolute_centers)
            kernel = symmetric_group_kernel(center, x, y, config)
            for entry_index in eachindex(entries)
                contribution[observation_index, center_index, entry_index] =
                    entries[entry_index].transfer * kernel
            end
        end
    end

    choices = ntuple(_ -> eachindex(entries), length(absolute_centers))
    best_score = -Inf
    best_indices = nothing
    best_focus = 0.0 + 0.0im
    for indices in Iterators.product(choices...)
        focus = sum(
            contribution[1, group_index, indices[group_index]]
            for group_index in eachindex(absolute_centers)
        )
        axial_mean = axial_count == 0 ? 0.0 : sum(
            abs2(sum(
                contribution[observation_index, group_index, indices[group_index]]
                for group_index in eachindex(absolute_centers)
            ))
            for observation_index in 2:(1 + axial_count)
        ) / axial_count
        transverse_count = length(observations) - 1 - axial_count
        transverse_mean = transverse_count == 0 ? 0.0 : sum(
            abs2(sum(
                contribution[observation_index, group_index, indices[group_index]]
                for group_index in eachindex(absolute_centers)
            ))
            for observation_index in (2 + axial_count):length(observations)
        ) / transverse_count
        focus_intensity = abs2(focus)
        penalty = objective.axial_weight * axial_mean +
                  objective.sidelobe_weight * transverse_mean
        score = if objective.minimum_focus_amplitude > 0
            abs(focus) >= objective.minimum_focus_amplitude ?
            -penalty + 1.0e-9 * focus_intensity : -Inf
        else
            focus_intensity - penalty
        end
        if score > best_score
            best_score = score
            best_indices = Tuple(indices)
            best_focus = focus
        end
    end

    isnothing(best_indices) && error("no lens satisfies the minimum focus amplitude")
    selected_entries = [
        entries[best_indices[findfirst(==(abs(center)), absolute_centers)]]
        for center in centers
    ]
    LensSelection(centers, selected_entries, best_focus)
end

function selection_objective_metrics(
    selection::LensSelection;
    config::LensConfig=LensConfig(),
    objective::LensObjectiveConfig=LensObjectiveConfig(),
)
    axial_points = [
        (config.focal_distance_mm + offset, 0.0)
        for offset in objective.axial_offsets_mm
        if config.focal_distance_mm + offset > 0
    ]
    transverse_points = [(config.focal_distance_mm, offset) for offset in objective.transverse_offsets_mm]
    focus_intensity = abs2(field_at(config.focal_distance_mm, 0.0, selection; config))
    axial_mean_intensity = isempty(axial_points) ? 0.0 : sum(
        abs2(field_at(x, y, selection; config)) for (x, y) in axial_points
    ) / length(axial_points)
    sidelobe_mean_intensity = isempty(transverse_points) ? 0.0 : sum(
        abs2(field_at(x, y, selection; config)) for (x, y) in transverse_points
    ) / length(transverse_points)
    (
        focus_intensity,
        axial_mean_intensity,
        sidelobe_mean_intensity,
        score=objective.minimum_focus_amplitude > 0 ?
              (sqrt(focus_intensity) >= objective.minimum_focus_amplitude ?
               -(objective.axial_weight * axial_mean_intensity +
                 objective.sidelobe_weight * sidelobe_mean_intensity) : -Inf) :
              focus_intensity - objective.axial_weight * axial_mean_intensity -
              objective.sidelobe_weight * sidelobe_mean_intensity,
    )
end

function field_at(
    x_mm::Real,
    y_mm::Real,
    selection::LensSelection;
    config::LensConfig=LensConfig(),
)
    sum(
        entry.transfer * element_kernel(center, x_mm, y_mm, config)
        for (center, entry) in zip(selection.centers_mm, selection.entries)
    )
end

function field_grid(
    x_mm::AbstractVector,
    y_mm::AbstractVector,
    selection::LensSelection;
    config::LensConfig=LensConfig(),
)
    [field_at(x, y, selection; config) for y in y_mm, x in x_mm]
end

function reference_field_grid(
    x_mm::AbstractVector,
    y_mm::AbstractVector;
    config::LensConfig=LensConfig(),
)
    centers = lens_centers_mm(config)
    [
        sum(element_kernel(center, x, y, config) for center in centers)
        for y in y_mm, x in x_mm
    ]
end

function phase_design_values(
    selection::LensSelection;
    config::LensConfig=LensConfig(),
)
    base_required = [
        -wavenumber_per_mm(config) *
        (hypot(config.focal_distance_mm, center) - config.focal_distance_mm)
        for center in selection.centers_mm
    ]
    actual = angle.(getproperty.(selection.entries, :transfer))
    weights = abs.(getproperty.(selection.entries, :transfer))
    offset = angle(sum(
        weight * cis(actual_phase - required_phase)
        for (weight, actual_phase, required_phase) in zip(weights, actual, base_required)
    ))
    target = base_required .+ offset
    selected = [
        target_phase + wrapped_phase_difference(actual_phase, target_phase)
        for (actual_phase, target_phase) in zip(actual, target)
    ]
    errors = selected .- target
    (target=target, selected=selected, errors=errors, offset=offset)
end

end
