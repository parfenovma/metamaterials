module MaterialLensMesher

using Gmsh: gmsh

if !isdefined(parentmodule(@__MODULE__), :SinusoidalMaterialLens)
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
end
using ..SinusoidalMaterialLens

export MaterialLensMeshConfig, build_material_lens_mesh

Base.@kwdef struct MaterialLensMeshConfig
    slot_width_mm::Float64 = 1.2
    output_length_mm::Float64 = 105.0
    lens_mesh_size_mm::Float64 = 0.30
    output_mesh_size_mm::Float64 = 1.25
    transition_distance_mm::Float64 = 12.0
    defect_radius_mm::Float64 = 0.0
    defect_center_x_mm::Float64 = 0.0
    defect_center_y_mm::Float64 = 0.0
    defect_mesh_size_mm::Float64 = 0.35
    defect_transition_distance_mm::Float64 = 6.0
    boundary_tolerance_mm::Float64 = 1.0e-5
end

function validate(cells, matching_thickness_mm, config)
    isempty(cells) && throw(ArgumentError("at least one lens element is required"))
    isodd(length(cells)) || throw(ArgumentError("a symmetric lens needs an odd element count"))
    length_value = first(cells).length_mm
    height_value = first(cells).height_mm
    all(cell -> cell.length_mm == length_value, cells) ||
        throw(ArgumentError("all element exits must share one x plane"))
    all(cell -> cell.height_mm == height_value, cells) ||
        throw(ArgumentError("all elements must have the same height"))
    matching_thickness_mm >= 0 || throw(ArgumentError("matching thickness must be non-negative"))
    config.slot_width_mm >= 0 || throw(ArgumentError("slot width must be non-negative"))
    config.output_length_mm > 0 || throw(ArgumentError("output length must be positive"))
    config.defect_radius_mm >= 0 || throw(ArgumentError("defect radius must be non-negative"))
    config.defect_mesh_size_mm > 0 || throw(ArgumentError("defect mesh size must be positive"))
    nothing
end

function polygon_surface(points_mm)
    clean_points = Tuple{Float64, Float64}[]
    for point in points_mm
        value = (Float64(point[1]), Float64(point[2]))
        (isempty(clean_points) || value != last(clean_points)) && push!(clean_points, value)
    end
    first(clean_points) == last(clean_points) && pop!(clean_points)
    point_tags = [gmsh.model.occ.addPoint(1.0e-3x, 1.0e-3y, 0.0) for (x, y) in clean_points]
    line_tags = [
        gmsh.model.occ.addLine(point_tags[index], point_tags[index == length(point_tags) ? 1 : index + 1])
        for index in eachindex(point_tags)
    ]
    loop = gmsh.model.occ.addCurveLoop(line_tags)
    gmsh.model.occ.addPlaneSurface([loop])
end

function rounded_notch_boundaries(
    cell::RoundedNotchCell,
    center_y_mm::Real;
    arc_samples::Integer=14,
)
    bottom_y = Float64(center_y_mm) - cell.height_mm / 2
    top_y = Float64(center_y_mm) + cell.height_mm / 2
    depth = cell.height_mm - cell.minimum_gap_mm
    radius = cell.notch_width_mm / 2
    stem = depth - radius
    centers = notch_centers_mm(cell)
    bottom_centers = [
        center for (index, center) in enumerate(centers)
        if isodd(index) == cell.start_from_bottom
    ]
    top_centers = [
        center for (index, center) in enumerate(centers)
        if isodd(index) != cell.start_from_bottom
    ]

    bottom = Tuple{Float64, Float64}[(0.0, bottom_y)]
    if depth > 0
        for center in bottom_centers
            left, right = center - radius, center + radius
            push!(bottom, (left, bottom_y), (left, bottom_y + stem))
            for theta in range(pi, 0.0; length=arc_samples + 1)
                push!(bottom, (
                    center + radius * cos(theta),
                    bottom_y + stem + radius * sin(theta),
                ))
            end
            push!(bottom, (right, bottom_y))
        end
    end
    push!(bottom, (cell.length_mm, bottom_y))

    top = Tuple{Float64, Float64}[(cell.length_mm, top_y)]
    if depth > 0
        for center in reverse(top_centers)
            right, left = center + radius, center - radius
            push!(top, (right, top_y), (right, top_y - stem))
            for theta in range(0.0, -pi; length=arc_samples + 1)
                push!(top, (
                    center + radius * cos(theta),
                    top_y - stem + radius * sin(theta),
                ))
            end
            push!(top, (left, top_y))
        end
    end
    push!(top, (0.0, top_y))
    bottom, top
end

function cell_boundary_points(cell::RoundedNotchCell, center_y_mm::Real; samples::Integer=160)
    bottom, top = rounded_notch_boundaries(cell, center_y_mm)
    vcat(bottom, top)
end

function cell_boundary_points(cell::SinusoidalCell, center_y_mm::Real; samples::Integer=160)
    amplitude = (cell.height_mm - cell.minimum_gap_mm) / 2
    bottom_y = center_y_mm - cell.height_mm / 2
    active = cell.length_mm - 2cell.end_margin_mm
    bottom = Tuple{Float64, Float64}[(0.0, bottom_y)]
    for index in 0:samples
        xi = index / samples
        x = cell.end_margin_mm + active * xi
        depth = amplitude * 0.5 * (1 - cos(2pi * cell.periods * xi))
        push!(bottom, (x, bottom_y + depth))
    end
    push!(bottom, (cell.length_mm, bottom_y))
    top = reverse([
        (x, 2Float64(center_y_mm) - y)
        for (x, y) in bottom
    ])
    vcat(bottom, top)
end

function boundary_entities()
    [
        tag for (_, tag) in gmsh.model.getEntities(1)
        if length(first(gmsh.model.getAdjacencies(1, tag))) == 1
    ]
end

function add_physical_group(dim, tags, id, name)
    isempty(tags) && return
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function build_material_lens_mesh(
    cells::AbstractVector{<:AbstractLensCell},
    output_path::AbstractString;
    matching_thickness_mm::Real=0.0,
    config::MaterialLensMeshConfig=MaterialLensMeshConfig(),
)
    validate(cells, matching_thickness_mm, config)
    mkpath(dirname(output_path))
    lens_length_mm = first(cells).length_mm
    element_height_mm = first(cells).height_mm
    pitch_mm = element_height_mm + config.slot_width_mm
    half_count = (length(cells) - 1) ÷ 2
    centers_mm = collect((-half_count):half_count) .* pitch_mm
    aperture_min_mm = first(centers_mm) - element_height_mm / 2
    aperture_height_mm = last(centers_mm) - first(centers_mm) + element_height_mm
    matching_end_mm = lens_length_mm + Float64(matching_thickness_mm)
    domain_end_mm = matching_end_mm + config.output_length_mm

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("profiled_material_lens")
        surfaces = Int[
            polygon_surface(cell_boundary_points(cell, center))
            for (cell, center) in zip(cells, centers_mm)
        ]
        if matching_thickness_mm > 0
            push!(surfaces, gmsh.model.occ.addRectangle(
                1.0e-3lens_length_mm,
                1.0e-3aperture_min_mm,
                0.0,
                1.0e-3matching_thickness_mm,
                1.0e-3aperture_height_mm,
            ))
        end
        output_surface = gmsh.model.occ.addRectangle(
            1.0e-3matching_end_mm,
            1.0e-3aperture_min_mm,
            0.0,
            1.0e-3config.output_length_mm,
            1.0e-3aperture_height_mm,
        )
        if config.defect_radius_mm > 0
            radius = config.defect_radius_mm
            center_x = config.defect_center_x_mm
            center_y = config.defect_center_y_mm
            center_x - radius > matching_end_mm ||
                throw(ArgumentError("defect intersects the lens/output interface"))
            center_x + radius < domain_end_mm ||
                throw(ArgumentError("defect intersects the downstream boundary"))
            center_y - radius > aperture_min_mm ||
                throw(ArgumentError("defect intersects the lower boundary"))
            center_y + radius < aperture_min_mm + aperture_height_mm ||
                throw(ArgumentError("defect intersects the upper boundary"))
            defect_disk = gmsh.model.occ.addDisk(
                1.0e-3center_x,
                1.0e-3center_y,
                0.0,
                1.0e-3radius,
                1.0e-3radius,
            )
            cut, _ = gmsh.model.occ.cut([(2, output_surface)], [(2, defect_disk)])
            output_surfaces = [tag for (dim, tag) in cut if dim == 2]
            length(output_surfaces) == 1 || error("defect cut did not preserve one output surface")
            push!(surfaces, only(output_surfaces))
        else
            push!(surfaces, output_surface)
        end

        gmsh.model.occ.fragment([(2, surfaces[1])], [(2, tag) for tag in surfaces[2:end]])
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()

        tolerance_m = config.boundary_tolerance_mm * 1.0e-3
        lens_surfaces, matching_surfaces, aluminium_surfaces = Int[], Int[], Int[]
        for (_, tag) in gmsh.model.getEntities(2)
            x = gmsh.model.occ.getCenterOfMass(2, tag)[1]
            if x < 1.0e-3lens_length_mm - tolerance_m
                push!(lens_surfaces, tag)
            elseif matching_thickness_mm > 0 && x < 1.0e-3matching_end_mm - tolerance_m
                push!(matching_surfaces, tag)
            else
                push!(aluminium_surfaces, tag)
            end
        end
        length(lens_surfaces) == length(cells) ||
            error("surface classification lost a lens strip")
        !isempty(aluminium_surfaces) || error("aluminium output domain was not found")

        source, radiation, free, defect_curves = Int[], Int[], Int[], Int[]
        y_min_m = aperture_min_mm * 1.0e-3
        y_max_m = (aperture_min_mm + aperture_height_mm) * 1.0e-3
        domain_end_m = domain_end_mm * 1.0e-3
        matching_end_m = matching_end_mm * 1.0e-3
        for tag in boundary_entities()
            center = gmsh.model.occ.getCenterOfMass(1, tag)
            x, y = center[1], center[2]
            if x < tolerance_m
                push!(source, tag)
            elseif abs(x - domain_end_m) < tolerance_m ||
                   (x > matching_end_m + tolerance_m &&
                    (abs(y - y_min_m) < tolerance_m || abs(y - y_max_m) < tolerance_m))
                push!(radiation, tag)
            else
                push!(free, tag)
            end
            if config.defect_radius_mm > 0
                box = gmsh.model.occ.getBoundingBox(1, tag)
                defect_box_tolerance_m = max(tolerance_m, 1.0e-6)
                target_xmin = (config.defect_center_x_mm - config.defect_radius_mm) * 1.0e-3
                target_xmax = (config.defect_center_x_mm + config.defect_radius_mm) * 1.0e-3
                target_ymin = (config.defect_center_y_mm - config.defect_radius_mm) * 1.0e-3
                target_ymax = (config.defect_center_y_mm + config.defect_radius_mm) * 1.0e-3
                if abs(box[1] - target_xmin) < defect_box_tolerance_m &&
                   abs(box[2] - target_ymin) < defect_box_tolerance_m &&
                   abs(box[4] - target_xmax) < defect_box_tolerance_m &&
                   abs(box[5] - target_ymax) < defect_box_tolerance_m
                    push!(defect_curves, tag)
                end
            end
        end
        !isempty(source) || error("source boundary was not found")
        !isempty(radiation) || error("radiation boundary was not found")

        source_by_strip = [Int[] for _ in cells]
        centers_m = centers_mm .* 1.0e-3
        for tag in source
            y = gmsh.model.occ.getCenterOfMass(1, tag)[2]
            strip_index = argmin(abs.(centers_m .- y))
            push!(source_by_strip[strip_index], tag)
        end
        all(tags -> !isempty(tags), source_by_strip) ||
            error("at least one strip source boundary was lost")

        add_physical_group(1, source, 101, "Source")
        add_physical_group(1, radiation, 102, "RadiationBoundary")
        add_physical_group(1, free, 103, "FreeSurface")
        add_physical_group(1, defect_curves, 104, "DefectBoundary")
        add_physical_group(2, lens_surfaces, 201, "Lens")
        add_physical_group(2, matching_surfaces, 202, "MatchingLayer")
        add_physical_group(2, aluminium_surfaces, 203, "Aluminium")
        for (index, tags) in enumerate(source_by_strip)
            add_physical_group(1, tags, 300 + index, "SourceStrip$(index)")
        end

        lens_curves = Int[]
        for tag in lens_surfaces
            append!(lens_curves, last(gmsh.model.getAdjacencies(2, tag)))
        end
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", unique(lens_curves))
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", config.lens_mesh_size_mm * 1.0e-3)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", config.output_mesh_size_mm * 1.0e-3)
        gmsh.model.mesh.field.setNumber(2, "DistMin", 1.0e-3)
        gmsh.model.mesh.field.setNumber(2, "DistMax", config.transition_distance_mm * 1.0e-3)
        background_field = 2
        if config.defect_radius_mm > 0
            !isempty(defect_curves) || error("defect boundary was not found")
            gmsh.model.mesh.field.add("Distance", 3)
            gmsh.model.mesh.field.setNumbers(3, "CurvesList", defect_curves)
            gmsh.model.mesh.field.add("Threshold", 4)
            gmsh.model.mesh.field.setNumber(4, "InField", 3)
            gmsh.model.mesh.field.setNumber(4, "SizeMin", config.defect_mesh_size_mm * 1.0e-3)
            gmsh.model.mesh.field.setNumber(4, "SizeMax", config.output_mesh_size_mm * 1.0e-3)
            gmsh.model.mesh.field.setNumber(4, "DistMin", 0.0)
            gmsh.model.mesh.field.setNumber(
                4,
                "DistMax",
                config.defect_transition_distance_mm * 1.0e-3,
            )
            gmsh.model.mesh.field.add("Min", 5)
            gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])
            background_field = 5
        end
        gmsh.model.mesh.field.setAsBackgroundMesh(background_field)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.model.mesh.generate(2)
        gmsh.write(output_path)
    finally
        gmsh.finalize()
    end
    (
        mesh_path=output_path,
        lens_length_mm,
        matching_thickness_mm=Float64(matching_thickness_mm),
        domain_end_mm,
        aperture_min_mm,
        aperture_height_mm,
        focus_x_mm=matching_end_mm + 60.0,
        defect_radius_mm=config.defect_radius_mm,
        defect_center_x_mm=config.defect_center_x_mm,
        defect_center_y_mm=config.defect_center_y_mm,
    )
end

end
