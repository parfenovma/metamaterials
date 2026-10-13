module OperatorO8BoundaryMesher

using Gmsh: gmsh

export O8BoundaryGeometry, validate, control_intervals, build_o8_boundary_mesh

Base.@kwdef struct O8BoundaryGeometry
    aperture_height_m::Float64
    block_length_m::Float64
    position_a::Float64
    width_a::Float64
    position_b::Float64
    width_b::Float64
    mesh_size_m::Float64
    output_bins::Int = 17
end

function validate(config::O8BoundaryGeometry)
    config.aperture_height_m > 0 || throw(ArgumentError("aperture must be positive"))
    config.block_length_m > 0 || throw(ArgumentError("length must be positive"))
    0 < config.position_a < 1 || throw(ArgumentError("position A must be normalized"))
    0 < config.position_b < 1 || throw(ArgumentError("position B must be normalized"))
    config.position_a < config.position_b || throw(ArgumentError("A must be inside B"))
    config.width_a > 0 || throw(ArgumentError("width A must be positive"))
    config.width_b > 0 || throw(ArgumentError("width B must be positive"))
    config.mesh_size_m > 0 || throw(ArgumentError("mesh size must be positive"))
    config.output_bins >= 3 || throw(ArgumentError("at least three output bins"))
    config
end

function control_intervals(config::O8BoundaryGeometry)
    validate(config)
    half = config.aperture_height_m / 2
    interval(position, width, sign) = begin
        center = sign * position * half
        radius = width * half
        (max(-half, center - radius), min(half, center + radius))
    end
    (
        a_minus=interval(config.position_a, config.width_a, -1),
        a_plus=interval(config.position_a, config.width_a, 1),
        b_minus=interval(config.position_b, config.width_b, -1),
        b_plus=interval(config.position_b, config.width_b, 1),
    )
end

function add_group(dimension, tags, id, name)
    isempty(tags) && error("physical group $name is empty")
    gmsh.model.addPhysicalGroup(dimension, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dimension, id, name)
end

function interval_lines(lines, coordinates, interval; tolerance=1.0e-10)
    lower, upper = interval
    [
        lines[index] for index in eachindex(lines)
        if coordinates[index] >= lower - tolerance &&
           coordinates[index + 1] <= upper + tolerance
    ]
end

function build_o8_boundary_mesh(
    output_path::AbstractString;
    config::O8BoundaryGeometry,
)
    validate(config)
    mkpath(dirname(output_path))
    intervals = control_intervals(config)
    half = config.aperture_height_m / 2
    endpoints = Float64[-half, half]
    for interval in values(intervals)
        append!(endpoints, interval)
    end
    left_y = sort(unique(endpoints))
    right_y = collect(range(-half, half; length=config.output_bins + 1))

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("operator_o8_stationary_boundary")
        left_points = [gmsh.model.occ.addPoint(0.0, y, 0.0) for y in left_y]
        right_points = [
            gmsh.model.occ.addPoint(config.block_length_m, y, 0.0) for y in right_y
        ]
        left_lines = [
            gmsh.model.occ.addLine(left_points[index], left_points[index + 1])
            for index in 1:(length(left_points) - 1)
        ]
        right_lines = [
            gmsh.model.occ.addLine(right_points[index], right_points[index + 1])
            for index in 1:config.output_bins
        ]
        top = gmsh.model.occ.addLine(last(left_points), last(right_points))
        bottom = gmsh.model.occ.addLine(first(right_points), first(left_points))
        loop = gmsh.model.occ.addCurveLoop(vcat(
            left_lines, [top], -reverse(right_lines), [bottom],
        ))
        surface = gmsh.model.occ.addPlaneSurface([loop])
        gmsh.model.occ.synchronize()

        add_group(2, [surface], 201, "Block")
        add_group(1, left_lines, 101, "Source")
        add_group(1, right_lines, 102, "RadiationBoundary")
        add_group(1, [top, bottom], 103, "FreeBoundary")
        add_group(1, interval_lines(left_lines, left_y, intervals.a_minus), 111, "ZoneAMinus")
        add_group(1, interval_lines(left_lines, left_y, intervals.a_plus), 112, "ZoneAPlus")
        add_group(1, interval_lines(left_lines, left_y, intervals.b_minus), 113, "ZoneBMinus")
        add_group(1, interval_lines(left_lines, left_y, intervals.b_plus), 114, "ZoneBPlus")
        for index in 1:config.output_bins
            add_group(1, [right_lines[index]], 300 + index, "OutputBin$(index)")
        end

        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.75config.mesh_size_m)
        gmsh.option.setNumber("Mesh.MeshSizeMax", config.mesh_size_m)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 1)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        gmsh.write(output_path)
        node_count = length(first(gmsh.model.mesh.getNodes()))
        _, tags, _ = gmsh.model.mesh.getElements(2)
        element_count = sum(length, tags)
        interval_rows = [
            (zone=name, lower_m=value[1], upper_m=value[2], length_m=value[2] - value[1])
            for (name, value) in pairs(intervals)
        ]
        (
            output_path=abspath(output_path), node_count, element_count,
            aperture_height_m=config.aperture_height_m,
            block_length_m=config.block_length_m,
            mesh_size_m=config.mesh_size_m, interval_rows,
        )
    finally
        gmsh.finalize()
    end
end

end

