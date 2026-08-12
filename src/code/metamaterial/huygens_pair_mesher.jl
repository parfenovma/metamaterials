module HuygensPairMesher

using Gmsh: gmsh

export HuygensPairConfig, validate_config, solid_rectangles, build_huygens_pair_mesh

Base.@kwdef struct HuygensPairConfig
    length_mm::Float64 = 24.0
    height_mm::Float64 = 4.2
    lead_length_mm::Float64 = 3.0
    spine_height_mm::Float64 = 1.2
    front_center_x_mm::Float64 = 8.0
    rear_center_x_mm::Float64 = 16.0
    front_mass_length_mm::Float64 = 3.0
    rear_mass_length_mm::Float64 = 3.0
    mass_height_mm::Float64 = 0.8
    wall_margin_mm::Float64 = 0.3
    front_neck_width_mm::Float64 = 0.6
    rear_neck_width_mm::Float64 = 0.6
end

function validate_config(config::HuygensPairConfig)
    config.length_mm > 0 || throw(ArgumentError("length must be positive"))
    config.height_mm > 0 || throw(ArgumentError("height must be positive"))
    0 < config.lead_length_mm < config.length_mm / 2 ||
        throw(ArgumentError("lead length must leave an internal cell"))
    0 < config.spine_height_mm < config.height_mm ||
        throw(ArgumentError("spine height must lie inside the strip"))
    config.front_mass_length_mm > 0 || throw(ArgumentError("front mass must be positive"))
    config.rear_mass_length_mm > 0 || throw(ArgumentError("rear mass must be positive"))
    config.mass_height_mm > 0 || throw(ArgumentError("mass height must be positive"))
    config.wall_margin_mm >= 0 || throw(ArgumentError("wall margin must be non-negative"))
    config.front_neck_width_mm > 0 || throw(ArgumentError("front neck must be positive"))
    config.rear_neck_width_mm > 0 || throw(ArgumentError("rear neck must be positive"))
    config.front_neck_width_mm <= config.front_mass_length_mm ||
        throw(ArgumentError("front neck cannot exceed its mass"))
    config.rear_neck_width_mm <= config.rear_mass_length_mm ||
        throw(ArgumentError("rear neck cannot exceed its mass"))

    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    top_spine = spine_bottom + config.spine_height_mm
    top_mass_bottom = config.height_mm - config.wall_margin_mm - config.mass_height_mm
    bottom_mass_top = config.wall_margin_mm + config.mass_height_mm
    top_mass_bottom > top_spine || throw(ArgumentError("top neck has no free length"))
    bottom_mass_top < spine_bottom || throw(ArgumentError("bottom neck has no free length"))

    left_limit = config.lead_length_mm
    right_limit = config.length_mm - config.lead_length_mm
    front_left = config.front_center_x_mm - config.front_mass_length_mm / 2
    front_right = config.front_center_x_mm + config.front_mass_length_mm / 2
    rear_left = config.rear_center_x_mm - config.rear_mass_length_mm / 2
    rear_right = config.rear_center_x_mm + config.rear_mass_length_mm / 2
    left_limit < front_left || throw(ArgumentError("front mass overlaps the left lead"))
    rear_right < right_limit || throw(ArgumentError("rear mass overlaps the right lead"))
    front_right < rear_left || throw(ArgumentError("front and rear masses overlap"))
    nothing
end

function resonator_rectangles(config, center_x_mm, mass_length_mm, neck_width_mm)
    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    top_spine = spine_bottom + config.spine_height_mm
    top_mass_bottom = config.height_mm - config.wall_margin_mm - config.mass_height_mm
    bottom_mass_top = config.wall_margin_mm + config.mass_height_mm
    mass_left = center_x_mm - mass_length_mm / 2
    neck_left = center_x_mm - neck_width_mm / 2
    [
        (mass_left, top_mass_bottom, mass_length_mm, config.mass_height_mm),
        (mass_left, config.wall_margin_mm, mass_length_mm, config.mass_height_mm),
        (neck_left, top_spine, neck_width_mm, top_mass_bottom - top_spine),
        (neck_left, bottom_mass_top, neck_width_mm, spine_bottom - bottom_mass_top),
    ]
end

function solid_rectangles(config::HuygensPairConfig; variant::Symbol=:pair)
    validate_config(config)
    variant in (:r0, :front, :rear, :pair) ||
        throw(ArgumentError("variant must be :r0, :front, :rear, or :pair"))
    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    rectangles = [
        (0.0, 0.0, config.lead_length_mm, config.height_mm),
        (
            config.lead_length_mm,
            spine_bottom,
            config.length_mm - 2config.lead_length_mm,
            config.spine_height_mm,
        ),
        (
            config.length_mm - config.lead_length_mm,
            0.0,
            config.lead_length_mm,
            config.height_mm,
        ),
    ]
    variant in (:front, :pair) && append!(rectangles, resonator_rectangles(
        config,
        config.front_center_x_mm,
        config.front_mass_length_mm,
        config.front_neck_width_mm,
    ))
    variant in (:rear, :pair) && append!(rectangles, resonator_rectangles(
        config,
        config.rear_center_x_mm,
        config.rear_mass_length_mm,
        config.rear_neck_width_mm,
    ))
    rectangles
end

function add_rectangle_mm(rectangle)
    x_mm, y_mm, width_mm, height_mm = rectangle
    gmsh.model.occ.addRectangle(
        x_mm * 1e-3,
        y_mm * 1e-3,
        0.0,
        width_mm * 1e-3,
        height_mm * 1e-3,
    )
end

function classify_boundaries(surface_tags, config)
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surface_tags], false, false, false)
    line_tags = unique(tag for (dimension, tag) in boundary if dimension == 1)
    source, microphone, free_surface = Int[], Int[], Int[]
    length_m = config.length_mm * 1e-3
    for tag in line_tags
        x = gmsh.model.occ.getCenterOfMass(1, tag)[1]
        if x < 1e-9
            push!(source, tag)
        elseif x > length_m - 1e-9
            push!(microphone, tag)
        else
            push!(free_surface, tag)
        end
    end
    source, microphone, free_surface
end

function build_huygens_pair_mesh(
    output_path::AbstractString;
    config::HuygensPairConfig=HuygensPairConfig(),
    variant::Symbol=:pair,
    size_min_mm::Real=0.06,
    size_max_mm::Real=0.28,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("huygens_pair_$(variant)")
    rectangle_tags = add_rectangle_mm.(solid_rectangles(config; variant))
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(rectangle_tags))],
        [(2, tag) for tag in Iterators.drop(rectangle_tags, 1)],
    )
    gmsh.model.occ.synchronize()
    surface_tags = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surface_tags) && error("OCC fuse produced no solid surface")
    source, microphone, free_surface = classify_boundaries(surface_tags, config)
    isempty(source) && error("source boundary was not detected")
    isempty(microphone) && error("microphone boundary was not detected")

    gmsh.model.addPhysicalGroup(1, source, 101)
    gmsh.model.setPhysicalName(1, 101, "Source")
    gmsh.model.addPhysicalGroup(1, microphone, 102)
    gmsh.model.setPhysicalName(1, 102, "Microphone")
    gmsh.model.addPhysicalGroup(1, free_surface, 103)
    gmsh.model.setPhysicalName(1, 103, "FreeSurface")
    gmsh.model.addPhysicalGroup(2, surface_tags, 201)
    gmsh.model.setPhysicalName(2, 201, "Domain")

    all_boundaries = vcat(source, microphone, free_surface)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", all_boundaries)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_min_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.10e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.0e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
