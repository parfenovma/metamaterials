module BalancedHuygensMesher

using Gmsh: gmsh

export BalancedHuygensConfig, validate_config, solid_rectangles, build_balanced_huygens_mesh

Base.@kwdef struct BalancedHuygensConfig
    length_mm::Float64 = 24.0
    height_mm::Float64 = 4.2
    lead_length_mm::Float64 = 3.0
    spine_height_mm::Float64 = 1.2
    series_center_x_mm::Float64 = 9.0
    series_length_mm::Float64 = 1.2
    series_height_mm::Float64 = 0.75
    shunt_center_x_mm::Float64 = 15.5
    shunt_mass_length_mm::Float64 = 3.0
    shunt_mass_height_mm::Float64 = 0.8
    wall_margin_mm::Float64 = 0.3
    shunt_neck_width_mm::Float64 = 0.6
end

function validate_config(config::BalancedHuygensConfig)
    config.length_mm > 0 || throw(ArgumentError("length must be positive"))
    config.height_mm > 0 || throw(ArgumentError("height must be positive"))
    0 < config.lead_length_mm < config.length_mm / 2 ||
        throw(ArgumentError("lead length must leave an internal cell"))
    0 < config.spine_height_mm < config.height_mm ||
        throw(ArgumentError("spine height must lie inside the strip"))
    0 < config.series_height_mm <= config.spine_height_mm ||
        throw(ArgumentError("series height must be positive and no larger than the spine"))
    config.series_length_mm > 0 || throw(ArgumentError("series length must be positive"))
    config.shunt_mass_length_mm > 0 || throw(ArgumentError("shunt mass must be positive"))
    config.shunt_mass_height_mm > 0 || throw(ArgumentError("shunt mass height must be positive"))
    config.wall_margin_mm >= 0 || throw(ArgumentError("wall margin must be non-negative"))
    0 < config.shunt_neck_width_mm <= config.shunt_mass_length_mm ||
        throw(ArgumentError("shunt neck must fit inside its mass"))

    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    spine_top = spine_bottom + config.spine_height_mm
    top_mass_bottom = config.height_mm - config.wall_margin_mm - config.shunt_mass_height_mm
    bottom_mass_top = config.wall_margin_mm + config.shunt_mass_height_mm
    top_mass_bottom > spine_top || throw(ArgumentError("top shunt neck has no free length"))
    bottom_mass_top < spine_bottom || throw(ArgumentError("bottom shunt neck has no free length"))

    series_left = config.series_center_x_mm - config.series_length_mm / 2
    series_right = config.series_center_x_mm + config.series_length_mm / 2
    shunt_left = config.shunt_center_x_mm - config.shunt_mass_length_mm / 2
    shunt_right = config.shunt_center_x_mm + config.shunt_mass_length_mm / 2
    config.lead_length_mm < series_left || throw(ArgumentError("series section overlaps left lead"))
    series_right < shunt_left || throw(ArgumentError("series section overlaps shunt mass"))
    shunt_right < config.length_mm - config.lead_length_mm ||
        throw(ArgumentError("shunt mass overlaps right lead"))
    nothing
end

function uniform_spine(config)
    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    [(
        config.lead_length_mm,
        spine_bottom,
        config.length_mm - 2config.lead_length_mm,
        config.spine_height_mm,
    )]
end

function series_spine(config)
    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    series_bottom = (config.height_mm - config.series_height_mm) / 2
    series_left = config.series_center_x_mm - config.series_length_mm / 2
    series_right = config.series_center_x_mm + config.series_length_mm / 2
    [
        (
            config.lead_length_mm,
            spine_bottom,
            series_left - config.lead_length_mm,
            config.spine_height_mm,
        ),
        (series_left, series_bottom, config.series_length_mm, config.series_height_mm),
        (
            series_right,
            spine_bottom,
            config.length_mm - config.lead_length_mm - series_right,
            config.spine_height_mm,
        ),
    ]
end

function shunt_rectangles(config)
    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    spine_top = spine_bottom + config.spine_height_mm
    top_mass_bottom = config.height_mm - config.wall_margin_mm - config.shunt_mass_height_mm
    bottom_mass_top = config.wall_margin_mm + config.shunt_mass_height_mm
    mass_left = config.shunt_center_x_mm - config.shunt_mass_length_mm / 2
    neck_left = config.shunt_center_x_mm - config.shunt_neck_width_mm / 2
    [
        (mass_left, top_mass_bottom, config.shunt_mass_length_mm, config.shunt_mass_height_mm),
        (mass_left, config.wall_margin_mm, config.shunt_mass_length_mm, config.shunt_mass_height_mm),
        (neck_left, spine_top, config.shunt_neck_width_mm, top_mass_bottom - spine_top),
        (neck_left, bottom_mass_top, config.shunt_neck_width_mm, spine_bottom - bottom_mass_top),
    ]
end

function solid_rectangles(config::BalancedHuygensConfig; variant::Symbol=:balanced)
    validate_config(config)
    variant in (:r0, :series, :shunt, :balanced) ||
        throw(ArgumentError("variant must be :r0, :series, :shunt, or :balanced"))
    rectangles = [
        (0.0, 0.0, config.lead_length_mm, config.height_mm),
        (
            config.length_mm - config.lead_length_mm,
            0.0,
            config.lead_length_mm,
            config.height_mm,
        ),
    ]
    append!(rectangles, variant in (:series, :balanced) ? series_spine(config) : uniform_spine(config))
    variant in (:shunt, :balanced) && append!(rectangles, shunt_rectangles(config))
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

function build_balanced_huygens_mesh(
    output_path::AbstractString;
    config::BalancedHuygensConfig=BalancedHuygensConfig(),
    variant::Symbol=:balanced,
    size_min_mm::Real=0.06,
    size_max_mm::Real=0.28,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("balanced_huygens_$(variant)")
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
