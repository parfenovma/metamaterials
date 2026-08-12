module PerforatedChannelMesher

using Gmsh: gmsh

export PerforatedChannelConfig,
       PERFORATED_VARIANTS,
       validate_config,
       active_length_mm,
       total_length_mm,
       hole_axes_mm,
       hole_centers_mm,
       hole_ellipses_mm,
       void_area_mm2,
       minimum_ligaments_mm,
       build_perforated_channel_mesh

const PERFORATED_VARIANTS = (:solid, :round, :ellipse_longitudinal, :ellipse_transverse)

"""Two aligned, symmetric rows of equal-area through-holes in a straight strip."""
Base.@kwdef struct PerforatedChannelConfig
    lead_length_mm::Float64 = 12.0
    active_length_mm::Float64 = 12.0
    port_height_mm::Float64 = 3.2
    holes_per_row::Int = 4
    pitch_mm::Float64 = 2.4
    row_center_mm::Float64 = 0.72
    round_radius_mm::Float64 = 0.40
    ellipse_aspect_ratio::Float64 = 1.5
    minimum_feature_mm::Float64 = 0.35
end

active_length_mm(config::PerforatedChannelConfig) = config.active_length_mm
total_length_mm(config::PerforatedChannelConfig) =
    2config.lead_length_mm + config.active_length_mm

function hole_axes_mm(config::PerforatedChannelConfig, variant::Symbol)
    variant in PERFORATED_VARIANTS || throw(ArgumentError("unknown perforated variant"))
    variant == :solid && return (0.0, 0.0)
    radius = config.round_radius_mm
    variant == :round && return (radius, radius)
    ratio = config.ellipse_aspect_ratio
    minor = radius / sqrt(ratio)
    major = radius * sqrt(ratio)
    variant == :ellipse_longitudinal ? (major, minor) : (minor, major)
end

function hole_centers_mm(config::PerforatedChannelConfig)
    occupied = (config.holes_per_row - 1) * config.pitch_mm
    first_x = config.lead_length_mm + (config.active_length_mm - occupied) / 2
    [
        (first_x + (index - 1) * config.pitch_mm, sign * config.row_center_mm)
        for index in 1:config.holes_per_row
        for sign in (-1.0, 1.0)
    ]
end

function hole_ellipses_mm(config::PerforatedChannelConfig, variant::Symbol)
    variant == :solid && return NTuple{4, Float64}[]
    radius_x, radius_y = hole_axes_mm(config, variant)
    [(x, y, radius_x, radius_y) for (x, y) in hole_centers_mm(config)]
end

function void_area_mm2(config::PerforatedChannelConfig, variant::Symbol)
    variant == :solid && return 0.0
    radius_x, radius_y = hole_axes_mm(config, variant)
    2config.holes_per_row * pi * radius_x * radius_y
end

function minimum_ligaments_mm(config::PerforatedChannelConfig, variant::Symbol)
    variant == :solid && return (
        between_rows=Inf,
        outer_wall=Inf,
        between_columns=Inf,
        active_end=Inf,
    )
    radius_x, radius_y = hole_axes_mm(config, variant)
    occupied = (config.holes_per_row - 1) * config.pitch_mm
    end_margin = (config.active_length_mm - occupied) / 2 - radius_x
    (
        between_rows=2(config.row_center_mm - radius_y),
        outer_wall=config.port_height_mm / 2 - config.row_center_mm - radius_y,
        between_columns=config.pitch_mm - 2radius_x,
        active_end=end_margin,
    )
end

function validate_config(config::PerforatedChannelConfig)
    config.lead_length_mm > 0 || throw(ArgumentError("lead length must be positive"))
    config.active_length_mm > 0 || throw(ArgumentError("active length must be positive"))
    config.port_height_mm > 0 || throw(ArgumentError("port height must be positive"))
    config.holes_per_row >= 2 || throw(ArgumentError("at least two holes per row are required"))
    config.pitch_mm > 0 || throw(ArgumentError("hole pitch must be positive"))
    config.row_center_mm > 0 || throw(ArgumentError("row center must be positive"))
    config.round_radius_mm > 0 || throw(ArgumentError("hole radius must be positive"))
    config.ellipse_aspect_ratio >= 1 || throw(ArgumentError("ellipse aspect ratio must be at least one"))
    for variant in Iterators.drop(PERFORATED_VARIANTS, 1)
        ligaments = minimum_ligaments_mm(config, variant)
        for (name, value) in pairs(ligaments)
            value + 1.0e-12 >= config.minimum_feature_mm || throw(ArgumentError(
                "$(variant) $(name) ligament $(value) mm violates the minimum feature",
            ))
        end
    end
    nothing
end

function build_domain(config, variant)
    base = gmsh.model.occ.addRectangle(
        0.0,
        -config.port_height_mm * 0.5e-3,
        0.0,
        total_length_mm(config) * 1e-3,
        config.port_height_mm * 1e-3,
    )
    variant == :solid && return [base]
    holes = Int[]
    for (x, y, radius_x, radius_y) in hole_ellipses_mm(config, variant)
        if radius_x >= radius_y
            push!(holes, gmsh.model.occ.addDisk(
                x * 1e-3,
                y * 1e-3,
                0.0,
                radius_x * 1e-3,
                radius_y * 1e-3,
            ))
        else
            tag = gmsh.model.occ.addDisk(
                x * 1e-3,
                y * 1e-3,
                0.0,
                radius_y * 1e-3,
                radius_x * 1e-3,
            )
            gmsh.model.occ.rotate(
                [(2, tag)],
                x * 1e-3,
                y * 1e-3,
                0.0,
                0.0,
                0.0,
                1.0,
                pi / 2,
            )
            push!(holes, tag)
        end
    end
    cut, _ = gmsh.model.occ.cut([(2, base)], [(2, tag) for tag in holes])
    surfaces = unique(tag for (dimension, tag) in cut if dimension == 2)
    isempty(surfaces) && error("perforated cut produced no solid domain")
    surfaces
end

function classify_boundaries(surfaces, config)
    gmsh.model.occ.synchronize()
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    curves = unique(tag for (dimension, tag) in boundary if dimension == 1)
    total_m = total_length_mm(config) * 1e-3
    tolerance = max(1.0e-10, total_m * 1.0e-8)
    source = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] < tolerance, curves)
    microphone = filter(
        tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] > total_m - tolerance,
        curves,
    )
    length(source) == 1 || error("perforated channel must have exactly one left port")
    length(microphone) == 1 || error("perforated channel must have exactly one right port")
    free = setdiff(curves, vcat(source, microphone))
    source, microphone, free
end

function tag_physical_groups(surfaces, source, microphone, free)
    for (id, tags, name) in (
        (101, source, "Source"),
        (102, microphone, "Microphone"),
        (103, free, "FreeSurface"),
    )
        gmsh.model.addPhysicalGroup(1, tags, id)
        gmsh.model.setPhysicalName(1, id, name)
    end
    gmsh.model.addPhysicalGroup(2, surfaces, 201)
    gmsh.model.setPhysicalName(2, 201, "Domain")
end

function configure_mesh(source, microphone, free, config; size_hole_mm, size_port_mm, size_lead_mm)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", vcat(source, microphone, free))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_hole_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_lead_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.04e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.45e-3)

    lead_m = config.lead_length_mm * 1e-3
    active_m = config.active_length_mm * 1e-3
    height_m = config.port_height_mm * 1e-3
    gmsh.model.mesh.field.add("Box", 3)
    gmsh.model.mesh.field.setNumber(3, "VIn", Float64(size_hole_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "VOut", Float64(size_lead_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "XMin", lead_m - 0.5e-3)
    gmsh.model.mesh.field.setNumber(3, "XMax", lead_m + active_m + 0.5e-3)
    gmsh.model.mesh.field.setNumber(3, "YMin", -height_m)
    gmsh.model.mesh.field.setNumber(3, "YMax", height_m)

    gmsh.model.mesh.field.add("Min", 4)
    gmsh.model.mesh.field.setNumbers(4, "FieldsList", [2, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(4)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
end

function build_perforated_channel_mesh(
    output_path::AbstractString;
    config::PerforatedChannelConfig=PerforatedChannelConfig(),
    variant::Symbol=:round,
    size_hole_mm::Real=0.06,
    size_port_mm::Real=0.08,
    size_lead_mm::Real=0.25,
)
    validate_config(config)
    variant in PERFORATED_VARIANTS || throw(ArgumentError("unknown perforated variant"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("perforated_channel_$(variant)")
    surfaces = build_domain(config, variant)
    source, microphone, free = classify_boundaries(surfaces, config)
    tag_physical_groups(surfaces, source, microphone, free)
    configure_mesh(source, microphone, free, config; size_hole_mm, size_port_mm, size_lead_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
