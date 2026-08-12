module StiffnessComponentMesher

using Gmsh: gmsh

export StiffnessConfig,
       validate_config,
       slot_height_mm,
       ligament_area_fraction,
       slot_boxes,
       build_stiffness_component_mesh

const SLOT_COUNT = 3

"""Three rounded longitudinal slots leaving four parallel axial ligaments."""
Base.@kwdef struct StiffnessConfig
    length_mm::Float64 = 1.6
    height_mm::Float64 = 4.2
    end_wall_mm::Float64 = 0.35
    ligament_height_mm::Float64 = 0.40
    minimum_feature_mm::Float64 = 0.35
    minimum_radius_mm::Float64 = 0.20
end

slot_height_mm(config::StiffnessConfig) =
    (config.height_mm - (SLOT_COUNT + 1) * config.ligament_height_mm) / SLOT_COUNT

ligament_area_fraction(config::StiffnessConfig) =
    (SLOT_COUNT + 1) * config.ligament_height_mm / config.height_mm

function validate_config(config::StiffnessConfig)
    config.length_mm > 0 || throw(ArgumentError("component length must be positive"))
    config.height_mm > 0 || throw(ArgumentError("component height must be positive"))
    config.minimum_feature_mm > 0 ||
        throw(ArgumentError("minimum feature must be positive"))
    config.minimum_radius_mm > 0 ||
        throw(ArgumentError("minimum radius must be positive"))
    config.end_wall_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("end wall violates the minimum printable feature"))
    config.ligament_height_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("ligament violates the minimum printable feature"))
    slot_length = config.length_mm - 2config.end_wall_mm
    slot_length > 0 || throw(ArgumentError("end walls leave no slot length"))
    height = slot_height_mm(config)
    height > 0 || throw(ArgumentError("ligaments leave no slot height"))
    height / 2 >= config.minimum_radius_mm ||
        throw(ArgumentError("rounded slot end violates the minimum radius"))
    slot_length > height ||
        throw(ArgumentError("slot must be longer than its rounded-end diameter"))
    nothing
end

"""Bounding boxes `(x, y, length, height)` of the three slots in millimetres."""
function slot_boxes(config::StiffnessConfig)
    validate_config(config)
    height = slot_height_mm(config)
    length = config.length_mm - 2config.end_wall_mm
    y = -config.height_mm / 2 + config.ligament_height_mm
    boxes = NTuple{4, Float64}[]
    for _ in 1:SLOT_COUNT
        push!(boxes, (config.end_wall_mm, y, length, height))
        y += height + config.ligament_height_mm
    end
    boxes
end

function add_stadium_mm(x_mm, y_mm, length_mm, height_mm)
    radius_mm = height_mm / 2
    center_y_mm = y_mm + radius_mm
    rectangle = gmsh.model.occ.addRectangle(
        (x_mm + radius_mm) * 1e-3,
        y_mm * 1e-3,
        0.0,
        (length_mm - 2radius_mm) * 1e-3,
        height_mm * 1e-3,
    )
    left_disk = gmsh.model.occ.addDisk(
        (x_mm + radius_mm) * 1e-3,
        center_y_mm * 1e-3,
        0.0,
        radius_mm * 1e-3,
        radius_mm * 1e-3,
    )
    right_disk = gmsh.model.occ.addDisk(
        (x_mm + length_mm - radius_mm) * 1e-3,
        center_y_mm * 1e-3,
        0.0,
        radius_mm * 1e-3,
        radius_mm * 1e-3,
    )
    fused, _ = gmsh.model.occ.fuse(
        [(2, rectangle)],
        [(2, left_disk), (2, right_disk)],
    )
    surfaces = unique(tag for (dimension, tag) in fused if dimension == 2)
    length(surfaces) == 1 || error("rounded slot did not fuse into one surface")
    only(surfaces)
end

function build_domain(config, variant)
    base = gmsh.model.occ.addRectangle(
        0.0,
        -config.height_mm * 0.5e-3,
        0.0,
        config.length_mm * 1e-3,
        config.height_mm * 1e-3,
    )
    variant == :solid && return [base]
    slot_tags = [add_stadium_mm(box...) for box in slot_boxes(config)]
    cut, _ = gmsh.model.occ.cut(
        [(2, base)],
        [(2, tag) for tag in slot_tags],
    )
    surfaces = unique(tag for (dimension, tag) in cut if dimension == 2)
    isempty(surfaces) && error("slot cut produced no solid domain")
    surfaces
end

function tag_boundaries(surfaces, length_mm)
    gmsh.model.occ.synchronize()
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    curves = unique(tag for (dimension, tag) in boundary if dimension == 1)
    length_m = length_mm * 1e-3
    tolerance = max(1.0e-10, length_m * 1.0e-8)
    left = filter(curves) do tag
        gmsh.model.occ.getCenterOfMass(1, tag)[1] < tolerance
    end
    right = filter(curves) do tag
        gmsh.model.occ.getCenterOfMass(1, tag)[1] > length_m - tolerance
    end
    length(left) == 1 || error("expected one left interface, found $(length(left))")
    length(right) == 1 || error("expected one right interface, found $(length(right))")
    fixed = vcat(left, right)
    free = setdiff(curves, fixed)
    for (id, tags, name) in (
        (101, left, "LeftInterface"),
        (102, right, "RightInterface"),
        (104, fixed, "FixedInterface"),
        (103, free, "FreeSurface"),
    )
        gmsh.model.addPhysicalGroup(1, tags, id)
        gmsh.model.setPhysicalName(1, id, name)
    end
    gmsh.model.addPhysicalGroup(2, surfaces, 201)
    gmsh.model.setPhysicalName(2, 201, "Domain")
    curves
end

function configure_mesh(curves; size_min_mm, size_max_mm)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", curves)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_min_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.08e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.7e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
end

function build_stiffness_component_mesh(
    output_path::AbstractString;
    config::StiffnessConfig=StiffnessConfig(),
    variant::Symbol=:slotted,
    size_min_mm::Real=0.045,
    size_max_mm::Real=0.15,
)
    validate_config(config)
    variant in (:solid, :slotted) ||
        throw(ArgumentError("variant must be :solid or :slotted"))
    size_min_mm > 0 || throw(ArgumentError("minimum mesh size must be positive"))
    size_max_mm >= size_min_mm ||
        throw(ArgumentError("maximum mesh size must not be smaller than minimum"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("stiffness_$(variant)")
    surfaces = build_domain(config, variant)
    curves = tag_boundaries(surfaces, config.length_mm)
    configure_mesh(curves; size_min_mm, size_max_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
