module SideMassComponentMesher

using Gmsh: gmsh

export build_bright_component_mesh, build_dark_component_mesh

function add_rectangle_mm(x_mm, y_mm, width_mm, height_mm)
    gmsh.model.occ.addRectangle(
        x_mm * 1e-3,
        y_mm * 1e-3,
        0.0,
        width_mm * 1e-3,
        height_mm * 1e-3,
    )
end

function fused_surface(rectangles)
    rectangle_tags = [add_rectangle_mm(rectangle...) for rectangle in rectangles]
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(rectangle_tags))],
        [(2, tag) for tag in Iterators.drop(rectangle_tags, 1)],
    )
    gmsh.model.occ.synchronize()
    surfaces = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surfaces) && error("component fuse produced no surface")
    surfaces
end

function physical_groups(surfaces, is_fixed)
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    curves = unique(tag for (dimension, tag) in boundary if dimension == 1)
    fixed = [tag for tag in curves if is_fixed(gmsh.model.occ.getCenterOfMass(1, tag))]
    free = setdiff(curves, fixed)
    isempty(fixed) && error("fixed interface was not detected")
    gmsh.model.addPhysicalGroup(1, fixed, 101)
    gmsh.model.setPhysicalName(1, 101, "FixedInterface")
    gmsh.model.addPhysicalGroup(1, free, 103)
    gmsh.model.setPhysicalName(1, 103, "FreeSurface")
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
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.8e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
end

function build_bright_component_mesh(
    output_path::AbstractString;
    config,
    size_min_mm::Real=0.04,
    size_max_mm::Real=0.16,
)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("bright_fixed_interface_component")
    spine_top = (config.height_mm + config.spine_height_mm) / 2
    mass_bottom = config.height_mm - config.wall_margin_mm - config.bright_mass_height_mm
    neck_height = mass_bottom - spine_top
    neck_height > 0 || error("bright component has non-positive neck height")
    neck_left = (config.bright_mass_length_mm - config.bright_neck_width_mm) / 2
    rectangles = [
        (0.0, neck_height, config.bright_mass_length_mm, config.bright_mass_height_mm),
        (neck_left, 0.0, config.bright_neck_width_mm, neck_height),
    ]
    surfaces = fused_surface(rectangles)
    tolerance = 1.0e-9
    curves = physical_groups(surfaces, center -> center[2] < tolerance)
    configure_mesh(curves; size_min_mm, size_max_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

function build_dark_component_mesh(
    output_path::AbstractString;
    config,
    size_min_mm::Real=0.04,
    size_max_mm::Real=0.16,
)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("dark_fixed_interface_component")
    bright_right = config.bright_center_x_mm + config.bright_mass_length_mm / 2
    dark_left = config.dark_center_x_mm - config.dark_mass_length_mm / 2
    bridge_length = dark_left - bright_right
    bridge_length > 0 || error("dark component has non-positive bridge length")
    bridge_bottom = (config.dark_mass_height_mm - config.dark_bridge_height_mm) / 2
    rectangles = [
        (0.0, bridge_bottom, bridge_length, config.dark_bridge_height_mm),
        (bridge_length, 0.0, config.dark_mass_length_mm, config.dark_mass_height_mm),
    ]
    surfaces = fused_surface(rectangles)
    tolerance = 1.0e-9
    curves = physical_groups(surfaces, center -> center[1] < tolerance)
    configure_mesh(curves; size_min_mm, size_max_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
