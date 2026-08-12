module InlineMassComponentMesher

using Gmsh: gmsh

export InlineMassConfig,
       validate_config,
       component_length_mm,
       solid_rectangles,
       build_inline_mass_component_mesh

"""Fixed-interface longitudinal mass between two identical axial necks."""
Base.@kwdef struct InlineMassConfig
    port_height_mm::Float64 = 4.2
    mass_length_mm::Float64 = 2.0
    mass_height_mm::Float64 = 3.0
    neck_length_mm::Float64 = 0.5
    neck_height_mm::Float64 = 0.6
    minimum_feature_mm::Float64 = 0.35
end

component_length_mm(config::InlineMassConfig) =
    config.mass_length_mm + 2config.neck_length_mm

function validate_config(config::InlineMassConfig)
    config.port_height_mm > 0 || throw(ArgumentError("port height must be positive"))
    config.mass_length_mm > 0 || throw(ArgumentError("mass length must be positive"))
    config.mass_height_mm > 0 || throw(ArgumentError("mass height must be positive"))
    config.neck_length_mm > 0 || throw(ArgumentError("neck length must be positive"))
    config.neck_height_mm > 0 || throw(ArgumentError("neck height must be positive"))
    config.minimum_feature_mm > 0 ||
        throw(ArgumentError("minimum feature must be positive"))
    config.mass_height_mm <= config.port_height_mm ||
        throw(ArgumentError("mass cannot be taller than the surrounding port"))
    config.neck_height_mm < config.mass_height_mm ||
        throw(ArgumentError("neck must be thinner than the central mass"))
    config.neck_height_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("neck violates the minimum printable feature"))
    config.neck_length_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("neck length violates the minimum printable feature"))
    nothing
end

"""Rectangles `(x, y, width, height)` in millimetres, centred about `y = 0`."""
function solid_rectangles(config::InlineMassConfig)
    validate_config(config)
    mass_bottom = -config.mass_height_mm / 2
    neck_bottom = -config.neck_height_mm / 2
    [
        (0.0, neck_bottom, config.neck_length_mm, config.neck_height_mm),
        (
            config.neck_length_mm,
            mass_bottom,
            config.mass_length_mm,
            config.mass_height_mm,
        ),
        (
            config.neck_length_mm + config.mass_length_mm,
            neck_bottom,
            config.neck_length_mm,
            config.neck_height_mm,
        ),
    ]
end

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
    tags = [add_rectangle_mm(rectangle...) for rectangle in rectangles]
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(tags))],
        [(2, tag) for tag in Iterators.drop(tags, 1)],
    )
    gmsh.model.occ.synchronize()
    surfaces = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surfaces) && error("inline-mass fuse produced no surface")
    surfaces
end

function tag_boundaries(surfaces, total_length_mm)
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    curves = unique(tag for (dimension, tag) in boundary if dimension == 1)
    total_length_m = total_length_mm * 1e-3
    tolerance = max(1.0e-10, total_length_m * 1.0e-8)
    fixed = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        center[1] < tolerance || center[1] > total_length_m - tolerance
    end
    length(fixed) == 2 || error("expected two fixed neck interfaces, found $(length(fixed))")
    free = setdiff(curves, fixed)
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
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.7e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
end

function build_inline_mass_component_mesh(
    output_path::AbstractString;
    config::InlineMassConfig=InlineMassConfig(),
    size_min_mm::Real=0.05,
    size_max_mm::Real=0.16,
)
    validate_config(config)
    size_min_mm > 0 || throw(ArgumentError("minimum mesh size must be positive"))
    size_max_mm >= size_min_mm ||
        throw(ArgumentError("maximum mesh size must not be smaller than minimum"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("inline_mass_fixed_interface")
    surfaces = fused_surface(solid_rectangles(config))
    curves = tag_boundaries(surfaces, component_length_mm(config))
    configure_mesh(curves; size_min_mm, size_max_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
