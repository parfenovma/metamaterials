module BaselineR0Mesher

using Gmsh: gmsh

export BaselineR0Config, total_length_mm, solid_rectangles, build_baseline_r0_mesh

Base.@kwdef struct BaselineR0Config
    narrow_length_mm::Float64 = 18.0
    lead_length_mm::Float64 = 12.0
    height_mm::Float64 = 4.2
    spine_height_mm::Float64 = 1.2
end

total_length_mm(config::BaselineR0Config) =
    config.narrow_length_mm + 2config.lead_length_mm

function validate_config(config::BaselineR0Config)
    config.narrow_length_mm > 0 || throw(ArgumentError("narrow length must be positive"))
    config.lead_length_mm > 0 || throw(ArgumentError("lead length must be positive"))
    config.height_mm > 0 || throw(ArgumentError("height must be positive"))
    0 < config.spine_height_mm < config.height_mm ||
        throw(ArgumentError("spine height must lie inside the full strip"))
    nothing
end

function solid_rectangles(config::BaselineR0Config)
    validate_config(config)
    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    [
        (0.0, 0.0, config.lead_length_mm, config.height_mm),
        (
            config.lead_length_mm,
            spine_bottom,
            config.narrow_length_mm,
            config.spine_height_mm,
        ),
        (
            config.lead_length_mm + config.narrow_length_mm,
            0.0,
            config.lead_length_mm,
            config.height_mm,
        ),
    ]
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
    length_m = total_length_mm(config) * 1e-3
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

function build_baseline_r0_mesh(
    output_path::AbstractString;
    config::BaselineR0Config=BaselineR0Config(),
    size_min_mm::Real=0.10,
    size_max_mm::Real=0.40,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("baseline_R0")
    rectangle_tags = add_rectangle_mm.(solid_rectangles(config))
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(rectangle_tags))],
        [(2, tag) for tag in Iterators.drop(rectangle_tags, 1)],
    )
    gmsh.model.occ.synchronize()
    surface_tags = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surface_tags) && error("OCC fuse produced no baseline_R0 surface")
    source, microphone, free_surface = classify_boundaries(surface_tags, config)
    length(source) == 1 || error("baseline_R0 must have exactly one left port")
    length(microphone) == 1 || error("baseline_R0 must have exactly one right port")

    gmsh.model.addPhysicalGroup(1, source, 101)
    gmsh.model.setPhysicalName(1, 101, "Source")
    gmsh.model.addPhysicalGroup(1, microphone, 102)
    gmsh.model.setPhysicalName(1, 102, "Microphone")
    gmsh.model.addPhysicalGroup(1, free_surface, 103)
    gmsh.model.setPhysicalName(1, 103, "FreeSurface")
    gmsh.model.addPhysicalGroup(2, surface_tags, 201)
    gmsh.model.setPhysicalName(2, 201, "Domain")

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", vcat(source, microphone, free_surface))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_min_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.10e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.2e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
