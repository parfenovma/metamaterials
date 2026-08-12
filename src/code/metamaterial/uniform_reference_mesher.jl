module UniformReferenceMesher

using Gmsh: gmsh

export UniformReferenceConfig, build_uniform_reference_mesh

Base.@kwdef struct UniformReferenceConfig
    length_mm::Float64 = 24.0
    height_mm::Float64 = 4.2
end

function validate_config(config::UniformReferenceConfig)
    config.length_mm > 0 || throw(ArgumentError("reference length must be positive"))
    config.height_mm > 0 || throw(ArgumentError("reference height must be positive"))
    nothing
end

function classify_boundaries(surface_tag, config)
    boundary = gmsh.model.getBoundary([(2, surface_tag)], false, false, false)
    source, microphone, free_surface = Int[], Int[], Int[]
    length_m = config.length_mm * 1e-3
    for (dimension, tag) in boundary
        dimension == 1 || continue
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

"""Build the full-height homogeneous strip used for absolute port calibration."""
function build_uniform_reference_mesh(
    output_path::AbstractString;
    config::UniformReferenceConfig=UniformReferenceConfig(),
    size_min_mm::Real=0.10,
    size_max_mm::Real=0.35,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("reference_uniform")
    surface_tag = gmsh.model.occ.addRectangle(
        0.0,
        0.0,
        0.0,
        config.length_mm * 1e-3,
        config.height_mm * 1e-3,
    )
    gmsh.model.occ.synchronize()
    source, microphone, free_surface = classify_boundaries(surface_tag, config)
    length(source) == 1 || error("uniform reference must have exactly one left port")
    length(microphone) == 1 || error("uniform reference must have exactly one right port")

    gmsh.model.addPhysicalGroup(1, source, 101)
    gmsh.model.setPhysicalName(1, 101, "Source")
    gmsh.model.addPhysicalGroup(1, microphone, 102)
    gmsh.model.setPhysicalName(1, 102, "Microphone")
    gmsh.model.addPhysicalGroup(1, free_surface, 103)
    gmsh.model.setPhysicalName(1, 103, "FreeSurface")
    gmsh.model.addPhysicalGroup(2, [surface_tag], 201)
    gmsh.model.setPhysicalName(2, 201, "Domain")

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", vcat(source, microphone, free_surface))
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
