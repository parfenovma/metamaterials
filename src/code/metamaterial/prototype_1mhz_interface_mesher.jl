module Prototype1MHzInterfaceMesher

using Gmsh: gmsh

include(joinpath(@__DIR__, "profiles.jl"))
using .MetamaterialProfiles

export PrototypeMeshConfig, build_meshes

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT = joinpath(ROOT, "tmp", "prototype_1mhz_interface")

Base.@kwdef struct PrototypeMeshConfig
    reference_frequency_hz::Float64 = 140.0e3
    target_frequency_hz::Float64 = 1.0e6
    aluminium_depth_mm::Float64 = 6.0
    aluminium_half_height_mm::Float64 = 3.0
    coarse_size_mm::Float64 = 0.10
    fine_size_mm::Float64 = 0.075
end

scale(config::PrototypeMeshConfig) =
    config.reference_frequency_hz / config.target_frequency_hz

function scaled_geometry(config::PrototypeMeshConfig)
    factor = scale(config)
    GeometryConfig(
        length_mm=17.0 * factor,
        height_mm=7.0 * factor,
        end_margin_mm=0.4 * factor,
        samples=200,
    )
end

function add_polygon(points_mm)
    points = [
        gmsh.model.occ.addPoint(x * 1.0e-3, y * 1.0e-3, 0.0)
        for (x, y) in points_mm
    ]
    lines = [
        gmsh.model.occ.addLine(points[index], points[mod1(index + 1, length(points))])
        for index in eachindex(points)
    ]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function add_group(dim, tags, id, name)
    unique_tags = sort(unique(Int.(tags)))
    isempty(unique_tags) && error("physical group $name is empty")
    gmsh.model.addPhysicalGroup(dim, unique_tags, id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function curve_length(tag)
    mass = gmsh.model.occ.getMass(1, tag)
    Float64(mass)
end

function build_case(role::Symbol, resolution::Symbol; config=PrototypeMeshConfig())
    role in (:straight, :stg) || throw(ArgumentError("role must be :straight or :stg"))
    resolution in (:coarse, :fine) ||
        throw(ArgumentError("resolution must be :coarse or :fine"))

    geometry = scaled_geometry(config)
    factor = scale(config)
    amplitude_mm = role == :straight ? 0.0 : 2.5 * factor
    profile = SinusoidalProfile(amplitude_mm; periods=2)
    bottom, top = generate_boundary_points(profile, geometry)
    y_shift = geometry.height_mm / 2
    polymer_points = [(x, y - y_shift) for (x, y) in vcat(bottom, top)]

    gmsh.clear()
    gmsh.model.add("prototype_1mhz_$(role)_$(resolution)")
    polymer = add_polygon(polymer_points)
    aluminium = gmsh.model.occ.addRectangle(
        geometry.length_mm * 1.0e-3,
        -config.aluminium_half_height_mm * 1.0e-3,
        0.0,
        config.aluminium_depth_mm * 1.0e-3,
        2 * config.aluminium_half_height_mm * 1.0e-3,
    )
    gmsh.model.occ.fragment([(2, polymer)], [(2, aluminium)])
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    surfaces = [tag for (dim, tag) in gmsh.model.getEntities(2) if dim == 2]
    length(surfaces) == 2 || error("expected two material surfaces, got $(length(surfaces))")
    polymer_surfaces = Int[]
    aluminium_surfaces = Int[]
    interface_x_m = geometry.length_mm * 1.0e-3
    for tag in surfaces
        center_x = gmsh.model.occ.getCenterOfMass(2, tag)[1]
        push!(center_x < interface_x_m ? polymer_surfaces : aluminium_surfaces, tag)
    end
    length(polymer_surfaces) == 1 || error("polymer surface classification failed")
    length(aluminium_surfaces) == 1 || error("aluminium surface classification failed")

    polymer_curves = unique(last.(gmsh.model.getBoundary(
        [(2, tag) for tag in polymer_surfaces], false, false, false,
    )))
    aluminium_curves = unique(last.(gmsh.model.getBoundary(
        [(2, tag) for tag in aluminium_surfaces], false, false, false,
    )))
    all_curves = sort(unique(vcat(polymer_curves, aluminium_curves)))

    tolerance_m = 1.0e-8
    source = Int[]
    interface = Int[]
    radiation = Int[]
    free = Int[]
    far_x_m = (geometry.length_mm + config.aluminium_depth_mm) * 1.0e-3
    half_height_m = config.aluminium_half_height_mm * 1.0e-3
    polymer_height_m = geometry.height_mm * 1.0e-3
    for tag in all_curves
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        x, y = center[1], center[2]
        length_m = curve_length(tag)
        if abs(x) < tolerance_m && tag in polymer_curves
            push!(source, tag)
        elseif abs(x - interface_x_m) < tolerance_m &&
               length_m <= 1.01 * polymer_height_m && tag in polymer_curves &&
               tag in aluminium_curves
            push!(interface, tag)
        elseif abs(x - far_x_m) < tolerance_m ||
               abs(abs(y) - half_height_m) < tolerance_m
            push!(radiation, tag)
        else
            push!(free, tag)
        end
    end

    length(source) == 1 || error("expected one source curve, got $(source)")
    length(interface) == 1 || error("expected one shared interface curve, got $(interface)")
    !isempty(radiation) || error("radiation boundary classification failed")

    add_group(1, source, 101, "Source")
    add_group(1, radiation, 102, "RadiationBoundary")
    add_group(1, free, 103, "FreeSurface")
    add_group(1, interface, 104, "Interface")
    add_group(2, polymer_surfaces, 201, "Polymer")
    add_group(2, aluminium_surfaces, 202, "Aluminium")

    size_mm = resolution == :coarse ? config.coarse_size_mm : config.fine_size_mm
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), size_mm * 1.0e-3)
    gmsh.option.setNumber("Mesh.MeshSizeMin", size_mm * 1.0e-3)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size_mm * 1.0e-3)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.mesh.generate(2)

    mesh_dir = joinpath(OUTPUT, "meshes")
    mkpath(mesh_dir)
    path = joinpath(mesh_dir, "mesh_$(role)_$(resolution).msh")
    gmsh.write(path)
    println("[+] $path")
    path
end

function build_meshes(; config=PrototypeMeshConfig())
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        [build_case(role, resolution; config) for
         role in (:straight, :stg), resolution in (:coarse, :fine)]
    finally
        gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    build_meshes()
end

end
