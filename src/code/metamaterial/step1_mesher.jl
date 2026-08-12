using Gmsh: gmsh

if !isdefined(@__MODULE__, :MetamaterialProfiles)
    include(joinpath(@__DIR__, "profiles.jl"))
end
using .MetamaterialProfiles

const DEFAULT_AMPLITUDES = [0.0, 0.6, 1.2, 1.6, 1.9, 2.2, 2.5, 2.6, 2.7, 2.8, 3.0, 3.2, 3.4, 3.6]

Base.@kwdef struct MeshConfig
    output_dir::String = "1_meshes"
    size_min_mm::Float64 = 0.15
    size_max_mm::Float64 = 0.60
    distance_min_mm::Float64 = 0.2
    distance_max_mm::Float64 = 2.0
end

function add_polygon_surface(points_mm)
    points_m = [(x * 1.0e-3, y * 1.0e-3) for (x, y) in points_mm]
    point_tags = Int[]
    previous = (-Inf, -Inf)

    for (x, y) in points_m
        if hypot(x - previous[1], y - previous[2]) > 1.0e-9
            push!(point_tags, gmsh.model.occ.addPoint(x, y, 0.0))
            previous = (x, y)
        end
    end

    line_tags = Int[]
    for i in eachindex(point_tags)
        next_i = i == length(point_tags) ? 1 : i + 1
        push!(line_tags, gmsh.model.occ.addLine(point_tags[i], point_tags[next_i]))
    end

    loop_tag = gmsh.model.occ.addCurveLoop(line_tags)
    surface_tag = gmsh.model.occ.addPlaneSurface([loop_tag])
    surface_tag, line_tags
end

function classify_boundaries(line_tags, geometry::GeometryConfig)
    source, microphone, free_surface = Int[], Int[], Int[]
    length_m = geometry.length_mm * 1.0e-3

    for tag in line_tags
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        x = center[1]
        if x < 1.0e-8
            push!(source, tag)
        elseif x > length_m - 1.0e-8
            push!(microphone, tag)
        else
            push!(free_surface, tag)
        end
    end

    source, microphone, free_surface
end

function build_mesh(
    profile::WallProfile;
    geometry::GeometryConfig=GeometryConfig(),
    mesh::MeshConfig=MeshConfig(),
)
    validate_geometry(profile, geometry)
    mkpath(mesh.output_dir)

    slug = profile_slug(profile)
    filename = joinpath(mesh.output_dir, "mesh_$(slug).msh")

    gmsh.clear()
    gmsh.model.add("Metamaterial_$(slug)")

    if iszero(amplitude_mm(profile))
        surface_tag = gmsh.model.occ.addRectangle(
            0.0,
            0.0,
            0.0,
            geometry.length_mm * 1.0e-3,
            geometry.height_mm * 1.0e-3,
        )
        gmsh.model.occ.synchronize()
        line_tags = [tag for (_, tag) in gmsh.model.getEntities(1)]
    else
        bottom, top = generate_boundary_points(profile, geometry)
        surface_tag, line_tags = add_polygon_surface(vcat(bottom, top))
        gmsh.model.occ.synchronize()
    end

    source, microphone, free_surface = classify_boundaries(line_tags, geometry)
    isempty(source) && error("source boundary was not detected")
    isempty(microphone) && error("microphone boundary was not detected")

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
    gmsh.model.mesh.field.setNumber(2, "SizeMin", mesh.size_min_mm * 1.0e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", mesh.size_max_mm * 1.0e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", mesh.distance_min_mm * 1.0e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", mesh.distance_max_mm * 1.0e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    println("  [+] Mesh saved: $filename")
    filename
end

function generate_meshes(
    profiles::AbstractVector{<:WallProfile};
    geometry::GeometryConfig=GeometryConfig(),
    mesh::MeshConfig=MeshConfig(),
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        [build_mesh(profile; geometry, mesh) for profile in profiles]
    finally
        gmsh.finalize()
    end
end

function main()
    # Legacy is the default so existing filenames and datasets remain usable.
    profiles = WallProfile[LegacySinusoidalProfile(A) for A in DEFAULT_AMPLITUDES]
    println("=== Generating meshes ===")
    generate_meshes(profiles)
    println("=== Generation completed ===")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
