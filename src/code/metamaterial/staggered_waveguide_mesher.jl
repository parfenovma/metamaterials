module StaggeredWaveguideMesher

using Gmsh: gmsh
include(joinpath(@__DIR__, "profiles.jl"))
using .MetamaterialProfiles

export extended_boundaries, build_staggered_waveguide_mesh, fixed_period_geometry

function fixed_period_geometry(periods; period_mm=6.48, samples_per_period=100)
    periods>0 || throw(ArgumentError("positive number of wall periods required"))
    GeometryConfig(length_mm=0.8+(periods+0.5)*period_mm,
        samples=samples_per_period*periods)
end

"""Translate the existing 17 mm STG core; never stretch its period or end transitions."""
function extended_boundaries(; lead_mm=0.0, straight=false,
    geometry=GeometryConfig(), amplitude_mm=2.5, periods=2)
    lead_mm >= 0 || throw(ArgumentError("lead length must be nonnegative"))
    profile = SinusoidalProfile(straight ? 0.0 : amplitude_mm; periods)
    bottom, top = generate_boundary_points(profile, geometry)
    bottom = [(x + lead_mm, y) for (x, y) in bottom]
    top = [(x + lead_mm, y) for (x, y) in top]
    if lead_mm > 0
        pushfirst!(bottom, (0.0, 0.0))
        push!(bottom, (geometry.length_mm + 2lead_mm, 0.0))
        pushfirst!(top, (geometry.length_mm + 2lead_mm, geometry.height_mm))
        push!(top, (0.0, geometry.height_mm))
    end
    bottom, top
end

function polygon(points)
    vertices = Int[]
    for p in points
        push!(vertices, gmsh.model.occ.addPoint(p[1]*1e-3, p[2]*1e-3, 0.0))
    end
    lines = [gmsh.model.occ.addLine(vertices[i], vertices[mod1(i+1,length(vertices))])
             for i in eachindex(vertices)]
    gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop(lines)])
end

function build_staggered_waveguide_mesh(path; lead_mm=0.0, straight=false,
    size_min_mm=0.15, size_max_mm=0.60, periods=2, geometry=GeometryConfig())
    lead_mm >= 0 || throw(ArgumentError("lead length must be nonnegative"))
    bottom, top = extended_boundaries(; straight, geometry, periods)
    core = polygon([(x+lead_mm,y) for (x,y) in vcat(bottom,top)])
    surfaces = [core]
    if lead_mm > 0
        left = gmsh.model.occ.addRectangle(0.,0.,0.,lead_mm*1e-3,geometry.height_mm*1e-3)
        right = gmsh.model.occ.addRectangle((lead_mm+geometry.length_mm)*1e-3,0.,0.,lead_mm*1e-3,geometry.height_mm*1e-3)
        fragments, _ = gmsh.model.occ.fragment([(2,core)],[(2,left),(2,right)])
        surfaces = [tag for (dim,tag) in fragments if dim==2]
    end
    gmsh.model.occ.synchronize()
    core_tags = [tag for tag in surfaces if
        lead_mm*1e-3-1e-9 < gmsh.model.occ.getCenterOfMass(2,tag)[1] <
        (lead_mm+geometry.length_mm)*1e-3+1e-9]
    boundary = gmsh.model.getBoundary([(2,tag) for tag in surfaces],true,false,false)
    total_m = (geometry.length_mm+2lead_mm)*1e-3
    source, output, free = Int[],Int[],Int[]
    for (dim,tag) in boundary
        dim==1 || continue
        x=gmsh.model.occ.getCenterOfMass(1,tag)[1]
        push!(abs(x)<1e-9 ? source : abs(x-total_m)<1e-9 ? output : free,tag)
    end
    for (tags,id,name) in [(source,101,"Source"),(output,102,"Microphone"),(free,103,"FreeSurface")]
        gmsh.model.addPhysicalGroup(1,tags,id);gmsh.model.setPhysicalName(1,id,name)
    end
    gmsh.model.addPhysicalGroup(2,surfaces,201);gmsh.model.setPhysicalName(2,201,"Domain")
    gmsh.model.addPhysicalGroup(2,core_tags,202);gmsh.model.setPhysicalName(2,202,"Core")
    gmsh.model.mesh.field.add("Distance",1)
    gmsh.model.mesh.field.setNumbers(1,"CurvesList",vcat(source,output,free))
    gmsh.model.mesh.field.add("Threshold",2)
    for (name,value) in [("InField",1.),("SizeMin",size_min_mm*1e-3),
        ("SizeMax",size_max_mm*1e-3),("DistMin",0.2e-3),("DistMax",2e-3)]
        gmsh.model.mesh.field.setNumber(2,name,value)
    end
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary",0)
    gmsh.model.mesh.generate(2)
    mkpath(dirname(path));gmsh.write(path)
    path
end

end
