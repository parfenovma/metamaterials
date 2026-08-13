module VirtualNDTStandMesher

using Gmsh: gmsh

if !isdefined(parentmodule(@__MODULE__), :NDTBenchmark)
    include(joinpath(@__DIR__, "ndt_benchmark.jl"))
end
using ..NDTBenchmark

export build_virtual_ndt_stand_mesh

function add_physical_group(dim, tags, id, name)
    isempty(tags) && return
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function boundary_entities()
    [
        tag for (_, tag) in gmsh.model.getEntities(1)
        if length(first(gmsh.model.getAdjacencies(1, tag))) == 1
    ]
end

function add_defect(defect::NDTDefectCase)
    defect.kind == :none && return nothing
    x = defect.center_x_mm * 1.0e-3
    y = defect.center_y_mm * 1.0e-3
    if defect.kind == :circular_void
        radius = defect.radius_mm * 1.0e-3
        return gmsh.model.occ.addDisk(x, y, 0.0, radius, radius)
    end
    major = defect.major_radius_mm * 1.0e-3
    minor = defect.minor_radius_mm * 1.0e-3
    tag = gmsh.model.occ.addDisk(x, y, 0.0, major, minor)
    if defect.angle_deg != 0
        gmsh.model.occ.rotate([(2, tag)], x, y, 0.0, 0.0, 0.0, 1.0, deg2rad(defect.angle_deg))
    end
    tag
end

function build_virtual_ndt_stand_mesh(
    output_path::AbstractString;
    config::NDTStandConfig=NDTStandConfig(),
    defect::NDTDefectCase=first(required_defects(config)),
)
    validate(config)
    validate(defect, config)
    mkpath(dirname(output_path))
    length_m = config.specimen_length_mm * 1.0e-3
    half_height_m = config.specimen_height_mm * 0.5e-3
    tolerance_m = 1.0e-8

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("virtual_ndt_stand_$(defect.id)")
        specimen = gmsh.model.occ.addRectangle(
            0.0,
            -half_height_m,
            0.0,
            length_m,
            2half_height_m,
        )
        defect_surface = add_defect(defect)
        specimen_surfaces = if isnothing(defect_surface)
            [specimen]
        else
            cut, _ = gmsh.model.occ.cut([(2, specimen)], [(2, defect_surface)])
            [tag for (dim, tag) in cut if dim == 2]
        end
        gmsh.model.occ.synchronize()
        isempty(specimen_surfaces) && error("specimen surface was lost")

        source, radiation, defect_curves = Int[], Int[], Int[]
        for tag in boundary_entities()
            center = gmsh.model.occ.getCenterOfMass(1, tag)
            x, y = center[1], center[2]
            if abs(x) < tolerance_m
                push!(source, tag)
            elseif abs(x - length_m) < tolerance_m || abs(abs(y) - half_height_m) < tolerance_m
                push!(radiation, tag)
            else
                push!(defect_curves, tag)
            end
        end
        !isempty(source) || error("source boundary was not found")
        !isempty(radiation) || error("radiation boundary was not found")
        defect.kind == :none || !isempty(defect_curves) || error("defect boundary was not found")

        add_physical_group(1, source, 101, "Source")
        add_physical_group(1, radiation, 102, "RadiationBoundary")
        add_physical_group(1, defect_curves, 103, "DefectBoundary")
        add_physical_group(2, specimen_surfaces, 201, "Specimen")

        gmsh.option.setNumber("Mesh.MeshSizeMin", config.defect_mesh_size_mm * 1.0e-3)
        gmsh.option.setNumber("Mesh.MeshSizeMax", config.mesh_size_mm * 1.0e-3)
        if !isempty(defect_curves)
            gmsh.model.mesh.field.add("Distance", 1)
            gmsh.model.mesh.field.setNumbers(1, "CurvesList", defect_curves)
            gmsh.model.mesh.field.add("Threshold", 2)
            gmsh.model.mesh.field.setNumber(2, "InField", 1)
            gmsh.model.mesh.field.setNumber(2, "SizeMin", config.defect_mesh_size_mm * 1.0e-3)
            gmsh.model.mesh.field.setNumber(2, "SizeMax", config.mesh_size_mm * 1.0e-3)
            gmsh.model.mesh.field.setNumber(2, "DistMin", 0.5e-3)
            gmsh.model.mesh.field.setNumber(
                2,
                "DistMax",
                config.defect_transition_distance_mm * 1.0e-3,
            )
            gmsh.model.mesh.field.setAsBackgroundMesh(2)
        end
        gmsh.model.mesh.generate(2)
        gmsh.write(output_path)
        (
            mesh_path=abspath(output_path),
            defect_id=defect.id,
            source_curve_count=length(source),
            radiation_curve_count=length(radiation),
            defect_curve_count=length(defect_curves),
            specimen_surface_count=length(specimen_surfaces),
        )
    finally
        gmsh.finalize()
    end
end

end
