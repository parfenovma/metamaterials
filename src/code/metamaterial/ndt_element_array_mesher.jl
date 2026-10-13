module NDTElementArrayMesher

using Gmsh: gmsh

if !isdefined(parentmodule(@__MODULE__), :NDTBenchmark)
    include(joinpath(@__DIR__, "ndt_benchmark.jl"))
end
using ..NDTBenchmark

export build_ndt_element_array_mesh, build_ndt_element_array_mesh_pair

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

function add_segmented_specimen(length_m, half_height_m, element_count)
    left_points = [
        gmsh.model.occ.addPoint(0.0, y, 0.0)
        for y in range(-half_height_m, half_height_m; length=element_count + 1)
    ]
    right_bottom = gmsh.model.occ.addPoint(length_m, -half_height_m, 0.0)
    right_top = gmsh.model.occ.addPoint(length_m, half_height_m, 0.0)
    left_lines = [
        gmsh.model.occ.addLine(left_points[index], left_points[index + 1])
        for index in 1:element_count
    ]
    top = gmsh.model.occ.addLine(last(left_points), right_top)
    backwall = gmsh.model.occ.addLine(right_top, right_bottom)
    bottom = gmsh.model.occ.addLine(right_bottom, first(left_points))
    loop = gmsh.model.occ.addCurveLoop(vcat(left_lines, top, backwall, bottom))
    gmsh.model.occ.addPlaneSurface([loop])
end

function add_defect(defect::NDTDefectCase)
    defect.kind == :none && return nothing
    x = defect.center_x_mm * 1.0e-3
    y = defect.center_y_mm * 1.0e-3
    if defect.kind == :circular_void
        radius = defect.radius_mm * 1.0e-3
        return gmsh.model.occ.addDisk(x, y, 0.0, radius, radius)
    end
    tag = gmsh.model.occ.addDisk(
        x, y, 0.0,
        defect.major_radius_mm * 1.0e-3,
        defect.minor_radius_mm * 1.0e-3,
    )
    if defect.angle_deg != 0
        gmsh.model.occ.rotate(
            [(2, tag)], x, y, 0.0, 0.0, 0.0, 1.0, deg2rad(defect.angle_deg),
        )
    end
    tag
end

"""
Build a paired free/void mesh from one conforming triangulation. The free
model contains both the outer specimen and the circular inclusion with the
same material. The void model removes only the inclusion surface elements,
so all exterior and surrounding elements are identical between hypotheses.
"""
function build_ndt_element_array_mesh_pair(
    free_output_path::AbstractString,
    void_output_path::AbstractString;
    config::NDTStandConfig=NDTStandConfig(),
    defect::NDTDefectCase=last(required_defects(config)),
    element_count::Integer=17,
)
    validate(config)
    validate(defect, config)
    defect.kind in (:circular_void, :slit) || throw(ArgumentError(
        "the common-mesh pair requires a void or elliptical slit",
    ))
    element_count >= 3 || throw(ArgumentError("at least three elements are required"))
    mkpath(dirname(free_output_path))
    mkpath(dirname(void_output_path))
    length_m = config.specimen_length_mm * 1.0e-3
    half_height_m = config.specimen_height_mm * 0.5e-3
    tolerance_m = 1.0e-8

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("ndt_element_array_common_pair")
        specimen = add_segmented_specimen(length_m, half_height_m, element_count)
        inclusion = add_defect(defect)
        gmsh.model.occ.fragment([(2, specimen)], [(2, inclusion)])
        gmsh.model.occ.synchronize()

        surfaces = [tag for (_, tag) in gmsh.model.getEntities(2)]
        length(surfaces) == 2 || error(
            "expected outer and inclusion surfaces, found $(length(surfaces))",
        )
        target_area = if defect.kind == :circular_void
            pi * (defect.radius_mm * 1.0e-3)^2
        else
            pi * defect.major_radius_mm * defect.minor_radius_mm * 1.0e-6
        end
        inclusion_surface = surfaces[argmin(abs.([
            gmsh.model.occ.getMass(2, tag) for tag in surfaces
        ] .- target_area))]
        outer_surfaces = filter(!=(inclusion_surface), surfaces)
        length(outer_surfaces) == 1 || error("outer specimen surface was not identified")

        source, radiation, backwall = Int[], Int[], Int[]
        for tag in boundary_entities()
            center = gmsh.model.occ.getCenterOfMass(1, tag)
            x, y = center[1], center[2]
            if abs(x) < tolerance_m
                push!(source, tag)
            elseif abs(x - length_m) < tolerance_m
                push!(backwall, tag)
            elseif abs(abs(y) - half_height_m) < tolerance_m
                push!(radiation, tag)
            end
        end
        defect_curves = [
            tag for (_, tag) in gmsh.model.getEntities(1)
            if length(first(gmsh.model.getAdjacencies(1, tag))) == 2
        ]
        sort!(source; by=tag -> gmsh.model.occ.getCenterOfMass(1, tag)[2])
        length(source) == element_count || error(
            "expected $element_count source segments, found $(length(source))",
        )
        !isempty(radiation) || error("lateral radiation boundary was not found")
        length(backwall) == 1 || error("expected one free backwall")
        length(defect_curves) == 1 || error("expected one circular interface")

        add_physical_group(1, source, 101, "Source")
        add_physical_group(1, radiation, 102, "RadiationBoundary")
        add_physical_group(1, backwall, 103, "Backwall")
        add_physical_group(1, defect_curves, 104, "DefectBoundary")
        add_physical_group(1, vcat(backwall, defect_curves), 105, "FreeSurface")
        add_physical_group(2, surfaces, 201, "Specimen")
        add_physical_group(2, [inclusion_surface], 202, "FilledInclusion")
        for (index, tag) in enumerate(source)
            add_physical_group(1, [tag], 300 + index, "SourceElement$(index)")
        end

        gmsh.option.setNumber("Mesh.MeshSizeMin", config.defect_mesh_size_mm * 1.0e-3)
        gmsh.option.setNumber("Mesh.MeshSizeMax", config.mesh_size_mm * 1.0e-3)
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", defect_curves)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", config.defect_mesh_size_mm * 1.0e-3)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", config.mesh_size_mm * 1.0e-3)
        gmsh.model.mesh.field.setNumber(2, "DistMin", 0.5e-3)
        gmsh.model.mesh.field.setNumber(
            2, "DistMax", config.defect_transition_distance_mm * 1.0e-3,
        )
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
        gmsh.model.mesh.generate(2)
        gmsh.write(free_output_path)

        # Preserve the already generated exterior mesh and remove only cells
        # that fill the circular inclusion. The shared interface nodes/edges
        # remain because they bound the outer specimen surface.
        gmsh.model.removePhysicalGroups([(2, 201), (2, 202)])
        add_physical_group(2, outer_surfaces, 201, "Specimen")
        gmsh.model.mesh.clear([(2, inclusion_surface)])
        gmsh.write(void_output_path)

        centers_mm = [
            gmsh.model.occ.getCenterOfMass(1, tag)[2] * 1.0e3
            for tag in source
        ]
        (
            free_mesh_path=abspath(free_output_path),
            void_mesh_path=abspath(void_output_path),
            defect_id=defect.id,
            element_count=Int(element_count),
            element_centers_mm=centers_mm,
            source_curve_count=length(source),
            radiation_curve_count=length(radiation),
            backwall_curve_count=length(backwall),
            defect_curve_count=length(defect_curves),
            outer_surface_count=length(outer_surfaces),
            common_triangulation=true,
        )
    finally
        gmsh.finalize()
    end
end

"""
Build a 2D specimen whose input face is split into individually tagged,
equal-height array elements. The downstream backwall is deliberately free and
has its own tag; only the lateral faces are absorbing.
"""
function build_ndt_element_array_mesh(
    output_path::AbstractString;
    config::NDTStandConfig=NDTStandConfig(),
    defect::NDTDefectCase=first(required_defects(config)),
    element_count::Integer=17,
)
    validate(config)
    validate(defect, config)
    element_count >= 3 || throw(ArgumentError("at least three elements are required"))
    mkpath(dirname(output_path))
    length_m = config.specimen_length_mm * 1.0e-3
    half_height_m = config.specimen_height_mm * 0.5e-3
    tolerance_m = 1.0e-8

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("ndt_element_array_$(defect.id)")
        specimen = add_segmented_specimen(length_m, half_height_m, element_count)
        defect_surface = add_defect(defect)
        specimen_surfaces = if isnothing(defect_surface)
            [specimen]
        else
            cut, _ = gmsh.model.occ.cut([(2, specimen)], [(2, defect_surface)])
            [tag for (dim, tag) in cut if dim == 2]
        end
        gmsh.model.occ.synchronize()
        length(specimen_surfaces) == 1 || error("expected one connected specimen surface")

        source, radiation, backwall, defect_curves = Int[], Int[], Int[], Int[]
        for tag in boundary_entities()
            center = gmsh.model.occ.getCenterOfMass(1, tag)
            x, y = center[1], center[2]
            if abs(x) < tolerance_m
                push!(source, tag)
            elseif abs(x - length_m) < tolerance_m
                push!(backwall, tag)
            elseif abs(abs(y) - half_height_m) < tolerance_m
                push!(radiation, tag)
            else
                push!(defect_curves, tag)
            end
        end
        sort!(source; by=tag -> gmsh.model.occ.getCenterOfMass(1, tag)[2])
        length(source) == element_count || error(
            "expected $element_count source segments, found $(length(source))",
        )
        !isempty(radiation) || error("lateral radiation boundary was not found")
        length(backwall) == 1 || error("expected one free backwall")
        defect.kind == :none || !isempty(defect_curves) ||
            error("defect boundary was not found")

        add_physical_group(1, source, 101, "Source")
        add_physical_group(1, radiation, 102, "RadiationBoundary")
        add_physical_group(1, backwall, 103, "Backwall")
        add_physical_group(1, defect_curves, 104, "DefectBoundary")
        add_physical_group(1, vcat(backwall, defect_curves), 105, "FreeSurface")
        add_physical_group(2, specimen_surfaces, 201, "Specimen")
        for (index, tag) in enumerate(source)
            add_physical_group(1, [tag], 300 + index, "SourceElement$(index)")
        end

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
                2, "DistMax", config.defect_transition_distance_mm * 1.0e-3,
            )
            gmsh.model.mesh.field.setAsBackgroundMesh(2)
        end
        gmsh.model.mesh.generate(2)
        gmsh.write(output_path)
        centers_mm = [
            gmsh.model.occ.getCenterOfMass(1, tag)[2] * 1.0e3
            for tag in source
        ]
        (
            mesh_path=abspath(output_path),
            defect_id=defect.id,
            element_count=Int(element_count),
            element_centers_mm=centers_mm,
            source_curve_count=length(source),
            radiation_curve_count=length(radiation),
            backwall_curve_count=length(backwall),
            defect_curve_count=length(defect_curves),
            specimen_surface_count=length(specimen_surfaces),
        )
    finally
        gmsh.finalize()
    end
end

end
