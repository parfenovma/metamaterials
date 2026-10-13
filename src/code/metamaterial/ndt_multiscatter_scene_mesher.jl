module NDTMultiscatterSceneMesher

using Gmsh: gmsh

if !isdefined(parentmodule(@__MODULE__), :NDTBenchmark)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "ndt_benchmark.jl"))
end
if !isdefined(parentmodule(@__MODULE__), :NDTElementArrayMesher)
    Base.include(parentmodule(@__MODULE__), joinpath(@__DIR__, "ndt_element_array_mesher.jl"))
end
using ..NDTBenchmark
using ..NDTElementArrayMesher

export build_multiscatter_scene_meshes

function add_physical_group(dim, tags, id, name)
    isempty(tags) && return
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function external_boundary_entities()
    [
        tag for (_, tag) in gmsh.model.getEntities(1)
        if length(first(gmsh.model.getAdjacencies(1, tag))) == 1
    ]
end

function nearest_surface(surfaces, center_x_m, center_y_m, area_m2)
    scores = map(surfaces) do tag
        center = gmsh.model.occ.getCenterOfMass(2, tag)
        mass = gmsh.model.occ.getMass(2, tag)
        hypot(center[1] - center_x_m, center[2] - center_y_m) /
            max(sqrt(area_m2), eps()) + abs(mass - area_m2) / area_m2
    end
    surfaces[argmin(scores)]
end

function nearest_interface(curves, center_x_m, center_y_m, radius_m)
    scores = map(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        mass = gmsh.model.occ.getMass(1, tag)
        hypot(center[1] - center_x_m, center[2] - center_y_m) /
            max(radius_m, eps()) + abs(mass - 2pi * radius_m) / (2pi * radius_m)
    end
    curves[argmin(scores)]
end

function set_specimen_group(kept_surfaces)
    gmsh.model.removePhysicalGroups([(2, 201)])
    add_physical_group(2, kept_surfaces, 201, "Specimen")
end

function write_variant(master_path, output_path, all_surfaces, removed_surfaces)
    gmsh.clear()
    gmsh.open(master_path)
    removed = Set(removed_surfaces)
    kept = filter(tag -> !(tag in removed), all_surfaces)
    set_specimen_group(kept)
    gmsh.model.mesh.clear([(2, tag) for tag in removed_surfaces])
    gmsh.write(output_path)
end

"""
Create one H0 and one H1 mesh per lateral target from a shared master mesh.
Background and selected target surface elements are removed without remeshing.
"""
function build_multiscatter_scene_meshes(
    output_directory::AbstractString,
    background_voids::AbstractVector{NDTDefectCase};
    target_y_mm=(-28.0, -14.0, 0.0, 14.0, 28.0),
    target_x_mm::Real=60.0,
    target_radius_mm::Real=1.0,
    config::NDTStandConfig=NDTStandConfig(mesh_size_mm=1.0, defect_mesh_size_mm=0.25),
    element_count::Integer=17,
)
    validate(config)
    all(defect -> defect.kind == :circular_void, background_voids) ||
        throw(ArgumentError("F1c backgrounds must be circular voids"))
    targets = [
        NDTDefectCase(
            id=Symbol("target_y$(replace(string(y), '.' => 'p', '-' => 'm'))"),
            center_x_mm=Float64(target_x_mm), center_y_mm=Float64(y),
            radius_mm=Float64(target_radius_mm),
        )
        for y in target_y_mm
    ]
    validate.(background_voids, Ref(config))
    validate.(targets, Ref(config))
    mkpath(output_directory)
    master_path = joinpath(output_directory, "master_filled.msh")
    h0_path = joinpath(output_directory, "h0_background.msh")
    h1_paths = Dict(target.id => joinpath(output_directory, "h1_$(target.id).msh") for target in targets)
    length_m = config.specimen_length_mm * 1.0e-3
    half_height_m = config.specimen_height_mm * 0.5e-3
    tolerance_m = 1.0e-8

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("ndt_f1c_multiscatter_scene")
        specimen = NDTElementArrayMesher.add_segmented_specimen(
            length_m, half_height_m, element_count,
        )
        inclusions = vcat(background_voids, targets)
        disks = [
            gmsh.model.occ.addDisk(
                defect.center_x_mm * 1.0e-3,
                defect.center_y_mm * 1.0e-3,
                0.0,
                defect.radius_mm * 1.0e-3,
                defect.radius_mm * 1.0e-3,
            )
            for defect in inclusions
        ]
        gmsh.model.occ.fragment([(2, specimen)], [(2, tag) for tag in disks])
        gmsh.model.occ.synchronize()

        surfaces = [tag for (_, tag) in gmsh.model.getEntities(2)]
        inclusion_surfaces = Dict{Symbol,Int}()
        used_surfaces = Set{Int}()
        for defect in inclusions
            radius_m = defect.radius_mm * 1.0e-3
            candidates = filter(tag -> !(tag in used_surfaces), surfaces)
            tag = nearest_surface(
                candidates,
                defect.center_x_mm * 1.0e-3,
                defect.center_y_mm * 1.0e-3,
                pi * radius_m^2,
            )
            inclusion_surfaces[defect.id] = tag
            push!(used_surfaces, tag)
        end
        outer_surfaces = filter(tag -> !(tag in used_surfaces), surfaces)
        length(outer_surfaces) == 1 || error("expected one outer specimen surface")

        source, radiation, backwall = Int[], Int[], Int[]
        for tag in external_boundary_entities()
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
        internal_curves = [
            tag for (_, tag) in gmsh.model.getEntities(1)
            if length(first(gmsh.model.getAdjacencies(1, tag))) == 2
        ]
        interface_curves = Dict{Symbol,Int}()
        used_curves = Set{Int}()
        for defect in inclusions
            candidates = filter(tag -> !(tag in used_curves), internal_curves)
            tag = nearest_interface(
                candidates,
                defect.center_x_mm * 1.0e-3,
                defect.center_y_mm * 1.0e-3,
                defect.radius_mm * 1.0e-3,
            )
            interface_curves[defect.id] = tag
            push!(used_curves, tag)
        end
        sort!(source; by=tag -> gmsh.model.occ.getCenterOfMass(1, tag)[2])
        length(source) == element_count || error("source segmentation was lost")
        length(backwall) == 1 || error("free backwall was not identified")

        add_physical_group(1, source, 101, "Source")
        add_physical_group(1, radiation, 102, "RadiationBoundary")
        add_physical_group(1, backwall, 103, "Backwall")
        add_physical_group(1, collect(values(interface_curves)), 104, "ScattererBoundary")
        add_physical_group(1, vcat(backwall, collect(values(interface_curves))), 105, "FreeSurface")
        add_physical_group(2, surfaces, 201, "Specimen")
        for (index, tag) in enumerate(source)
            add_physical_group(1, [tag], 300 + index, "SourceElement$(index)")
        end
        for (index, target) in enumerate(targets)
            add_physical_group(1, [interface_curves[target.id]], 400 + index, "TargetBoundary$(index)")
        end
        for (index, background) in enumerate(background_voids)
            add_physical_group(1, [interface_curves[background.id]], 500 + index, "BackgroundBoundary$(index)")
        end

        all_interfaces = collect(values(interface_curves))
        gmsh.option.setNumber("Mesh.MeshSizeMin", config.defect_mesh_size_mm * 1.0e-3)
        gmsh.option.setNumber("Mesh.MeshSizeMax", config.mesh_size_mm * 1.0e-3)
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", all_interfaces)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", config.defect_mesh_size_mm * 1.0e-3)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", config.mesh_size_mm * 1.0e-3)
        gmsh.model.mesh.field.setNumber(2, "DistMin", 0.25e-3)
        gmsh.model.mesh.field.setNumber(2, "DistMax", config.defect_transition_distance_mm * 1.0e-3)
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
        gmsh.model.mesh.generate(2)
        gmsh.write(master_path)

        background_surfaces = [inclusion_surfaces[defect.id] for defect in background_voids]
        write_variant(master_path, h0_path, surfaces, background_surfaces)
        for target in targets
            write_variant(
                master_path,
                h1_paths[target.id],
                surfaces,
                vcat(background_surfaces, inclusion_surfaces[target.id]),
            )
        end
        (
            master_path=abspath(master_path),
            h0_path=abspath(h0_path),
            h1_paths=Dict(id => abspath(path) for (id, path) in h1_paths),
            target_ids=getproperty.(targets, :id),
            target_y_mm=Float64.(collect(target_y_mm)),
            background_ids=getproperty.(background_voids, :id),
            common_master_mesh=true,
        )
    finally
        gmsh.finalize()
    end
end

end
