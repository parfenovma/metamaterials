module PointRadiatorApertureMesher

using Gmsh: gmsh

export PointRadiatorApertureConfig,
       radiator_centers_mm,
       build_point_radiator_aperture_mesh

Base.@kwdef struct PointRadiatorApertureConfig
    element_count::Int = 15
    radiator_width_mm::Float64 = 1.6
    pitch_mm::Float64 = 4.8
    receiver_length_mm::Float64 = 70.0
    receiver_half_height_mm::Float64 = 45.0
    focal_distance_mm::Float64 = 35.0
end

function radiator_centers_mm(config::PointRadiatorApertureConfig)
    half = (config.element_count - 1) ÷ 2
    collect((-half):half) .* config.pitch_mm
end

function validate_config(config)
    isodd(config.element_count) || throw(ArgumentError("element count must be odd"))
    config.element_count >= 3 || throw(ArgumentError("at least three radiators are required"))
    0 < config.radiator_width_mm < config.pitch_mm ||
        throw(ArgumentError("radiator width must lie inside one pitch"))
    config.receiver_length_mm > config.focal_distance_mm > 0 ||
        throw(ArgumentError("focus must lie inside the receiver"))
    maximum(abs, radiator_centers_mm(config)) + config.radiator_width_mm / 2 <
        config.receiver_half_height_mm || throw(ArgumentError("aperture does not fit receiver"))
    nothing
end

function boundary_entities()
    [
        tag for (_, tag) in gmsh.model.getEntities(1)
        if length(first(gmsh.model.getAdjacencies(1, tag))) == 1
    ]
end

function add_group(dim, tags, id, name)
    isempty(tags) && return
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function build_partitioned_receiver(config)
    half_width = config.radiator_width_mm / 2
    boundaries = Float64[-config.receiver_half_height_mm, config.receiver_half_height_mm]
    for center in radiator_centers_mm(config)
        push!(boundaries, center - half_width, center + half_width)
    end
    sort!(unique!(boundaries))
    surfaces = Int[]
    for (lower, upper) in zip(boundaries[1:end-1], boundaries[2:end])
        push!(surfaces, gmsh.model.occ.addRectangle(
            0.0,
            lower * 1e-3,
            0.0,
            config.receiver_length_mm * 1e-3,
            (upper - lower) * 1e-3,
        ))
    end
    fragmented, _ = gmsh.model.occ.fragment(
        [(2, first(surfaces))],
        [(2, tag) for tag in Iterators.drop(surfaces, 1)],
    )
    unique(tag for (dimension, tag) in fragmented if dimension == 2)
end

function classify_boundaries(config)
    tolerance = 1e-9
    end_x = config.receiver_length_mm * 1e-3
    half_height = config.receiver_half_height_mm * 1e-3
    curves = boundary_entities()
    left = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] < tolerance, curves)
    radiation = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1] - end_x) < tolerance || abs(abs(center[2]) - half_height) < tolerance
    end
    radiator_groups = Vector{Vector{Int}}()
    half_width_m = config.radiator_width_mm * 0.5e-3
    for center_mm in radiator_centers_mm(config)
        center_m = center_mm * 1e-3
        group = filter(left) do tag
            y = gmsh.model.occ.getCenterOfMass(1, tag)[2]
            abs(y - center_m) < half_width_m + tolerance
        end
        length(group) == 1 || error("radiator at y=$center_mm mm was not isolated")
        push!(radiator_groups, group)
    end
    radiator_curves = reduce(vcat, radiator_groups)
    free = setdiff(curves, vcat(radiation, radiator_curves))
    radiator_groups, radiation, free
end

function configure_mesh(config, radiator_curves; size_radiator_mm, size_focus_mm, size_max_mm)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", radiator_curves)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_radiator_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.2e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 10.0e-3)

    gmsh.model.mesh.field.add("Box", 3)
    gmsh.model.mesh.field.setNumber(3, "VIn", Float64(size_focus_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "VOut", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "XMin", (config.focal_distance_mm - 12) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "XMax", (config.focal_distance_mm + 12) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "YMin", -18e-3)
    gmsh.model.mesh.field.setNumber(3, "YMax", 18e-3)

    gmsh.model.mesh.field.add("Min", 4)
    gmsh.model.mesh.field.setNumbers(4, "FieldsList", [2, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(4)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
end

function build_point_radiator_aperture_mesh(
    output_path::AbstractString;
    config::PointRadiatorApertureConfig=PointRadiatorApertureConfig(),
    size_radiator_mm::Real=0.25,
    size_focus_mm::Real=0.55,
    size_max_mm::Real=0.90,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("measured_point_radiator_aperture")
        surfaces = build_partitioned_receiver(config)
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        surfaces = last.(gmsh.model.getEntities(2))
        radiator_groups, radiation, free = classify_boundaries(config)
        for (index, group) in enumerate(radiator_groups)
            add_group(1, group, 100 + index, "Radiator_$index")
        end
        add_group(1, radiation, 201, "RadiationBoundary")
        add_group(1, free, 202, "FreeSurface")
        add_group(2, surfaces, 301, "Domain")
        configure_mesh(config, reduce(vcat, radiator_groups);
                       size_radiator_mm, size_focus_mm, size_max_mm)
        gmsh.model.mesh.generate(2)
        gmsh.write(output_path)
    finally
        gmsh.finalize()
    end
    println("[+] $output_path")
    output_path
end

end
