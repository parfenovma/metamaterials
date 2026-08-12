module PassiveFeedManifoldMesher

using Gmsh: gmsh

export PassiveFeedManifoldConfig,
       aperture_width_mm,
       required_transformer_length_mm,
       output_center_y_mm,
       validate_config,
       build_passive_feed_manifold_mesh

Base.@kwdef struct PassiveFeedManifoldConfig
    input_width_mm::Float64 = 7.0
    channel_count::Int = 15
    channel_width_mm::Float64 = 7.0
    channel_gap_mm::Float64 = 1.2
    input_straight_length_mm::Float64 = 20.0
    transformer_length_mm::Float64 = 194.0
    output_straight_length_mm::Float64 = 20.0
    lower_band_frequency_hz::Float64 = 226.04e3
    pressure_wave_speed_m_s::Float64 = 6122.102437409232
    maximum_adiabatic_parameter::Float64 = 0.05
    profile_samples::Int = 480
end

aperture_width_mm(config) =
    config.channel_count * config.channel_width_mm +
    (config.channel_count - 1) * config.channel_gap_mm

function required_transformer_length_mm(config)
    ratio = aperture_width_mm(config) / config.input_width_mm
    k_min_per_mm = 2pi * config.lower_band_frequency_hz /
                   config.pressure_wave_speed_m_s * 1e-3
    pi * log(ratio) / (4k_min_per_mm * config.maximum_adiabatic_parameter)
end

function output_center_y_mm(config, index)
    1 <= index <= config.channel_count || throw(BoundsError())
    -aperture_width_mm(config) / 2 + config.channel_width_mm / 2 +
    (index - 1) * (config.channel_width_mm + config.channel_gap_mm)
end

function validate_config(config::PassiveFeedManifoldConfig)
    config.input_width_mm > 0 || throw(ArgumentError("input width must be positive"))
    config.channel_count >= 3 || throw(ArgumentError("at least three channels are required"))
    isodd(config.channel_count) || throw(ArgumentError("channel count must be odd"))
    config.channel_width_mm > 0 || throw(ArgumentError("channel width must be positive"))
    config.channel_gap_mm > 0 || throw(ArgumentError("channel gap must be positive"))
    config.input_straight_length_mm > 0 || throw(ArgumentError("input straight must be positive"))
    config.output_straight_length_mm > config.channel_gap_mm / 2 ||
        throw(ArgumentError("output straight must exceed the rounded slot radius"))
    config.transformer_length_mm + 1e-9 >= required_transformer_length_mm(config) ||
        throw(ArgumentError("transformer violates the adiabatic lower-band gate"))
    config.profile_samples >= 200 || throw(ArgumentError("horn profile is undersampled"))
    nothing
end

transition_q(local_x_mm, length_mm) =
    (1 - cospi(Float64(local_x_mm) / Float64(length_mm))) / 2

function horn_width_mm(config, local_x_mm)
    ratio = aperture_width_mm(config) / config.input_width_mm
    config.input_width_mm * exp(
        transition_q(local_x_mm, config.transformer_length_mm) * log(ratio),
    )
end

function outer_polygon_mm(config)
    upper = NTuple{2, Float64}[
        (0.0, config.input_width_mm / 2),
        (config.input_straight_length_mm, config.input_width_mm / 2),
    ]
    lower = NTuple{2, Float64}[
        (0.0, -config.input_width_mm / 2),
        (config.input_straight_length_mm, -config.input_width_mm / 2),
    ]
    for local_x in range(
        0.0, config.transformer_length_mm; length=config.profile_samples,
    )[2:end]
        x = config.input_straight_length_mm + local_x
        half_width = horn_width_mm(config, local_x) / 2
        push!(upper, (x, half_width))
        push!(lower, (x, -half_width))
    end
    end_x = config.input_straight_length_mm + config.transformer_length_mm +
            config.output_straight_length_mm
    push!(upper, (end_x, aperture_width_mm(config) / 2))
    push!(lower, (end_x, -aperture_width_mm(config) / 2))
    vcat(upper, reverse(lower))
end

function add_polygon_surface_mm(points_mm)
    points = [gmsh.model.occ.addPoint(x * 1e-3, y * 1e-3, 0.0) for (x, y) in points_mm]
    lines = [gmsh.model.occ.addLine(
        points[index], points[mod1(index + 1, length(points))],
    ) for index in eachindex(points)]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function slot_center_y_mm(config, index)
    lower = -aperture_width_mm(config) / 2
    lower + index * config.channel_width_mm +
    (index - 0.5) * config.channel_gap_mm
end

function add_rounded_slot(config, index)
    radius_m = config.channel_gap_mm * 0.5e-3
    mouth_x_m = (config.input_straight_length_mm + config.transformer_length_mm) * 1e-3
    end_x_m = mouth_x_m + config.output_straight_length_mm * 1e-3
    center_x_m = mouth_x_m + radius_m
    center_y_m = slot_center_y_mm(config, index) * 1e-3
    disk = gmsh.model.occ.addDisk(center_x_m, center_y_m, 0.0, radius_m, radius_m)
    rectangle = gmsh.model.occ.addRectangle(
        center_x_m, center_y_m - radius_m, 0.0,
        end_x_m - center_x_m, 2radius_m,
    )
    fused, _ = gmsh.model.occ.fuse([(2, disk)], [(2, rectangle)])
    only((dimension, tag) for (dimension, tag) in fused if dimension == 2)
end

function external_curves()
    gmsh.model.occ.synchronize()
    [tag for (dimension, tag) in gmsh.model.getEntities(1)
     if dimension == 1 && length(first(gmsh.model.getAdjacencies(1, tag))) == 1]
end

function add_physical_group(dim, tags, id, name)
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function build_passive_feed_manifold_mesh(
    output_path::AbstractString;
    config::PassiveFeedManifoldConfig=PassiveFeedManifoldConfig(),
    mesh_size_mm::Real=0.48,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("passive_feed_manifold")
    outer = add_polygon_surface_mm(outer_polygon_mm(config))
    slots = [add_rounded_slot(config, index) for index in 1:(config.channel_count - 1)]
    cut, _ = gmsh.model.occ.cut([(2, outer)], slots)
    surfaces = unique(tag for (dimension, tag) in cut if dimension == 2)
    isempty(surfaces) && error("manifold cut produced no domain")
    gmsh.model.occ.synchronize()
    curves = external_curves()
    tolerance = 1e-8
    end_x_m = (
        config.input_straight_length_mm + config.transformer_length_mm +
        config.output_straight_length_mm
    ) * 1e-3
    source = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1]) < tolerance
    end
    outputs = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1] - end_x_m) < tolerance
    end
    sort!(outputs; by=tag -> gmsh.model.occ.getCenterOfMass(1, tag)[2])
    length(source) == 1 || error("manifold source classification failed")
    length(outputs) == config.channel_count || error(
        "expected $(config.channel_count) output ports, found $(length(outputs))",
    )
    free = setdiff(curves, vcat(source, outputs))
    add_physical_group(1, source, 101, "Source")
    for (index, tag) in enumerate(outputs)
        add_physical_group(1, [tag], 110 + index, "Output$(lpad(index, 2, '0'))")
    end
    add_physical_group(1, free, 104, "FreeSurface")
    add_physical_group(2, surfaces, 201, "Domain")
    gmsh.option.setNumber("Mesh.MeshSizeMin", Float64(mesh_size_mm) * 1e-3)
    gmsh.option.setNumber("Mesh.MeshSizeMax", Float64(mesh_size_mm) * 1e-3)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    node_count = length(first(gmsh.model.mesh.getNodes()))
    element_count = sum(length, gmsh.model.mesh.getElements()[2])
    println("[+] passive feed manifold mesh: $node_count nodes, $element_count elements")
    (; node_count, element_count, output_path)
end

end
