module TopologyHybridTTDGRINMesher

using Gmsh: gmsh

if !isdefined(parentmodule(@__MODULE__), :TopologyHybridTTDGRINDesign)
    include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
    using .TopologyHybridTTDGRINDesign
else
    using ..TopologyHybridTTDGRINDesign
end

export HybridChannelMeshConfig,
       build_hybrid_channel_mesh,
       variant_geometry,
       material_tag_name

Base.@kwdef struct HybridChannelMeshConfig
    design::HybridTTDGRINConfig = HybridTTDGRINConfig()
    centerline_samples::Int = 1201
    block_samples::Int = 28
    transition_steps::Int = 5
    mesh_size_path_mm::Float64 = 0.16
    mesh_size_port_mm::Float64 = 0.16
    mesh_size_max_mm::Float64 = 0.30
end

material_key(fraction) = round(Int, 1000Float64(fraction))
material_tag_name(key::Integer) = "MaterialF$(lpad(string(Int(key)), 4, '0'))"

function variant_geometry(variant::Symbol, config::HybridTTDGRINConfig)
    summary = design_summary(config)
    active = config.device_depth_mm - 2config.lead_length_mm
    if variant == :straight
        return (amplitude_mm=0.0, extra_path_mm=0.0, delay_s=0.0, use_insert=false)
    elseif variant == :path_only
        return (
            amplitude_mm=summary.hybrid_amplitude_mm,
            extra_path_mm=summary.hybrid_extra_path_mm,
            delay_s=maximum_delay_s(config),
            use_insert=false,
        )
    elseif variant == :hybrid
        return (
            amplitude_mm=summary.hybrid_amplitude_mm,
            extra_path_mm=summary.hybrid_extra_path_mm,
            delay_s=maximum_delay_s(config),
            use_insert=true,
        )
    elseif variant == :pure_same_depth
        return (
            amplitude_mm=summary.pure_same_depth_amplitude_mm,
            extra_path_mm=summary.pure_ttd_extra_path_mm,
            delay_s=maximum_delay_s(config),
            use_insert=false,
        )
    end
    throw(ArgumentError("unknown hybrid-channel variant: $variant"))
end

function centerline_table(mesh_config, variant)
    config = mesh_config.design
    geometry = variant_geometry(variant, config)
    active = config.device_depth_mm - 2config.lead_length_mm
    x = collect(range(0.0, config.device_depth_mm; length=mesh_config.centerline_samples))
    y = zeros(length(x))
    slope = zeros(length(x))
    for index in eachindex(x)
        local_x = x[index] - config.lead_length_mm
        if 0 <= local_x <= active
            y[index] = TopologyHybridTTDGRINDesign.smooth_centerline_mm(
                local_x,
                active,
                geometry.amplitude_mm,
            )
            slope[index] = TopologyHybridTTDGRINDesign.smooth_slope(
                local_x,
                active,
                geometry.amplitude_mm,
            )
        end
    end
    s = zeros(length(x))
    for index in 2:length(x)
        dx = x[index] - x[index - 1]
        s[index] = s[index - 1] + dx *
                   (hypot(1.0, slope[index - 1]) + hypot(1.0, slope[index])) / 2
    end
    (; x, y, slope, s, geometry)
end

function linear_interpolate(x, values, query)
    query <= first(x) && return first(values)
    query >= last(x) && return last(values)
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function schedule_blocks(table, mesh_config, variant)
    total_s = last(table.s)
    if variant != :hybrid
        return [(s0=0.0, s1=total_s, copper_fraction=0.0,
                 width_mm=mesh_config.design.guide_width_mm)]
    end
    insert = slow_insert_profile(
        table.geometry.delay_s,
        mesh_config.design;
        transition_steps=mesh_config.transition_steps,
    )
    start_s = (total_s - insert.total_length_mm) / 2
    blocks = NamedTuple[]
    start_s > 1e-10 && push!(blocks, (
        s0=0.0,
        s1=start_s,
        copper_fraction=0.0,
        width_mm=mesh_config.design.guide_width_mm,
    ))
    cursor = start_s
    for segment in insert.segments
        push!(blocks, (
            s0=cursor,
            s1=cursor + segment.length_mm,
            copper_fraction=segment.copper_fraction,
            width_mm=segment.width_mm,
        ))
        cursor += segment.length_mm
    end
    cursor < total_s - 1e-10 && push!(blocks, (
        s0=cursor,
        s1=total_s,
        copper_fraction=0.0,
        width_mm=mesh_config.design.guide_width_mm,
    ))
    blocks
end

function point_at_s(table, s_value)
    x = linear_interpolate(table.s, table.x, s_value)
    y = linear_interpolate(table.x, table.y, x)
    slope = linear_interpolate(table.x, table.slope, x)
    (x=x, y=y, slope=slope)
end

function add_strip_surface(table, block, sample_count)
    sample_count >= 2 || throw(ArgumentError("block needs at least two samples"))
    coordinates = [
        point_at_s(table, s)
        for s in range(block.s0, block.s1; length=sample_count)
    ]
    half = block.width_mm / 2
    lower = Tuple{Float64, Float64}[]
    upper = Tuple{Float64, Float64}[]
    for point in coordinates
        norm = hypot(1.0, point.slope)
        nx, ny = -point.slope / norm, 1 / norm
        push!(lower, (point.x - half * nx, point.y - half * ny))
        push!(upper, (point.x + half * nx, point.y + half * ny))
    end
    polygon = vcat(lower, reverse(upper))
    points = [gmsh.model.occ.addPoint(x * 1e-3, y * 1e-3, 0.0) for (x, y) in polygon]
    lines = [
        gmsh.model.occ.addLine(points[index], points[mod1(index + 1, length(points))])
        for index in eachindex(points)
    ]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function fragment_surfaces(tags)
    length(tags) == 1 && return tags, [[(2, only(tags))]]
    output, mapping = gmsh.model.occ.fragment(
        [(2, first(tags))],
        [(2, tag) for tag in Iterators.drop(tags, 1)],
    )
    surfaces = unique(tag for (dimension, tag) in output if dimension == 2)
    surfaces, mapping
end

function boundary_entities(surfaces)
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    unique(tag for (dimension, tag) in boundary if dimension == 1 &&
           length(first(gmsh.model.getAdjacencies(1, tag))) == 1)
end

function build_hybrid_channel_mesh(
    output_path::AbstractString;
    variant::Symbol=:hybrid,
    config::HybridChannelMeshConfig=HybridChannelMeshConfig(),
)
    validate(config.design)
    config.centerline_samples >= 201 || throw(ArgumentError("centerline is undersampled"))
    config.block_samples >= 4 || throw(ArgumentError("material block is undersampled"))
    mkpath(dirname(output_path))
    table = centerline_table(config, variant)
    blocks = schedule_blocks(table, config, variant)
    gmsh.clear()
    gmsh.model.add("hybrid_ttd_grin_$(variant)")
    original_tags = [
        add_strip_surface(table, block, config.block_samples)
        for block in blocks
    ]
    surfaces, mapping = fragment_surfaces(original_tags)
    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    groups = Dict{Int, Vector{Int}}()
    if length(original_tags) == 1
        groups[0] = surfaces
    else
        for (block, mapped) in zip(blocks, mapping)
            key = material_key(block.copper_fraction)
            tags = get!(groups, key, Int[])
            append!(tags, [tag for (dimension, tag) in mapped if dimension == 2])
        end
        for key in keys(groups)
            groups[key] = unique(filter(in(surfaces), groups[key]))
        end
    end
    isempty(groups) && error("material surface classification failed")

    curves = boundary_entities(surfaces)
    tolerance = 2e-7
    source = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] < tolerance, curves)
    microphone = filter(
        tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] >
               config.design.device_depth_mm * 1e-3 - tolerance,
        curves,
    )
    length(source) == 1 || error("hybrid channel needs one source port; found $(length(source))")
    length(microphone) == 1 || error("hybrid channel needs one output port; found $(length(microphone))")
    free = setdiff(curves, vcat(source, microphone))
    for (id, tags, name) in (
        (101, source, "Source"),
        (102, microphone, "Microphone"),
        (103, free, "FreeSurface"),
    )
        gmsh.model.addPhysicalGroup(1, tags, id)
        gmsh.model.setPhysicalName(1, id, name)
    end
    for (offset, key) in enumerate(sort(collect(keys(groups))))
        tags = groups[key]
        isempty(tags) && continue
        gmsh.model.addPhysicalGroup(2, tags, 200 + offset)
        gmsh.model.setPhysicalName(2, 200 + offset, material_tag_name(key))
    end

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", curves)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", config.mesh_size_path_mm * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", config.mesh_size_max_mm * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.05e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.0e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    (
        mesh_path=abspath(output_path),
        variant,
        blocks,
        material_keys=sort(collect(keys(groups))),
        centreline_length_mm=last(table.s),
        amplitude_mm=table.geometry.amplitude_mm,
        extra_path_mm=table.geometry.extra_path_mm,
    )
end

end
