module TopologyHybridMonostaticMesher

using Gmsh: gmsh

include(joinpath(@__DIR__, "ndt_benchmark.jl"))
using .NDTBenchmark
include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
using .TopologyHybridTTDGRINDesign

export HybridMonostaticMeshConfig,
       build_hybrid_monostatic_mesh,
       topology_channel_geometry,
       source_pressure_scale

Base.@kwdef struct HybridMonostaticMeshConfig
    stand::NDTStandConfig = NDTStandConfig()
    design::HybridTTDGRINConfig = HybridTTDGRINConfig()
    topology::Symbol = :hybrid
    defect::NDTDefectCase = required_defects(stand)[1]
    specimen_mesh_mm::Float64 = 1.35
    guide_mesh_mm::Float64 = 0.45
    defect_mesh_mm::Float64 = 0.30
    common_extra_path_mm::Float64 = 1.50
    centreline_samples::Int = 1201
    polygon_samples_per_block::Int = 16
    transition_steps::Int = 5
end

const DEVICE_TOPOLOGIES = (:straight_footprint, :pure_ttd, :hybrid)
const ALL_TOPOLOGIES = (:bare, DEVICE_TOPOLOGIES...)

function validate(config::HybridMonostaticMeshConfig)
    config.topology in ALL_TOPOLOGIES ||
        throw(ArgumentError("unsupported topology: $(config.topology)"))
    NDTBenchmark.validate(config.stand)
    NDTBenchmark.validate(config.defect, config.stand)
    TopologyHybridTTDGRINDesign.validate(config.design)
    config.stand.aperture_height_mm == config.design.aperture_height_mm ||
        throw(ArgumentError("stand and device apertures differ"))
    config.design.device_depth_mm <= config.stand.maximum_device_depth_mm ||
        throw(ArgumentError("device exceeds the common depth gate"))
    config.specimen_mesh_mm > 0 || throw(ArgumentError("invalid specimen mesh"))
    config.guide_mesh_mm > 0 || throw(ArgumentError("invalid guide mesh"))
    config.defect_mesh_mm > 0 || throw(ArgumentError("invalid defect mesh"))
    config.common_extra_path_mm >= 0 || throw(ArgumentError("invalid common path"))
    config.centreline_samples >= 401 || throw(ArgumentError("centreline undersampled"))
    config.polygon_samples_per_block >= 6 || throw(ArgumentError("strip polygon undersampled"))
    config.transition_steps >= 2 || throw(ArgumentError("insert transition undersampled"))
    config
end

function delay_fraction(y_mm, config)
    delay = geometric_delay_s(y_mm, config.design)
    delay / maximum_delay_s(config.design)
end

"""Geometry actually used by the coupled G3 stand for one aperture channel."""
function topology_channel_geometry(y_mm::Real, config::HybridMonostaticMeshConfig)
    validate(config)
    topology = config.topology
    fraction = delay_fraction(y_mm, config)
    active_length = config.design.device_depth_mm - 2config.design.lead_length_mm
    if topology == :straight_footprint
        return (
            delay_fraction=fraction,
            relative_extra_path_mm=0.0,
            total_extra_path_mm=0.0,
            amplitude_mm=0.0,
            insert=slow_insert_profile(0.0, config.design),
        )
    end
    topology in (:pure_ttd, :hybrid) || throw(ArgumentError("bare has no channels"))
    maximum_relative = topology == :pure_ttd ?
        design_summary(config.design).pure_ttd_extra_path_mm :
        design_summary(config.design).hybrid_extra_path_mm
    relative_extra = fraction * maximum_relative
    total_extra = config.common_extra_path_mm + relative_extra
    amplitude = smooth_amplitude_mm(
        active_length,
        total_extra;
        samples=config.design.profile_samples,
    )
    insert = topology == :hybrid ?
        slow_insert_profile(
            fraction * maximum_delay_s(config.design),
            config.design;
            transition_steps=config.transition_steps,
        ) : slow_insert_profile(0.0, config.design)
    (
        delay_fraction=fraction,
        relative_extra_path_mm=relative_extra,
        total_extra_path_mm=total_extra,
        amplitude_mm=amplitude,
        insert,
    )
end

function centreline_table(y0_mm, geometry, config)
    depth = config.design.device_depth_mm
    lead = config.design.lead_length_mm
    active = depth - 2lead
    xlocal = collect(range(0.0, depth; length=config.centreline_samples))
    y = fill(Float64(y0_mm), length(xlocal))
    slope = zeros(length(xlocal))
    for index in eachindex(xlocal)
        active_x = xlocal[index] - lead
        if 0 <= active_x <= active
            y[index] += TopologyHybridTTDGRINDesign.smooth_centerline_mm(
                active_x, active, geometry.amplitude_mm,
            )
            slope[index] = TopologyHybridTTDGRINDesign.smooth_slope(
                active_x, active, geometry.amplitude_mm,
            )
        end
    end
    s = zeros(length(xlocal))
    for index in 2:length(xlocal)
        dx = xlocal[index] - xlocal[index - 1]
        s[index] = s[index - 1] + dx *
            (hypot(1.0, slope[index - 1]) + hypot(1.0, slope[index])) / 2
    end
    (; xlocal, y, slope, s)
end

function interpolate_linear(x, values, query)
    query <= first(x) && return first(values)
    query >= last(x) && return last(values)
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function point_at_s(table, s_query)
    xlocal = interpolate_linear(table.s, table.xlocal, s_query)
    (
        xlocal,
        y=interpolate_linear(table.xlocal, table.y, xlocal),
        slope=interpolate_linear(table.xlocal, table.slope, xlocal),
    )
end

function material_blocks(table, geometry, config)
    total = last(table.s)
    insert = geometry.insert
    isempty(insert.segments) && return [(
        s0=0.0, s1=total, fraction=0.0, width_mm=config.design.guide_width_mm,
    )]
    start = (total - insert.total_length_mm) / 2
    start >= 0 || error("insert does not fit the channel")
    blocks = NamedTuple[]
    start > 1e-9 && push!(blocks, (
        s0=0.0, s1=start, fraction=0.0, width_mm=config.design.guide_width_mm,
    ))
    cursor = start
    for segment in insert.segments
        push!(blocks, (
            s0=cursor,
            s1=cursor + segment.length_mm,
            fraction=segment.copper_fraction,
            width_mm=segment.width_mm,
        ))
        cursor += segment.length_mm
    end
    cursor < total - 1e-9 && push!(blocks, (
        s0=cursor, s1=total, fraction=0.0, width_mm=config.design.guide_width_mm,
    ))
    blocks
end

function add_strip_block(table, block, config)
    points_on_line = [
        point_at_s(table, s)
        for s in range(block.s0, block.s1; length=config.polygon_samples_per_block)
    ]
    half = block.width_mm / 2
    lower, upper = Tuple{Float64, Float64}[], Tuple{Float64, Float64}[]
    for point in points_on_line
        scale = hypot(1.0, point.slope)
        nx, ny = -point.slope / scale, 1 / scale
        # The device occupies x=-depth..0; xlocal runs input -> specimen.
        x = point.xlocal - config.design.device_depth_mm
        push!(lower, (x - half * nx, point.y - half * ny))
        push!(upper, (x + half * nx, point.y + half * ny))
    end
    polygon = vcat(lower, reverse(upper))
    point_tags = [gmsh.model.occ.addPoint(x * 1e-3, y * 1e-3, 0.0) for (x, y) in polygon]
    line_tags = [
        gmsh.model.occ.addLine(point_tags[index], point_tags[mod1(index + 1, length(point_tags))])
        for index in eachindex(point_tags)
    ]
    loop = gmsh.model.occ.addCurveLoop(line_tags)
    gmsh.model.occ.addPlaneSurface([loop])
end

function add_defect_surface(defect)
    defect.kind == :none && return nothing
    x = defect.center_x_mm * 1e-3
    y = defect.center_y_mm * 1e-3
    if defect.kind == :circular_void
        r = defect.radius_mm * 1e-3
        return gmsh.model.occ.addDisk(x, y, 0.0, r, r)
    end
    major = defect.major_radius_mm * 1e-3
    minor = defect.minor_radius_mm * 1e-3
    tag = gmsh.model.occ.addDisk(x, y, 0.0, major, minor)
    defect.angle_deg == 0 || gmsh.model.occ.rotate(
        [(2, tag)], x, y, 0.0, 0.0, 0.0, 1.0, deg2rad(defect.angle_deg),
    )
    tag
end

function external_boundary_curves(surface_tags)
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surface_tags], false, false, false)
    unique(tag for (dimension, tag) in boundary if dimension == 1 &&
           length(first(gmsh.model.getAdjacencies(1, tag))) == 1)
end

function add_group(dimension, tags, id, name)
    tags = sort(unique(tags))
    isempty(tags) && return
    gmsh.model.addPhysicalGroup(dimension, tags, id)
    gmsh.model.setPhysicalName(dimension, id, name)
end

material_key(fraction) = round(Int, 1000Float64(fraction))
material_name(key) = "MaterialF$(lpad(string(key), 4, '0'))"

"""Pressure multiplier that equalises nominal incident energy across source lengths."""
source_pressure_scale(source_length_mm::Real, stand::NDTStandConfig=NDTStandConfig()) =
    sqrt(stand.aperture_height_mm / Float64(source_length_mm))

function minimum_channel_clearance_mm(config)
    config.topology in (:pure_ttd, :hybrid) || return Inf
    centers = channel_centers_mm(config.design)
    amplitudes = [topology_channel_geometry(y, config).amplitude_mm for y in centers]
    minimum(diff(centers) .- abs.(diff(amplitudes))) - config.design.guide_width_mm
end

function build_hybrid_monostatic_mesh(
    output_path::AbstractString;
    config::HybridMonostaticMeshConfig=HybridMonostaticMeshConfig(),
)
    validate(config)
    config.topology in (:pure_ttd, :hybrid) && minimum_channel_clearance_mm(config) <= 0 &&
        error("2D guide routes intersect; increase common_extra_path_mm")
    mkpath(dirname(output_path))
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("hybrid_monostatic_$(config.topology)_$(config.defect.id)")
        half_height = config.stand.specimen_height_mm * 0.5e-3
        specimen = gmsh.model.occ.addRectangle(
            0.0, -half_height, 0.0,
            config.stand.specimen_length_mm * 1e-3,
            config.stand.specimen_height_mm * 1e-3,
        )
        defect_surface = add_defect_surface(config.defect)
        if !isnothing(defect_surface)
            cut, _ = gmsh.model.occ.cut([(2, specimen)], [(2, defect_surface)])
            specimen_tags = [tag for (dimension, tag) in cut if dimension == 2]
        else
            specimen_tags = [specimen]
        end

        original_tags = copy(specimen_tags)
        original_keys = fill(0, length(specimen_tags))
        channel_summaries = NamedTuple[]
        if config.topology != :bare
            for center in channel_centers_mm(config.design)
                geometry = topology_channel_geometry(center, config)
                table = centreline_table(center, geometry, config)
                blocks = material_blocks(table, geometry, config)
                for block in blocks
                    push!(original_tags, add_strip_block(table, block, config))
                    push!(original_keys, material_key(block.fraction))
                end
                push!(channel_summaries, (
                    center_y_mm=center,
                    delay_fraction=geometry.delay_fraction,
                    amplitude_mm=geometry.amplitude_mm,
                    relative_extra_path_mm=geometry.relative_extra_path_mm,
                    total_extra_path_mm=geometry.total_extra_path_mm,
                    insert_length_mm=geometry.insert.total_length_mm,
                ))
            end
        end

        # Fragmenting every surface together creates a conforming channel/specimen interface.
        if length(original_tags) == 1
            final_tags = copy(original_tags)
            mapping = [[(2, only(original_tags))]]
        else
            output, mapping = gmsh.model.occ.fragment(
                [(2, first(original_tags))],
                [(2, tag) for tag in Iterators.drop(original_tags, 1)],
            )
            final_tags = unique(tag for (dimension, tag) in output if dimension == 2)
        end
        gmsh.model.occ.removeAllDuplicates()
        gmsh.model.occ.synchronize()
        final_tags = [tag for (dimension, tag) in gmsh.model.getEntities(2)]

        groups = Dict{Int, Vector{Int}}()
        for (key, mapped) in zip(original_keys, mapping)
            append!(get!(groups, key, Int[]), [
                tag for (dimension, tag) in mapped if dimension == 2 && tag in final_tags
            ])
        end
        for key in keys(groups)
            groups[key] = unique(groups[key])
        end
        isempty(get(groups, 0, Int[])) && error("Al domain classification failed")

        curves = collect(external_boundary_curves(final_tags))
        tolerance = 5e-7
        source_x = config.topology == :bare ? 0.0 : -config.design.device_depth_mm * 1e-3
        source = filter(curves) do tag
            center = gmsh.model.occ.getCenterOfMass(1, tag)
            abs(center[1] - source_x) < tolerance
        end
        radiation = filter(curves) do tag
            center = gmsh.model.occ.getCenterOfMass(1, tag)
            x, y = center[1], center[2]
            abs(x - config.stand.specimen_length_mm * 1e-3) < tolerance ||
                (x > -tolerance && abs(abs(y) - half_height) < tolerance)
        end
        defect_curves = filter(curves) do tag
            config.defect.kind == :none && return false
            center = gmsh.model.occ.getCenterOfMass(1, tag)
            hypot(
                center[1] - config.defect.center_x_mm * 1e-3,
                center[2] - config.defect.center_y_mm * 1e-3,
            ) < 1.1max(config.defect.radius_mm, config.defect.major_radius_mm) * 1e-3
        end
        free = setdiff(curves, vcat(source, radiation, defect_curves))
        isempty(source) && error("source aperture was not found")
        isempty(radiation) && error("radiation boundary was not found")
        config.defect.kind == :none || !isempty(defect_curves) || error("defect boundary missing")

        add_group(1, source, 101, "Source")
        add_group(1, radiation, 102, "RadiationBoundary")
        add_group(1, defect_curves, 103, "DefectBoundary")
        add_group(1, free, 104, "FreeSurface")
        for (offset, key) in enumerate(sort(collect(keys(groups))))
            add_group(2, groups[key], 200 + offset, material_name(key))
        end

        fine_curves = setdiff(curves, radiation)
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", fine_curves)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", config.guide_mesh_mm * 1e-3)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", config.specimen_mesh_mm * 1e-3)
        gmsh.model.mesh.field.setNumber(2, "DistMin", 0.8e-3)
        gmsh.model.mesh.field.setNumber(2, "DistMax", 5.0e-3)
        field_id = 2
        if !isempty(defect_curves)
            gmsh.model.mesh.field.add("Distance", 3)
            gmsh.model.mesh.field.setNumbers(3, "CurvesList", defect_curves)
            gmsh.model.mesh.field.add("Threshold", 4)
            gmsh.model.mesh.field.setNumber(4, "InField", 3)
            gmsh.model.mesh.field.setNumber(4, "SizeMin", config.defect_mesh_mm * 1e-3)
            gmsh.model.mesh.field.setNumber(4, "SizeMax", config.specimen_mesh_mm * 1e-3)
            gmsh.model.mesh.field.setNumber(4, "DistMin", 0.5e-3)
            gmsh.model.mesh.field.setNumber(4, "DistMax", 7.0e-3)
            gmsh.model.mesh.field.add("Min", 5)
            gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])
            field_id = 5
        end
        gmsh.model.mesh.field.setAsBackgroundMesh(field_id)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        element_types, element_tags, _ = gmsh.model.mesh.getElements(2)
        element_count = sum(length, element_tags)
        gmsh.write(output_path)
        source_length_mm = sum(source) do tag
            gmsh.model.occ.getMass(1, tag) * 1e3
        end
        (
            mesh_path=abspath(output_path),
            topology=config.topology,
            defect_id=config.defect.id,
            node_count=length(node_tags),
            element_count,
            material_keys=sort(collect(keys(groups))),
            source_curve_count=length(source),
            source_length_mm,
            pressure_scale=source_pressure_scale(source_length_mm, config.stand),
            minimum_clearance_mm=minimum_channel_clearance_mm(config),
            channel_summaries,
        )
    finally
        gmsh.finalize()
    end
end

end
