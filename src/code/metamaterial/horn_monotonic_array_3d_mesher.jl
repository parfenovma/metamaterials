module HornMonotonicArray3DMesher

using Gmsh: gmsh

if !isdefined(parentmodule(@__MODULE__), :HornPointRadiatorMesher)
    include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
end
if !isdefined(parentmodule(@__MODULE__), :MonotonicHornLens)
    include(joinpath(@__DIR__, "monotonic_horn_lens.jl"))
end
using ..HornPointRadiatorMesher: HornPointRadiatorConfig, horn_height_mm
using ..MonotonicHornLens: smooth_centerline_mm, smooth_slope

export HornMonotonicArray3DConfig,
       validate_config,
       channel_centers_mm,
       outlet_x_mm,
       focus_x_mm,
       channel_center_z_mm,
       channel_polygon_mm,
       probe_points_mm,
       build_horn_monotonic_array_3d_mesh,
       build_horn_monotonic_array_3d_half_mesh

Base.@kwdef struct HornMonotonicArray3DConfig
    input_height_mm::Float64 = 3.2
    throat_height_mm::Float64 = 1.6
    channel_depth_mm::Float64 = 3.2
    channel_pitch_mm::Float64 = 4.8
    horn_length_mm::Float64 = 25.0
    guide_axial_length_mm::Float64 = 50.663944208799734
    bend_amplitude_mm::Vector{Float64} = zeros(15)
    diffuser_output_height_mm::Float64 = 1.6
    diffuser_length_mm::Float64 = 0.0
    receiver_length_mm::Float64 = 55.0
    receiver_half_width_mm::Float64 = 40.0
    receiver_half_height_mm::Float64 = 12.0
    focal_distance_mm::Float64 = 35.0
    profile_samples::Int = 72
end

function channel_centers_mm(config::HornMonotonicArray3DConfig)
    count = length(config.bend_amplitude_mm)
    isodd(count) || throw(ArgumentError("the monotonic array needs an odd channel count"))
    half = (count - 1) ÷ 2
    collect((-half):half) .* config.channel_pitch_mm
end

guide_outlet_x_mm(config::HornMonotonicArray3DConfig) =
    config.horn_length_mm + config.guide_axial_length_mm

outlet_x_mm(config::HornMonotonicArray3DConfig) =
    guide_outlet_x_mm(config) + config.diffuser_length_mm

focus_x_mm(config::HornMonotonicArray3DConfig) =
    outlet_x_mm(config) + config.focal_distance_mm

function validate_config(config::HornMonotonicArray3DConfig)
    count = length(config.bend_amplitude_mm)
    count >= 3 && isodd(count) || throw(ArgumentError("use an odd array with at least 3 channels"))
    all(>=(0.0), config.bend_amplitude_mm) ||
        throw(ArgumentError("bend amplitudes must be non-negative"))
    config.bend_amplitude_mm == reverse(config.bend_amplitude_mm) ||
        throw(ArgumentError("the v0 lens must be mirror symmetric"))
    config.input_height_mm > config.throat_height_mm > 0 ||
        throw(ArgumentError("the horn must narrow to a positive throat"))
    config.channel_pitch_mm > config.channel_depth_mm > 0 ||
        throw(ArgumentError("neighbouring channel prisms overlap"))
    config.horn_length_mm > 0 || throw(ArgumentError("horn length must be positive"))
    config.guide_axial_length_mm > 0 ||
        throw(ArgumentError("guide axial length must be positive"))
    config.diffuser_length_mm >= 0 ||
        throw(ArgumentError("diffuser length must be non-negative"))
    if config.diffuser_length_mm > 0
        config.diffuser_output_height_mm >= config.throat_height_mm ||
            throw(ArgumentError("output diffuser must not contract below the throat"))
    end
    config.receiver_length_mm > config.focal_distance_mm > 0 ||
        throw(ArgumentError("focus must lie inside the receiver"))
    maximum(abs, channel_centers_mm(config)) + config.channel_depth_mm / 2 <
        config.receiver_half_width_mm || throw(ArgumentError("array does not fit the receiver width"))
    config.receiver_half_height_mm > config.throat_height_mm ||
        throw(ArgumentError("receiver is too thin at the output plane"))
    config.profile_samples >= 48 || throw(ArgumentError("profiles need at least 48 samples"))
    nothing
end

function horn_profile_config(config)
    HornPointRadiatorConfig(
        input_height_mm=config.input_height_mm,
        throat_height_mm=config.throat_height_mm,
        horn_length_mm=config.horn_length_mm,
        straight_guide_length_mm=config.guide_axial_length_mm,
        rounded_axial_length_mm=config.guide_axial_length_mm - 0.1,
        receiver_length_mm=config.receiver_length_mm,
        receiver_half_height_mm=config.receiver_half_height_mm,
        profile_samples=config.profile_samples,
    )
end

channel_center_z_mm(config, channel_index, local_x_mm) = smooth_centerline_mm(
    local_x_mm,
    config.guide_axial_length_mm,
    config.bend_amplitude_mm[channel_index],
)

function channel_polygon_mm(config, channel_index)
    profile = horn_profile_config(config)
    horn_x_mm = collect(range(0.0, config.horn_length_mm; length=config.profile_samples))
    horn_half_mm = horn_height_mm.(Ref(profile), horn_x_mm) ./ 2
    horn_upper = collect(zip(horn_x_mm, horn_half_mm))
    horn_lower = collect(zip(horn_x_mm, -horn_half_mm))
    amplitude_mm = config.bend_amplitude_mm[channel_index]
    half_width_mm = config.throat_height_mm / 2
    guide_lower = NTuple{2, Float64}[]
    guide_upper = NTuple{2, Float64}[]
    for local_x_mm in range(0.0, config.guide_axial_length_mm; length=config.profile_samples)
        center_z_mm = smooth_centerline_mm(
            local_x_mm, config.guide_axial_length_mm, amplitude_mm,
        )
        slope = smooth_slope(local_x_mm, config.guide_axial_length_mm, amplitude_mm)
        normalization = hypot(1.0, slope)
        normal_x, normal_z = -slope / normalization, 1 / normalization
        global_x_mm = config.horn_length_mm + local_x_mm
        push!(guide_lower, (
            global_x_mm - half_width_mm * normal_x,
            center_z_mm - half_width_mm * normal_z,
        ))
        push!(guide_upper, (
            global_x_mm + half_width_mm * normal_x,
            center_z_mm + half_width_mm * normal_z,
        ))
    end
    if config.diffuser_length_mm > 0
        local_x_mm = collect(range(
            0.0, config.diffuser_length_mm; length=config.profile_samples,
        ))
        q = (1 .- cospi.(local_x_mm ./ config.diffuser_length_mm)) ./ 2
        height_mm = config.throat_height_mm .* exp.(
            q .* log(config.diffuser_output_height_mm / config.throat_height_mm),
        )
        diffuser_x_mm = guide_outlet_x_mm(config) .+ local_x_mm
        append!(guide_upper, collect(zip(diffuser_x_mm[2:end], height_mm[2:end] ./ 2)))
        append!(guide_lower, collect(zip(diffuser_x_mm[2:end], -height_mm[2:end] ./ 2)))
    end
    vcat(horn_upper, guide_upper, reverse(guide_lower), reverse(horn_lower))
end

function add_polygon_surface_xz_mm(points_mm, y_mm)
    clean = NTuple{2, Float64}[]
    for (x_mm, z_mm) in points_mm
        point = (Float64(x_mm), Float64(z_mm))
        (isempty(clean) || point != last(clean)) && push!(clean, point)
    end
    first(clean) == last(clean) && pop!(clean)
    points = [
        gmsh.model.occ.addPoint(x_mm * 1e-3, y_mm * 1e-3, z_mm * 1e-3)
        for (x_mm, z_mm) in clean
    ]
    lines = [
        gmsh.model.occ.addLine(points[index], points[mod1(index + 1, length(points))])
        for index in eachindex(points)
    ]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function add_channel_volume(config, channel_index, center_y_mm)
    y_min_mm = center_y_mm - config.channel_depth_mm / 2
    surface = add_polygon_surface_xz_mm(
        channel_polygon_mm(config, channel_index),
        y_min_mm,
    )
    extruded = gmsh.model.occ.extrude(
        [(2, surface)],
        0.0,
        config.channel_depth_mm * 1e-3,
        0.0,
    )
    volumes = unique(tag for (dimension, tag) in extruded if dimension == 3)
    length(volumes) == 1 || error("channel extrusion did not produce one volume")
    only(volumes)
end

function build_domain(config)
    centers_mm = channel_centers_mm(config)
    channels = [
        add_channel_volume(config, index, center_y_mm)
        for (index, center_y_mm) in enumerate(centers_mm)
    ]
    receiver = gmsh.model.occ.addBox(
        outlet_x_mm(config) * 1e-3,
        -config.receiver_half_width_mm * 1e-3,
        -config.receiver_half_height_mm * 1e-3,
        config.receiver_length_mm * 1e-3,
        2config.receiver_half_width_mm * 1e-3,
        2config.receiver_half_height_mm * 1e-3,
    )
    fused, _ = gmsh.model.occ.fuse([(3, receiver)], [(3, tag) for tag in channels])
    volumes = unique(tag for (dimension, tag) in fused if dimension == 3)
    isempty(volumes) && error("monotonic-array fuse produced no volume")
    volumes, centers_mm
end

function build_half_domain(config)
    volumes, centers_mm = build_domain(config)
    maximum_bend_mm = maximum(config.bend_amplitude_mm)
    x_min_mm = -1.0
    z_min_mm = -config.receiver_half_height_mm - 2.0
    half_box = gmsh.model.occ.addBox(
        x_min_mm * 1e-3,
        0.0,
        z_min_mm * 1e-3,
        (outlet_x_mm(config) + config.receiver_length_mm - x_min_mm + 2.0) * 1e-3,
        (config.receiver_half_width_mm + 1.0) * 1e-3,
        (2config.receiver_half_height_mm + maximum_bend_mm + 4.0) * 1e-3,
    )
    clipped, _ = gmsh.model.occ.intersect(
        [(3, tag) for tag in volumes],
        [(3, half_box)],
    )
    half_volumes = unique(tag for (dimension, tag) in clipped if dimension == 3)
    isempty(half_volumes) && error("half-domain intersection produced no volume")
    half_centers_mm = centers_mm[centers_mm .>= 0.0]
    half_volumes, half_centers_mm
end

function external_surfaces()
    gmsh.model.occ.synchronize()
    [
        tag for (dimension, tag) in gmsh.model.getEntities(2)
        if dimension == 2 && length(first(gmsh.model.getAdjacencies(2, tag))) == 1
    ]
end

function classify_boundaries(config, centers_mm)
    surfaces = external_surfaces()
    tolerance = 1.0e-8
    outlet_m = outlet_x_mm(config) * 1e-3
    end_x_m = (outlet_x_mm(config) + config.receiver_length_mm) * 1e-3
    half_y_m = config.receiver_half_width_mm * 1e-3
    half_z_m = config.receiver_half_height_mm * 1e-3
    sources = Vector{Vector{Int}}()
    for center_y_mm in centers_mm
        source = filter(surfaces) do tag
            center = gmsh.model.occ.getCenterOfMass(2, tag)
            abs(center[1]) < tolerance &&
                abs(center[2] - center_y_mm * 1e-3) < config.channel_depth_mm * 0.26e-3
        end
        length(source) == 1 || error("channel at y=$center_y_mm mm must have one source face")
        push!(sources, source)
    end
    all_sources = reduce(vcat, sources)
    radiation = filter(surfaces) do tag
        center = gmsh.model.occ.getCenterOfMass(2, tag)
        abs(center[1] - end_x_m) < tolerance ||
        (center[1] > outlet_m + tolerance && abs(abs(center[2]) - half_y_m) < tolerance) ||
        (center[1] > outlet_m + tolerance && abs(abs(center[3]) - half_z_m) < tolerance)
    end
    isempty(radiation) && error("receiver radiation faces were not found")
    sources, radiation, setdiff(surfaces, vcat(all_sources, radiation))
end

function classify_half_boundaries(config, centers_mm)
    surfaces = external_surfaces()
    tolerance = 1.0e-8
    outlet_m = outlet_x_mm(config) * 1e-3
    end_x_m = (outlet_x_mm(config) + config.receiver_length_mm) * 1e-3
    half_y_m = config.receiver_half_width_mm * 1e-3
    half_z_m = config.receiver_half_height_mm * 1e-3
    source_tags = filter(surfaces) do tag
        center = gmsh.model.occ.getCenterOfMass(2, tag)
        abs(center[1]) < tolerance && center[2] >= -tolerance
    end
    sort!(source_tags; by=tag -> gmsh.model.occ.getCenterOfMass(2, tag)[2])
    length(source_tags) == length(centers_mm) || error(
        "half-domain needs $(length(centers_mm)) source faces, found $(length(source_tags))",
    )
    sources = [[tag] for tag in source_tags]
    symmetry = filter(surfaces) do tag
        center = gmsh.model.occ.getCenterOfMass(2, tag)
        abs(center[2]) < tolerance
    end
    isempty(symmetry) && error("half-domain symmetry faces were not found")
    radiation = filter(surfaces) do tag
        center = gmsh.model.occ.getCenterOfMass(2, tag)
        abs(center[1] - end_x_m) < tolerance ||
        (center[1] > outlet_m + tolerance && abs(center[2] - half_y_m) < tolerance) ||
        (center[1] > outlet_m + tolerance && abs(abs(center[3]) - half_z_m) < tolerance)
    end
    isempty(radiation) && error("half-domain receiver radiation faces were not found")
    excluded = vcat(reduce(vcat, sources), radiation, symmetry)
    sources, radiation, symmetry, setdiff(surfaces, excluded)
end

function add_physical_group(dim, tags, id, name)
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function configure_mesh(
    config;
    size_path_mm,
    size_focus_mm,
    size_receiver_mm,
)
    outlet_mm = outlet_x_mm(config)
    maximum_bend_mm = maximum(config.bend_amplitude_mm)
    aperture_edge_mm = maximum(abs, channel_centers_mm(config)) + config.channel_depth_mm

    gmsh.model.mesh.field.add("Box", 1)
    gmsh.model.mesh.field.setNumber(1, "VIn", Float64(size_path_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "VOut", Float64(size_receiver_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "XMin", -0.2e-3)
    gmsh.model.mesh.field.setNumber(1, "XMax", (outlet_mm + 1.5) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "YMin", -aperture_edge_mm * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "YMax", aperture_edge_mm * 1e-3)
    gmsh.model.mesh.field.setNumber(
        1, "ZMin", -(config.input_height_mm / 2 + 2.0) * 1e-3,
    )
    gmsh.model.mesh.field.setNumber(1, "ZMax", (maximum_bend_mm + 2.0) * 1e-3)

    gmsh.model.mesh.field.add("Box", 2)
    gmsh.model.mesh.field.setNumber(2, "VIn", Float64(size_focus_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "VOut", Float64(size_receiver_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "XMin", (focus_x_mm(config) - 13.0) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "XMax", (focus_x_mm(config) + 13.0) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "YMin", -18.0e-3)
    gmsh.model.mesh.field.setNumber(2, "YMax", 18.0e-3)
    gmsh.model.mesh.field.setNumber(2, "ZMin", -8.0e-3)
    gmsh.model.mesh.field.setNumber(2, "ZMax", 8.0e-3)

    gmsh.model.mesh.field.add("Min", 3)
    gmsh.model.mesh.field.setNumbers(3, "FieldsList", [1, 2])
    gmsh.model.mesh.field.setAsBackgroundMesh(3)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
end

function probe_points_mm(config::HornMonotonicArray3DConfig)
    centers_mm = channel_centers_mm(config)
    names = String[]
    x_mm = Float64[]
    y_mm = Float64[]
    z_mm = Float64[]
    local_x_mm = config.guide_axial_length_mm - 1.0
    for (index, center_y_mm) in enumerate(centers_mm)
        push!(names, "preout_$index")
        push!(x_mm, config.horn_length_mm + local_x_mm)
        push!(y_mm, center_y_mm)
        push!(z_mm, channel_center_z_mm(config, index, local_x_mm))
        push!(names, "receiver_$index")
        push!(x_mm, outlet_x_mm(config) + 2.0)
        push!(y_mm, center_y_mm)
        push!(z_mm, 0.0)
    end
    push!(names, "focus")
    push!(x_mm, focus_x_mm(config))
    push!(y_mm, 0.0)
    push!(z_mm, 0.0)
    (; names, x_mm, y_mm, z_mm)
end

function build_horn_monotonic_array_3d_mesh(
    output_path::AbstractString;
    config::HornMonotonicArray3DConfig,
    size_path_mm::Real=0.90,
    size_focus_mm::Real=1.80,
    size_receiver_mm::Real=2.60,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("horn_monotonic_array_3d")
    volumes, centers_mm = build_domain(config)
    sources, radiation, free = classify_boundaries(config, centers_mm)
    for (index, source) in enumerate(sources)
        add_physical_group(2, source, 100 + index, "Source_$index")
    end
    add_physical_group(2, radiation, 130, "RadiationBoundary")
    add_physical_group(2, free, 131, "FreeSurface")
    add_physical_group(3, volumes, 201, "Domain")
    configure_mesh(config; size_path_mm, size_focus_mm, size_receiver_mm)
    gmsh.model.mesh.generate(3)
    node_count = length(first(gmsh.model.mesh.getNodes()))
    element_count = sum(length, gmsh.model.mesh.getElements()[2])
    gmsh.write(output_path)
    println("[+] monotonic 3D mesh: $node_count nodes, $element_count elements")
    println("[+] $output_path")
    (; output_path, node_count, element_count)
end

function build_horn_monotonic_array_3d_half_mesh(
    output_path::AbstractString;
    config::HornMonotonicArray3DConfig,
    size_path_mm::Real=0.90,
    size_focus_mm::Real=1.80,
    size_receiver_mm::Real=2.60,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("horn_monotonic_array_3d_half")
    volumes, centers_mm = build_half_domain(config)
    sources, radiation, symmetry, free = classify_half_boundaries(config, centers_mm)
    for (index, source) in enumerate(sources)
        add_physical_group(2, source, 100 + index, "Source_$index")
    end
    add_physical_group(2, radiation, 130, "RadiationBoundary")
    add_physical_group(2, free, 131, "FreeSurface")
    add_physical_group(2, symmetry, 132, "SymmetryBoundary")
    add_physical_group(3, volumes, 201, "Domain")
    configure_mesh(config; size_path_mm, size_focus_mm, size_receiver_mm)
    gmsh.model.mesh.generate(3)
    node_count = length(first(gmsh.model.mesh.getNodes()))
    element_count = sum(length, gmsh.model.mesh.getElements()[2])
    gmsh.write(output_path)
    println("[+] monotonic half-domain 3D mesh: $node_count nodes, $element_count elements")
    println("[+] $output_path")
    (; output_path, node_count, element_count, source_count=length(sources))
end

end
