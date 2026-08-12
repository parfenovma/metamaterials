module HornNeighbor3DMesher

using Gmsh: gmsh

if !isdefined(parentmodule(@__MODULE__), :HornPointRadiatorMesher)
    include(joinpath(@__DIR__, "horn_point_radiator_mesher.jl"))
end
using ..HornPointRadiatorMesher

export HornNeighbor3DConfig,
       validate_config,
       guide_amplitude_mm,
       smooth_guide_amplitude_mm,
       guide_center_z_mm,
       outlet_x_mm,
       channel_centers_mm,
       probe_points_mm,
       build_horn_neighbor_3d_mesh

Base.@kwdef struct HornNeighbor3DConfig
    input_height_mm::Float64 = 3.2
    throat_height_mm::Float64 = 1.6
    channel_depth_mm::Float64 = 3.2
    channel_pitch_mm::Float64 = 4.8
    horn_length_mm::Float64 = 25.0
    device_guide_length_mm::Float64 = 26.0
    gentle_guide_axial_length_mm::Float64 = 26.0
    gentle_guide_path_length_mm::Float64 = 30.42
    receiver_length_mm::Float64 = 15.0
    receiver_half_width_mm::Float64 = 10.0
    receiver_half_height_mm::Float64 = 12.0
    profile_samples::Int = 96
end

function profile_config(config::HornNeighbor3DConfig)
    HornPointRadiatorConfig(
        input_height_mm=config.input_height_mm,
        throat_height_mm=config.throat_height_mm,
        horn_length_mm=config.horn_length_mm,
        straight_guide_length_mm=config.gentle_guide_path_length_mm,
        rounded_axial_length_mm=config.gentle_guide_axial_length_mm,
        receiver_length_mm=config.receiver_length_mm,
        receiver_half_height_mm=config.receiver_half_height_mm,
        profile_samples=config.profile_samples,
    )
end

guide_amplitude_mm(config::HornNeighbor3DConfig) =
    rounded_guide_amplitude_mm(profile_config(config))

function smooth_centerline_z_mm(config, local_x_mm, amplitude_mm)
    t = clamp(Float64(local_x_mm) / config.gentle_guide_axial_length_mm, 0.0, 1.0)
    amplitude_mm * sinpi(t)^4
end

function smooth_centerline_slope(config, local_x_mm, amplitude_mm)
    t = clamp(Float64(local_x_mm) / config.gentle_guide_axial_length_mm, 0.0, 1.0)
    4amplitude_mm * pi / config.gentle_guide_axial_length_mm * sinpi(t)^3 * cospi(t)
end

function smooth_guide_length_mm(config, amplitude_mm; samples=4001)
    x = range(0.0, config.gentle_guide_axial_length_mm; length=samples)
    integrand = hypot.(1.0, smooth_centerline_slope.(Ref(config), x, amplitude_mm))
    step = config.gentle_guide_axial_length_mm / (length(x) - 1)
    step * (sum(integrand) - (first(integrand) + last(integrand)) / 2)
end

function smooth_guide_amplitude_mm(config::HornNeighbor3DConfig)
    lower, upper = 0.0, config.gentle_guide_path_length_mm
    while smooth_guide_length_mm(config, upper) < config.gentle_guide_path_length_mm
        upper *= 2
    end
    for _ in 1:70
        middle = (lower + upper) / 2
        if smooth_guide_length_mm(config, middle) < config.gentle_guide_path_length_mm
            lower = middle
        else
            upper = middle
        end
    end
    (lower + upper) / 2
end

function guide_center_z_mm(config::HornNeighbor3DConfig, state::Symbol, local_x_mm::Real)
    state == :device && return 0.0
    state == :gentle &&
        return guide_center_y_mm(profile_config(config), :collector_gentle, local_x_mm)
    state == :smooth &&
        return smooth_centerline_z_mm(config, local_x_mm, smooth_guide_amplitude_mm(config))
    throw(ArgumentError("unknown channel state: $state"))
end

outlet_x_mm(config::HornNeighbor3DConfig) =
    config.horn_length_mm + config.device_guide_length_mm

function channel_centers_mm(config::HornNeighbor3DConfig, count::Integer)
    count in (1, 3) || throw(ArgumentError("the neighbour gate supports one or three channels"))
    count == 1 ? [0.0] : [-config.channel_pitch_mm, 0.0, config.channel_pitch_mm]
end

function validate_config(config::HornNeighbor3DConfig, states)
    length(states) in (1, 3) || throw(ArgumentError("use one or three channels"))
    all(state -> state in (:device, :gentle, :smooth), states) ||
        throw(ArgumentError("channel states must be :device, :gentle or :smooth"))
    config.input_height_mm > config.throat_height_mm > 0 ||
        throw(ArgumentError("the horn must narrow to a positive throat"))
    config.channel_depth_mm > 0 || throw(ArgumentError("channel depth must be positive"))
    config.channel_pitch_mm > config.channel_depth_mm ||
        throw(ArgumentError("neighbouring channels overlap"))
    config.device_guide_length_mm == config.gentle_guide_axial_length_mm ||
        throw(ArgumentError("device and gentle outlets must be coplanar"))
    config.gentle_guide_path_length_mm > config.gentle_guide_axial_length_mm ||
        throw(ArgumentError("gentle guide must add path length"))
    config.receiver_half_width_mm > config.channel_pitch_mm + config.channel_depth_mm / 2 ||
        throw(ArgumentError("receiver is too narrow for the three-channel gate"))
    config.receiver_half_height_mm >
        max(guide_amplitude_mm(config), smooth_guide_amplitude_mm(config)) + config.throat_height_mm ||
        throw(ArgumentError("receiver is too short in z for the gentle guide"))
    config.profile_samples >= 48 || throw(ArgumentError("profiles need at least 48 samples"))
    nothing
end

function add_polygon_surface_xz_mm(points_mm, y_mm)
    clean = NTuple{2, Float64}[]
    for (x, z) in points_mm
        point = (Float64(x), Float64(z))
        (isempty(clean) || point != last(clean)) && push!(clean, point)
    end
    first(clean) == last(clean) && pop!(clean)
    points = [
        gmsh.model.occ.addPoint(x * 1e-3, y_mm * 1e-3, z * 1e-3)
        for (x, z) in clean
    ]
    lines = [
        gmsh.model.occ.addLine(points[index], points[mod1(index + 1, length(points))])
        for index in eachindex(points)
    ]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function add_horn_surface(config, y_min_mm)
    profile = profile_config(config)
    x = range(0.0, config.horn_length_mm; length=config.profile_samples)
    half_height = horn_height_mm.(Ref(profile), x) ./ 2
    points = vcat(collect(zip(x, half_height)), collect(zip(reverse(x), -reverse(half_height))))
    add_polygon_surface_xz_mm(points, y_min_mm)
end

function add_device_surface(config, y_min_mm)
    x0 = config.horn_length_mm
    half_height = config.throat_height_mm / 2
    add_polygon_surface_xz_mm(
        [(x0, -half_height), (outlet_x_mm(config), -half_height),
         (outlet_x_mm(config), half_height), (x0, half_height)],
        y_min_mm,
    )
end

function gentle_centerline_slope(config, local_x_mm)
    amplitude = guide_amplitude_mm(config)
    amplitude * pi / config.gentle_guide_axial_length_mm *
    sinpi(2Float64(local_x_mm) / config.gentle_guide_axial_length_mm)
end

function guide_centerline_slope(config, state, local_x_mm)
    state == :gentle && return gentle_centerline_slope(config, local_x_mm)
    state == :smooth &&
        return smooth_centerline_slope(config, local_x_mm, smooth_guide_amplitude_mm(config))
    state == :device && return 0.0
    throw(ArgumentError("unknown channel state: $state"))
end

function add_gentle_surface(config, y_min_mm)
    half_width = config.throat_height_mm / 2
    lower = NTuple{2, Float64}[]
    upper = NTuple{2, Float64}[]
    for local_x in range(0.0, config.gentle_guide_axial_length_mm; length=config.profile_samples)
        center_z = guide_center_z_mm(config, :gentle, local_x)
        slope = gentle_centerline_slope(config, local_x)
        normalization = hypot(1.0, slope)
        normal_x, normal_z = -slope / normalization, 1 / normalization
        global_x = config.horn_length_mm + local_x
        push!(lower, (global_x - half_width * normal_x, center_z - half_width * normal_z))
        push!(upper, (global_x + half_width * normal_x, center_z + half_width * normal_z))
    end
    add_polygon_surface_xz_mm(vcat(lower, reverse(upper)), y_min_mm)
end

function channel_polygon_mm(config, state)
    profile = profile_config(config)
    horn_x = collect(range(0.0, config.horn_length_mm; length=config.profile_samples))
    horn_half = horn_height_mm.(Ref(profile), horn_x) ./ 2
    horn_upper = collect(zip(horn_x, horn_half))
    horn_lower = collect(zip(horn_x, -horn_half))
    if state == :device
        half_width = config.throat_height_mm / 2
        guide_upper = [(config.horn_length_mm, half_width), (outlet_x_mm(config), half_width)]
        guide_lower = [(config.horn_length_mm, -half_width), (outlet_x_mm(config), -half_width)]
    else
        half_width = config.throat_height_mm / 2
        guide_lower = NTuple{2, Float64}[]
        guide_upper = NTuple{2, Float64}[]
        for local_x in range(0.0, config.gentle_guide_axial_length_mm; length=config.profile_samples)
            center_z = guide_center_z_mm(config, state, local_x)
            slope = guide_centerline_slope(config, state, local_x)
            normalization = hypot(1.0, slope)
            normal_x, normal_z = -slope / normalization, 1 / normalization
            global_x = config.horn_length_mm + local_x
            push!(guide_lower, (global_x - half_width * normal_x, center_z - half_width * normal_z))
            push!(guide_upper, (global_x + half_width * normal_x, center_z + half_width * normal_z))
        end
    end
    vcat(horn_upper, guide_upper, reverse(guide_lower), reverse(horn_lower))
end

function add_channel_volume(config, state, center_y_mm)
    y_min = center_y_mm - config.channel_depth_mm / 2
    surface = add_polygon_surface_xz_mm(channel_polygon_mm(config, state), y_min)
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

function build_domain(config, states)
    centers = channel_centers_mm(config, length(states))
    channels = [add_channel_volume(config, state, center) for (state, center) in zip(states, centers)]
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
    isempty(volumes) && error("3D horn-neighbour fuse produced no volume")
    volumes, centers
end

function external_surfaces()
    gmsh.model.occ.synchronize()
    [
        tag for (dimension, tag) in gmsh.model.getEntities(2)
        if dimension == 2 && length(first(gmsh.model.getAdjacencies(2, tag))) == 1
    ]
end

function classify_boundaries(config, centers)
    surfaces = external_surfaces()
    tolerance = 1.0e-8
    outlet = outlet_x_mm(config) * 1e-3
    end_x = outlet + config.receiver_length_mm * 1e-3
    half_y = config.receiver_half_width_mm * 1e-3
    half_z = config.receiver_half_height_mm * 1e-3
    sources = Vector{Vector{Int}}()
    for center_y in centers
        source = filter(surfaces) do tag
            center = gmsh.model.occ.getCenterOfMass(2, tag)
            abs(center[1]) < tolerance &&
                abs(center[2] - center_y * 1e-3) < config.channel_depth_mm * 0.26e-3
        end
        length(source) == 1 || error("channel at y=$center_y mm must have one source face")
        push!(sources, source)
    end
    all_sources = reduce(vcat, sources)
    radiation = filter(surfaces) do tag
        center = gmsh.model.occ.getCenterOfMass(2, tag)
        abs(center[1] - end_x) < tolerance ||
        (center[1] > outlet + tolerance && abs(abs(center[2]) - half_y) < tolerance) ||
        (center[1] > outlet + tolerance && abs(abs(center[3]) - half_z) < tolerance)
    end
    isempty(radiation) && error("receiver radiation faces were not found")
    sources, radiation, setdiff(surfaces, vcat(all_sources, radiation))
end

function add_physical_group(dim, tags, id, name)
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function configure_mesh(config; size_path_mm, size_receiver_mm)
    outlet = outlet_x_mm(config)
    amplitude = max(guide_amplitude_mm(config), smooth_guide_amplitude_mm(config))
    gmsh.model.mesh.field.add("Box", 1)
    gmsh.model.mesh.field.setNumber(1, "VIn", Float64(size_path_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "VOut", Float64(size_receiver_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "XMin", -0.2e-3)
    gmsh.model.mesh.field.setNumber(1, "XMax", (outlet + 2.0) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "YMin", -8.0e-3)
    gmsh.model.mesh.field.setNumber(1, "YMax", 8.0e-3)
    gmsh.model.mesh.field.setNumber(1, "ZMin", -(amplitude + 2.0) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "ZMax", (amplitude + 2.0) * 1e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
end

function probe_points_mm(config::HornNeighbor3DConfig, states)
    centers = channel_centers_mm(config, length(states))
    names = String[]
    x = Float64[]
    y = Float64[]
    z = Float64[]
    for (index, (state, center_y)) in enumerate(zip(states, centers))
        local_x = config.device_guide_length_mm - 1.0
        push!(names, "preout_$index")
        push!(x, config.horn_length_mm + local_x)
        push!(y, center_y)
        push!(z, guide_center_z_mm(config, state, local_x))
        push!(names, "receiver_$index")
        push!(x, outlet_x_mm(config) + 2.0)
        push!(y, center_y)
        push!(z, 0.0)
    end
    push!(names, "receiver_axis_8mm")
    push!(x, outlet_x_mm(config) + 8.0)
    push!(y, 0.0)
    push!(z, 0.0)
    (names=names, x_mm=x, y_mm=y, z_mm=z)
end

function build_horn_neighbor_3d_mesh(
    output_path::AbstractString;
    config::HornNeighbor3DConfig=HornNeighbor3DConfig(),
    states::AbstractVector{Symbol}=[:device, :gentle, :device],
    size_path_mm::Real=0.65,
    size_receiver_mm::Real=1.40,
)
    validate_config(config, states)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("horn_neighbor_3d_$(join(String.(states), '_'))")
    volumes, centers = build_domain(config, states)
    sources, radiation, free = classify_boundaries(config, centers)
    for (index, source) in enumerate(sources)
        add_physical_group(2, source, 100 + index, "Source_$index")
    end
    add_physical_group(2, radiation, 120, "RadiationBoundary")
    add_physical_group(2, free, 121, "FreeSurface")
    add_physical_group(3, volumes, 201, "Domain")
    configure_mesh(config; size_path_mm, size_receiver_mm)
    gmsh.model.mesh.generate(3)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
