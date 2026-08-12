module TaperedShuntTrimMesher

using Gmsh: gmsh

export TaperedShuntTrimConfig,
       validate_config,
       total_length_mm,
       taper_half_height_mm,
       shunt_rectangles,
       build_tapered_shunt_trim_mesh

"""Smooth direct spine with a symmetric pair of shunt masses."""
Base.@kwdef struct TaperedShuntTrimConfig
    lead_length_mm::Float64 = 12.0
    port_height_mm::Float64 = 3.2
    active_length_mm::Float64 = 12.0
    taper_length_mm::Float64 = 3.0
    spine_height_mm::Float64 = 1.6
    mid_mass_length_mm::Float64 = 2.4
    max_mass_length_mm::Float64 = 4.0
    mass_height_mm::Float64 = 0.4
    neck_width_mm::Float64 = 0.4
    neck_height_mm::Float64 = 0.4
    profile_segments_per_taper::Int = 12
    minimum_feature_mm::Float64 = 0.35
end

total_length_mm(config::TaperedShuntTrimConfig) =
    2config.lead_length_mm + config.active_length_mm

function validate_config(config::TaperedShuntTrimConfig)
    config.lead_length_mm > 0 || throw(ArgumentError("lead length must be positive"))
    config.port_height_mm > 0 || throw(ArgumentError("port height must be positive"))
    0 < config.spine_height_mm < config.port_height_mm ||
        throw(ArgumentError("spine height must lie inside the port"))
    0 < 2config.taper_length_mm < config.active_length_mm ||
        throw(ArgumentError("tapers must leave a central constant-spine region"))
    config.mid_mass_length_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("mid mass is too short"))
    config.max_mass_length_mm >= config.mid_mass_length_mm ||
        throw(ArgumentError("max mass must not be shorter than mid mass"))
    config.mass_height_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("mass height violates the minimum feature"))
    config.neck_width_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("neck width violates the minimum feature"))
    config.neck_height_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("neck height violates the minimum feature"))
    config.profile_segments_per_taper >= 4 ||
        throw(ArgumentError("at least four taper segments are required"))
    available_half_height = (config.port_height_mm - config.spine_height_mm) / 2
    isapprox(
        config.mass_height_mm + config.neck_height_mm,
        available_half_height;
        atol=1.0e-12,
    ) || throw(ArgumentError("mass and neck must exactly fill the half-height above the spine"))
    config.max_mass_length_mm < config.active_length_mm - 2config.taper_length_mm ||
        throw(ArgumentError("largest shunt mass must fit inside the constant-spine region"))
    nothing
end

"""Half-height of the cosine-tapered direct path at local active coordinate x."""
function taper_half_height_mm(config::TaperedShuntTrimConfig, local_x_mm::Real)
    validate_config(config)
    0 <= local_x_mm <= config.active_length_mm ||
        throw(ArgumentError("active coordinate is outside the cell"))
    full = config.port_height_mm / 2
    spine = config.spine_height_mm / 2
    if local_x_mm <= config.taper_length_mm
        fraction = Float64(local_x_mm) / config.taper_length_mm
        return spine + (full - spine) * (1 + cospi(fraction)) / 2
    elseif local_x_mm >= config.active_length_mm - config.taper_length_mm
        fraction = (config.active_length_mm - Float64(local_x_mm)) / config.taper_length_mm
        return spine + (full - spine) * (1 + cospi(fraction)) / 2
    end
    spine
end

function profile_samples(config)
    left = collect(range(0.0, config.taper_length_mm; length=config.profile_segments_per_taper + 1))
    right = collect(range(
        config.active_length_mm - config.taper_length_mm,
        config.active_length_mm;
        length=config.profile_segments_per_taper + 1,
    ))
    sort(unique(vcat(left, config.active_length_mm / 2, right)))
end

function add_spine_surface(config)
    local_x = profile_samples(config)
    global_x = config.lead_length_mm .+ local_x
    half_height = taper_half_height_mm.(Ref(config), local_x)
    coordinates = vcat(
        collect(zip(global_x, half_height)),
        collect(zip(reverse(global_x), -reverse(half_height))),
    )
    points = [gmsh.model.occ.addPoint(x * 1e-3, y * 1e-3, 0.0) for (x, y) in coordinates]
    lines = [
        gmsh.model.occ.addLine(points[index], points[mod1(index + 1, length(points))])
        for index in eachindex(points)
    ]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function mass_length_mm(config, variant)
    variant == :mid && return config.mid_mass_length_mm
    variant == :max && return config.max_mass_length_mm
    throw(ArgumentError("shunt variant must be :mid or :max"))
end

"""Upper/lower masses and necks as `(x,y,width,height)` rectangles in mm."""
function shunt_rectangles(config::TaperedShuntTrimConfig, variant::Symbol)
    validate_config(config)
    length_mm = mass_length_mm(config, variant)
    center_x = config.lead_length_mm + config.active_length_mm / 2
    mass_left = center_x - length_mm / 2
    neck_left = center_x - config.neck_width_mm / 2
    spine_half = config.spine_height_mm / 2
    upper_mass_bottom = spine_half + config.neck_height_mm
    lower_mass_bottom = -config.port_height_mm / 2
    [
        (mass_left, upper_mass_bottom, length_mm, config.mass_height_mm),
        (mass_left, lower_mass_bottom, length_mm, config.mass_height_mm),
        (neck_left, spine_half, config.neck_width_mm, config.neck_height_mm),
        (neck_left, -spine_half - config.neck_height_mm, config.neck_width_mm, config.neck_height_mm),
    ]
end

function add_rectangle_mm(rectangle)
    x, y, width, height = rectangle
    gmsh.model.occ.addRectangle(x * 1e-3, y * 1e-3, 0.0, width * 1e-3, height * 1e-3)
end

function fused_surfaces(config, variant)
    if variant == :solid
        return [gmsh.model.occ.addRectangle(
            0.0,
            -config.port_height_mm * 0.5e-3,
            0.0,
            total_length_mm(config) * 1e-3,
            config.port_height_mm * 1e-3,
        )]
    end
    tags = [
        gmsh.model.occ.addRectangle(
            0.0,
            -config.port_height_mm * 0.5e-3,
            0.0,
            config.lead_length_mm * 1e-3,
            config.port_height_mm * 1e-3,
        ),
        add_spine_surface(config),
        gmsh.model.occ.addRectangle(
            (config.lead_length_mm + config.active_length_mm) * 1e-3,
            -config.port_height_mm * 0.5e-3,
            0.0,
            config.lead_length_mm * 1e-3,
            config.port_height_mm * 1e-3,
        ),
    ]
    append!(tags, add_rectangle_mm.(shunt_rectangles(config, variant)))
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(tags))],
        [(2, tag) for tag in Iterators.drop(tags, 1)],
    )
    surfaces = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surfaces) && error("tapered shunt fuse produced no domain")
    surfaces
end

function classify_boundaries(surfaces, config)
    gmsh.model.occ.synchronize()
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    curves = unique(tag for (dimension, tag) in boundary if dimension == 1)
    total_m = total_length_mm(config) * 1e-3
    tolerance = max(1.0e-10, total_m * 1.0e-8)
    source = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] < tolerance, curves)
    microphone = filter(
        tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] > total_m - tolerance,
        curves,
    )
    length(source) == 1 || error("tapered shunt cell must have exactly one left port")
    length(microphone) == 1 || error("tapered shunt cell must have exactly one right port")
    free = setdiff(curves, vcat(source, microphone))
    source, microphone, free
end

function tag_physical_groups(surfaces, source, microphone, free)
    for (id, tags, name) in (
        (101, source, "Source"),
        (102, microphone, "Microphone"),
        (103, free, "FreeSurface"),
    )
        gmsh.model.addPhysicalGroup(1, tags, id)
        gmsh.model.setPhysicalName(1, id, name)
    end
    gmsh.model.addPhysicalGroup(2, surfaces, 201)
    gmsh.model.setPhysicalName(2, 201, "Domain")
end

function configure_mesh(source, microphone, config; size_cell_mm, size_port_mm, size_lead_mm)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", vcat(source, microphone))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_port_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_lead_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.05e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.2e-3)

    lead_m = config.lead_length_mm * 1e-3
    active_m = config.active_length_mm * 1e-3
    height_m = config.port_height_mm * 1e-3
    gmsh.model.mesh.field.add("Box", 3)
    gmsh.model.mesh.field.setNumber(3, "VIn", Float64(size_cell_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "VOut", Float64(size_lead_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "XMin", lead_m - 0.5e-3)
    gmsh.model.mesh.field.setNumber(3, "XMax", lead_m + active_m + 0.5e-3)
    gmsh.model.mesh.field.setNumber(3, "YMin", -height_m)
    gmsh.model.mesh.field.setNumber(3, "YMax", height_m)

    gmsh.model.mesh.field.add("Min", 4)
    gmsh.model.mesh.field.setNumbers(4, "FieldsList", [2, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(4)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
end

function build_tapered_shunt_trim_mesh(
    output_path::AbstractString;
    config::TaperedShuntTrimConfig=TaperedShuntTrimConfig(),
    variant::Symbol=:max,
    size_cell_mm::Real=0.055,
    size_port_mm::Real=0.08,
    size_lead_mm::Real=0.25,
)
    validate_config(config)
    variant in (:solid, :mid, :max) || throw(ArgumentError("unknown tapered shunt variant"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("tapered_shunt_trim_$(variant)")
    surfaces = fused_surfaces(config, variant)
    source, microphone, free = classify_boundaries(surfaces, config)
    tag_physical_groups(surfaces, source, microphone, free)
    configure_mesh(source, microphone, config; size_cell_mm, size_port_mm, size_lead_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
