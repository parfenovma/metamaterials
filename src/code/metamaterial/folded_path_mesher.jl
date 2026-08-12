module FoldedPathMesher

using Gmsh: gmsh

export FoldedPathConfig,
       FOLDED_PATH_VARIANTS,
       validate_config,
       unfolded_active_length_mm,
       total_axial_length_mm,
       sharp_fold_amplitude_mm,
       rounded_fold_amplitude_mm,
       rounded_centerline_mm,
       rounded_centerline_length_mm,
       minimum_rounded_inner_radius_mm,
       build_folded_path_mesh

const FOLDED_PATH_VARIANTS = (
    :straight_device,
    :straight_unfolded,
    :sharp_fold,
    :rounded_fold,
)

"""Binary elastic-fold pilot with equal unfolded path lengths."""
Base.@kwdef struct FoldedPathConfig
    lead_length_mm::Float64 = 12.0
    active_axial_length_mm::Float64 = 20.0
    extra_path_length_mm::Float64 = 10.42
    guide_width_mm::Float64 = 3.2
    rounded_samples::Int = 180
    arc_integration_samples::Int = 4001
    minimum_inner_radius_mm::Float64 = 0.20
end

unfolded_active_length_mm(config::FoldedPathConfig) =
    config.active_axial_length_mm + config.extra_path_length_mm

function total_axial_length_mm(config::FoldedPathConfig, variant::Symbol)
    variant in FOLDED_PATH_VARIANTS || throw(ArgumentError("unknown folded-path variant"))
    active = variant == :straight_unfolded ?
             unfolded_active_length_mm(config) : config.active_axial_length_mm
    2config.lead_length_mm + active
end

function sharp_fold_amplitude_mm(config::FoldedPathConfig)
    half_path = unfolded_active_length_mm(config) / 2
    half_axial = config.active_axial_length_mm / 2
    sqrt(half_path^2 - half_axial^2)
end

rounded_centerline_y_mm(config::FoldedPathConfig, x_mm::Real, amplitude_mm::Real) =
    Float64(amplitude_mm) / 2 * (1 - cospi(2Float64(x_mm) / config.active_axial_length_mm))

rounded_centerline_slope(config::FoldedPathConfig, x_mm::Real, amplitude_mm::Real) =
    Float64(amplitude_mm) * pi / config.active_axial_length_mm *
    sinpi(2Float64(x_mm) / config.active_axial_length_mm)

function rounded_centerline_length_mm(
    config::FoldedPathConfig,
    amplitude_mm::Real;
    samples::Integer=config.arc_integration_samples,
)
    samples >= 3 || throw(ArgumentError("at least three arc-integration samples are required"))
    x = range(0.0, config.active_axial_length_mm; length=Int(samples))
    integrand = [
        hypot(1.0, rounded_centerline_slope(config, coordinate, amplitude_mm))
        for coordinate in x
    ]
    step = config.active_axial_length_mm / (length(x) - 1)
    step * (sum(integrand) - (first(integrand) + last(integrand)) / 2)
end

function rounded_fold_amplitude_mm(config::FoldedPathConfig)
    target = unfolded_active_length_mm(config)
    lower = 0.0
    upper = max(config.active_axial_length_mm, config.extra_path_length_mm)
    while rounded_centerline_length_mm(config, upper) < target
        upper *= 2
    end
    for _ in 1:70
        middle = (lower + upper) / 2
        if rounded_centerline_length_mm(config, middle) < target
            lower = middle
        else
            upper = middle
        end
    end
    (lower + upper) / 2
end

function rounded_centerline_mm(
    config::FoldedPathConfig;
    samples::Integer=config.rounded_samples,
)
    samples >= 8 || throw(ArgumentError("at least eight rounded-path samples are required"))
    amplitude = rounded_fold_amplitude_mm(config)
    [
        (
            x,
            rounded_centerline_y_mm(config, x, amplitude),
            rounded_centerline_slope(config, x, amplitude),
        )
        for x in range(0.0, config.active_axial_length_mm; length=Int(samples))
    ]
end

function minimum_rounded_inner_radius_mm(config::FoldedPathConfig)
    amplitude = rounded_fold_amplitude_mm(config)
    maximum_curvature = 2pi^2 * amplitude / config.active_axial_length_mm^2
    inv(maximum_curvature) - config.guide_width_mm / 2
end

function validate_config(config::FoldedPathConfig)
    config.lead_length_mm > 0 || throw(ArgumentError("lead length must be positive"))
    config.active_axial_length_mm > 0 || throw(ArgumentError("active axial length must be positive"))
    config.extra_path_length_mm > 0 || throw(ArgumentError("extra path length must be positive"))
    config.guide_width_mm > 0 || throw(ArgumentError("guide width must be positive"))
    config.rounded_samples >= 32 || throw(ArgumentError("rounded path needs at least 32 samples"))
    config.arc_integration_samples >= 101 ||
        throw(ArgumentError("arc integration needs at least 101 samples"))
    minimum_rounded_inner_radius_mm(config) + 1.0e-12 >= config.minimum_inner_radius_mm ||
        throw(ArgumentError("rounded fold violates the minimum inner radius"))
    nothing
end

function add_polygon_surface_mm(points_mm)
    points = [gmsh.model.occ.addPoint(x * 1e-3, y * 1e-3, 0.0) for (x, y) in points_mm]
    lines = [
        gmsh.model.occ.addLine(points[index], points[mod1(index + 1, length(points))])
        for index in eachindex(points)
    ]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function add_rectangle_mm(x, y, width, height)
    gmsh.model.occ.addRectangle(x * 1e-3, y * 1e-3, 0.0, width * 1e-3, height * 1e-3)
end

function add_segment_strip_mm(first_point, second_point, width_mm)
    x0, y0 = first_point
    x1, y1 = second_point
    dx = x1 - x0
    dy = y1 - y0
    length_mm = hypot(dx, dy)
    normal = (-dy / length_mm, dx / length_mm)
    offset = (normal[1] * width_mm / 2, normal[2] * width_mm / 2)
    add_polygon_surface_mm([
        (x0 + offset[1], y0 + offset[2]),
        (x1 + offset[1], y1 + offset[2]),
        (x1 - offset[1], y1 - offset[2]),
        (x0 - offset[1], y0 - offset[2]),
    ])
end

function add_rounded_active_surface(config)
    half_width = config.guide_width_mm / 2
    centerline = rounded_centerline_mm(config)
    lower = NTuple{2, Float64}[]
    upper = NTuple{2, Float64}[]
    for (local_x, center_y, slope) in centerline
        normalization = hypot(1.0, slope)
        normal_x = -slope / normalization
        normal_y = 1 / normalization
        global_x = config.lead_length_mm + local_x
        push!(lower, (global_x - half_width * normal_x, center_y - half_width * normal_y))
        push!(upper, (global_x + half_width * normal_x, center_y + half_width * normal_y))
    end
    add_polygon_surface_mm(vcat(lower, reverse(upper)))
end

function fuse_surfaces(tags)
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(tags))],
        [(2, tag) for tag in Iterators.drop(tags, 1)],
    )
    surfaces = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surfaces) && error("folded-path fuse produced no solid domain")
    surfaces
end

function build_domain(config, variant)
    half_width = config.guide_width_mm / 2
    if variant in (:straight_device, :straight_unfolded)
        return [add_rectangle_mm(
            0.0,
            -half_width,
            total_axial_length_mm(config, variant),
            config.guide_width_mm,
        )]
    end

    left = add_rectangle_mm(0.0, -half_width, config.lead_length_mm, config.guide_width_mm)
    right_start = config.lead_length_mm + config.active_axial_length_mm
    right = add_rectangle_mm(
        right_start,
        -half_width,
        config.lead_length_mm,
        config.guide_width_mm,
    )
    active_tags = if variant == :rounded_fold
        [add_rounded_active_surface(config)]
    elseif variant == :sharp_fold
        amplitude = sharp_fold_amplitude_mm(config)
        first_point = (config.lead_length_mm, 0.0)
        apex = (config.lead_length_mm + config.active_axial_length_mm / 2, amplitude)
        last_point = (right_start, 0.0)
        [
            add_segment_strip_mm(first_point, apex, config.guide_width_mm),
            add_segment_strip_mm(apex, last_point, config.guide_width_mm),
        ]
    else
        throw(ArgumentError("unknown folded-path variant"))
    end
    fuse_surfaces(vcat(left, active_tags, right))
end

function classify_boundaries(surfaces, config, variant)
    gmsh.model.occ.synchronize()
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    curves = unique(tag for (dimension, tag) in boundary if dimension == 1)
    total_m = total_axial_length_mm(config, variant) * 1e-3
    tolerance = max(1.0e-10, total_m * 1.0e-8)
    source = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] < tolerance, curves)
    microphone = filter(
        tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] > total_m - tolerance,
        curves,
    )
    length(source) == 1 || error("folded path must have exactly one left port")
    length(microphone) == 1 || error("folded path must have exactly one right port")
    source, microphone, setdiff(curves, vcat(source, microphone))
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

function configure_mesh(source, microphone, free; size_path_mm, size_port_mm, size_max_mm)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", vcat(source, microphone))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_port_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.05e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.0e-3)

    gmsh.model.mesh.field.add("Distance", 3)
    gmsh.model.mesh.field.setNumbers(3, "CurvesList", free)
    gmsh.model.mesh.field.add("Threshold", 4)
    gmsh.model.mesh.field.setNumber(4, "InField", 3)
    gmsh.model.mesh.field.setNumber(4, "SizeMin", Float64(size_path_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(4, "SizeMax", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(4, "DistMin", 0.04e-3)
    gmsh.model.mesh.field.setNumber(4, "DistMax", 0.8e-3)

    gmsh.model.mesh.field.add("Min", 5)
    gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])
    gmsh.model.mesh.field.setAsBackgroundMesh(5)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
end

function build_folded_path_mesh(
    output_path::AbstractString;
    config::FoldedPathConfig=FoldedPathConfig(),
    variant::Symbol=:rounded_fold,
    size_path_mm::Real=0.10,
    size_port_mm::Real=0.08,
    size_max_mm::Real=0.22,
)
    validate_config(config)
    variant in FOLDED_PATH_VARIANTS || throw(ArgumentError("unknown folded-path variant"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("folded_path_$(variant)")
    surfaces = build_domain(config, variant)
    source, microphone, free = classify_boundaries(surfaces, config, variant)
    tag_physical_groups(surfaces, source, microphone, free)
    configure_mesh(source, microphone, free; size_path_mm, size_port_mm, size_max_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
