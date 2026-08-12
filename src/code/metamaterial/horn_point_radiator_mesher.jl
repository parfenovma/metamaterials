module HornPointRadiatorMesher

using Gmsh: gmsh

export HornPointRadiatorConfig,
       HORN_POINT_VARIANTS,
       validate_config,
       horn_height_mm,
       adiabatic_parameter,
       diffuser_height_mm,
       diffuser_adiabatic_parameter,
       rounded_guide_amplitude_mm,
       rounded_guide_length_mm,
       smooth_guide_amplitude_mm,
       smooth_guide_length_mm,
       minimum_smooth_inner_radius_mm,
       outlet_x_mm,
       guide_center_y_mm,
       probe_points_mm,
       build_horn_point_radiator_mesh

const HORN_POINT_VARIANTS = (
    :abrupt_straight,
    :collector_device,
    :collector_straight,
    :collector_rounded,
    :collector_gentle,
    :collector_smooth,
    :collector_straight_diffuser,
    :collector_smooth_diffuser,
)

is_rounded_variant(variant) = variant in (
    :collector_rounded,
    :collector_gentle,
    :collector_smooth,
    :collector_smooth_diffuser,
)
has_smooth_guide(variant) = variant in (:collector_smooth, :collector_smooth_diffuser)
has_diffuser(variant) = variant in (:collector_straight_diffuser, :collector_smooth_diffuser)

"""Deterministic collector → narrow guide → point-radiator pilot from research-plan B1."""
Base.@kwdef struct HornPointRadiatorConfig
    input_height_mm::Float64 = 3.2
    throat_height_mm::Float64 = 1.6
    horn_length_mm::Float64 = 25.0
    lower_band_frequency_hz::Float64 = 162.2e3
    pressure_wave_speed_m_s::Float64 = 2340.0
    straight_guide_length_mm::Float64 = 30.42
    rounded_axial_length_mm::Float64 = 20.0
    receiver_length_mm::Float64 = 55.0
    receiver_half_height_mm::Float64 = 15.0
    target_distance_mm::Float64 = 35.0
    diffuser_output_height_mm::Float64 = 3.2
    diffuser_length_mm::Float64 = 25.0
    profile_samples::Int = 120
    arc_integration_samples::Int = 4001
    minimum_inner_radius_mm::Float64 = 0.35
end

function horn_height_mm(config::HornPointRadiatorConfig, x_mm::Real)
    0 <= x_mm <= config.horn_length_mm ||
        throw(ArgumentError("horn coordinate is outside the collector"))
    q = (1 - cospi(Float64(x_mm) / config.horn_length_mm)) / 2
    config.input_height_mm * exp(q * log(config.throat_height_mm / config.input_height_mm))
end

function adiabatic_parameter(config::HornPointRadiatorConfig)
    k_min_per_mm = 2pi * config.lower_band_frequency_hz /
                   config.pressure_wave_speed_m_s * 1e-3
    maximum_log_impedance_slope = pi * abs(log(config.input_height_mm / config.throat_height_mm)) /
                                  (2config.horn_length_mm)
    maximum_log_impedance_slope / (2k_min_per_mm)
end

function diffuser_height_mm(config::HornPointRadiatorConfig, x_mm::Real)
    0 <= x_mm <= config.diffuser_length_mm ||
        throw(ArgumentError("diffuser coordinate is outside the output flare"))
    q = (1 - cospi(Float64(x_mm) / config.diffuser_length_mm)) / 2
    config.throat_height_mm *
    exp(q * log(config.diffuser_output_height_mm / config.throat_height_mm))
end

function diffuser_adiabatic_parameter(config::HornPointRadiatorConfig)
    k_min_per_mm = 2pi * config.lower_band_frequency_hz /
                   config.pressure_wave_speed_m_s * 1e-3
    maximum_log_impedance_slope =
        pi * abs(log(config.diffuser_output_height_mm / config.throat_height_mm)) /
        (2config.diffuser_length_mm)
    maximum_log_impedance_slope / (2k_min_per_mm)
end

function rounded_centerline_y_mm(config, x_mm, amplitude_mm)
    amplitude_mm / 2 * (1 - cospi(2Float64(x_mm) / config.rounded_axial_length_mm))
end


function rounded_centerline_slope(config, x_mm, amplitude_mm)
    amplitude_mm * pi / config.rounded_axial_length_mm *
    sinpi(2Float64(x_mm) / config.rounded_axial_length_mm)
end

function rounded_guide_length_mm(config, amplitude_mm; samples=config.arc_integration_samples)
    x = range(0.0, config.rounded_axial_length_mm; length=samples)
    integrand = hypot.(1.0, rounded_centerline_slope.(Ref(config), x, amplitude_mm))
    step = config.rounded_axial_length_mm / (length(x) - 1)
    step * (sum(integrand) - (first(integrand) + last(integrand)) / 2)
end

function rounded_guide_amplitude_mm(config::HornPointRadiatorConfig)
    lower, upper = 0.0, config.straight_guide_length_mm
    while rounded_guide_length_mm(config, upper) < config.straight_guide_length_mm
        upper *= 2
    end
    for _ in 1:70
        middle = (lower + upper) / 2
        if rounded_guide_length_mm(config, middle) < config.straight_guide_length_mm
            lower = middle
        else
            upper = middle
        end
    end
    (lower + upper) / 2
end

function smooth_centerline_y_mm(config, x_mm, amplitude_mm)
    t = clamp(Float64(x_mm) / config.rounded_axial_length_mm, 0.0, 1.0)
    amplitude_mm * sinpi(t)^4
end


function smooth_centerline_slope(config, x_mm, amplitude_mm)
    t = clamp(Float64(x_mm) / config.rounded_axial_length_mm, 0.0, 1.0)
    4amplitude_mm * pi / config.rounded_axial_length_mm * sinpi(t)^3 * cospi(t)
end


function smooth_guide_length_mm(config, amplitude_mm; samples=config.arc_integration_samples)
    x = range(0.0, config.rounded_axial_length_mm; length=samples)
    integrand = hypot.(1.0, smooth_centerline_slope.(Ref(config), x, amplitude_mm))
    step = config.rounded_axial_length_mm / (length(x) - 1)
    step * (sum(integrand) - (first(integrand) + last(integrand)) / 2)
end


function smooth_centerline_curvature_per_mm(config, x_mm, amplitude_mm)
    t = clamp(Float64(x_mm) / config.rounded_axial_length_mm, 0.0, 1.0)
    sine = sinpi(t)
    cosine = cospi(t)
    slope = 4amplitude_mm * pi / config.rounded_axial_length_mm * sine^3 * cosine
    second = 4amplitude_mm * pi^2 / config.rounded_axial_length_mm^2 *
             (3sine^2 * cosine^2 - sine^4)
    abs(second) / (1 + slope^2)^(3 / 2)
end


function minimum_smooth_inner_radius_mm(config::HornPointRadiatorConfig)
    amplitude = smooth_guide_amplitude_mm(config)
    coordinate = range(0.0, config.rounded_axial_length_mm;
                       length=config.arc_integration_samples)
    maximum_curvature = maximum(
        smooth_centerline_curvature_per_mm.(Ref(config), coordinate, amplitude),
    )
    iszero(maximum_curvature) && return Inf
    inv(maximum_curvature) - config.throat_height_mm / 2
end


function smooth_guide_amplitude_mm(config::HornPointRadiatorConfig)
    lower, upper = 0.0, config.straight_guide_length_mm
    while smooth_guide_length_mm(config, upper) < config.straight_guide_length_mm
        upper *= 2
    end
    for _ in 1:70
        middle = (lower + upper) / 2
        if smooth_guide_length_mm(config, middle) < config.straight_guide_length_mm
            lower = middle
        else
            upper = middle
        end
    end
    (lower + upper) / 2
end

function minimum_rounded_inner_radius_mm(config)
    amplitude = rounded_guide_amplitude_mm(config)
    maximum_curvature = 2pi^2 * amplitude / config.rounded_axial_length_mm^2
    inv(maximum_curvature) - config.throat_height_mm / 2
end

function validate_config(config::HornPointRadiatorConfig; variant::Symbol=:collector_rounded)
    config.input_height_mm > 0 && config.throat_height_mm > 0 ||
        throw(ArgumentError("input transition dimensions must be positive"))
    !isapprox(config.input_height_mm, config.throat_height_mm; atol=1e-12) ||
        throw(ArgumentError("input transition needs unequal endpoint heights"))
    config.horn_length_mm > 0 || throw(ArgumentError("horn length must be positive"))
    config.lower_band_frequency_hz > 0 || throw(ArgumentError("lower band frequency must be positive"))
    config.pressure_wave_speed_m_s > 0 || throw(ArgumentError("P-wave speed must be positive"))
    config.straight_guide_length_mm > config.rounded_axial_length_mm > 0 ||
        throw(ArgumentError("rounded guide must pack a longer path into a shorter axial length"))
    config.receiver_length_mm > config.target_distance_mm > 0 ||
        throw(ArgumentError("target must lie inside the receiver"))
    config.receiver_half_height_mm > config.input_height_mm ||
        throw(ArgumentError("receiver is too narrow"))
    config.profile_samples >= 32 || throw(ArgumentError("profiles need at least 32 samples"))
    config.arc_integration_samples >= 101 ||
        throw(ArgumentError("arc integration needs at least 101 samples"))
    adiabatic_parameter(config) <= 0.051 ||
        throw(ArgumentError("collector violates the epsilon_ad <= 0.05 design target"))
    if has_diffuser(variant)
        config.diffuser_output_height_mm > 0 ||
            throw(ArgumentError("output transition height must be positive"))
        !isapprox(config.diffuser_output_height_mm, config.throat_height_mm; atol=1e-12) ||
            throw(ArgumentError("output transition needs unequal endpoint heights"))
        config.diffuser_length_mm > 0 ||
            throw(ArgumentError("diffuser length must be positive"))
        diffuser_adiabatic_parameter(config) <= 0.051 ||
            throw(ArgumentError("diffuser violates the epsilon_ad <= 0.05 design target"))
    end
    if variant == :collector_smooth
        minimum_smooth_inner_radius_mm(config) >= config.minimum_inner_radius_mm ||
            throw(ArgumentError("smooth delay guide violates the minimum inner radius"))
    elseif is_rounded_variant(variant)
        minimum_rounded_inner_radius_mm(config) >= config.minimum_inner_radius_mm ||
            throw(ArgumentError("rounded delay guide violates the minimum inner radius"))
    end
    nothing
end

function outlet_x_mm(config::HornPointRadiatorConfig, variant::Symbol)
    variant in HORN_POINT_VARIANTS || throw(ArgumentError("unknown horn-point variant"))
    guide_axial = is_rounded_variant(variant) ?
                  config.rounded_axial_length_mm : config.straight_guide_length_mm
    config.horn_length_mm + guide_axial +
    (has_diffuser(variant) ? config.diffuser_length_mm : 0.0)
end

function guide_center_y_mm(config::HornPointRadiatorConfig, variant::Symbol, local_x_mm::Real)
    is_rounded_variant(variant) || return 0.0
    if has_smooth_guide(variant)
        return smooth_centerline_y_mm(
            config,
            clamp(Float64(local_x_mm), 0.0, config.rounded_axial_length_mm),
            smooth_guide_amplitude_mm(config),
        )
    end
    rounded_centerline_y_mm(
        config,
        clamp(Float64(local_x_mm), 0.0, config.rounded_axial_length_mm),
        rounded_guide_amplitude_mm(config),
    )
end

function probe_points_mm(config::HornPointRadiatorConfig, variant::Symbol)
    outlet = outlet_x_mm(config, variant)
    throat_local_x = 2.0
    names = String["throat", "near_radiator", "target_axis"]
    x = Float64[
        config.horn_length_mm + throat_local_x,
        outlet + 2.0,
        outlet + config.target_distance_mm,
    ]
    y = Float64[
        guide_center_y_mm(config, variant, throat_local_x),
        0.0,
        0.0,
    ]
    for transverse_y in -12.0:1.0:12.0
        push!(names, "target_y_$(replace(string(transverse_y), "." => "p", "-" => "m"))")
        push!(x, outlet + config.target_distance_mm)
        push!(y, transverse_y)
    end
    (names=names, x_mm=x, y_mm=y)
end

function add_polygon_surface_mm(points_mm)
    clean = NTuple{2, Float64}[]
    for (x, y) in points_mm
        point = (Float64(x), Float64(y))
        (isempty(clean) || point != last(clean)) && push!(clean, point)
    end
    first(clean) == last(clean) && pop!(clean)
    points = [gmsh.model.occ.addPoint(x * 1e-3, y * 1e-3, 0.0) for (x, y) in clean]
    lines = [
        gmsh.model.occ.addLine(points[index], points[mod1(index + 1, length(points))])
        for index in eachindex(points)
    ]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function add_horn_surface(config)
    x = range(0.0, config.horn_length_mm; length=config.profile_samples)
    half_height = horn_height_mm.(Ref(config), x) ./ 2
    points = vcat(collect(zip(x, half_height)), collect(zip(reverse(x), -reverse(half_height))))
    add_polygon_surface_mm(points)
end

function add_straight_guide(config)
    gmsh.model.occ.addRectangle(
        config.horn_length_mm * 1e-3,
        -config.throat_height_mm * 0.5e-3,
        0.0,
        config.straight_guide_length_mm * 1e-3,
        config.throat_height_mm * 1e-3,
    )
end

function add_rounded_guide(config)
    amplitude = rounded_guide_amplitude_mm(config)
    half_width = config.throat_height_mm / 2
    lower = NTuple{2, Float64}[]
    upper = NTuple{2, Float64}[]
    for local_x in range(0.0, config.rounded_axial_length_mm; length=config.profile_samples)
        center_y = rounded_centerline_y_mm(config, local_x, amplitude)
        slope = rounded_centerline_slope(config, local_x, amplitude)
        normalization = hypot(1.0, slope)
        normal_x, normal_y = -slope / normalization, 1 / normalization
        global_x = config.horn_length_mm + local_x
        push!(lower, (global_x - half_width * normal_x, center_y - half_width * normal_y))
        push!(upper, (global_x + half_width * normal_x, center_y + half_width * normal_y))
    end
    add_polygon_surface_mm(vcat(lower, reverse(upper)))
end


function add_smooth_guide(config)
    amplitude = smooth_guide_amplitude_mm(config)
    half_width = config.throat_height_mm / 2
    lower = NTuple{2, Float64}[]
    upper = NTuple{2, Float64}[]
    for local_x in range(0.0, config.rounded_axial_length_mm; length=config.profile_samples)
        center_y = smooth_centerline_y_mm(config, local_x, amplitude)
        slope = smooth_centerline_slope(config, local_x, amplitude)
        normalization = hypot(1.0, slope)
        normal_x, normal_y = -slope / normalization, 1 / normalization
        global_x = config.horn_length_mm + local_x
        push!(lower, (global_x - half_width * normal_x, center_y - half_width * normal_y))
        push!(upper, (global_x + half_width * normal_x, center_y + half_width * normal_y))
    end
    add_polygon_surface_mm(vcat(lower, reverse(upper)))
end

function add_diffuser(config, start_x_mm)
    local_x = range(0.0, config.diffuser_length_mm; length=config.profile_samples)
    half_height = diffuser_height_mm.(Ref(config), local_x) ./ 2
    global_x = start_x_mm .+ local_x
    points = vcat(
        collect(zip(global_x, half_height)),
        collect(zip(reverse(global_x), -reverse(half_height))),
    )
    add_polygon_surface_mm(points)
end

function add_input_surface(config, variant)
    if variant == :abrupt_straight
        return gmsh.model.occ.addRectangle(
            0.0,
            -config.input_height_mm * 0.5e-3,
            0.0,
            config.horn_length_mm * 1e-3,
            config.input_height_mm * 1e-3,
        )
    end
    add_horn_surface(config)
end

function build_domain(config, variant)
    outlet = outlet_x_mm(config, variant)
    input = add_input_surface(config, variant)
    guide = has_smooth_guide(variant) ? add_smooth_guide(config) :
            is_rounded_variant(variant) ? add_rounded_guide(config) : add_straight_guide(config)
    guide_axial = is_rounded_variant(variant) ?
                  config.rounded_axial_length_mm : config.straight_guide_length_mm
    diffuser = has_diffuser(variant) ?
               add_diffuser(config, config.horn_length_mm + guide_axial) : nothing
    receiver = gmsh.model.occ.addRectangle(
        outlet * 1e-3,
        -config.receiver_half_height_mm * 1e-3,
        0.0,
        config.receiver_length_mm * 1e-3,
        2config.receiver_half_height_mm * 1e-3,
    )
    tools = isnothing(diffuser) ? [(2, guide), (2, receiver)] :
            [(2, guide), (2, diffuser), (2, receiver)]
    fused, _ = gmsh.model.occ.fuse([(2, input)], tools)
    surfaces = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surfaces) && error("horn-point fuse produced no solid domain")
    surfaces
end

function classify_boundaries(surfaces, config, variant)
    gmsh.model.occ.synchronize()
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    curves = unique(tag for (dimension, tag) in boundary if dimension == 1)
    outlet = outlet_x_mm(config, variant) * 1e-3
    end_x = outlet + config.receiver_length_mm * 1e-3
    receiver_half = config.receiver_half_height_mm * 1e-3
    tolerance = 1.0e-8
    source = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] < tolerance, curves)
    radiation = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1] - end_x) < tolerance ||
        (center[1] > outlet + tolerance && abs(abs(center[2]) - receiver_half) < tolerance)
    end
    length(source) == 1 || error("horn-point geometry must have one input boundary")
    !isempty(radiation) || error("receiver radiation boundary was not found")
    source, radiation, setdiff(curves, vcat(source, radiation))
end

function add_physical_group(dim, tags, id, name)
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function configure_mesh(config, variant, free; size_path_mm, size_receiver_mm, size_radiator_mm)
    outlet = outlet_x_mm(config, variant)
    guide_amplitude_mm = has_smooth_guide(variant) ? smooth_guide_amplitude_mm(config) :
                         is_rounded_variant(variant) ? rounded_guide_amplitude_mm(config) : 0.0
    path_lower_mm = min(-config.input_height_mm / 2, -config.throat_height_mm / 2) - 2.0
    path_upper_mm = max(config.input_height_mm / 2,
                        guide_amplitude_mm + config.throat_height_mm / 2) + 2.0
    gmsh.model.mesh.field.add("Box", 1)
    gmsh.model.mesh.field.setNumber(1, "VIn", Float64(size_path_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "VOut", Float64(size_receiver_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "XMin", -0.2e-3)
    gmsh.model.mesh.field.setNumber(1, "XMax", (outlet + 2.0) * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "YMin", path_lower_mm * 1e-3)
    gmsh.model.mesh.field.setNumber(1, "YMax", path_upper_mm * 1e-3)

    radiator_curves = filter(free) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1] - outlet * 1e-3) < 0.2e-3
    end
    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "CurvesList", radiator_curves)
    gmsh.model.mesh.field.add("Threshold", 3)
    gmsh.model.mesh.field.setNumber(3, "InField", 2)
    gmsh.model.mesh.field.setNumber(3, "SizeMin", Float64(size_radiator_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "SizeMax", Float64(size_receiver_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "DistMin", 0.1e-3)
    gmsh.model.mesh.field.setNumber(3, "DistMax", 8.0e-3)

    gmsh.model.mesh.field.add("Min", 4)
    gmsh.model.mesh.field.setNumbers(4, "FieldsList", [1, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(4)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
end

function build_horn_point_radiator_mesh(
    output_path::AbstractString;
    config::HornPointRadiatorConfig=HornPointRadiatorConfig(),
    variant::Symbol=:collector_straight,
    size_path_mm::Real=0.20,
    size_receiver_mm::Real=0.65,
    size_radiator_mm::Real=0.22,
)
    validate_config(config; variant)
    variant in HORN_POINT_VARIANTS || throw(ArgumentError("unknown horn-point variant"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("horn_point_radiator_$(variant)")
    surfaces = build_domain(config, variant)
    source, radiation, free = classify_boundaries(surfaces, config, variant)
    add_physical_group(1, source, 101, "Source")
    add_physical_group(1, radiation, 102, "RadiationBoundary")
    add_physical_group(1, free, 103, "FreeSurface")
    add_physical_group(2, surfaces, 201, "Domain")
    configure_mesh(config, variant, free; size_path_mm, size_receiver_mm, size_radiator_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
