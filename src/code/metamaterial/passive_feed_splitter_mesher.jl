module PassiveFeedSplitterMesher

using Gmsh: gmsh

export PassiveFeedSplitterConfig,
       validate_config,
       small_junction_width_mm,
       large_junction_width_mm,
       target_output_amplitude_ratio,
       adiabatic_parameter,
       output_x_mm,
       build_passive_feed_splitter_mesh,
       build_passive_feed_straight_control_mesh

Base.@kwdef struct PassiveFeedSplitterConfig
    input_width_mm::Float64 = 7.0
    small_power_fraction::Float64 = 0.14446029370996585
    input_straight_length_mm::Float64 = 20.0
    transformer_length_mm::Float64 = 131.00271326349065
    output_straight_length_mm::Float64 = 20.0
    output_width_mm::Float64 = 7.0
    output_center_offset_mm::Float64 = 6.0
    splitter_tip_radius_mm::Float64 = 0.0
    lower_band_frequency_hz::Float64 = 226.04e3
    pressure_wave_speed_m_s::Float64 = 6122.102437409232
    profile_samples::Int = 240
end

small_junction_width_mm(config) = config.input_width_mm * config.small_power_fraction
large_junction_width_mm(config) =
    config.input_width_mm * (1 - config.small_power_fraction)
target_output_amplitude_ratio(config) = sqrt(
    config.small_power_fraction / (1 - config.small_power_fraction),
)
output_x_mm(config) = config.input_straight_length_mm +
                      config.transformer_length_mm +
                      config.output_straight_length_mm

function adiabatic_parameter(config, initial_width_mm)
    k_min_per_mm = 2pi * config.lower_band_frequency_hz /
                   config.pressure_wave_speed_m_s * 1e-3
    maximum_log_slope = pi * abs(log(config.output_width_mm / initial_width_mm)) /
                        (2config.transformer_length_mm)
    maximum_log_slope / (2k_min_per_mm)
end

function validate_config(config::PassiveFeedSplitterConfig)
    config.input_width_mm > 0 || throw(ArgumentError("input width must be positive"))
    0 < config.small_power_fraction < 0.5 ||
        throw(ArgumentError("small power fraction must lie between zero and one half"))
    config.output_width_mm > 0 || throw(ArgumentError("output width must be positive"))
    config.input_straight_length_mm > 0 ||
        throw(ArgumentError("input straight length must be positive"))
    config.transformer_length_mm > 0 ||
        throw(ArgumentError("transformer length must be positive"))
    config.output_straight_length_mm > 0 ||
        throw(ArgumentError("output straight length must be positive"))
    config.output_center_offset_mm > config.output_width_mm / 2 ||
        throw(ArgumentError("output branches overlap"))
    maximum_tip_radius = config.output_center_offset_mm - config.output_width_mm / 2
    0 <= config.splitter_tip_radius_mm <= maximum_tip_radius ||
        throw(ArgumentError("splitter tip radius must fit inside the final branch gap"))
    config.profile_samples >= 80 || throw(ArgumentError("splitter profile is undersampled"))
    adiabatic_parameter(config, small_junction_width_mm(config)) <= 0.051 ||
        throw(ArgumentError("small-branch transformer violates epsilon_ad <= 0.05"))
    adiabatic_parameter(config, large_junction_width_mm(config)) <= 0.051 ||
        throw(ArgumentError("large-branch transformer violates epsilon_ad <= 0.05"))
    nothing
end

transition_q(local_x_mm, length_mm) =
    (1 - cospi(Float64(local_x_mm) / Float64(length_mm))) / 2

function branch_polygon_mm(config, branch::Symbol)
    branch in (:small, :large) || throw(ArgumentError("branch must be small or large"))
    small_width = small_junction_width_mm(config)
    large_width = large_junction_width_mm(config)
    if branch == :small
        initial_width = small_width
        initial_center = config.input_width_mm / 2 - small_width / 2
        final_center = config.output_center_offset_mm
    else
        initial_width = large_width
        initial_center = -config.input_width_mm / 2 + large_width / 2
        final_center = -config.output_center_offset_mm
    end
    upper = NTuple{2, Float64}[]
    lower = NTuple{2, Float64}[]
    for local_x in range(
        0.0, config.transformer_length_mm; length=config.profile_samples,
    )
        q = transition_q(local_x, config.transformer_length_mm)
        center = initial_center + q * (final_center - initial_center)
        width = initial_width * exp(q * log(config.output_width_mm / initial_width))
        x = config.input_straight_length_mm + local_x
        branch_upper = center + width / 2
        branch_lower = center - width / 2
        radius = config.splitter_tip_radius_mm
        if radius > 0
            final_half_gap = config.output_center_offset_mm - config.output_width_mm / 2
            half_gap = if local_x <= radius
                sqrt(max(radius^2 - (local_x - radius)^2, 0.0))
            else
                local_q = transition_q(
                    local_x - radius, config.transformer_length_mm - radius,
                )
                radius + local_q * (final_half_gap - radius)
            end
            if branch == :small
                branch_lower = half_gap
            else
                branch_upper = -half_gap
            end
        end
        push!(upper, (x, branch_upper))
        push!(lower, (x, branch_lower))
    end
    x_end = output_x_mm(config)
    push!(upper, (x_end, final_center + config.output_width_mm / 2))
    push!(lower, (x_end, final_center - config.output_width_mm / 2))
    vcat(upper, reverse(lower))
end

function add_polygon_surface_mm(points_mm)
    clean = NTuple{2, Float64}[]
    for point in points_mm
        (isempty(clean) || point != last(clean)) && push!(clean, point)
    end
    first(clean) == last(clean) && pop!(clean)
    points = [gmsh.model.occ.addPoint(x * 1e-3, y * 1e-3, 0.0) for (x, y) in clean]
    lines = [gmsh.model.occ.addLine(
        points[index], points[mod1(index + 1, length(points))],
    ) for index in eachindex(points)]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
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

function build_passive_feed_splitter_mesh(
    output_path::AbstractString;
    config::PassiveFeedSplitterConfig=PassiveFeedSplitterConfig(),
    mesh_size_mm::Real=0.38,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("passive_feed_splitter")
    input = gmsh.model.occ.addRectangle(
        0.0, -config.input_width_mm * 0.5e-3, 0.0,
        config.input_straight_length_mm * 1e-3,
        config.input_width_mm * 1e-3,
    )
    small = add_polygon_surface_mm(branch_polygon_mm(config, :small))
    large = add_polygon_surface_mm(branch_polygon_mm(config, :large))
    fused, _ = gmsh.model.occ.fuse([(2, input)], [(2, small), (2, large)])
    surfaces = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surfaces) && error("splitter fuse produced no domain")
    gmsh.model.occ.synchronize()
    curves = external_curves()
    tolerance = 1e-8
    end_x_m = output_x_mm(config) * 1e-3
    source = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1]) < tolerance
    end
    output_small = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1] - end_x_m) < tolerance && center[2] > 0
    end
    output_large = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1] - end_x_m) < tolerance && center[2] < 0
    end
    length(source) == length(output_small) == length(output_large) == 1 ||
        error("splitter source/output boundary classification failed")
    free = setdiff(curves, vcat(source, output_small, output_large))
    add_physical_group(1, source, 101, "Source")
    add_physical_group(1, output_small, 102, "OutputSmall")
    add_physical_group(1, output_large, 103, "OutputLarge")
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
    println("[+] passive splitter mesh: $node_count nodes, $element_count elements")
    (; node_count, element_count, output_path)
end

"""Straight, unsplit 7 mm guide with the outlet divided into two audit tags."""
function build_passive_feed_straight_control_mesh(
    output_path::AbstractString;
    config::PassiveFeedSplitterConfig=PassiveFeedSplitterConfig(),
    mesh_size_mm::Real=0.38,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("passive_feed_straight_control")
    length_m = output_x_mm(config) * 1e-3
    half_width_m = config.input_width_mm * 0.5e-3
    lower = gmsh.model.occ.addRectangle(
        0.0, -half_width_m, 0.0, length_m, half_width_m,
    )
    upper = gmsh.model.occ.addRectangle(
        0.0, 0.0, 0.0, length_m, half_width_m,
    )
    fragmented, _ = gmsh.model.occ.fragment([(2, lower)], [(2, upper)])
    surfaces = unique(tag for (dimension, tag) in fragmented if dimension == 2)
    gmsh.model.occ.synchronize()
    curves = external_curves()
    tolerance = 1e-8
    source = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1]) < tolerance
    end
    output_small = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1] - length_m) < tolerance && center[2] > 0
    end
    output_large = filter(curves) do tag
        center = gmsh.model.occ.getCenterOfMass(1, tag)
        abs(center[1] - length_m) < tolerance && center[2] < 0
    end
    length(source) == 2 || error("straight-control source classification failed")
    length(output_small) == length(output_large) == 1 ||
        error("straight-control output classification failed")
    free = setdiff(curves, vcat(source, output_small, output_large))
    add_physical_group(1, source, 101, "Source")
    add_physical_group(1, output_small, 102, "OutputSmall")
    add_physical_group(1, output_large, 103, "OutputLarge")
    add_physical_group(1, vcat(output_small, output_large), 105, "Output")
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
    println("[+] passive straight control mesh: $node_count nodes, $element_count elements")
    (; node_count, element_count, output_path)
end

end
