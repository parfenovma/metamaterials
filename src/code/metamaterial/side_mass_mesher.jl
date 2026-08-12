module SideMassMesher

using Gmsh: gmsh

export SideMassConfig, validate_config, build_side_mass_mesh, build_side_mass_array_mesh

Base.@kwdef struct SideMassConfig
    length_mm::Float64 = 30.0
    height_mm::Float64 = 7.0
    lead_length_mm::Float64 = 4.0
    spine_height_mm::Float64 = 3.0
    bright_center_x_mm::Float64 = 10.0
    dark_center_x_mm::Float64 = 16.0
    bright_mass_length_mm::Float64 = 3.0
    bright_mass_height_mm::Float64 = 1.1
    dark_mass_length_mm::Float64 = 3.0
    dark_mass_height_mm::Float64 = 1.1
    wall_margin_mm::Float64 = 0.35
    bright_neck_width_mm::Float64 = 0.6
    dark_bridge_height_mm::Float64 = 0.5
end

function validate_config(config::SideMassConfig)
    config.length_mm > 0 || throw(ArgumentError("length must be positive"))
    config.height_mm > 0 || throw(ArgumentError("height must be positive"))
    0 < config.lead_length_mm < config.length_mm / 2 ||
        throw(ArgumentError("lead length must leave a central resonant section"))
    0 < config.spine_height_mm < config.height_mm ||
        throw(ArgumentError("spine height must lie inside the element"))
    config.bright_mass_length_mm > 0 ||
        throw(ArgumentError("bright mass length must be positive"))
    config.bright_mass_height_mm > 0 ||
        throw(ArgumentError("bright mass height must be positive"))
    config.dark_mass_length_mm > 0 ||
        throw(ArgumentError("dark mass length must be positive"))
    config.dark_mass_height_mm > 0 ||
        throw(ArgumentError("dark mass height must be positive"))
    config.bright_neck_width_mm > 0 || throw(ArgumentError("bright neck must be positive"))
    config.dark_bridge_height_mm > 0 || throw(ArgumentError("dark bridge must be positive"))
    config.bright_neck_width_mm <= config.bright_mass_length_mm ||
        throw(ArgumentError("bright neck cannot be wider than the bright mass"))
    config.dark_bridge_height_mm <=
        min(config.bright_mass_height_mm, config.dark_mass_height_mm) ||
        throw(ArgumentError("dark bridge cannot be taller than the connected masses"))
    config.wall_margin_mm >= 0 || throw(ArgumentError("wall margin must be non-negative"))

    spine_top = (config.height_mm + config.spine_height_mm) / 2
    bright_mass_bottom =
        config.height_mm - config.wall_margin_mm - config.bright_mass_height_mm
    dark_mass_bottom =
        config.height_mm - config.wall_margin_mm - config.dark_mass_height_mm
    bright_mass_bottom > spine_top ||
        throw(ArgumentError("bright masses leave no neck clearance"))
    dark_mass_bottom > spine_top ||
        throw(ArgumentError("dark masses overlap the central spine"))

    left_limit = config.lead_length_mm
    right_limit = config.length_mm - config.lead_length_mm
    bright_half_length = config.bright_mass_length_mm / 2
    dark_half_length = config.dark_mass_length_mm / 2
    left_limit < config.bright_center_x_mm - bright_half_length ||
        throw(ArgumentError("bright mass overlaps the left lead"))
    config.dark_center_x_mm + dark_half_length < right_limit ||
        throw(ArgumentError("dark mass overlaps the right lead"))
    bright_right = config.bright_center_x_mm + bright_half_length
    dark_left = config.dark_center_x_mm - dark_half_length
    bright_right < dark_left || throw(ArgumentError("bright and dark masses overlap"))
    nothing
end

function add_rectangle_mm(x_mm, y_mm, width_mm, height_mm)
    gmsh.model.occ.addRectangle(
        x_mm * 1e-3,
        y_mm * 1e-3,
        0.0,
        width_mm * 1e-3,
        height_mm * 1e-3,
    )
end

function solid_rectangles(
    config::SideMassConfig;
    resonators::Bool=true,
    variant::Union{Nothing, Symbol}=nothing,
)
    selected_variant = isnothing(variant) ? (resonators ? :bd : :r0) : variant
    selected_variant in (:r0, :b, :bd) ||
        throw(ArgumentError("variant must be :r0, :b, or :bd"))
    spine_bottom = (config.height_mm - config.spine_height_mm) / 2
    rectangles = [
        (0.0, 0.0, config.lead_length_mm, config.height_mm),
        (
            config.lead_length_mm,
            spine_bottom,
            config.length_mm - 2config.lead_length_mm,
            config.spine_height_mm,
        ),
        (
            config.length_mm - config.lead_length_mm,
            0.0,
            config.lead_length_mm,
            config.height_mm,
        ),
    ]
    selected_variant == :r0 && return rectangles

    bright_mass_bottom =
        config.height_mm - config.wall_margin_mm - config.bright_mass_height_mm
    top_spine = spine_bottom + config.spine_height_mm
    neck_height = bright_mass_bottom - top_spine
    bright_mass_left = config.bright_center_x_mm - config.bright_mass_length_mm / 2
    push!(rectangles, (
        bright_mass_left,
        bright_mass_bottom,
        config.bright_mass_length_mm,
        config.bright_mass_height_mm,
    ))
    push!(rectangles, (
        bright_mass_left,
        config.wall_margin_mm,
        config.bright_mass_length_mm,
        config.bright_mass_height_mm,
    ))

    # Only the bright pair is connected directly to the spine.
    bright_neck_left = config.bright_center_x_mm - config.bright_neck_width_mm / 2
    push!(rectangles, (bright_neck_left, top_spine, config.bright_neck_width_mm, neck_height))
    push!(rectangles, (
        bright_neck_left,
        config.wall_margin_mm + config.bright_mass_height_mm,
        config.bright_neck_width_mm,
        neck_height,
    ))

    selected_variant == :b && return rectangles

    dark_mass_bottom =
        config.height_mm - config.wall_margin_mm - config.dark_mass_height_mm
    dark_mass_left = config.dark_center_x_mm - config.dark_mass_length_mm / 2
    push!(rectangles, (
        dark_mass_left,
        dark_mass_bottom,
        config.dark_mass_length_mm,
        config.dark_mass_height_mm,
    ))
    push!(rectangles, (
        dark_mass_left,
        config.wall_margin_mm,
        config.dark_mass_length_mm,
        config.dark_mass_height_mm,
    ))

    bright_right = config.bright_center_x_mm + config.bright_mass_length_mm / 2
    dark_left = config.dark_center_x_mm - config.dark_mass_length_mm / 2
    bridge_length = dark_left - bright_right
    top_bridge_y = max(bright_mass_bottom, dark_mass_bottom) +
                   (min(config.bright_mass_height_mm, config.dark_mass_height_mm) -
                    config.dark_bridge_height_mm) / 2
    bottom_bridge_y = config.wall_margin_mm +
                      (min(config.bright_mass_height_mm, config.dark_mass_height_mm) -
                       config.dark_bridge_height_mm) / 2
    push!(rectangles, (bright_right, top_bridge_y, bridge_length, config.dark_bridge_height_mm))
    push!(rectangles, (bright_right, bottom_bridge_y, bridge_length, config.dark_bridge_height_mm))
    rectangles
end

function classify_boundaries(surface_tags, config)
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surface_tags], false, false, false)
    line_tags = unique(tag for (dimension, tag) in boundary if dimension == 1)
    source, microphone, free_surface = Int[], Int[], Int[]
    length_m = config.length_mm * 1e-3
    for tag in line_tags
        x = gmsh.model.occ.getCenterOfMass(1, tag)[1]
        if x < 1e-9
            push!(source, tag)
        elseif x > length_m - 1e-9
            push!(microphone, tag)
        else
            push!(free_surface, tag)
        end
    end
    source, microphone, free_surface
end

function build_side_mass_mesh(
    output_path::AbstractString;
    config::SideMassConfig=SideMassConfig(),
    resonators::Bool=true,
    variant::Union{Nothing, Symbol}=nothing,
    size_min_mm::Real=0.10,
    size_max_mm::Real=0.45,
)
    validate_config(config)
    mkpath(dirname(output_path))
    gmsh.clear()
    selected_variant = isnothing(variant) ? (resonators ? :bd : :r0) : variant
    gmsh.model.add("side_mass_$(selected_variant)")
    rectangle_tags = [
        add_rectangle_mm(rectangle...)
        for rectangle in solid_rectangles(config; resonators, variant=selected_variant)
    ]
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(rectangle_tags))],
        [(2, tag) for tag in Iterators.drop(rectangle_tags, 1)],
    )
    gmsh.model.occ.synchronize()
    surface_tags = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surface_tags) && error("OCC fuse produced no solid surface")
    source, microphone, free_surface = classify_boundaries(surface_tags, config)
    isempty(source) && error("source boundary was not detected")
    isempty(microphone) && error("microphone boundary was not detected")

    gmsh.model.addPhysicalGroup(1, source, 101)
    gmsh.model.setPhysicalName(1, 101, "Source")
    gmsh.model.addPhysicalGroup(1, microphone, 102)
    gmsh.model.setPhysicalName(1, 102, "Microphone")
    gmsh.model.addPhysicalGroup(1, free_surface, 103)
    gmsh.model.setPhysicalName(1, 103, "FreeSurface")
    gmsh.model.addPhysicalGroup(2, surface_tags, 201)
    gmsh.model.setPhysicalName(2, 201, "Domain")

    all_boundaries = vcat(source, microphone, free_surface)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", all_boundaries)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_min_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.15e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.5e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

function build_side_mass_array_mesh(
    output_path::AbstractString;
    config::SideMassConfig=SideMassConfig(),
    element_count::Integer=2,
    pitch_mm::Real=8.2,
    resonators::Bool=true,
    variant::Union{Nothing, Symbol}=nothing,
    size_min_mm::Real=0.10,
    size_max_mm::Real=0.45,
)
    validate_config(config)
    element_count >= 2 || throw(ArgumentError("array mesh requires at least two elements"))
    pitch = Float64(pitch_mm)
    pitch >= config.height_mm || throw(ArgumentError("pitch must not overlap neighbouring cells"))
    total_height = config.height_mm + (element_count - 1) * pitch
    rectangles = [
        (0.0, 0.0, config.lead_length_mm, total_height),
        (
            config.length_mm - config.lead_length_mm,
            0.0,
            config.lead_length_mm,
            total_height,
        ),
    ]
    selected_variant = isnothing(variant) ? (resonators ? :bd : :r0) : variant
    cell_rectangles = solid_rectangles(config; resonators, variant=selected_variant)
    internal_rectangles = vcat(cell_rectangles[2:2], cell_rectangles[4:end])
    for index in 0:(element_count - 1)
        y_offset = index * pitch
        append!(rectangles, [(x, y + y_offset, width, height) for (x, y, width, height) in internal_rectangles])
    end

    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("side_mass_array_$(selected_variant)")
    rectangle_tags = [add_rectangle_mm(rectangle...) for rectangle in rectangles]
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(rectangle_tags))],
        [(2, tag) for tag in Iterators.drop(rectangle_tags, 1)],
    )
    gmsh.model.occ.synchronize()
    surface_tags = unique(tag for (dimension, tag) in fused if dimension == 2)
    source, microphone, free_surface = classify_boundaries(surface_tags, config)
    isempty(source) && error("array source boundary was not detected")
    isempty(microphone) && error("array microphone boundary was not detected")
    gmsh.model.addPhysicalGroup(1, source, 101)
    gmsh.model.setPhysicalName(1, 101, "Source")
    gmsh.model.addPhysicalGroup(1, microphone, 102)
    gmsh.model.setPhysicalName(1, 102, "Microphone")
    gmsh.model.addPhysicalGroup(1, free_surface, 103)
    gmsh.model.setPhysicalName(1, 103, "FreeSurface")
    gmsh.model.addPhysicalGroup(2, surface_tags, 201)
    gmsh.model.setPhysicalName(2, 201, "Domain")

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", vcat(source, microphone, free_surface))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_min_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_max_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.15e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.5e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
