module CombinedKMKMesher

using Gmsh: gmsh

include(joinpath(@__DIR__, "inline_mass_component_mesher.jl"))
using .InlineMassComponentMesher
include(joinpath(@__DIR__, "stiffness_component_mesher.jl"))
using .StiffnessComponentMesher

export CombinedKMKConfig,
       validate_config,
       cell_length_mm,
       total_length_mm,
       solid_rectangles,
       build_combined_kmk_mesh

Base.@kwdef struct CombinedKMKConfig
    lead_length_mm::Float64 = 12.0
    mass::InlineMassConfig = InlineMassConfig(
        mass_length_mm=1.4,
        mass_height_mm=2.4,
        neck_length_mm=0.35,
        neck_height_mm=0.95,
    )
    stiffness::StiffnessConfig = StiffnessConfig(ligament_height_mm=0.60)
end

cell_length_mm(config::CombinedKMKConfig) =
    2config.stiffness.length_mm + component_length_mm(config.mass)
total_length_mm(config::CombinedKMKConfig) =
    2config.lead_length_mm + cell_length_mm(config)

function validate_config(config::CombinedKMKConfig)
    config.lead_length_mm > 0 || throw(ArgumentError("lead length must be positive"))
    InlineMassComponentMesher.validate_config(config.mass)
    StiffnessComponentMesher.validate_config(config.stiffness)
    isapprox(config.mass.port_height_mm, config.stiffness.height_mm; atol=1.0e-12) ||
        throw(ArgumentError("M and K must use the same strip height"))
    nothing
end

function offset_rectangle(rectangle, offset_x_mm)
    x, y, width, height = rectangle
    (x + offset_x_mm, y, width, height)
end

"""Solid rectangles before subtracting K-slots."""
function solid_rectangles(config::CombinedKMKConfig; variant::Symbol=:combined)
    validate_config(config)
    variant in (:solid, :k_only, :m_only, :combined) ||
        throw(ArgumentError("variant must be :solid, :k_only, :m_only, or :combined"))
    height = config.stiffness.height_mm
    y_bottom = -height / 2
    total = total_length_mm(config)
    variant == :solid && return [(0.0, y_bottom, total, height)]

    lead = config.lead_length_mm
    k_length = config.stiffness.length_mm
    m_length = component_length_mm(config.mass)
    if variant == :k_only
        return [(0.0, y_bottom, total, height)]
    end

    rectangles = [
        (0.0, y_bottom, lead + k_length, height),
    ]
    m_start = lead + k_length
    append!(
        rectangles,
        [offset_rectangle(rectangle, m_start) for rectangle in
         InlineMassComponentMesher.solid_rectangles(config.mass)],
    )
    push!(rectangles, (m_start + m_length, y_bottom, k_length + lead, height))
    rectangles
end

function add_rectangle_mm(rectangle)
    x, y, width, height = rectangle
    gmsh.model.occ.addRectangle(x * 1e-3, y * 1e-3, 0.0, width * 1e-3, height * 1e-3)
end

function fused_surfaces(rectangles)
    tags = add_rectangle_mm.(rectangles)
    length(tags) == 1 && return tags
    fused, _ = gmsh.model.occ.fuse(
        [(2, first(tags))],
        [(2, tag) for tag in Iterators.drop(tags, 1)],
    )
    surfaces = unique(tag for (dimension, tag) in fused if dimension == 2)
    isempty(surfaces) && error("KMK fuse produced no solid domain")
    surfaces
end

function cut_k_slots(surfaces, config)
    lead = config.lead_length_mm
    k_length = config.stiffness.length_mm
    m_length = component_length_mm(config.mass)
    offsets = (lead, lead + k_length + m_length)
    slots = Int[]
    for offset in offsets, box in StiffnessComponentMesher.slot_boxes(config.stiffness)
        x, y, length, height = box
        push!(slots, StiffnessComponentMesher.add_stadium_mm(x + offset, y, length, height))
    end
    cut, _ = gmsh.model.occ.cut(
        [(2, tag) for tag in surfaces],
        [(2, tag) for tag in slots],
    )
    result = unique(tag for (dimension, tag) in cut if dimension == 2)
    isempty(result) && error("KMK slot cut produced no solid domain")
    result
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
    length(source) == 1 || error("KMK must have exactly one left port")
    length(microphone) == 1 || error("KMK must have exactly one right port")
    free = setdiff(curves, vcat(source, microphone))
    source, microphone, free, curves
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
    # Resolve the actual guided-mode profile at both modal ports.
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", vcat(source, microphone))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_port_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_lead_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.05e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 1.2e-3)

    # Uniform refinement around the whole KMK cell, including rounded slots.
    lead_m = config.lead_length_mm * 1e-3
    cell_m = cell_length_mm(config) * 1e-3
    height_m = config.stiffness.height_mm * 1e-3
    gmsh.model.mesh.field.add("Box", 3)
    gmsh.model.mesh.field.setNumber(3, "VIn", Float64(size_cell_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "VOut", Float64(size_lead_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(3, "XMin", lead_m - 0.5e-3)
    gmsh.model.mesh.field.setNumber(3, "XMax", lead_m + cell_m + 0.5e-3)
    gmsh.model.mesh.field.setNumber(3, "YMin", -height_m)
    gmsh.model.mesh.field.setNumber(3, "YMax", height_m)

    gmsh.model.mesh.field.add("Min", 4)
    gmsh.model.mesh.field.setNumbers(4, "FieldsList", [2, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(4)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
end

function build_combined_kmk_mesh(
    output_path::AbstractString;
    config::CombinedKMKConfig=CombinedKMKConfig(),
    variant::Symbol=:combined,
    size_cell_mm::Real=0.055,
    size_port_mm::Real=0.10,
    size_lead_mm::Real=0.30,
)
    validate_config(config)
    variant in (:solid, :k_only, :m_only, :combined) ||
        throw(ArgumentError("unknown KMK variant"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("KMK_$(variant)")
    surfaces = fused_surfaces(solid_rectangles(config; variant))
    variant in (:k_only, :combined) && (surfaces = cut_k_slots(surfaces, config))
    source, microphone, free, _ = classify_boundaries(surfaces, config)
    tag_physical_groups(surfaces, source, microphone, free)
    configure_mesh(
        source,
        microphone,
        config;
        size_cell_mm,
        size_port_mm,
        size_lead_mm,
    )
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    println("[+] $output_path")
    output_path
end

end
