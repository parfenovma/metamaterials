module DistributedSlowWaveMesher

using Gmsh: gmsh

include(joinpath(@__DIR__, "stiffness_component_mesher.jl"))
using .StiffnessComponentMesher

export DistributedSlowWaveConfig,
       validate_config,
       active_length_mm,
       total_length_mm,
       section_starts_mm,
       active_section_indices,
       slot_boxes,
       build_distributed_slow_wave_mesh

"""Symmetric series of weak, longitudinally slotted slow-wave sections."""
Base.@kwdef struct DistributedSlowWaveConfig
    lead_length_mm::Float64 = 12.0
    port_height_mm::Float64 = 3.2
    section_count::Int = 4
    section_length_mm::Float64 = 2.0
    spacer_length_mm::Float64 = 0.5
    end_wall_mm::Float64 = 0.35
    ligament_height_mm::Float64 = 0.45
    minimum_feature_mm::Float64 = 0.35
    minimum_radius_mm::Float64 = 0.20
end

active_length_mm(config::DistributedSlowWaveConfig) =
    config.section_count * config.section_length_mm +
    (config.section_count - 1) * config.spacer_length_mm

total_length_mm(config::DistributedSlowWaveConfig) =
    2config.lead_length_mm + active_length_mm(config)

section_starts_mm(config::DistributedSlowWaveConfig) = [
    config.lead_length_mm +
    (index - 1) * (config.section_length_mm + config.spacer_length_mm)
    for index in 1:config.section_count
]

function local_stiffness_config(config::DistributedSlowWaveConfig)
    StiffnessConfig(
        length_mm=config.section_length_mm,
        height_mm=config.port_height_mm,
        end_wall_mm=config.end_wall_mm,
        ligament_height_mm=config.ligament_height_mm,
        minimum_feature_mm=config.minimum_feature_mm,
        minimum_radius_mm=config.minimum_radius_mm,
    )
end

function validate_config(config::DistributedSlowWaveConfig)
    config.lead_length_mm > 0 || throw(ArgumentError("lead length must be positive"))
    config.section_count >= 2 || throw(ArgumentError("at least two sections are required"))
    iseven(config.section_count) ||
        throw(ArgumentError("section count must be even for the three-state pilot"))
    config.spacer_length_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("solid spacer violates the minimum printable feature"))
    StiffnessComponentMesher.validate_config(local_stiffness_config(config))
    nothing
end

function active_section_indices(config::DistributedSlowWaveConfig, variant::Symbol)
    validate_config(config)
    variant in (:solid, :mid, :max) ||
        throw(ArgumentError("variant must be :solid, :mid, or :max"))
    variant == :solid && return Int[]
    variant == :max && return collect(1:config.section_count)
    middle = config.section_count ÷ 2
    [middle, middle + 1]
end

"""Global slot boxes `(x, y, length, height)` for a pilot state, in mm."""
function slot_boxes(config::DistributedSlowWaveConfig, variant::Symbol)
    local_boxes = StiffnessComponentMesher.slot_boxes(local_stiffness_config(config))
    starts = section_starts_mm(config)
    [
        (starts[index] + x, y, length, height)
        for index in active_section_indices(config, variant)
        for (x, y, length, height) in local_boxes
    ]
end

function build_domain(config, variant)
    base = gmsh.model.occ.addRectangle(
        0.0,
        -config.port_height_mm * 0.5e-3,
        0.0,
        total_length_mm(config) * 1e-3,
        config.port_height_mm * 1e-3,
    )
    variant == :solid && return [base]
    slots = [StiffnessComponentMesher.add_stadium_mm(box...) for box in slot_boxes(config, variant)]
    cut, _ = gmsh.model.occ.cut([(2, base)], [(2, tag) for tag in slots])
    surfaces = unique(tag for (dimension, tag) in cut if dimension == 2)
    isempty(surfaces) && error("slow-wave slot cut produced no solid domain")
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
    length(source) == 1 || error("slow-wave cell must have exactly one left port")
    length(microphone) == 1 || error("slow-wave cell must have exactly one right port")
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
    active_m = active_length_mm(config) * 1e-3
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

function build_distributed_slow_wave_mesh(
    output_path::AbstractString;
    config::DistributedSlowWaveConfig=DistributedSlowWaveConfig(),
    variant::Symbol=:max,
    size_cell_mm::Real=0.055,
    size_port_mm::Real=0.08,
    size_lead_mm::Real=0.25,
)
    validate_config(config)
    variant in (:solid, :mid, :max) || throw(ArgumentError("unknown slow-wave variant"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("distributed_slow_wave_$(variant)")
    surfaces = build_domain(config, variant)
    source, microphone, free = classify_boundaries(surfaces, config)
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
