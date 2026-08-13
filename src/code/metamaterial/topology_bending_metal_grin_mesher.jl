module TopologyBendingMetalGRINMesher

using Gmsh: gmsh

export BendingMetalGRINConfig,
       BENDING_METAL_GRIN_VARIANTS,
       validate_config,
       variant_parameters,
       total_length_mm,
       minimum_feature_mm,
       relative_areal_density,
       build_bending_metal_grin_mesh

"""
2-D plane-strain surrogate of a stacked, wire-EDM aluminium cell.

The propagating direction is `x`.  Every cell contains a full-height boundary
node, a retained central inertia block, and four re-entrant two-segment
ligaments with corner masses.  Axial compliance is controlled primarily by
ligament thickness and re-entrant offset; retained inertia is controlled by
the central/corner blocks.  It is therefore not the periodically waisted rod
connectivity rejected by the earlier metal-GRIN branch.
"""
Base.@kwdef struct BendingMetalGRINConfig
    port_height_mm::Float64 = 7.0
    period_mm::Float64 = 5.0
    period_count::Int = 9
    lead_length_mm::Float64 = 5.0
    boundary_node_width_mm::Float64 = 0.9
    plate_stack_thickness_mm::Float64 = 1.5
    nominal_minimum_feature_mm::Float64 = 0.8
end

const BENDING_METAL_GRIN_VARIANTS = (
    :solid,
    :dense_a,
    :dense_b,
    :dense_c,
    :bd_a,
    :bd_b,
    :bd_c,
    :bd_d,
    :bd_e,
    :bd_f,
    :bd_g,
)

# Two independent design controls are intentionally not collapsed into one
# porosity parameter: ligament thickness controls flexural stiffness, while
# block width/height and corner size control inertia.
const PARAMETERS = Dict(
    :dense_a => (ligament_mm=2.00, mass_width_mm=2.30, mass_height_mm=4.20,
                 corner_mass_mm=1.20, elbow_y_mm=1.80),
    :dense_b => (ligament_mm=2.35, mass_width_mm=2.60, mass_height_mm=4.80,
                 corner_mass_mm=1.45, elbow_y_mm=1.95),
    :dense_c => (ligament_mm=2.70, mass_width_mm=2.90, mass_height_mm=5.40,
                 corner_mass_mm=1.70, elbow_y_mm=2.10),
    :bd_a => (ligament_mm=1.70, mass_width_mm=1.70, mass_height_mm=2.80,
              corner_mass_mm=0.90, elbow_y_mm=1.70),
    :bd_b => (ligament_mm=1.50, mass_width_mm=1.90, mass_height_mm=3.10,
              corner_mass_mm=1.00, elbow_y_mm=1.55),
    :bd_c => (ligament_mm=1.30, mass_width_mm=2.10, mass_height_mm=3.40,
              corner_mass_mm=1.10, elbow_y_mm=1.40),
    :bd_d => (ligament_mm=1.15, mass_width_mm=2.30, mass_height_mm=3.70,
              corner_mass_mm=1.20, elbow_y_mm=1.25),
    :bd_e => (ligament_mm=1.00, mass_width_mm=2.50, mass_height_mm=4.00,
              corner_mass_mm=1.30, elbow_y_mm=1.10),
    :bd_f => (ligament_mm=0.90, mass_width_mm=2.65, mass_height_mm=4.20,
              corner_mass_mm=1.40, elbow_y_mm=0.95),
    :bd_g => (ligament_mm=0.80, mass_width_mm=2.80, mass_height_mm=4.40,
              corner_mass_mm=1.50, elbow_y_mm=0.80),
)

variant_parameters(::BendingMetalGRINConfig, variant::Symbol) =
    variant == :solid ?
    (ligament_mm=7.0, mass_width_mm=5.0, mass_height_mm=7.0,
     corner_mass_mm=7.0, elbow_y_mm=0.0) :
    get(PARAMETERS, variant) do
        throw(ArgumentError("unknown bending-metal-GRIN variant: $variant"))
    end

total_length_mm(config::BendingMetalGRINConfig) =
    2config.lead_length_mm + config.period_count * config.period_mm

minimum_feature_mm(config::BendingMetalGRINConfig, variant::Symbol) =
    variant == :solid ? config.port_height_mm : min(
        config.boundary_node_width_mm,
        variant_parameters(config, variant).ligament_mm,
        variant_parameters(config, variant).corner_mass_mm,
    )

function validate_config(config::BendingMetalGRINConfig)
    config.port_height_mm > 0 || throw(ArgumentError("port height must be positive"))
    config.period_mm > 0 || throw(ArgumentError("period must be positive"))
    config.period_count >= 1 || throw(ArgumentError("at least one cell is required"))
    config.lead_length_mm > 0 || throw(ArgumentError("lead length must be positive"))
    config.boundary_node_width_mm >= config.nominal_minimum_feature_mm ||
        throw(ArgumentError("boundary node violates minimum feature"))
    for variant in BENDING_METAL_GRIN_VARIANTS
        minimum_feature_mm(config, variant) + 1e-12 >= config.nominal_minimum_feature_mm ||
            throw(ArgumentError("$variant violates the nominal minimum feature"))
        variant == :solid && continue
        p = variant_parameters(config, variant)
        p.mass_width_mm < config.period_mm - config.boundary_node_width_mm ||
            throw(ArgumentError("$variant mass overlaps the boundary node"))
        p.mass_height_mm < config.port_height_mm ||
            throw(ArgumentError("$variant mass height exceeds the port"))
        0 < p.elbow_y_mm < config.port_height_mm / 2 ||
            throw(ArgumentError("$variant has invalid elbow position"))
    end
    nothing
end

function add_centered_rectangle(cx_mm, cy_mm, width_mm, height_mm)
    gmsh.model.occ.addRectangle(
        (cx_mm - width_mm / 2) * 1e-3,
        (cy_mm - height_mm / 2) * 1e-3,
        0.0,
        width_mm * 1e-3,
        height_mm * 1e-3,
    )
end

function add_beam(p1, p2, thickness_mm)
    dx = p2[1] - p1[1]
    dy = p2[2] - p1[2]
    length_mm = hypot(dx, dy)
    angle = atan(dy, dx)
    # A small axial overlap prevents tolerance-sized disconnected joints.
    overlap = 0.12 * thickness_mm
    tag = gmsh.model.occ.addRectangle(
        (p1[1] - overlap) * 1e-3,
        (p1[2] - thickness_mm / 2) * 1e-3,
        0.0,
        (length_mm + 2overlap) * 1e-3,
        thickness_mm * 1e-3,
    )
    gmsh.model.occ.rotate([(2, tag)], p1[1] * 1e-3, p1[2] * 1e-3, 0.0, 0.0, 0.0, 1.0, angle)
    tag
end

function cell_surfaces(config, variant, cell_index)
    p = variant_parameters(config, variant)
    a = config.period_mm
    active_start = config.lead_length_mm
    xl = active_start + (cell_index - 1) * a
    xr = xl + a
    xc = (xl + xr) / 2
    nw = config.boundary_node_width_mm
    surfaces = Int[]

    # Full-height boundary nodes are inertia/matching nodes, not axial rails.
    push!(surfaces, add_centered_rectangle(xl, 0.0, nw, config.port_height_mm))
    push!(surfaces, add_centered_rectangle(xr, 0.0, nw, config.port_height_mm))
    push!(surfaces, add_centered_rectangle(xc, 0.0, p.mass_width_mm, p.mass_height_mm))

    x_elbow_left = xl + 0.29a
    x_elbow_right = xr - 0.29a
    y_port = config.port_height_mm / 2 - p.ligament_mm / 2
    y_mass = p.mass_height_mm / 2 - p.ligament_mm / 2
    for sign in (-1.0, 1.0)
        left_port = (xl + nw / 2 - 0.05p.ligament_mm, sign * y_port)
        left_elbow = (x_elbow_left, sign * p.elbow_y_mm)
        left_mass = (xc - p.mass_width_mm / 2 + 0.05p.ligament_mm, sign * y_mass)
        right_mass = (xc + p.mass_width_mm / 2 - 0.05p.ligament_mm, sign * y_mass)
        right_elbow = (x_elbow_right, sign * p.elbow_y_mm)
        right_port = (xr - nw / 2 + 0.05p.ligament_mm, sign * y_port)
        append!(surfaces, (
            add_beam(left_port, left_elbow, p.ligament_mm),
            add_beam(left_elbow, left_mass, p.ligament_mm),
            add_beam(right_mass, right_elbow, p.ligament_mm),
            add_beam(right_elbow, right_port, p.ligament_mm),
            add_centered_rectangle(left_elbow[1], left_elbow[2], p.corner_mass_mm, p.corner_mass_mm),
            add_centered_rectangle(right_elbow[1], right_elbow[2], p.corner_mass_mm, p.corner_mass_mm),
        ))
    end
    surfaces
end

function fuse_surfaces(tags)
    result, _ = gmsh.model.occ.fuse([(2, first(tags))], [(2, tag) for tag in Iterators.drop(tags, 1)])
    surfaces = unique(tag for (dim, tag) in result if dim == 2)
    isempty(surfaces) && error("boolean union returned no surface")
    surfaces
end

function construct_geometry(config, variant)
    if variant == :solid
        return [gmsh.model.occ.addRectangle(
            0.0, -config.port_height_mm / 2 * 1e-3, 0.0,
            total_length_mm(config) * 1e-3, config.port_height_mm * 1e-3,
        )]
    end
    tags = Int[]
    # Solid launch sections overlap the first and final lattice nodes.
    push!(tags, gmsh.model.occ.addRectangle(
        0.0, -config.port_height_mm / 2 * 1e-3, 0.0,
        (config.lead_length_mm + config.boundary_node_width_mm / 2) * 1e-3,
        config.port_height_mm * 1e-3,
    ))
    push!(tags, gmsh.model.occ.addRectangle(
        (total_length_mm(config) - config.lead_length_mm - config.boundary_node_width_mm / 2) * 1e-3,
        -config.port_height_mm / 2 * 1e-3, 0.0,
        (config.lead_length_mm + config.boundary_node_width_mm / 2) * 1e-3,
        config.port_height_mm * 1e-3,
    ))
    for cell_index in 1:config.period_count
        append!(tags, cell_surfaces(config, variant, cell_index))
    end
    fuse_surfaces(tags)
end

function classify_boundaries(surfaces, config)
    gmsh.model.occ.synchronize()
    boundary = gmsh.model.getBoundary([(2, tag) for tag in surfaces], false, false, false)
    curves = unique(tag for (dim, tag) in boundary if dim == 1)
    total_m = total_length_mm(config) * 1e-3
    tol = 1e-10
    source = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] < tol, curves)
    microphone = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] > total_m - tol, curves)
    length(source) == 1 || error("expected one left port, found $(length(source))")
    length(microphone) == 1 || error("expected one right port, found $(length(microphone))")
    source, microphone, setdiff(curves, vcat(source, microphone))
end

function add_group(dim, tags, id, name)
    gmsh.model.addPhysicalGroup(dim, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dim, id, name)
end

function build_bending_metal_grin_mesh(
    output_path::AbstractString;
    config::BendingMetalGRINConfig=BendingMetalGRINConfig(),
    variant::Symbol=:bd_d,
    size_ligament_mm::Real=0.16,
    size_bulk_mm::Real=0.34,
)
    validate_config(config)
    variant in BENDING_METAL_GRIN_VARIANTS || throw(ArgumentError("unknown variant"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("topology_bending_metal_grin_$(variant)")
    surfaces = construct_geometry(config, variant)
    source, microphone, free = classify_boundaries(surfaces, config)
    add_group(1, source, 101, "Source")
    add_group(1, microphone, 102, "Microphone")
    add_group(1, free, 103, "FreeSurface")
    add_group(2, surfaces, 201, "Domain")

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", free)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_ligament_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_bulk_mm) * 1e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.10e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.75e-3)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    output_path
end

function relative_areal_density(config::BendingMetalGRINConfig, variant::Symbol)
    variant == :solid && return 1.0
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.clear()
        gmsh.model.add("density_$(variant)")
        surfaces = construct_geometry(config, variant)
        gmsh.model.occ.synchronize()
        area_m2 = sum(gmsh.model.occ.getMass(2, tag) for tag in surfaces)
        area_m2 / (total_length_mm(config) * config.port_height_mm * 1e-6)
    finally
        gmsh.finalize()
    end
end

end
