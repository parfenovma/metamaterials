module TopologyMetalGRINMesher

using Gmsh: gmsh

export MetalGRINCellConfig,
       METAL_GRIN_VARIANTS,
       validate_config,
       neck_height_mm,
       local_height_mm,
       total_length_mm,
       minimum_feature_mm,
       minimum_curvature_radius_mm,
       build_metal_grin_cell_mesh

"""
Axisymmetric-ring surrogate cut in a 2-D meridional section.

Each radial zone is an aluminium annulus with a smooth, periodically waisted
axial load path.  In 3-D the two free sides become the inner and outer
cylindrical surfaces of a separately machined ring.  The aligned waists are
not holes: they are circumferential grooves machined from both radial sides.
"""
Base.@kwdef struct MetalGRINCellConfig
    lead_length_mm::Float64 = 5.5
    active_length_mm::Float64 = 44.0
    port_height_mm::Float64 = 3.2
    period_mm::Float64 = 4.0
    period_count::Int = 11
    transition_fraction::Float64 = 0.25
    samples_per_period::Int = 20
    nominal_minimum_feature_mm::Float64 = 0.8
end

const METAL_GRIN_VARIANTS = (
    :solid,
    :neck_2p8,
    :neck_2p4,
    :neck_2p0,
    :neck_1p6,
    :neck_1p2,
    :neck_0p8,
)

const VARIANT_NECK_MM = Dict(
    :solid => 3.2,
    :neck_2p8 => 2.8,
    :neck_2p4 => 2.4,
    :neck_2p0 => 2.0,
    :neck_1p6 => 1.6,
    :neck_1p2 => 1.2,
    :neck_0p8 => 0.8,
)

total_length_mm(config::MetalGRINCellConfig) =
    2config.lead_length_mm + config.active_length_mm

function neck_height_mm(config::MetalGRINCellConfig, variant::Symbol)
    variant in METAL_GRIN_VARIANTS || throw(ArgumentError("unknown metal-GRIN variant: $variant"))
    variant == :solid ? config.port_height_mm : VARIANT_NECK_MM[variant]
end

minimum_feature_mm(config::MetalGRINCellConfig, variant::Symbol) =
    neck_height_mm(config, variant)

"""Smooth periodic height with full-height joins and a waist at cell centre."""
function local_height_mm(config::MetalGRINCellConfig, variant::Symbol, x_mm::Real)
    neck = neck_height_mm(config, variant)
    active_start = config.lead_length_mm
    active_end = active_start + config.active_length_mm
    (variant == :solid || x_mm <= active_start || x_mm >= active_end) &&
        return config.port_height_mm
    xi = mod((Float64(x_mm) - active_start) / config.period_mm, 1.0)
    transition = config.transition_fraction
    if xi < transition
        neck + (config.port_height_mm - neck) * (1 + cospi(xi / transition)) / 2
    elseif xi <= 1 - transition
        neck
    else
        eta = (xi - (1 - transition)) / transition
        neck + (config.port_height_mm - neck) * (1 - cospi(eta)) / 2
    end
end

"""
Analytic lower bound on the radius of curvature at the smooth extrema.
It is a tool-access descriptor, not a fatigue/stress-concentration rating.
"""
function minimum_curvature_radius_mm(config::MetalGRINCellConfig, variant::Symbol)
    reduction = config.port_height_mm - neck_height_mm(config, variant)
    reduction <= 0 && return Inf
    transition_length = config.transition_fraction * config.period_mm
    4transition_length^2 / (pi^2 * reduction)
end

function validate_config(config::MetalGRINCellConfig)
    config.lead_length_mm > 0 || throw(ArgumentError("lead length must be positive"))
    config.active_length_mm > 0 || throw(ArgumentError("active length must be positive"))
    config.port_height_mm > 0 || throw(ArgumentError("port height must be positive"))
    config.period_mm > 0 || throw(ArgumentError("period must be positive"))
    config.period_count >= 2 || throw(ArgumentError("at least two periods are required"))
    0 < config.transition_fraction < 0.5 ||
        throw(ArgumentError("transition fraction must lie in (0, 0.5)"))
    isapprox(
        config.active_length_mm,
        config.period_count * config.period_mm;
        atol=1.0e-12,
        rtol=0.0,
    ) || throw(ArgumentError("active length must equal period_count * period"))
    config.samples_per_period >= 8 ||
        throw(ArgumentError("at least eight boundary samples per period are required"))
    config.nominal_minimum_feature_mm > 0 ||
        throw(ArgumentError("minimum feature must be positive"))
    for variant in METAL_GRIN_VARIANTS
        0 < neck_height_mm(config, variant) <= config.port_height_mm ||
            throw(ArgumentError("invalid neck height for $variant"))
        minimum_feature_mm(config, variant) + 1.0e-12 >=
            config.nominal_minimum_feature_mm ||
            throw(ArgumentError("$variant violates the nominal minimum feature"))
    end
    nothing
end

function boundary_points_mm(config::MetalGRINCellConfig, variant::Symbol)
    count = config.period_count * config.samples_per_period
    active_x = range(
        config.lead_length_mm,
        config.lead_length_mm + config.active_length_mm;
        length=count + 1,
    )
    top = Tuple{Float64, Float64}[(0.0, config.port_height_mm / 2)]
    append!(top, [(Float64(x), local_height_mm(config, variant, x) / 2) for x in active_x])
    push!(top, (total_length_mm(config), config.port_height_mm / 2))
    bottom = [(x, -y) for (x, y) in reverse(top)]
    vcat(top, bottom)
end

function polygon_surface(points_mm)
    tags = [
        gmsh.model.occ.addPoint(x * 1.0e-3, y * 1.0e-3, 0.0)
        for (x, y) in points_mm
    ]
    lines = [
        gmsh.model.occ.addLine(tags[index], tags[index == length(tags) ? 1 : index + 1])
        for index in eachindex(tags)
    ]
    loop = gmsh.model.occ.addCurveLoop(lines)
    gmsh.model.occ.addPlaneSurface([loop])
end

function classify_boundaries(surface, config)
    gmsh.model.occ.synchronize()
    boundary = gmsh.model.getBoundary([(2, surface)], false, false, false)
    curves = unique(tag for (dimension, tag) in boundary if dimension == 1)
    total_m = total_length_mm(config) * 1.0e-3
    tolerance = 1.0e-10
    source = filter(tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] < tolerance, curves)
    microphone = filter(
        tag -> gmsh.model.occ.getCenterOfMass(1, tag)[1] > total_m - tolerance,
        curves,
    )
    length(source) == 1 || error("metal-GRIN cell must have exactly one left port")
    length(microphone) == 1 || error("metal-GRIN cell must have exactly one right port")
    source, microphone, setdiff(curves, vcat(source, microphone))
end

function add_physical_group(dimension, tags, id, name)
    gmsh.model.addPhysicalGroup(dimension, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dimension, id, name)
end

function configure_mesh(source, microphone, free, config; size_waist_mm, size_lead_mm)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", free)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", Float64(size_waist_mm) * 1.0e-3)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", Float64(size_lead_mm) * 1.0e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.05e-3)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.8e-3)

    gmsh.model.mesh.field.add("Box", 3)
    gmsh.model.mesh.field.setNumber(3, "VIn", Float64(size_waist_mm) * 1.0e-3)
    gmsh.model.mesh.field.setNumber(3, "VOut", Float64(size_lead_mm) * 1.0e-3)
    gmsh.model.mesh.field.setNumber(3, "XMin", (config.lead_length_mm - 0.5) * 1.0e-3)
    gmsh.model.mesh.field.setNumber(
        3,
        "XMax",
        (config.lead_length_mm + config.active_length_mm + 0.5) * 1.0e-3,
    )
    gmsh.model.mesh.field.setNumber(3, "YMin", -config.port_height_mm * 1.0e-3)
    gmsh.model.mesh.field.setNumber(3, "YMax", config.port_height_mm * 1.0e-3)

    gmsh.model.mesh.field.add("Min", 4)
    gmsh.model.mesh.field.setNumbers(4, "FieldsList", [2, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(4)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
end

function build_metal_grin_cell_mesh(
    output_path::AbstractString;
    config::MetalGRINCellConfig=MetalGRINCellConfig(),
    variant::Symbol=:neck_1p6,
    size_waist_mm::Real=0.10,
    size_lead_mm::Real=0.30,
)
    validate_config(config)
    variant in METAL_GRIN_VARIANTS || throw(ArgumentError("unknown metal-GRIN variant"))
    mkpath(dirname(output_path))
    gmsh.clear()
    gmsh.model.add("topology_metal_grin_$(variant)")
    surface = polygon_surface(boundary_points_mm(config, variant))
    source, microphone, free = classify_boundaries(surface, config)
    add_physical_group(1, source, 101, "Source")
    add_physical_group(1, microphone, 102, "Microphone")
    add_physical_group(1, free, 103, "FreeSurface")
    add_physical_group(2, [surface], 201, "Domain")
    configure_mesh(source, microphone, free, config; size_waist_mm, size_lead_mm)
    gmsh.model.mesh.generate(2)
    gmsh.write(output_path)
    output_path
end

end
