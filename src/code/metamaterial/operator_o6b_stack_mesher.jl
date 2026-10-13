module OperatorO6BStackMesher

using Gmsh: gmsh

export O6BStackGeometry, build_o6b_stack_mesh

Base.@kwdef struct O6BStackGeometry
    aperture_height_mm::Float64 = 71.6
    pzt_thickness_mm::Float64 = 6.0
    bond_thickness_mm::Float64 = 0.10
    aluminium_length_mm::Float64 = 35.0
    pzt_mesh_mm::Float64 = 0.60
    bond_mesh_mm::Float64 = 0.10
    aluminium_mesh_mm::Float64 = 1.00
end

function validate(config::O6BStackGeometry)
    config.aperture_height_mm > 0 || throw(ArgumentError("aperture must be positive"))
    config.pzt_thickness_mm > 0 || throw(ArgumentError("PZT thickness must be positive"))
    config.bond_thickness_mm > 0 || throw(ArgumentError("bond thickness must be positive"))
    config.aluminium_length_mm > 0 || throw(ArgumentError("Al load must be positive"))
    config.pzt_mesh_mm > 0 || throw(ArgumentError("PZT mesh must be positive"))
    config.bond_mesh_mm > 0 || throw(ArgumentError("bond mesh must be positive"))
    config.aluminium_mesh_mm > 0 || throw(ArgumentError("Al mesh must be positive"))
    nothing
end

function add_group(dimension, tags, id, name)
    isempty(tags) && return
    gmsh.model.addPhysicalGroup(dimension, sort(unique(tags)), id)
    gmsh.model.setPhysicalName(dimension, id, name)
end

function build_o6b_stack_mesh(
    output_path::AbstractString;
    config::O6BStackGeometry=O6BStackGeometry(),
)
    validate(config)
    mkpath(dirname(output_path))
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.model.add("operator_o6b_pzt_bond_al")
        half_height = 0.5config.aperture_height_mm * 1e-3
        pzt_back_x = -(config.pzt_thickness_mm + config.bond_thickness_mm) * 1e-3
        bond_back_x = -config.bond_thickness_mm * 1e-3
        aluminium_end_x = config.aluminium_length_mm * 1e-3
        pzt = gmsh.model.occ.addRectangle(
            pzt_back_x, -half_height, 0.0,
            config.pzt_thickness_mm * 1e-3, 2half_height,
        )
        bond = gmsh.model.occ.addRectangle(
            bond_back_x, -half_height, 0.0,
            config.bond_thickness_mm * 1e-3, 2half_height,
        )
        aluminium = gmsh.model.occ.addRectangle(
            0.0, -half_height, 0.0, aluminium_end_x, 2half_height,
        )
        gmsh.model.occ.fragment([(2, pzt), (2, bond)], [(2, aluminium)])
        gmsh.model.occ.synchronize()

        pzt_surfaces, bond_surfaces, aluminium_surfaces = Int[], Int[], Int[]
        for (_, tag) in gmsh.model.getEntities(2)
            x = gmsh.model.occ.getCenterOfMass(2, tag)[1]
            if x < bond_back_x - 1e-9
                push!(pzt_surfaces, tag)
            elseif x < -1e-9
                push!(bond_surfaces, tag)
            else
                push!(aluminium_surfaces, tag)
            end
        end
        length(pzt_surfaces) == 1 || error("expected one PZT surface")
        length(bond_surfaces) == 1 || error("expected one bond surface")
        length(aluminium_surfaces) == 1 || error("expected one Al surface")

        source, backing, radiation, free, contact = Int[], Int[], Int[], Int[], Int[]
        tolerance = 2e-8
        for (_, tag) in gmsh.model.getEntities(1)
            upward, downward = gmsh.model.getAdjacencies(1, tag)
            center = gmsh.model.occ.getCenterOfMass(1, tag)
            x, y = center[1], center[2]
            if length(upward) == 2 && abs(x) < tolerance
                push!(contact, tag)
            elseif length(upward) == 1
                if abs(x - pzt_back_x) < tolerance
                    push!(source, tag)
                    push!(backing, tag)
                elseif abs(x - aluminium_end_x) < tolerance
                    push!(radiation, tag)
                elseif abs(abs(y) - half_height) < tolerance
                    if x > tolerance
                        push!(radiation, tag)
                    else
                        push!(free, tag)
                    end
                else
                    push!(free, tag)
                end
            end
        end
        length(source) == 1 || error("expected one PZT back source")
        length(contact) == 1 || error("expected one bond/Al contact interface")
        length(radiation) == 3 || error("expected Al far/top/bottom radiation boundaries")

        add_group(2, pzt_surfaces, 201, "PZT")
        add_group(2, bond_surfaces, 202, "Bond")
        add_group(2, aluminium_surfaces, 203, "Aluminium")
        add_group(1, source, 101, "Source")
        add_group(1, backing, 102, "BackingBoundary")
        add_group(1, radiation, 103, "RadiationBoundary")
        add_group(1, free, 104, "FreeSurface")
        add_group(1, contact, 105, "ContactInterface")

        gmsh.model.mesh.field.add("Box", 1)
        gmsh.model.mesh.field.setNumber(1, "VIn", config.pzt_mesh_mm * 1e-3)
        gmsh.model.mesh.field.setNumber(1, "VOut", config.aluminium_mesh_mm * 1e-3)
        gmsh.model.mesh.field.setNumber(1, "XMin", pzt_back_x - 1e-6)
        gmsh.model.mesh.field.setNumber(1, "XMax", bond_back_x + 1e-6)
        gmsh.model.mesh.field.setNumber(1, "YMin", -half_height - 1e-6)
        gmsh.model.mesh.field.setNumber(1, "YMax", half_height + 1e-6)

        gmsh.model.mesh.field.add("Box", 2)
        gmsh.model.mesh.field.setNumber(2, "VIn", config.bond_mesh_mm * 1e-3)
        gmsh.model.mesh.field.setNumber(2, "VOut", config.aluminium_mesh_mm * 1e-3)
        gmsh.model.mesh.field.setNumber(2, "XMin", bond_back_x - 1e-6)
        gmsh.model.mesh.field.setNumber(2, "XMax", 1e-6)
        gmsh.model.mesh.field.setNumber(2, "YMin", -half_height - 1e-6)
        gmsh.model.mesh.field.setNumber(2, "YMax", half_height + 1e-6)

        gmsh.model.mesh.field.add("Box", 3)
        gmsh.model.mesh.field.setNumber(3, "VIn", 0.75e-3)
        gmsh.model.mesh.field.setNumber(3, "VOut", config.aluminium_mesh_mm * 1e-3)
        gmsh.model.mesh.field.setNumber(3, "XMin", -1e-6)
        gmsh.model.mesh.field.setNumber(3, "XMax", 5e-3)
        gmsh.model.mesh.field.setNumber(3, "YMin", -half_height - 1e-6)
        gmsh.model.mesh.field.setNumber(3, "YMax", half_height + 1e-6)

        gmsh.model.mesh.field.add("Min", 4)
        gmsh.model.mesh.field.setNumbers(4, "FieldsList", [1, 2, 3])
        gmsh.model.mesh.field.setAsBackgroundMesh(4)
        gmsh.option.setNumber("Mesh.MeshSizeMin", config.bond_mesh_mm * 1e-3)
        gmsh.option.setNumber("Mesh.MeshSizeMax", config.aluminium_mesh_mm * 1e-3)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        gmsh.write(output_path)
        node_count = length(first(gmsh.model.mesh.getNodes()))
        _, tags, _ = gmsh.model.mesh.getElements(2)
        element_count = sum(length, tags)
        (
            output_path=abspath(output_path),
            node_count,
            element_count,
            pzt_back_x_m=pzt_back_x,
            contact_x_m=0.0,
            aluminium_end_x_m=aluminium_end_x,
        )
    finally
        gmsh.finalize()
    end
end

end
