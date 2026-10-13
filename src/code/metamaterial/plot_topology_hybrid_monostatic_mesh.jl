module PlotTopologyHybridMonostaticMesh

using Gmsh: gmsh
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_MESH = joinpath(
    PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_monostatic_g3",
    "meshes", "mesh_refined_hybrid_void_r2_focus.msh",
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HYBRID_ROBUSTNESS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_robustness_sweep"),
)

function triangle_segments(entity_tags, node_index, x, y; stride=1)
    xs, ys = Float64[], Float64[]
    counter = 0
    for entity in entity_tags
        element_types, _, connectivity = gmsh.model.mesh.getElements(2, entity)
        for (element_type, nodes) in zip(element_types, connectivity)
            name, _, _, node_count, _, _ = gmsh.model.mesh.getElementProperties(element_type)
            startswith(name, "Triangle") || continue
            for offset in 1:node_count:length(nodes)
                counter += 1
                counter % stride == 0 || continue
                corners = nodes[offset:(offset + 2)]
                indices = getindex.(Ref(node_index), Int.(corners))
                append!(xs, (x[indices[1]], x[indices[2]], x[indices[3]], x[indices[1]], NaN))
                append!(ys, (y[indices[1]], y[indices[2]], y[indices[3]], y[indices[1]], NaN))
            end
        end
    end
    xs, ys
end

function material_groups()
    groups = NamedTuple[]
    for (dimension, physical_tag) in gmsh.model.getPhysicalGroups(2)
        name = gmsh.model.getPhysicalName(dimension, physical_tag)
        startswith(name, "MaterialF") || continue
        key = parse(Int, replace(name, "MaterialF" => ""))
        entities = Int.(gmsh.model.getEntitiesForPhysicalGroup(dimension, physical_tag))
        push!(groups, (; key, name, entities))
    end
    sort(groups; by=group -> group.key)
end

function mesh_panel(groups, node_index, x, y; xlims, ylims, title, stride=1)
    panel = plot(
        ; aspect_ratio=:equal, xlims, ylims, title, xlabel="x, mm", ylabel="y, mm",
        grid=false, legend=false, background_color=:white,
    )
    colors = cgrad(:thermal, max(2, length(groups)), categorical=true)
    for (index, group) in enumerate(groups)
        xs, ys = triangle_segments(group.entities, node_index, x, y; stride)
        plot!(panel, xs, ys; color=colors[index], linewidth=0.12, alpha=0.82, label=false)
    end
    panel
end

function plot_mesh(mesh_path=DEFAULT_MESH)
    isfile(mesh_path) || error("mesh missing: $mesh_path")
    mkpath(OUTPUT_ROOT)
    ENV["GKSwstype"] = "100"
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.open(mesh_path)
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        node_index = Dict(Int(tag) => index for (index, tag) in enumerate(node_tags))
        x = coordinates[1:3:end] .* 1e3
        y = coordinates[2:3:end] .* 1e3
        groups = material_groups()
        panels = [
            mesh_panel(groups, node_index, x, y; xlims=(-57, 126), ylims=(-38, 38), title="Refined hybrid + Al specimen + r=2 mm void", stride=3),
            mesh_panel(groups, node_index, x, y; xlims=(-56, 4), ylims=(-37, 37), title="15-channel hybrid device", stride=1),
            mesh_panel(groups, node_index, x, y; xlims=(-37, -18), ylims=(-13, 13), title="Central Al-Cu graded inserts", stride=1),
            mesh_panel(groups, node_index, x, y; xlims=(54, 66), ylims=(-6, 6), title="Refined mesh around explicit void", stride=1),
        ]
        figure = plot(
            panels...; layout=(2, 2), size=(1700, 1000), margin=3Plots.mm,
            plot_title="Actual Gmsh mesh used for refined monostatic G3",
        )
        path = joinpath(OUTPUT_ROOT, "hybrid_refined_g3_mesh.png")
        savefig(figure, path)
        println("[+] $path")
        path
    finally
        gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_mesh()
end

end
