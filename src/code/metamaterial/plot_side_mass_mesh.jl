ENV["GKSwstype"] = "100"

using Gmsh: gmsh
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_MESH_PATH = joinpath(
    PROJECT_ROOT,
    "tmp",
    "side_mass_refined_lossless",
    "meshes",
    "side_mass_Tb_0p7",
    "mesh_cell.msh",
)
const DEFAULT_OUTPUT_PATH = joinpath(
    PROJECT_ROOT,
    "tmp",
    "side_mass_refined_lossless",
    "side_mass_actual_mesh.png",
)

function triangle_segments(mesh_path)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        gmsh.open(mesh_path)
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        nodes = Dict(
            Int(tag) => (coordinates[3index - 2] * 1e3, coordinates[3index - 1] * 1e3)
            for (index, tag) in enumerate(node_tags)
        )
        x_segments = Float64[]
        y_segments = Float64[]
        element_types, _, element_nodes = gmsh.model.mesh.getElements(2)
        for (element_type, flattened_nodes) in zip(element_types, element_nodes)
            _, _, _, nodes_per_element, _, primary_nodes =
                gmsh.model.mesh.getElementProperties(element_type)
            primary_nodes >= 3 || continue
            for offset in 1:nodes_per_element:length(flattened_nodes)
                corners = Int.(flattened_nodes[offset:(offset + 2)])
                for (left, right) in ((1, 2), (2, 3), (3, 1))
                    x1, y1 = nodes[corners[left]]
                    x2, y2 = nodes[corners[right]]
                    append!(x_segments, (x1, x2, NaN))
                    append!(y_segments, (y1, y2, NaN))
                end
            end
        end
        x_segments, y_segments
    finally
        gmsh.finalize()
    end
end

function mesh_panel(x, y; title, xlims)
    plot(
        x,
        y;
        color=:steelblue,
        linewidth=0.22,
        alpha=0.62,
        label=false,
        aspect_ratio=:equal,
        xlabel="x, mm",
        ylabel="y, mm",
        title,
        xlims,
        ylims=(-0.25, 7.25),
        grid=false,
        background_color=:white,
        foreground_color=:black,
    )
end

function main(args=ARGS)
    mesh_path = isempty(args) ? DEFAULT_MESH_PATH : args[1]
    output_path = length(args) >= 2 ? args[2] : DEFAULT_OUTPUT_PATH
    x, y = triangle_segments(mesh_path)
    full = mesh_panel(x, y; title="Actual FEM mesh: full cell", xlims=(-0.5, 30.5))
    zoom = mesh_panel(x, y; title="Mesh near side resonators", xlims=(7.5, 18.5))
    figure = plot(full, zoom; layout=(2, 1), size=(1200, 800), margin=5Plots.mm)
    mkpath(dirname(output_path))
    savefig(figure, output_path)
    println("[+] $output_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
