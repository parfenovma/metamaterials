ENV["GKSwstype"] = "100"

using Gmsh: gmsh
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PILOT_MESH_ROOT = joinpath(PROJECT_ROOT, "tmp", "pilot_220khz", "meshes")
const DEFAULT_OUTPUT_ROOT = joinpath(
    PROJECT_ROOT,
    "tmp",
    "pilot_220khz",
    "reference_meshes",
)

const REFERENCE_MESH_CASES = [
    (
        id="rectangular",
        title="Rectangular reference",
        filename="mesh_sin_A_0.0_N_2.msh",
    ),
    (
        id="sinusoidal",
        title="Sinusoidal, A=2.5 mm, g=4.5 mm",
        filename="mesh_sin_A_2.5_N_2_STG.msh",
    ),
    (
        id="exponential_k0p5",
        title="Soft exponential, κ=0.5, g=4.5 mm",
        filename="mesh_exp_A_2.5_N_2_K_0.5_STG.msh",
    ),
    (
        id="exponential_k0p5_n1",
        title="Stretched exponential, κ=0.5, N=1, g=4.5 mm",
        filename="mesh_exp_A_2.5_N_1_K_0.5_STG.msh",
    ),
    (
        id="exponential_k1",
        title="Exponential, κ=1, g=4.5 mm",
        filename="mesh_exp_A_2.5_N_2_K_1.0_STG.msh",
    ),
    (
        id="exponential_k3",
        title="Exponential, κ=3, g=4.5 mm",
        filename="mesh_exp_A_2.5_N_2_K_3.0_STG.msh",
    ),
    (
        id="power_p2",
        title="Power, p=2, g=4.5 mm",
        filename="mesh_pow_A_2.5_N_2_P_2.0_STG.msh",
    ),
    (
        id="power_p4",
        title="Power, p=4, g=4.5 mm",
        filename="mesh_pow_A_2.5_N_2_P_4.0_STG.msh",
    ),
    (
        id="rounded_notch",
        title="Rounded U-notch, N=4, g=4.5 mm",
        filename="mesh_notch_A_2.5_N_4_W_1.4_B_STG.msh",
    ),
    (
        id="legacy",
        title="Legacy sinusoidal, A=2.5 mm",
        filename="mesh_A_2.5.msh",
    ),
]

function triangle_segments(mesh_path)
    gmsh.clear()
    gmsh.open(mesh_path)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    nodes = Dict(
        Int(tag) => (
            coordinates[3index - 2] * 1.0e3,
            coordinates[3index - 1] * 1.0e3,
        )
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
end

function mesh_panel(x, y; title)
    plot(
        x,
        y;
        color=:steelblue,
        linewidth=0.20,
        alpha=0.68,
        label=false,
        aspect_ratio=:equal,
        xlabel="x, mm",
        ylabel="y, mm",
        title,
        xlims=(-0.25, 17.25),
        ylims=(-0.25, 7.25),
        grid=false,
        background_color=:white,
        foreground_color=:black,
        titlefontsize=10,
        guidefontsize=9,
        tickfontsize=7,
        margin=3Plots.mm,
    )
end

function plot_reference_meshes(; output_root::AbstractString=DEFAULT_OUTPUT_ROOT)
    mkpath(output_root)
    panels = Any[]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for case in REFERENCE_MESH_CASES
            mesh_path = joinpath(PILOT_MESH_ROOT, case.filename)
            isfile(mesh_path) || error("reference mesh not found: $mesh_path")
            x, y = triangle_segments(mesh_path)
            panel = mesh_panel(x, y; title=case.title)
            push!(panels, panel)
            individual_path = joinpath(output_root, "reference_mesh_$(case.id).png")
            savefig(plot(panel; size=(1200, 520)), individual_path)
            println("[+] $individual_path")
        end
    finally
        gmsh.finalize()
    end

    gallery = plot(
        panels...;
        layout=(ceil(Int, length(panels) / 2), 2),
        size=(1600, 310 * ceil(Int, length(panels) / 2)),
        plot_title="Reference FEM meshes: staggered pilot wall topologies",
    )
    gallery_path = joinpath(output_root, "reference_meshes_all_topologies.png")
    savefig(gallery, gallery_path)
    println("[+] $gallery_path")
    gallery_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_root = isempty(ARGS) ? DEFAULT_OUTPUT_ROOT : ARGS[1]
    plot_reference_meshes(; output_root)
end
