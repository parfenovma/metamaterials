using Gridap
using GridapGmsh

function model_filename(mesh_path::AbstractString)
    stem = splitext(basename(mesh_path))[1]
    startswith(stem, "mesh_") || error("expected a mesh_*.msh file, got $mesh_path")
    "model_$(stem[6:end]).json"
end

"""
Convert a Gmsh mesh to Gridap's native JSON representation.

This command intentionally runs separately from the FEM solver. On macOS,
`gmsh_jll` and Homebrew Julia can load different LLVM OpenMP runtimes; keeping
Gmsh out of the solver process avoids that conflict without unsafe overrides.
"""
function convert_mesh(
    mesh_path::AbstractString;
    output_dir::AbstractString="1_models",
)
    isfile(mesh_path) || error("mesh not found: $mesh_path")
    mkpath(output_dir)
    model = GmshDiscreteModel(mesh_path)
    output_path = joinpath(output_dir, model_filename(mesh_path))
    Gridap.Io.to_json_file(model, output_path)
    println("  [+] Gridap model saved: $output_path")
    output_path
end

function convert_directory(
    mesh_dir::AbstractString="1_meshes";
    output_dir::AbstractString="1_models",
)
    meshes = sort(collect(
        joinpath(mesh_dir, name)
        for name in readdir(mesh_dir)
        if startswith(name, "mesh_") && endswith(name, ".msh")
    ))
    isempty(meshes) && error("no mesh_*.msh files found in $mesh_dir")
    [convert_mesh(path; output_dir) for path in meshes]
end

function main(args=ARGS)
    if isempty(args)
        convert_directory()
    elseif isdir(args[1])
        output_dir = length(args) >= 2 ? args[2] : "1_models"
        convert_directory(args[1]; output_dir)
    else
        output_dir = length(args) >= 2 ? args[2] : "1_models"
        convert_mesh(args[1]; output_dir)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
