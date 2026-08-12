module RunRoundedNotchLensFEM

using JLD2

include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
using .SinusoidalMaterialLens
include(joinpath(@__DIR__, "material_lens_mesher.jl"))
using .MaterialLensMesher
include(joinpath(@__DIR__, "material_lens_harmonic_solver.jl"))
using .MaterialLensHarmonicSolver

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_ROOT = get(
    ENV,
    "METAMATERIALS_ROUNDED_NOTCH_DESIGN",
    joinpath(PROJECT_ROOT, "tmp", "rounded_notch_material_lens_242khz"),
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_ROUNDED_NOTCH_FEM_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "rounded_notch_lens_fem_242khz"),
)
const CONFIGURATIONS = (:polymer_direct, :polymer_matched, :aluminium)
const FEM_ORDER = parse(Int, get(ENV, "METAMATERIALS_ROUNDED_NOTCH_FEM_ORDER", "1"))
const FORCE = lowercase(get(ENV, "METAMATERIALS_ROUNDED_NOTCH_FEM_FORCE", "false")) == "true"

function csv_rows(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1) if !isempty(strip(line))]
end

function selected_cells(configuration)
    rows = csv_rows(joinpath(DESIGN_ROOT, "selected_$(configuration).csv"))
    [
        RoundedNotchCell(
            length_mm=parse(Float64, row["length_mm"]),
            minimum_gap_mm=parse(Float64, row["minimum_gap_mm"]),
            notch_count=parse(Int, row["notch_count"]),
            notch_width_mm=parse(Float64, row["notch_width_mm"]),
            start_from_bottom=parse(Bool, row["start_from_bottom"]),
        )
        for row in rows
    ]
end

function matching_thickness_mm(configuration)
    configuration == :polymer_matched || return 0.0
    override = get(ENV, "METAMATERIALS_ROUNDED_NOTCH_MATCHING_THICKNESS_MM", "")
    !isempty(override) && return parse(Float64, override)
    rows = csv_rows(joinpath(DESIGN_ROOT, "rounded_notch_summary.csv"))
    row = only(filter(row -> row["configuration"] == "polymer_matched", rows))
    parse(Float64, row["matching_thickness_mm"])
end

function materials(configuration)
    polymer = photopolymer()
    aluminium = aluminium_6061()
    configuration == :polymer_direct && return polymer, aluminium, nothing
    configuration == :polymer_matched &&
        return polymer, aluminium, geometric_matching_material(polymer, aluminium)
    configuration == :aluminium && return aluminium, aluminium, nothing
    throw(ArgumentError("unknown configuration: $configuration"))
end

function case_cells(configuration, role)
    selected = selected_cells(configuration)
    role == :selected && return selected
    role == :uniform || throw(ArgumentError("role must be selected or uniform"))
    reference = first(selected)
    straight = RoundedNotchCell(
        length_mm=reference.length_mm,
        height_mm=reference.height_mm,
        minimum_gap_mm=reference.height_mm,
        notch_count=1,
        notch_width_mm=min(reference.notch_width_mm, reference.length_mm / 4),
        start_from_bottom=true,
    )
    fill(straight, length(selected))
end

function case_paths(configuration, role)
    stem = "$(configuration)_$(role)"
    (
        mesh=joinpath(OUTPUT_ROOT, "meshes", "$(stem).msh"),
        result=joinpath(OUTPUT_ROOT, "results", "$(stem).jld2"),
    )
end

function build_case_mesh(configuration, role)
    paths = case_paths(configuration, role)
    if isfile(paths.mesh)
        println("[=] Existing mesh: $(paths.mesh)")
        return paths.mesh
    end
    cells = case_cells(configuration, role)
    info = build_material_lens_mesh(
        cells,
        paths.mesh;
        matching_thickness_mm=matching_thickness_mm(configuration),
        config=MaterialLensMeshConfig(
            output_length_mm=105.0,
            lens_mesh_size_mm=FEM_ORDER == 1 ? 0.36 : 0.62,
            output_mesh_size_mm=FEM_ORDER == 1 ? 1.55 : 2.60,
            transition_distance_mm=14.0,
        ),
    )
    println("[+] Mesh: $(paths.mesh)")
    info
end

function solve_case(configuration, role)
    paths = case_paths(configuration, role)
    isfile(paths.result) && !FORCE && return paths.result
    mesh_info = build_case_mesh(configuration, role)
    mesh_path = mesh_info isa AbstractString ? mesh_info : mesh_info.mesh_path
    cells = case_cells(configuration, role)
    thickness_mm = matching_thickness_mm(configuration)
    focus_x_mm = first(cells).length_mm + thickness_mm + 60.0
    lens_material, output_material, matching_material = materials(configuration)

    result = solve_material_lens_harmonic(
        mesh_path,
        lens_material,
        output_material;
        matching_material,
        config=MaterialLensHarmonicConfig(
            frequency_hz=242.0e3,
            focus_x_mm=focus_x_mm,
            element_order=FEM_ORDER,
            quadrature_degree=2FEM_ORDER,
            scan_before_focus_mm=40.0,
            scan_after_focus_mm=30.0,
            scan_transverse_half_width_mm=30.0,
            scan_step_mm=2.0,
        ),
    )
    mkpath(dirname(paths.result))
    jldsave(
        paths.result;
        format_version=1,
        geometry="alternating_rounded_notch",
        configuration=String(configuration),
        role=String(role),
        frequency_hz=242.0e3,
        focus_x_mm,
        matching_thickness_mm=thickness_mm,
        focus_displacement=result.focus_displacement,
        focus_longitudinal_amplitude_m=result.focus_longitudinal_amplitude_m,
        focus_total_amplitude_m=result.focus_total_amplitude_m,
        local_peak_amplitude_m=result.local_peak_amplitude_m,
        local_peak_x_mm=result.local_peak_x_mm,
        local_peak_y_mm=result.local_peak_y_mm,
        local_peak_displacement=result.local_peak_displacement,
        scan_x_mm=result.scan_x_mm,
        scan_y_mm=result.scan_y_mm,
        scan_amplitude_m=result.scan_amplitude_m,
        element_order=FEM_ORDER,
        quadrature_degree=2FEM_ORDER,
        mesh_path=abspath(mesh_path),
    )
    println("[+] $(configuration)/$(role): |ux(focus)|=$(result.focus_longitudinal_amplitude_m) m")
    println("[+] Result: $(paths.result)")
    paths.result
end

function write_summary()
    rows = NamedTuple[]
    for configuration in CONFIGURATIONS
        selected_path = case_paths(configuration, :selected).result
        uniform_path = case_paths(configuration, :uniform).result
        isfile(selected_path) && isfile(uniform_path) || continue
        selected = load(selected_path)
        uniform = load(uniform_path)
        push!(rows, (
            configuration=String(configuration),
            selected_focus_ux_m=selected["focus_longitudinal_amplitude_m"],
            uniform_focus_ux_m=uniform["focus_longitudinal_amplitude_m"],
            fem_focus_gain=selected["focus_longitudinal_amplitude_m"] /
                           uniform["focus_longitudinal_amplitude_m"],
            selected_local_peak_m=selected["local_peak_amplitude_m"],
            local_peak_x_mm=selected["local_peak_x_mm"],
            local_peak_y_mm=selected["local_peak_y_mm"],
            element_order=selected["element_order"],
        ))
    end
    isempty(rows) && return nothing
    path = joinpath(OUTPUT_ROOT, "rounded_notch_fem_summary.csv")
    open(path, "w") do io
        columns = propertynames(first(rows))
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
    println("[+] FEM summary: $path")
    path
end

function parse_case(value)
    parts = split(value, ':')
    length(parts) == 2 || error("case must be CONFIGURATION:ROLE")
    configuration = Symbol(parts[1])
    role = Symbol(parts[2])
    configuration in CONFIGURATIONS || error("unknown configuration: $configuration")
    role in (:selected, :uniform) || error("unknown role: $role")
    configuration, role
end

function main(args=ARGS)
    mkpath(OUTPUT_ROOT)
    if length(args) == 2 && args[1] == "--mesh"
        configuration, role = parse_case(args[2])
        build_case_mesh(configuration, role)
    elseif length(args) == 2 && args[1] == "--solve"
        configuration, role = parse_case(args[2])
        solve_case(configuration, role)
    elseif args == ["--summary"]
        write_summary()
    else
        error("usage: run_rounded_notch_lens_fem.jl --mesh|--solve CONFIGURATION:ROLE, or --summary")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
