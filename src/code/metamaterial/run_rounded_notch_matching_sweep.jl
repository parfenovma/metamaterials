module RunRoundedNotchMatchingSweep

using JLD2

include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
using .SinusoidalMaterialLens
include(joinpath(@__DIR__, "material_lens_mesher.jl"))
using .MaterialLensMesher
include(joinpath(@__DIR__, "material_lens_harmonic_solver.jl"))
using .MaterialLensHarmonicSolver

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_ROOT = joinpath(PROJECT_ROOT, "tmp", "rounded_notch_material_lens_242khz")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "rounded_notch_matching_sweep_242khz")

label(value) = replace(string(Float64(value)), "." => "p")

function selected_cells()
    lines = readlines(joinpath(DESIGN_ROOT, "selected_polymer_matched.csv"))
    header = split(first(lines), ',')
    rows = [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1)]
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

function paths(thickness_mm)
    stem = "matching_t$(label(thickness_mm))"
    (
        mesh=joinpath(OUTPUT_ROOT, "meshes", "$(stem).msh"),
        result=joinpath(OUTPUT_ROOT, "results", "$(stem).jld2"),
    )
end

function run_thickness(thickness_mm::Real)
    thickness_mm >= 0 || throw(ArgumentError("matching thickness must be non-negative"))
    case_paths = paths(thickness_mm)
    cells = selected_cells()
    if !isfile(case_paths.mesh)
        build_material_lens_mesh(
            cells,
            case_paths.mesh;
            matching_thickness_mm=thickness_mm,
            config=MaterialLensMeshConfig(
                output_length_mm=105.0,
                lens_mesh_size_mm=0.36,
                output_mesh_size_mm=1.55,
                transition_distance_mm=14.0,
            ),
        )
        println("[+] Mesh: $(case_paths.mesh)")
    end
    isfile(case_paths.result) && return case_paths.result

    polymer = photopolymer()
    aluminium = aluminium_6061()
    matching = thickness_mm > 0 ? geometric_matching_material(polymer, aluminium) : nothing
    focus_x_mm = first(cells).length_mm + Float64(thickness_mm) + 60.0
    result = solve_material_lens_harmonic(
        case_paths.mesh,
        polymer,
        aluminium;
        matching_material=matching,
        config=MaterialLensHarmonicConfig(
            frequency_hz=242.0e3,
            focus_x_mm=focus_x_mm,
            element_order=1,
            quadrature_degree=2,
            scan_before_focus_mm=0.0,
            scan_after_focus_mm=0.0,
            scan_transverse_half_width_mm=0.0,
            scan_step_mm=2.0,
        ),
    )
    mkpath(dirname(case_paths.result))
    jldsave(
        case_paths.result;
        format_version=1,
        matching_thickness_mm=Float64(thickness_mm),
        focus_x_mm,
        focus_longitudinal_amplitude_m=result.focus_longitudinal_amplitude_m,
        focus_total_amplitude_m=result.focus_total_amplitude_m,
    )
    println("[+] t=$(thickness_mm) mm: |ux(focus)|=$(result.focus_longitudinal_amplitude_m) m")
    case_paths.result
end

function write_summary()
    result_dir = joinpath(OUTPUT_ROOT, "results")
    files = isdir(result_dir) ? filter(path -> endswith(path, ".jld2"), readdir(result_dir; join=true)) : String[]
    rows = [
        let data=load(path)
            (
                matching_thickness_mm=data["matching_thickness_mm"],
                focus_longitudinal_amplitude_m=data["focus_longitudinal_amplitude_m"],
                focus_total_amplitude_m=data["focus_total_amplitude_m"],
            )
        end
        for path in files
    ]
    sort!(rows; by=row -> row.matching_thickness_mm)
    isempty(rows) && return nothing
    path = joinpath(OUTPUT_ROOT, "rounded_notch_matching_sweep.csv")
    open(path, "w") do io
        columns = propertynames(first(rows))
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
    println("[+] Summary: $path")
    path
end

function main(args=ARGS)
    if length(args) == 1 && startswith(args[1], "--thickness=")
        thickness = parse(Float64, split(args[1], '='; limit=2)[2])
        run_thickness(thickness)
    elseif args == ["--summary"]
        write_summary()
    else
        error("usage: run_rounded_notch_matching_sweep.jl --thickness=MM or --summary")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
