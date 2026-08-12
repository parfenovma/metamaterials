module AnalyzeSinusoidalMaterialTransfer

using JLD2
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN4 = joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_242khz")
const DESIGN5 = joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_5cycle_242khz")
const FEM_LINEAR = joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_fem_242khz")
const FEM_MATCHED5 = joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_fem_5cycle_242khz")
const FEM_Q2 = joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_fem_q2_242khz")
const FEM_Q2_REFINED = joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_fem_q2_refined_242khz")
const MATCHING_SWEEP = joinpath(PROJECT_ROOT, "tmp", "matching_layer_fem_sweep_242khz")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_TRANSFER_ANALYSIS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_transfer_analysis_242khz"),
)

const FREQUENCY_HZ = 242.0e3
const ALUMINIUM_CP_M_S = 6122.102437409232
const ALUMINIUM_CS_M_S = 3083.810277185563

result(root, name) = load(joinpath(root, "results", "$name.jld2"))

function csv_rows(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1) if !isempty(strip(line))]
end

function write_rows(path, rows)
    open(path, "w") do io
        columns = propertynames(first(rows))
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
    path
end

function identical_design(configuration)
    read(joinpath(DESIGN4, "selected_$(configuration).csv"), String) ==
    read(joinpath(DESIGN5, "selected_$(configuration).csv"), String)
end

function material_rows()
    identical_design(:polymer_direct) || error("5-cycle direct-polymer geometry changed")
    identical_design(:aluminium) || error("5-cycle aluminium geometry changed")

    polymer_direct_selected = result(FEM_LINEAR, "polymer_direct_selected")
    polymer_direct_uniform = result(FEM_LINEAR, "polymer_direct_uniform")
    polymer_matched_selected = result(FEM_MATCHED5, "polymer_matched_selected")
    polymer_matched_uniform = result(FEM_LINEAR, "polymer_matched_uniform")
    aluminium_selected = result(FEM_LINEAR, "aluminium_selected")
    aluminium_uniform = result(FEM_LINEAR, "aluminium_uniform")
    polymer_matched_selected["focus_x_mm"] == polymer_matched_uniform["focus_x_mm"] ||
        error("matched-polymer reference has a different focus")
    polymer_matched_selected["matching_thickness_mm"] ==
        polymer_matched_uniform["matching_thickness_mm"] ||
        error("matched-polymer reference has a different layer")

    configurations = (
        (
            "photopolymer_direct", polymer_direct_selected, polymer_direct_uniform,
            "5-cycle profile equals 4-cycle profile exactly",
        ),
        (
            "photopolymer_match_3p5mm", polymer_matched_selected, polymer_matched_uniform,
            "5-cycle profile; geometric-mean matching layer",
        ),
        (
            "aluminium_to_aluminium", aluminium_selected, aluminium_uniform,
            "5-cycle profile equals 4-cycle profile exactly",
        ),
    )
    [
        (
            configuration=name,
            selected_focus_ux_m=selected["focus_longitudinal_amplitude_m"],
            uniform_focus_ux_m=uniform["focus_longitudinal_amplitude_m"],
            matching_thickness_mm=selected["matching_thickness_mm"],
            frequency_hz=selected["frequency_hz"],
            source_traction_pa=1.0e6,
            fem_element_order=selected["element_order"],
            focus_gain=selected["focus_longitudinal_amplitude_m"] /
                       uniform["focus_longitudinal_amplitude_m"],
            note,
        )
        for (name, selected, uniform, note) in configurations
    ]
end

function convergence_rows()
    levels = (
        ("linear_fine", FEM_LINEAR, 1, 0.42, 1.55),
        ("quadratic_coarse", FEM_Q2, 2, 0.75, 2.60),
        ("quadratic_refined", FEM_Q2_REFINED, 2, 0.60, 2.00),
    )
    [
        let selected=result(root, "aluminium_selected"), uniform=result(root, "aluminium_uniform")
            (
                level,
                selected_focus_ux_m=selected["focus_longitudinal_amplitude_m"],
                uniform_focus_ux_m=uniform["focus_longitudinal_amplitude_m"],
                lens_mesh_size_mm=h_lens,
                output_mesh_size_mm=h_output,
                element_order=order,
                focus_gain=selected["focus_longitudinal_amplitude_m"] /
                           uniform["focus_longitudinal_amplitude_m"],
                lens_h_over_lambda_s=h_lens / (1.0e3ALUMINIUM_CS_M_S / FREQUENCY_HZ),
                output_h_over_lambda_s=h_output / (1.0e3ALUMINIUM_CS_M_S / FREQUENCY_HZ),
            )
        end
        for (level, root, order, h_lens, h_output) in levels
    ]
end

function matching_rows()
    rows = csv_rows(joinpath(MATCHING_SWEEP, "matching_layer_fem_sweep.csv"))
    sort([
        (
            thickness_mm=parse(Float64, row["matching_thickness_mm"]),
            focus_ux_m=parse(Float64, row["focus_longitudinal_amplitude_m"]),
        )
        for row in rows
    ]; by=row -> row.thickness_mm)
end

function make_plot(materials, convergence, matching)
    labels = ["polymer\ndirect", "polymer +\n3.5 mm layer", "Al → Al"]
    amplitudes_nm = 1.0e9 .* getproperty.(materials, :selected_focus_ux_m)
    p1 = bar(
        labels, amplitudes_nm;
        ylabel="|uₓ(focus)|, nm", title="Fine linear FEM, same 1 MPa traction",
        label=false, color=[:darkorange, :mediumpurple, :steelblue],
        ylim=(0, 1.12maximum(amplitudes_nm)),
    )
    for (index, value) in enumerate(amplitudes_nm)
        annotate!(p1, index, value + 0.7, text(string(round(value, digits=2)), 9))
    end

    thickness = getproperty.(matching, :thickness_mm)
    matching_nm = 1.0e9 .* getproperty.(matching, :focus_ux_m)
    p2 = plot(
        thickness, matching_nm;
        xlabel="matching-layer thickness, mm", ylabel="|uₓ(focus)|, nm",
        title="Fixed direct-polymer profile", lw=2.5, marker=:circle,
        label="full FEM",
    )
    best = argmax(matching_nm)
    scatter!(p2, [thickness[best]], [matching_nm[best]]; ms=7, label="best sampled")

    levels = ["P1 fine", "P2 coarse", "P2 refined"]
    p3 = plot(
        levels, 1.0e9 .* getproperty.(convergence, :selected_focus_ux_m);
        ylabel="|uₓ(focus)|, nm", title="Al → Al mesh/order check",
        marker=:circle, lw=2.5, label="lens",
    )
    plot!(p3, levels, 1.0e9 .* getproperty.(convergence, :uniform_focus_ux_m);
          marker=:circle, lw=2.5, ls=:dash, label="straight")

    path = joinpath(OUTPUT_ROOT, "sinusoidal_material_transfer.png")
    savefig(plot(
        p1, p2, p3;
        layout=(1, 3), size=(1800, 600),
        bottom_margin=8Plots.mm, top_margin=5Plots.mm,
        left_margin=5Plots.mm, right_margin=5Plots.mm,
    ), path)
    path
end

function main()
    mkpath(OUTPUT_ROOT)
    materials = material_rows()
    convergence = convergence_rows()
    matching = matching_rows()
    material_path = write_rows(joinpath(OUTPUT_ROOT, "material_comparison_5cycle.csv"), materials)
    convergence_path = write_rows(joinpath(OUTPUT_ROOT, "aluminium_convergence.csv"), convergence)
    matching_path = write_rows(joinpath(OUTPUT_ROOT, "matching_layer_sweep.csv"), matching)
    plot_path = make_plot(materials, convergence, matching)
    best_match = matching[argmax(getproperty.(matching, :focus_ux_m))]
    refined = last(convergence)
    println("[+] 5-cycle fine-FEM material comparison: $material_path")
    println("[+] Al convergence: $convergence_path")
    println("[+] Matching sweep: $matching_path")
    println("[+] Plot: $plot_path")
    println("[+] Best fixed-profile layer: $(best_match.thickness_mm) mm, $(best_match.focus_ux_m) m")
    println("[+] Refined Al gain: $(refined.focus_gain)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
