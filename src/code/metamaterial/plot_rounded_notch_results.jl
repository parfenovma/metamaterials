ENV["GKSwstype"] = "100"

using JLD2
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_ROOT = joinpath(PROJECT_ROOT, "tmp", "rounded_notch_material_lens_242khz")
const BASE_FEM_ROOT = joinpath(PROJECT_ROOT, "tmp", "rounded_notch_lens_fem_242khz")
const OPTIMAL_Q1_ROOT = joinpath(PROJECT_ROOT, "tmp", "rounded_notch_lens_fem_t4p2_242khz")
const OPTIMAL_Q2_ROOT = joinpath(PROJECT_ROOT, "tmp", "rounded_notch_lens_fem_t4p2_q2_242khz")
const SWEEP_ROOT = joinpath(PROJECT_ROOT, "tmp", "rounded_notch_matching_sweep_242khz")
const SINUSOIDAL_Q2_ROOT = joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_fem_q2_242khz")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "rounded_notch_results_242khz")

function csv_rows(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1) if !isempty(strip(line))]
end

function result(root, configuration, role)
    load(joinpath(root, "results", "$(configuration)_$(role).jld2"))
end

function save_base_focus_maps()
    configurations = (:polymer_direct, :polymer_matched, :aluminium)
    titles = Dict(
        :polymer_direct => "photopolymer → Al",
        :polymer_matched => "photopolymer → 2.2 mm match → Al",
        :aluminium => "Al → Al",
    )
    panels = Plots.Plot[]
    for configuration in configurations
        data = result(BASE_FEM_ROOT, configuration, :selected)
        panel = heatmap(
            data["scan_x_mm"],
            data["scan_y_mm"],
            data["scan_amplitude_m"] .* 1.0e9;
            xlabel="x, mm",
            ylabel="y, mm",
            title=titles[configuration],
            color=:viridis,
            colorbar_title="|u|, nm",
            aspect_ratio=:equal,
            grid=false,
        )
        scatter!(panel, [data["focus_x_mm"]], [0.0]; marker=:xcross, color=:white,
                 markerstrokewidth=2, markersize=7, label=false)
        push!(panels, panel)
    end
    path = joinpath(OUTPUT_ROOT, "rounded_notch_base_focus_maps.png")
    savefig(plot(panels...; layout=(3, 1), size=(1000, 1250), left_margin=7Plots.mm), path)
    path
end

function save_optimal_focus_maps()
    selected = result(OPTIMAL_Q2_ROOT, :polymer_matched, :selected)
    uniform = result(OPTIMAL_Q2_ROOT, :polymer_matched, :uniform)
    panels = Plots.Plot[]
    for (data, title) in ((selected, "rounded-notch lens"), (uniform, "straight control"))
        panel = heatmap(
            data["scan_x_mm"],
            data["scan_y_mm"],
            data["scan_amplitude_m"] .* 1.0e9;
            xlabel="x, mm",
            ylabel="y, mm",
            title=title,
            color=:viridis,
            colorbar_title="|u|, nm",
            aspect_ratio=:equal,
            grid=false,
        )
        scatter!(panel, [data["focus_x_mm"]], [0.0]; marker=:xcross, color=:white,
                 markerstrokewidth=2, markersize=7, label=false)
        push!(panels, panel)
    end
    path = joinpath(OUTPUT_ROOT, "rounded_notch_optimal_q2_focus_maps.png")
    savefig(plot(panels...; layout=(2, 1), size=(1000, 900), left_margin=7Plots.mm), path)
    path
end

function save_matching_sweep()
    rows = csv_rows(joinpath(SWEEP_ROOT, "rounded_notch_matching_sweep.csv"))
    thickness = parse.(Float64, getindex.(rows, "matching_thickness_mm"))
    amplitude_nm = 1.0e9 .* parse.(Float64, getindex.(rows, "focus_longitudinal_amplitude_m"))
    best = argmax(amplitude_nm)
    figure = plot(
        thickness,
        amplitude_nm;
        marker=:circle,
        linewidth=2.5,
        markersize=5,
        label="linear FEM",
        xlabel="matching-layer thickness, mm",
        ylabel="|ux| at target, nm",
        title="Rounded-notch photopolymer lens at 242 kHz",
        gridalpha=0.25,
        size=(950, 600),
        left_margin=8Plots.mm,
        bottom_margin=7Plots.mm,
    )
    scatter!(figure, [thickness[best]], [amplitude_nm[best]]; color=:red,
             markersize=8, label="sampled optimum: $(thickness[best]) mm")
    path = joinpath(OUTPUT_ROOT, "rounded_notch_matching_sweep.png")
    savefig(figure, path)
    path
end

function save_verified_comparison()
    rounded_selected = result(OPTIMAL_Q2_ROOT, :polymer_matched, :selected)
    rounded_uniform = result(OPTIMAL_Q2_ROOT, :polymer_matched, :uniform)
    sinusoidal_selected = result(SINUSOIDAL_Q2_ROOT, :aluminium, :selected)
    sinusoidal_uniform = result(SINUSOIDAL_Q2_ROOT, :aluminium, :uniform)
    names = ["rounded notch\npolymer + 4.2 mm match", "sinusoidal\naluminium"]
    selected_nm = 1.0e9 .* [
        rounded_selected["focus_longitudinal_amplitude_m"],
        sinusoidal_selected["focus_longitudinal_amplitude_m"],
    ]
    uniform_nm = 1.0e9 .* [
        rounded_uniform["focus_longitudinal_amplitude_m"],
        sinusoidal_uniform["focus_longitudinal_amplitude_m"],
    ]
    x = collect(eachindex(names))
    figure = bar(
        x .- 0.16,
        uniform_nm;
        bar_width=0.30,
        label="straight control",
        ylabel="|ux| at target, nm",
        title="Quadratic FEM, 242 kHz, 1 MPa incident traction",
        xticks=(x, names),
        gridalpha=0.25,
        size=(950, 620),
        left_margin=8Plots.mm,
        bottom_margin=14Plots.mm,
    )
    bar!(figure, x .+ 0.16, selected_nm; bar_width=0.30, label="profiled lens")
    path = joinpath(OUTPUT_ROOT, "verified_profile_comparison_q2.png")
    savefig(figure, path)
    path
end

mkpath(OUTPUT_ROOT)
println("[+] Base focus maps: $(save_base_focus_maps())")
println("[+] Optimized focus maps: $(save_optimal_focus_maps())")
println("[+] Matching sweep: $(save_matching_sweep())")
println("[+] Verified comparison: $(save_verified_comparison())")
