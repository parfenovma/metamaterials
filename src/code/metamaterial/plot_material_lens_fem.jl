ENV["GKSwstype"] = "100"

using JLD2
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_FEM_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_fem_242khz"),
)
const CONFIGURATIONS = (:polymer_direct, :polymer_matched, :aluminium)

function selected_result(configuration)
    load(joinpath(OUTPUT_ROOT, "results", "$(configuration)_selected.jld2"))
end

function save_focus_maps()
    panels = Plots.Plot[]
    titles = Dict(
        :polymer_direct => "photopolymer → Al",
        :polymer_matched => "photopolymer → match → Al",
        :aluminium => "Al → Al",
    )
    for configuration in CONFIGURATIONS
        data = selected_result(configuration)
        amplitude_um = data["scan_amplitude_m"] .* 1.0e6
        panel = heatmap(
            data["scan_x_mm"],
            data["scan_y_mm"],
            amplitude_um;
            xlabel="x, mm",
            ylabel="y, mm",
            title=titles[configuration],
            color=:viridis,
            colorbar_title="|u|, µm",
            aspect_ratio=:equal,
            grid=false,
        )
        scatter!(
            panel,
            [data["focus_x_mm"]],
            [0.0];
            marker=:xcross,
            color=:white,
            markersize=7,
            markerstrokewidth=2,
            label=false,
        )
        scatter!(
            panel,
            [data["local_peak_x_mm"]],
            [data["local_peak_y_mm"]];
            marker=:circle,
            color=:red,
            markersize=4,
            label=false,
        )
        push!(panels, panel)
    end
    path = joinpath(OUTPUT_ROOT, "material_lens_fem_focus_maps.png")
    savefig(plot(panels...; layout=(3, 1), size=(1000, 1250), left_margin=7Plots.mm), path)
    path
end

function save_comparison()
    summary_path = joinpath(OUTPUT_ROOT, "material_lens_fem_summary.csv")
    lines = readlines(summary_path)
    header = split(first(lines), ',')
    rows = [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1)]
    names = replace.(getindex.(rows, "configuration"), '_' => ' ')
    selected = parse.(Float64, getindex.(rows, "selected_focus_ux_m")) .* 1.0e6
    uniform = parse.(Float64, getindex.(rows, "uniform_focus_ux_m")) .* 1.0e6
    peak = parse.(Float64, getindex.(rows, "selected_local_peak_m")) .* 1.0e6
    x = collect(eachindex(names))
    figure = bar(
        x .- 0.24,
        uniform;
        bar_width=0.22,
        label="uniform at target",
        ylabel="displacement amplitude, µm",
        title="Full-lens FEM at 242 kHz, 1 MPa incident traction",
        size=(1050, 650),
        xticks=(x, names),
        xrotation=12,
        gridalpha=0.25,
        left_margin=8Plots.mm,
        bottom_margin=10Plots.mm,
    )
    bar!(figure, x, selected; bar_width=0.22, label="lens at target")
    bar!(figure, x .+ 0.24, peak; bar_width=0.22, label="lens local peak")
    path = joinpath(OUTPUT_ROOT, "material_lens_fem_comparison.png")
    savefig(figure, path)
    path
end

println("[+] Focus maps: $(save_focus_maps())")
println("[+] Comparison: $(save_comparison())")
