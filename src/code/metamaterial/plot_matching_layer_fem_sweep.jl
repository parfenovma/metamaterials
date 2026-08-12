ENV["GKSwstype"] = "100"

using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "matching_layer_fem_sweep_242khz")

function rows(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(zip(header, split(line, ','))) for line in Iterators.drop(lines, 1)]
end

summary_path = joinpath(OUTPUT_ROOT, "matching_layer_fem_sweep.csv")
data = rows(summary_path)
thickness = parse.(Float64, getindex.(data, "matching_thickness_mm"))
target = parse.(Float64, getindex.(data, "focus_longitudinal_amplitude_m")) .* 1.0e6
local_peak = parse.(Float64, getindex.(data, "local_peak_amplitude_m")) .* 1.0e6
best_index = argmax(target)

figure = plot(
    thickness,
    target;
    marker=:circle,
    linewidth=2,
    label="|ux| at target",
    xlabel="matching-layer thickness, mm",
    ylabel="displacement amplitude, µm",
    title="Photopolymer lens → matching layer → Al, 242 kHz",
    gridalpha=0.25,
    size=(1000, 620),
    left_margin=8Plots.mm,
    bottom_margin=8Plots.mm,
)
plot!(figure, thickness, local_peak; marker=:diamond, linewidth=2, label="local |u| peak")
vline!(figure, [3.9100550098641302]; linestyle=:dash, color=:gray, label="P quarter-wave")
scatter!(
    figure,
    [thickness[best_index]],
    [target[best_index]];
    color=:red,
    markersize=7,
    label="best target",
)
path = joinpath(OUTPUT_ROOT, "matching_layer_fem_sweep.png")
savefig(figure, path)
println("[+] Matching-layer plot: $path")
