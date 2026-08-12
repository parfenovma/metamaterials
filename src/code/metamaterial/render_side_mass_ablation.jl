module RenderSideMassAblation

ENV["GKSwstype"] = "100"

using JLD2
using Plots
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_SIDE_MASS_ABLATION_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "side_mass_ablation"),
)
const VARIANTS = (:r0, :b, :bd)
const FIELD_LABEL = Symbol(get(ENV, "METAMATERIALS_ABLATION_FIELD_LABEL", "transparency"))
const FRAME_COUNT = parse(Int, get(ENV, "METAMATERIALS_ABLATION_ANIMATION_FRAMES", "32"))
const FRAME_RATE = parse(Int, get(ENV, "METAMATERIALS_ABLATION_ANIMATION_FPS", "16"))

struct CompactField
    variant::Symbol
    frequency_hz::Float64
    x_mm::Vector{Float64}
    y_mm::Vector{Float64}
    ux_mm::Vector{ComplexF64}
    uy_mm::Vector{ComplexF64}
end

field_path(variant) = joinpath(OUTPUT_ROOT, "fields", "field_$(FIELD_LABEL)_$(variant).jld2")

function compact_field(variant)
    data = JLD2.load(field_path(variant))
    x = Float64.(data["node_x_m"]) .* 1e3
    y = Float64.(data["node_y_m"]) .* 1e3
    ux = ComplexF64.(data["displacement_x_m"]) .* 1e3
    uy = ComplexF64.(data["displacement_y_m"]) .* 1e3
    seen = Set{Tuple{Int, Int}}()
    indices = Int[]
    for index in eachindex(x)
        key = (round(Int, x[index] * 1e9), round(Int, y[index] * 1e9))
        if key ∉ seen
            push!(seen, key)
            push!(indices, index)
        end
    end
    CompactField(
        variant,
        Float64(data["frequency_hz"]),
        x[indices],
        y[indices],
        ux[indices],
        uy[indices],
    )
end

function rectangle_mask(field, center_x, length_x, y_min, y_max)
    x_min = center_x - length_x / 2
    x_max = center_x + length_x / 2
    (field.x_mm .>= x_min) .& (field.x_mm .<= x_max) .&
        (field.y_mm .>= y_min) .& (field.y_mm .<= y_max)
end

bright_mask(field) = rectangle_mask(field, 10.0, 3.0, 5.55, 6.65)
dark_mask(field) = rectangle_mask(field, 16.0, 3.0, 5.55, 6.65)

function region_metrics(field, mask)
    count(mask) > 0 || return (rms_mm=NaN, mean_ux_mm=ComplexF64(NaN, NaN))
    rms = sqrt(mean(abs2.(field.ux_mm[mask]) .+ abs2.(field.uy_mm[mask])))
    (rms_mm=rms, mean_ux_mm=mean(field.ux_mm[mask]))
end

function phase_difference_deg(a, b)
    (abs(a) > eps() && abs(b) > eps()) || return NaN
    rad2deg(mod(angle(a / b) + pi, 2pi) - pi)
end

function variant_title(field)
    if field.variant == :r0
        return "R0: direct channel"
    end
    bright = region_metrics(field, bright_mask(field))
    if field.variant == :b
        return "B: bright only  |u_b|=$(round(bright.rms_mm * 1e6; digits=2)) nm"
    end
    dark = region_metrics(field, dark_mask(field))
    ratio = dark.rms_mm / bright.rms_mm
    delta = phase_difference_deg(dark.mean_ux_mm, bright.mean_ux_mm)
    "BD: |u_d|/|u_b|=$(round(ratio; digits=2)),  Δφ=$(round(delta; digits=1))°"
end

function add_region_guides!(plot_object, variant)
    variant == :r0 && return plot_object
    plot!(
        plot_object,
        Shape([8.5, 11.5, 11.5, 8.5], [5.55, 5.55, 6.65, 6.65]);
        fillalpha=0.0,
        linecolor=:darkorange,
        linestyle=:dash,
        linewidth=1.5,
        label=false,
    )
    annotate!(plot_object, 10.0, 6.9, text("bright", 8, :darkorange))
    if variant == :bd
        plot!(
            plot_object,
            Shape([14.5, 17.5, 17.5, 14.5], [5.55, 5.55, 6.65, 6.65]);
            fillalpha=0.0,
            linecolor=:royalblue,
            linestyle=:dash,
            linewidth=1.5,
            label=false,
        )
        annotate!(plot_object, 16.0, 6.9, text("dark", 8, :royalblue))
    end
    plot_object
end

function frame_panel(field, phase, deformation_scale, color_scale)
    factor = exp(-im * phase)
    ux = real.(field.ux_mm .* factor)
    uy = real.(field.uy_mm .* factor)
    x_deformed = field.x_mm .+ deformation_scale .* ux
    y_deformed = field.y_mm .+ deformation_scale .* uy
    panel = scatter(
        x_deformed,
        y_deformed;
        marker_z=ux ./ color_scale,
        color=:balance,
        colorbar=false,
        clims=(-1, 1),
        markersize=1.8,
        markerstrokewidth=0,
        label=false,
        aspect_ratio=:equal,
        xlims=(-0.5, 30.5),
        ylims=(-0.7, 7.3),
        xlabel="x, mm",
        ylabel="y, mm",
        title=variant_title(field),
        grid=false,
        framestyle=:box,
    )
    add_region_guides!(panel, field.variant)
end

function save_metrics(fields)
    path = joinpath(OUTPUT_ROOT, "mass_motion_metrics_$(FIELD_LABEL).csv")
    open(path, "w") do io
        println(io, "variant,frequency_hz,bright_rms_mm,dark_rms_mm,dark_to_bright,phase_dark_minus_bright_deg")
        for field in fields
            bright = field.variant == :r0 ? (rms_mm=NaN, mean_ux_mm=ComplexF64(NaN, NaN)) :
                     region_metrics(field, bright_mask(field))
            dark = field.variant == :bd ? region_metrics(field, dark_mask(field)) :
                   (rms_mm=NaN, mean_ux_mm=ComplexF64(NaN, NaN))
            ratio = dark.rms_mm / bright.rms_mm
            delta = phase_difference_deg(dark.mean_ux_mm, bright.mean_ux_mm)
            println(io, join((
                uppercase(String(field.variant)),
                field.frequency_hz,
                bright.rms_mm,
                dark.rms_mm,
                ratio,
                delta,
            ), ','))
        end
    end
    println("[+] $path")
end

function encode_animation(frame_dir)
    ffmpeg = get(ENV, "METAMATERIALS_FFMPEG", "ffmpeg")
    input_pattern = joinpath(frame_dir, "frame_%03d.png")
    mp4_path = joinpath(OUTPUT_ROOT, "bright_dark_motion_$(FIELD_LABEL).mp4")
    gif_path = joinpath(OUTPUT_ROOT, "bright_dark_motion_$(FIELD_LABEL).gif")
    run(Cmd(String[
        ffmpeg,
        "-y",
        "-framerate",
        string(FRAME_RATE),
        "-i",
        input_pattern,
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        mp4_path,
    ]))
    run(Cmd(String[
        ffmpeg,
        "-y",
        "-framerate",
        string(FRAME_RATE),
        "-i",
        input_pattern,
        "-vf",
        "fps=$(FRAME_RATE),scale=1200:-1:flags=lanczos",
        "-loop",
        "0",
        gif_path,
    ]))
    println("[+] $mp4_path")
    println("[+] $gif_path")
end


function main()
    fields = compact_field.(VARIANTS)
    all_amplitudes = reduce(vcat, [sqrt.(abs2.(field.ux_mm) .+ abs2.(field.uy_mm)) for field in fields])
    maximum_displacement_mm = maximum(all_amplitudes)
    maximum_displacement_mm > 0 || error("harmonic fields are identically zero")
    deformation_scale = 0.35 / maximum_displacement_mm
    color_scale = maximum(reduce(vcat, abs.(field.ux_mm) for field in fields))
    save_metrics(fields)

    frame_dir = joinpath(OUTPUT_ROOT, "animation_frames_$(FIELD_LABEL)")
    mkpath(frame_dir)
    frequency_khz = first(fields).frequency_hz / 1e3
    for frame_index in 1:FRAME_COUNT
        phase = 2pi * (frame_index - 1) / FRAME_COUNT
        panels = [frame_panel(field, phase, deformation_scale, color_scale) for field in fields]
        figure = plot(
            panels...;
            layout=(3, 1),
            size=(1500, 980),
            margin=4Plots.mm,
            plot_title="$(FIELD_LABEL) at $(round(frequency_khz; digits=3)) kHz — deformation ×$(round(deformation_scale; digits=0))",
        )
        path = joinpath(frame_dir, "frame_$(lpad(frame_index, 3, '0')).png")
        savefig(figure, path)
        println("[frame $frame_index/$FRAME_COUNT] $path")
    end
    encode_animation(frame_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
