module ProfileSmoothnessDiagnostic

ENV["GKSwstype"] = "100"

using Plots

include(joinpath(@__DIR__, "profiles.jl"))
using .MetamaterialProfiles

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "profile_smoothness")
const GEOMETRY = GeometryConfig()
const AMPLITUDE_MM = 2.5
const CASES = (
    (label="κ=0; N=2 (sinusoidal limit)", sharpness=0.0, periods=2),
    (label="κ=0.5; N=2", sharpness=0.5, periods=2),
    (label="κ=0.5; N=1 (stretched)", sharpness=0.5, periods=1),
    (label="κ=1; N=2", sharpness=1.0, periods=2),
    (label="κ=3; N=2", sharpness=3.0, periods=2),
)

function exponential_warp_derivatives(sharpness, raised_cosine_value)
    if abs(sharpness) < 1.0e-8
        return (value=raised_cosine_value, first=1.0, second=0.0)
    end
    denominator = expm1(sharpness)
    exponential = exp(sharpness * raised_cosine_value)
    (
        value=expm1(sharpness * raised_cosine_value) / denominator,
        first=sharpness * exponential / denominator,
        second=sharpness^2 * exponential / denominator,
    )
end

function profile_derivatives(sharpness; periods=2, label="", sample_count=4001)
    profile = ExponentialProfile(AMPLITUDE_MM; periods, sharpness)
    active_length = MetamaterialProfiles.active_span_mm(profile, GEOMETRY)
    xi = collect(range(0.0, 1.0; length=sample_count))
    phase = 2pi * periods .* xi
    seed = 0.5 .* (1 .- cos.(phase))
    seed_first = pi * periods .* sin.(phase)
    seed_second = 2pi^2 * periods^2 .* cos.(phase)
    warp = exponential_warp_derivatives.(Ref(Float64(sharpness)), seed)
    depth = AMPLITUDE_MM .* getproperty.(warp, :value)
    slope = AMPLITUDE_MM .* getproperty.(warp, :first) .* seed_first ./ active_length
    curvature = AMPLITUDE_MM .* (
        getproperty.(warp, :second) .* seed_first.^2 .+
        getproperty.(warp, :first) .* seed_second
    ) ./ active_length^2
    (
        label=String(label),
        sharpness=Float64(sharpness),
        periods=Int(periods),
        active_length_mm=active_length,
        x_mm=active_length .* xi,
        depth_mm=depth,
        slope,
        curvature_per_mm=curvature,
        max_abs_slope=maximum(abs, slope),
        max_abs_curvature_per_mm=maximum(abs, curvature),
        join_abs_curvature_per_mm=abs(first(curvature)),
        peak_abs_curvature_per_mm=abs(curvature[argmax(depth)]),
    )
end

function write_summary(rows)
    path = joinpath(OUTPUT_ROOT, "profile_smoothness_summary.csv")
    open(path, "w") do io
        println(io, "label,sharpness,periods,active_length_mm,max_abs_slope,max_abs_curvature_per_mm,join_abs_curvature_per_mm,peak_abs_curvature_per_mm")
        for row in rows
            println(io, join((
                row.label,
                row.sharpness,
                row.periods,
                row.active_length_mm,
                row.max_abs_slope,
                row.max_abs_curvature_per_mm,
                row.join_abs_curvature_per_mm,
                row.peak_abs_curvature_per_mm,
            ), ','))
        end
    end
    println("[+] $path")
    path
end

function save_figure(rows)
    depth_panel = plot(
        xlabel="local active coordinate, mm",
        ylabel="indentation, mm",
        title="Profile",
        gridalpha=0.25,
    )
    slope_panel = plot(
        xlabel="local active coordinate, mm",
        ylabel="dy/dx",
        title="First derivative",
        gridalpha=0.25,
    )
    curvature_panel = plot(
        xlabel="local active coordinate, mm",
        ylabel="d²y/dx², 1/mm",
        title="Second derivative",
        gridalpha=0.25,
    )
    colors = (:black, :deepskyblue3, :navy, :darkorange, :firebrick)
    for (index, row) in enumerate(rows)
        plot!(depth_panel, row.x_mm, row.depth_mm;
              color=colors[index], linewidth=2.2, label=row.label)
        plot!(slope_panel, row.x_mm, row.slope;
              color=colors[index], linewidth=2.2, label=row.label)
        plot!(curvature_panel, row.x_mm, row.curvature_per_mm;
              color=colors[index], linewidth=2.2, label=row.label)
    end
    figure = plot(
        depth_panel,
        slope_panel,
        curvature_panel;
        layout=(3, 1),
        size=(1200, 1100),
        margin=5Plots.mm,
        plot_title="Exponential wall smoothness, A=2.5 mm",
    )
    path = joinpath(OUTPUT_ROOT, "profile_smoothness.png")
    savefig(figure, path)
    println("[+] $path")
    path
end

function main()
    mkpath(OUTPUT_ROOT)
    rows = [
        profile_derivatives(
            case.sharpness;
            periods=case.periods,
            label=case.label,
        )
        for case in CASES
    ]
    write_summary(rows)
    save_figure(rows)
    for row in rows
        println(
            "$(row.label): max|slope|=$(round(row.max_abs_slope; digits=4)), ",
            "join |curvature|=$(round(row.join_abs_curvature_per_mm; digits=4)), ",
            "max |curvature|=$(round(row.max_abs_curvature_per_mm; digits=4))",
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
