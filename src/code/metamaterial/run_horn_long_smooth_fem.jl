module HornLongSmoothFEM

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const APERTURE_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_binary_lens_242khz")
const DESIGN_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_long_smooth_aperture")
const MESH_PATH = joinpath(APERTURE_ROOT, "mesh_point_radiator_aperture.msh")
const UNIFORM_PATH = joinpath(APERTURE_ROOT, "uniform_fem_order2.jld2")
const OUTPUT_ROOT = DESIGN_ROOT

ENV["GKSwstype"] = "100"
using JLD2
using Plots

include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
include(joinpath(@__DIR__, "measured_aperture_harmonic_solver.jl"))
include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .SinusoidalMaterialLens
using .MeasuredApertureHarmonicSolver
using .ImpulseRiskAnalysis: contiguous_width

result_path() = joinpath(OUTPUT_ROOT, "long_smooth_lens_fem_order2.jld2")

function solve_lens()
    design = JLD2.load(joinpath(DESIGN_ROOT, "long_smooth_aperture_impulse.jld2"))
    weights = ComplexF64.(design["context_coefficients"])
    result = solve_measured_aperture_harmonic(
        MESH_PATH,
        weights,
        photopolymer();
        config=MeasuredApertureHarmonicConfig(
            frequency_hz=242.0e3,
            focal_distance_mm=35.0,
            element_order=2,
            quadrature_degree=4,
            scan_step_mm=0.5,
        ),
    )
    jldsave(
        result_path();
        format_version=1,
        mask=design["mask"],
        weights,
        focus_displacement=result.focus_displacement,
        focus_longitudinal_amplitude_m=result.focus_longitudinal_amplitude_m,
        focus_total_amplitude_m=result.focus_total_amplitude_m,
        scan_x_mm=result.scan_x_mm,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=result.scan_ux_m,
        scan_uy_m=result.scan_uy_m,
        scan_total_amplitude_m=result.scan_total_amplitude_m,
    )
    println("[+] long-smooth q2 |ux(focus)|=$(result.focus_longitudinal_amplitude_m) m")
end

function fem_metrics(data)
    x_mm = data["scan_x_mm"]
    y_mm = data["scan_y_mm"]
    ux = abs.(data["scan_ux_m"])
    focus_x = argmin(abs.(x_mm .- 35.0))
    axis_y = argmin(abs.(y_mm))
    transverse = ux[:, focus_x]
    axial = ux[axis_y, :]
    transverse_width = contiguous_width(y_mm, transverse)
    target_index = argmin(abs.(x_mm .- 35.0))
    threshold = axial[target_index] / sqrt(2.0)
    left, right = target_index, target_index
    while left > firstindex(axial) && axial[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(axial) && axial[right + 1] >= threshold
        right += 1
    end
    local_peak_index = left - 1 + argmax(axial[left:right])
    center_peak = argmax(transverse)
    left_minimum = center_peak
    while left_minimum > firstindex(transverse) + 1
        left_minimum -= 1
        transverse[left_minimum] <= transverse[left_minimum - 1] &&
            transverse[left_minimum] <= transverse[left_minimum + 1] && break
    end
    right_minimum = center_peak
    while right_minimum < lastindex(transverse) - 1
        right_minimum += 1
        transverse[right_minimum] <= transverse[right_minimum - 1] &&
            transverse[right_minimum] <= transverse[right_minimum + 1] && break
    end
    outside = vcat(collect(firstindex(transverse):(left_minimum - 1)),
                   collect((right_minimum + 1):lastindex(transverse)))
    sidelobe = isempty(outside) ? 0.0 : maximum(transverse[outside]) / transverse_width.peak_amplitude
    (
        transverse,
        axial,
        transverse_fwhm_mm=transverse_width.width,
        axial_dof_mm=x_mm[right] - x_mm[left],
        peak_x_mm=x_mm[local_peak_index],
        peak_y_mm=transverse_width.peak_coordinate,
        sidelobe_amplitude_ratio=sidelobe,
    )
end

function report()
    lens = JLD2.load(result_path())
    uniform = JLD2.load(UNIFORM_PATH)
    lens_metrics = fem_metrics(lens)
    uniform_metrics = fem_metrics(uniform)
    gain = lens["focus_longitudinal_amplitude_m"] / uniform["focus_longitudinal_amplitude_m"]
    summary_path = joinpath(OUTPUT_ROOT, "long_smooth_fem_summary_order2.csv")
    open(summary_path, "w") do io
        println(io, "mask,focus_ux_m,uniform_focus_ux_m,focus_gain,transverse_fwhm_mm,uniform_transverse_fwhm_mm,axial_dof_mm,peak_x_mm,peak_y_mm,sidelobe_amplitude_ratio")
        println(io, join((
            lens["mask"],
            lens["focus_longitudinal_amplitude_m"],
            uniform["focus_longitudinal_amplitude_m"],
            gain,
            lens_metrics.transverse_fwhm_mm,
            uniform_metrics.transverse_fwhm_mm,
            lens_metrics.axial_dof_mm,
            lens_metrics.peak_x_mm,
            lens_metrics.peak_y_mm,
            lens_metrics.sidelobe_amplitude_ratio,
        ), ','))
    end

    maps = Any[]
    for (data, title) in ((lens, "context-calibrated long-smooth lens"),
                          (uniform, "uniform point aperture"))
        panel = heatmap(
            data["scan_x_mm"], data["scan_y_mm"], data["scan_total_amplitude_m"] .* 1e9;
            xlabel="x, mm", ylabel="y, mm", title,
            color=:viridis, colorbar_title="|u|, nm", aspect_ratio=:equal, grid=false,
        )
        scatter!(panel, [35.0], [0.0]; marker=:xcross, color=:white,
                 markerstrokewidth=2, label=false)
        push!(maps, panel)
    end
    map_path = joinpath(OUTPUT_ROOT, "long_smooth_fem_maps_order2.png")
    savefig(plot(maps...; layout=(2, 1), size=(1050, 920), left_margin=7Plots.mm), map_path)
    profile = plot(
        lens["scan_y_mm"], lens_metrics.transverse .* 1e9;
        linewidth=2.5, label="long-smooth binary", xlabel="y at x=35 mm, mm",
        ylabel="|ux|, nm", title="Quadratic FEM transverse focus", gridalpha=0.25,
    )
    plot!(profile, uniform["scan_y_mm"], uniform_metrics.transverse .* 1e9;
          linewidth=2, label="uniform aperture")
    profile_path = joinpath(OUTPUT_ROOT, "long_smooth_fem_profiles_order2.png")
    savefig(profile, profile_path)
    println("[+] q2 gain=$gain, FWHM=$(lens_metrics.transverse_fwhm_mm) mm")
    println("[+] peak x=$(lens_metrics.peak_x_mm) mm, sidelobe=$(lens_metrics.sidelobe_amplitude_ratio)")
    println("[+] $summary_path")
    println("[+] $map_path")
    println("[+] $profile_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if isempty(ARGS) || first(ARGS) == "solve"
        solve_lens()
    elseif first(ARGS) == "report"
        report()
    else
        error("usage: julia run_horn_long_smooth_fem.jl [solve|report]")
    end
end

end
