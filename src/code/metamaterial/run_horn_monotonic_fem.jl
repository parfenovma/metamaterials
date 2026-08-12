module HornMonotonicFEM

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const APERTURE_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_binary_lens_242khz")
const DESIGN_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_monotonic_aperture")
const BINARY_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_long_smooth_aperture")
const MESH_PATH = joinpath(APERTURE_ROOT, "mesh_point_radiator_aperture.msh")
const UNIFORM_PATH = joinpath(APERTURE_ROOT, "uniform_fem_order2.jld2")
const BINARY_PATH = joinpath(BINARY_ROOT, "long_smooth_lens_fem_order2.jld2")
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

result_path() = joinpath(OUTPUT_ROOT, "monotonic_lens_fem_order2.jld2")

function solve_lens()
    design = JLD2.load(joinpath(DESIGN_ROOT, "monotonic_aperture_impulse.jld2"))
    weights = ComplexF64.(design["weights"])
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
        transfer_source=design["transfer_source"],
        maximum_extra_path_mm=design["maximum_extra_path_mm"],
        common_axial_length_mm=design["common_axial_length_mm"],
        centers_mm=design["centers_mm"],
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
    println("[+] monotonic q2 |ux(focus)|=$(result.focus_longitudinal_amplitude_m) m")
end

function sidelobe_ratio(profile, peak_index)
    left_minimum = peak_index
    while left_minimum > firstindex(profile) + 1
        left_minimum -= 1
        profile[left_minimum] <= profile[left_minimum - 1] &&
            profile[left_minimum] <= profile[left_minimum + 1] && break
    end
    right_minimum = peak_index
    while right_minimum < lastindex(profile) - 1
        right_minimum += 1
        profile[right_minimum] <= profile[right_minimum - 1] &&
            profile[right_minimum] <= profile[right_minimum + 1] && break
    end
    outside = vcat(
        collect(firstindex(profile):(left_minimum - 1)),
        collect((right_minimum + 1):lastindex(profile)),
    )
    isempty(outside) ? 0.0 : maximum(profile[outside]) / profile[peak_index]
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
    peak_index = argmax(transverse)
    (
        transverse,
        axial,
        transverse_fwhm_mm=transverse_width.width,
        axial_dof_mm=x_mm[right] - x_mm[left],
        peak_x_mm=x_mm[local_peak_index],
        peak_y_mm=transverse_width.peak_coordinate,
        sidelobe_amplitude_ratio=sidelobe_ratio(transverse, peak_index),
    )
end

function report()
    lens = JLD2.load(result_path())
    uniform = JLD2.load(UNIFORM_PATH)
    binary = JLD2.load(BINARY_PATH)
    lens_metrics = fem_metrics(lens)
    uniform_metrics = fem_metrics(uniform)
    binary_metrics = fem_metrics(binary)
    gain = lens["focus_longitudinal_amplitude_m"] /
           uniform["focus_longitudinal_amplitude_m"]
    binary_gain = binary["focus_longitudinal_amplitude_m"] /
                  uniform["focus_longitudinal_amplitude_m"]
    gain_over_binary = lens["focus_longitudinal_amplitude_m"] /
                       binary["focus_longitudinal_amplitude_m"]
    summary_path = joinpath(OUTPUT_ROOT, "monotonic_fem_summary_order2.csv")
    open(summary_path, "w") do io
        println(io, "transfer_source,focus_ux_m,uniform_focus_ux_m,binary_focus_ux_m,focus_gain,binary_gain,gain_over_binary,transverse_fwhm_mm,binary_transverse_fwhm_mm,uniform_transverse_fwhm_mm,axial_dof_mm,peak_x_mm,peak_y_mm,sidelobe_amplitude_ratio")
        println(io, join((
            lens["transfer_source"],
            lens["focus_longitudinal_amplitude_m"],
            uniform["focus_longitudinal_amplitude_m"],
            binary["focus_longitudinal_amplitude_m"],
            gain,
            binary_gain,
            gain_over_binary,
            lens_metrics.transverse_fwhm_mm,
            binary_metrics.transverse_fwhm_mm,
            uniform_metrics.transverse_fwhm_mm,
            lens_metrics.axial_dof_mm,
            lens_metrics.peak_x_mm,
            lens_metrics.peak_y_mm,
            lens_metrics.sidelobe_amplitude_ratio,
        ), ','))
    end

    common_limit_nm = maximum(
        maximum(data["scan_total_amplitude_m"]) for data in (lens, binary, uniform)
    ) * 1e9
    maps = Any[]
    for (data, title) in (
        (lens, "monotonic true-time-delay lens"),
        (binary, "best binary lens"),
        (uniform, "uniform point aperture"),
    )
        panel = heatmap(
            data["scan_x_mm"], data["scan_y_mm"], data["scan_total_amplitude_m"] .* 1e9;
            xlabel="x, mm", ylabel="y, mm", title,
            color=:viridis, colorbar_title="|u|, nm", aspect_ratio=:equal, grid=false,
            clims=(0.0, common_limit_nm),
        )
        scatter!(panel, [35.0], [0.0]; marker=:xcross, color=:white,
                 markerstrokewidth=2, label=false)
        push!(maps, panel)
    end
    map_path = joinpath(OUTPUT_ROOT, "monotonic_fem_maps_order2.png")
    savefig(plot(maps...; layout=(3, 1), size=(1050, 1320), left_margin=7Plots.mm), map_path)

    profile = plot(
        lens["scan_y_mm"], lens_metrics.transverse .* 1e9;
        linewidth=2.8, label="monotonic TTD", xlabel="y at x=35 mm, mm",
        ylabel="|ux|, nm", title="Quadratic FEM transverse focus", gridalpha=0.25,
    )
    plot!(profile, binary["scan_y_mm"], binary_metrics.transverse .* 1e9;
          linewidth=2.2, label="best binary")
    plot!(profile, uniform["scan_y_mm"], uniform_metrics.transverse .* 1e9;
          linewidth=2, label="uniform aperture")
    profile_path = joinpath(OUTPUT_ROOT, "monotonic_fem_profiles_order2.png")
    savefig(profile, profile_path)
    println("[+] q2 gain=$gain, gain/binary=$gain_over_binary, FWHM=$(lens_metrics.transverse_fwhm_mm) mm")
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
        error("usage: julia run_horn_monotonic_fem.jl [solve|report]")
    end
end

end
