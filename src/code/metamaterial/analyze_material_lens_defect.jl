module AnalyzeMaterialLensDefect

using JLD2
using Plots
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASE_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_DEFECT_BASE_ROOT",
    joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_fem_q2_refined_242khz"),
)
const DEFECT_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_DEFECT_ROOT",
    joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_defect_q2_refined_242khz"),
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_DEFECT_ANALYSIS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_defect_analysis_242khz"),
)

const FREQUENCY_HZ = 242.0e3
const ALUMINIUM_DENSITY_KG_M3 = 2700.0
const ALUMINIUM_CP_M_S = 6122.102437409232
const ALUMINIUM_CS_M_S = 3083.810277185563
const DEFECT_RADIUS_MM = 2.0
const DEFECT_X_MM = 88.0
const DEFECT_Y_MM = 0.0
const RECEIVER_X_MM = 68.0
const DEFECT_SUFFIX = "_defect_r2p0_x88p0_y0p0"

result_path(root, role; defect=false) = joinpath(
    root,
    "results",
    "aluminium_$(role)$(defect ? DEFECT_SUFFIX : "").jld2",
)

finite_complex(value) = isfinite(real(value)) && isfinite(imag(value))

function kinetic_energy_density(displacement)
    omega = 2pi * FREQUENCY_HZ
    0.25 * ALUMINIUM_DENSITY_KG_M3 * omega^2 * sum(abs2, displacement)
end

function rms_finite(values)
    clean = filter(isfinite, vec(values))
    isempty(clean) && return NaN
    sqrt(mean(abs2, clean))
end

function load_pair(role)
    baseline = load(result_path(BASE_ROOT, role))
    defect = load(result_path(DEFECT_ROOT, role; defect=true))
    baseline["scan_x_mm"] == defect["scan_x_mm"] || error("x scan grids differ")
    baseline["scan_y_mm"] == defect["scan_y_mm"] || error("y scan grids differ")
    baseline, defect
end

function scattered_fields(baseline, defect)
    ux = defect["scan_ux_m"] .- baseline["scan_ux_m"]
    uy = defect["scan_uy_m"] .- baseline["scan_uy_m"]
    total = sqrt.(abs2.(ux) .+ abs2.(uy))
    ux, uy, total
end

function metrics(role, baseline, defect)
    x_mm = baseline["scan_x_mm"]
    receiver_index = findfirst(==(RECEIVER_X_MM), x_mm)
    isnothing(receiver_index) && error("receiver x=$RECEIVER_X_MM mm is outside scan")
    scattered_ux, scattered_uy, scattered_total = scattered_fields(baseline, defect)
    receiver_ux = abs.(scattered_ux[:, receiver_index])
    receiver_total = scattered_total[:, receiver_index]
    focus_displacement = baseline["focus_displacement"]
    (
        role=String(role),
        focus_ux_m=baseline["focus_longitudinal_amplitude_m"],
        focus_total_m=baseline["focus_total_amplitude_m"],
        local_kinetic_energy_density_j_m3=kinetic_energy_density(focus_displacement),
        receiver_x_mm=RECEIVER_X_MM,
        receiver_scattered_ux_rms_m=rms_finite(receiver_ux),
        receiver_scattered_ux_peak_m=maximum(filter(isfinite, receiver_ux)),
        receiver_scattered_total_rms_m=rms_finite(receiver_total),
        scan_scattered_total_rms_m=rms_finite(scattered_total),
        scattered_ux,
        scattered_uy,
        scattered_total,
        receiver_ux,
        receiver_total,
    )
end

function write_summary(selected, uniform)
    amplitude_gain = selected.focus_ux_m / uniform.focus_ux_m
    local_energy_gain = selected.local_kinetic_energy_density_j_m3 /
                        uniform.local_kinetic_energy_density_j_m3
    scattering_gain = selected.receiver_scattered_ux_rms_m /
                      uniform.receiver_scattered_ux_rms_m
    snr_improvement_db = 20log10(scattering_gain)
    path = joinpath(OUTPUT_ROOT, "material_lens_defect_summary.csv")
    open(path, "w") do io
        println(io, join((
            "frequency_hz", "source_traction_pa", "focus_x_mm", "defect_radius_mm",
            "receiver_x_mm", "selected_focus_ux_m", "uniform_focus_ux_m",
            "selected_local_kinetic_energy_density_j_m3",
            "uniform_local_kinetic_energy_density_j_m3",
            "selected_receiver_scattered_ux_rms_m",
            "uniform_receiver_scattered_ux_rms_m",
            "selected_receiver_scattered_ux_peak_m",
            "uniform_receiver_scattered_ux_peak_m",
            "lambda_p_mm", "lambda_s_mm", "focus_amplitude_gain",
            "local_energy_gain", "scattering_gain_equal_noise_snr_gain",
            "snr_improvement_db", "defect_diameter_over_lambda_s",
        ), ','))
        println(io, join((
            FREQUENCY_HZ, 1.0e6, DEFECT_X_MM, DEFECT_RADIUS_MM, RECEIVER_X_MM,
            selected.focus_ux_m, uniform.focus_ux_m,
            selected.local_kinetic_energy_density_j_m3,
            uniform.local_kinetic_energy_density_j_m3,
            selected.receiver_scattered_ux_rms_m,
            uniform.receiver_scattered_ux_rms_m,
            selected.receiver_scattered_ux_peak_m,
            uniform.receiver_scattered_ux_peak_m,
            1.0e3ALUMINIUM_CP_M_S / FREQUENCY_HZ,
            1.0e3ALUMINIUM_CS_M_S / FREQUENCY_HZ,
            amplitude_gain, local_energy_gain, scattering_gain,
            snr_improvement_db,
            2DEFECT_RADIUS_MM / (1.0e3ALUMINIUM_CS_M_S / FREQUENCY_HZ),
        ), ','))
    end
    (
        path,
        amplitude_gain,
        local_energy_gain,
        scattering_gain,
        snr_improvement_db,
    )
end

function write_receiver_profiles(y_mm, selected, uniform)
    path = joinpath(OUTPUT_ROOT, "receiver_scattered_profiles.csv")
    open(path, "w") do io
        println(io, "y_mm,selected_scattered_ux_m,uniform_scattered_ux_m,selected_scattered_total_m,uniform_scattered_total_m")
        for index in eachindex(y_mm)
            println(io, join((
                y_mm[index], selected.receiver_ux[index], uniform.receiver_ux[index],
                selected.receiver_total[index], uniform.receiver_total[index],
            ), ','))
        end
    end
    path
end

function make_plot(baseline_selected, selected, uniform, summary)
    x_mm = baseline_selected["scan_x_mm"]
    y_mm = baseline_selected["scan_y_mm"]
    focus_nm = 1.0e9 .* baseline_selected["scan_amplitude_m"]
    scattered_nm = 1.0e9 .* selected.scattered_total

    p1 = heatmap(
        x_mm, y_mm, focus_nm;
        xlabel="x, mm", ylabel="y, mm", colorbar_title="|u|, nm",
        title="Defect-free Al lens", aspect_ratio=:equal, c=:viridis,
    )
    scatter!(p1, [DEFECT_X_MM], [DEFECT_Y_MM]; marker=:circle, ms=5, mc=:white,
             markerstrokecolor=:black, label="future defect")

    p2 = heatmap(
        x_mm, y_mm, scattered_nm;
        xlabel="x, mm", ylabel="y, mm", colorbar_title="|u_sc|, nm",
        title="Scattered field: lens", aspect_ratio=:equal, c=:magma,
    )
    vline!(p2, [RECEIVER_X_MM]; color=:cyan, lw=2, ls=:dash, label="receiver")
    scatter!(p2, [DEFECT_X_MM], [DEFECT_Y_MM]; marker=:circle, ms=5, mc=:white,
             markerstrokecolor=:black, label="hole r=2 mm")

    p3 = plot(
        y_mm, 1.0e9 .* selected.receiver_ux;
        xlabel="y, mm", ylabel="|u_sc,x|, nm", lw=2.5,
        label="lens", title="Receiver line, x=68 mm",
    )
    plot!(p3, y_mm, 1.0e9 .* uniform.receiver_ux; lw=2.5, ls=:dash, label="straight")

    energy_values = [
        selected.local_kinetic_energy_density_j_m3,
        uniform.local_kinetic_energy_density_j_m3,
    ]
    p4 = bar(
        ["lens", "straight"], energy_values;
        ylabel="mean kinetic energy density, J/m³", title="At future defect centre",
        label=false, color=[:steelblue, :gray55],
    )
    annotate!(p4, 1.5, maximum(energy_values) * 0.86,
              text("energy gain = $(round(summary.local_energy_gain, digits=3))\nSNR proxy = $(round(summary.snr_improvement_db, digits=2)) dB", 9))

    path = joinpath(OUTPUT_ROOT, "material_lens_defect_analysis.png")
    savefig(plot(p1, p2, p3, p4; layout=(2, 2), size=(1400, 980)), path)
    path
end

function main()
    mkpath(OUTPUT_ROOT)
    baseline_selected, defect_selected = load_pair(:selected)
    baseline_uniform, defect_uniform = load_pair(:uniform)
    selected = metrics(:selected, baseline_selected, defect_selected)
    uniform = metrics(:uniform, baseline_uniform, defect_uniform)
    summary_path, amplitude_gain, local_energy_gain, scattering_gain, snr_improvement_db =
        write_summary(selected, uniform)
    receiver_path = write_receiver_profiles(baseline_selected["scan_y_mm"], selected, uniform)
    plot_path = make_plot(
        baseline_selected,
        selected,
        uniform,
        (; local_energy_gain, snr_improvement_db),
    )
    jld_path = joinpath(OUTPUT_ROOT, "material_lens_defect_analysis.jld2")
    jldsave(
        jld_path;
        format_version=1,
        selected,
        uniform,
        amplitude_gain,
        local_energy_gain,
        scattering_gain_equal_noise_snr_gain=scattering_gain,
        snr_improvement_db,
        receiver_x_mm=RECEIVER_X_MM,
        defect_radius_mm=DEFECT_RADIUS_MM,
        defect_x_mm=DEFECT_X_MM,
        defect_y_mm=DEFECT_Y_MM,
    )
    println("[+] Focus amplitude gain: $amplitude_gain")
    println("[+] Local kinetic-energy gain: $local_energy_gain")
    println("[+] Receiver scattered-field gain: $scattering_gain")
    println("[+] Equal-noise SNR improvement: $snr_improvement_db dB")
    println("[+] Summary: $summary_path")
    println("[+] Receiver profiles: $receiver_path")
    println("[+] Plot: $plot_path")
    println("[+] JLD2: $jld_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
