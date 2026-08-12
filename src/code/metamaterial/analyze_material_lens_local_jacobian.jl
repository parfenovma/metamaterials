module AnalyzeMaterialLensLocalJacobian

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BASE_ROOT = joinpath(PROJECT_ROOT, "tmp", "aluminium_material_lens_response_matrix_242khz")
const PERTURBED_ROOT = joinpath(PROJECT_ROOT, "tmp", "aluminium_material_lens_jacobian_g7plus_242khz")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_JACOBIAN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_material_lens_local_jacobian_242khz"),
)
const DELTA_GAP_MM = 0.25
const PERTURBED_GROUP = 7
const CENTERS_MM = collect(0.0:8.2:57.4)

matrix_path(root) = joinpath(root, "aluminium_strip_group_response_matrix.jld2")

wrap_phase(value) = mod(value + pi, 2pi) - pi

function amplitude_directional_derivative(value, derivative)
    abs(value) > 0 || return NaN
    real(conj(value) * derivative) / abs(value)
end

function write_group_csv(base, perturbed, derivative)
    path = joinpath(OUTPUT_ROOT, "local_gap_jacobian_by_group.csv")
    open(path, "w") do io
        println(io, "group,abs_y_mm,base_focus_ux_abs_m,perturbed_focus_ux_abs_m,amplitude_derivative_m_per_mm,phase_derivative_deg_per_mm,complex_derivative_real_m_per_mm,complex_derivative_imag_m_per_mm")
        for group in eachindex(base)
            phase_derivative = rad2deg(wrap_phase(angle(perturbed[group] / base[group]))) /
                               DELTA_GAP_MM
            println(io, join((
                group, CENTERS_MM[group], abs(base[group]), abs(perturbed[group]),
                (abs(perturbed[group]) - abs(base[group])) / DELTA_GAP_MM,
                phase_derivative, real(derivative[group]), imag(derivative[group]),
            ), ','))
        end
    end
    path
end

function make_plot(base, perturbed, derivative, metrics)
    amplitude_derivative_nm_mm = 1.0e9 .* (abs.(perturbed) .- abs.(base)) ./ DELTA_GAP_MM
    phase_derivative_deg_mm = rad2deg.([
        wrap_phase(angle(perturbed[index] / base[index])) / DELTA_GAP_MM
        for index in eachindex(base)
    ])
    p1 = bar(
        CENTERS_MM, amplitude_derivative_nm_mm;
        xlabel="|y|, mm", ylabel="d|C_i|/dg₇, nm/mm",
        title="Coupled amplitude Jacobian", label=false,
        color=[index == PERTURBED_GROUP ? :darkorange : :steelblue for index in eachindex(base)],
        bar_width=5.5,
    )
    p2 = bar(
        CENTERS_MM, phase_derivative_deg_mm;
        xlabel="|y|, mm", ylabel="d phase(C_i)/dg₇, deg/mm",
        title="Coupled phase Jacobian", label=false,
        color=[index == PERTURBED_GROUP ? :darkorange : :seagreen for index in eachindex(base)],
        bar_width=5.5,
    )
    p3 = scatter(
        real.(derivative) .* 1.0e9,
        imag.(derivative) .* 1.0e9;
        xlabel="Re(dC/dg₇), nm/mm", ylabel="Im(dC/dg₇), nm/mm",
        title="Complex focal Jacobian", marker_z=1:length(base), ms=8,
        color=:viridis, colorbar_title="group", label=false,
    )
    hline!(p3, [0.0]; color=:gray60, ls=:dot, label=false)
    vline!(p3, [0.0]; color=:gray60, ls=:dot, label=false)
    p4 = bar(
        ["base", "g₇ + 0.25 mm", "target"],
        [metrics.base_gain, metrics.perturbed_gain, 2.0];
        ylabel="carrier gain", title="One-sided physical step",
        color=[:steelblue, :darkorange, :gray60], label=false,
    )
    path = joinpath(OUTPUT_ROOT, "local_gap_jacobian.png")
    savefig(plot(p1, p2, p3, p4; layout=(2, 2), size=(1400, 950), margin=5Plots.mm), path)
    path
end

function main()
    mkpath(OUTPUT_ROOT)
    base_data = load(matrix_path(BASE_ROOT))
    perturbed_data = load(matrix_path(PERTURBED_ROOT))
    base = ComplexF64.(base_data["selected"].focus)
    perturbed = ComplexF64.(perturbed_data["selected"].focus)
    uniform = ComplexF64.(base_data["uniform"].focus)
    derivative = (perturbed .- base) ./ DELTA_GAP_MM
    base_total = sum(base)
    perturbed_total = sum(perturbed)
    uniform_total = sum(uniform)
    base_gain = abs(base_total) / abs(uniform_total)
    perturbed_gain = abs(perturbed_total) / abs(uniform_total)
    focus_amplitude_derivative_m_per_mm =
        amplitude_directional_derivative(base_total, sum(derivative))
    gain_derivative_per_mm = focus_amplitude_derivative_m_per_mm / abs(uniform_total)
    perturbed_response_fraction = abs2(derivative[PERTURBED_GROUP]) / sum(abs2, derivative)
    cross_coupled_response_fraction = 1 - perturbed_response_fraction
    other_groups = filter(!=(PERTURBED_GROUP), collect(eachindex(derivative)))
    cross_to_direct_l2_ratio = norm(derivative[other_groups]) /
                               abs(derivative[PERTURBED_GROUP])
    base_phase_error = wrap_phase(angle(base[PERTURBED_GROUP]) - angle(base_total))
    perturbed_phase_error = wrap_phase(angle(perturbed[PERTURBED_GROUP]) - angle(perturbed_total))
    metrics = (
        perturbed_group=PERTURBED_GROUP,
        base_gap_mm=1.5,
        perturbed_gap_mm=1.75,
        delta_gap_mm=DELTA_GAP_MM,
        base_focus_ux_m=abs(base_total),
        perturbed_focus_ux_m=abs(perturbed_total),
        focus_amplitude_relative_change=(abs(perturbed_total) / abs(base_total) - 1),
        focus_amplitude_derivative_m_per_mm,
        base_gain,
        perturbed_gain,
        gain_derivative_per_mm,
        base_group7_phase_error_deg=rad2deg(base_phase_error),
        perturbed_group7_phase_error_deg=rad2deg(perturbed_phase_error),
        group7_phase_error_derivative_deg_per_mm=
            rad2deg(wrap_phase(perturbed_phase_error - base_phase_error)) / DELTA_GAP_MM,
        cross_coupled_response_fraction,
        cross_to_direct_l2_ratio,
    )
    summary_path = joinpath(OUTPUT_ROOT, "local_gap_jacobian_summary.csv")
    open(summary_path, "w") do io
        columns = propertynames(metrics)
        println(io, join(string.(columns), ','))
        println(io, join((getproperty(metrics, column) for column in columns), ','))
    end
    group_path = write_group_csv(base, perturbed, derivative)
    plot_path = make_plot(base, perturbed, derivative, metrics)
    data_path = joinpath(OUTPUT_ROOT, "local_gap_jacobian.jld2")
    jldsave(data_path; format_version=1, base, perturbed, uniform, derivative, metrics)
    println("[+] Focus relative change: $(100 * metrics.focus_amplitude_relative_change)%")
    println("[+] Gain: $(metrics.base_gain) -> $(metrics.perturbed_gain)")
    println("[+] d(gain)/dg7: $(metrics.gain_derivative_per_mm) per mm")
    println("[+] Group-7 phase error: $(metrics.base_group7_phase_error_deg) -> $(metrics.perturbed_group7_phase_error_deg) deg")
    println("[+] Cross-coupled derivative fraction: $(metrics.cross_coupled_response_fraction)")
    println("[+] Summary: $summary_path")
    println("[+] Per-group Jacobian: $group_path")
    println("[+] Plot: $plot_path")
    println("[+] Data: $data_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
