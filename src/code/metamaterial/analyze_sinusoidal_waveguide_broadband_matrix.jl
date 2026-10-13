module AnalyzeSinusoidalWaveguideBroadbandMatrix

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Statistics

include(joinpath(@__DIR__, "operator_time_reversal_audit.jl"))
using .OperatorTimeReversalAudit

export analyze_broadband_matrix,
       reflection_project,
       parity_bases,
       trapezoid_pulse_weights,
       power_whiten

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PACKAGE_ROOT = joinpath(
    PROJECT_ROOT, "results", "sinusoidal_material_transfer_p6_242khz",
)
const INPUT_ROOT = joinpath(
    PACKAGE_ROOT, "time_reversal_broadband_matrix", "frequency_matrices",
)
const OUTPUT_ROOT = joinpath(
    PACKAGE_ROOT, "time_reversal_broadband_matrix", "analysis",
)
const REFINEMENT_ROOT = joinpath(
    PACKAGE_ROOT, "time_reversal_broadband_matrix", "refinement_321700hz",
)
const FREQUENCIES_HZ = Float64[
    162_200, 178_160, 194_120, 202_100, 210_080, 218_060, 226_040,
    234_020, 242_000, 249_970, 257_940, 265_910, 273_880, 281_850,
    289_820, 305_760, 321_700,
]
const CARRIER_HZ = 242.0e3
const PRESSURE_SPEED_M_S = 6122.102437409232
const PRESSURE_WAVELENGTH_MM = 1.0e3 * PRESSURE_SPEED_M_S / CARRIER_HZ
const ROLES = (:sinusoidal, :straight)
const DYNAMIC_RANGES_DB = (20.0, 30.0, 40.0, 60.0)

result_path(role, frequency_hz) = joinpath(
    INPUT_ROOT, "$(role)_f$(round(Int, frequency_hz))hz.jld2",
)

reflection_project(matrix) = 0.5 .* (
    matrix + reverse(reverse(matrix; dims=1); dims=2)
)

function parity_bases(count::Integer)
    isodd(count) || throw(ArgumentError("parity basis requires an odd count"))
    centre = (count + 1) ÷ 2
    even = zeros(Float64, count, centre)
    odd = zeros(Float64, count, centre - 1)
    even[centre, 1] = 1.0
    for offset in 1:(centre - 1)
        lower = centre - offset
        upper = centre + offset
        even[lower, offset + 1] = inv(sqrt(2))
        even[upper, offset + 1] = inv(sqrt(2))
        odd[lower, offset] = inv(sqrt(2))
        odd[upper, offset] = -inv(sqrt(2))
    end
    even, odd
end

function trapezoid_pulse_weights(frequencies_hz; carrier_hz=CARRIER_HZ, cycles=5.0)
    issorted(frequencies_hz) || throw(ArgumentError("frequencies must be sorted"))
    length(frequencies_hz) >= 2 || throw(ArgumentError("at least two frequencies required"))
    all(>(0), diff(frequencies_hz)) || throw(ArgumentError("frequencies must be unique"))
    quadrature = similar(Float64.(frequencies_hz))
    quadrature[1] = (frequencies_hz[2] - frequencies_hz[1]) / 2
    for index in 2:(length(frequencies_hz) - 1)
        quadrature[index] = (frequencies_hz[index + 1] - frequencies_hz[index - 1]) / 2
    end
    quadrature[end] = (frequencies_hz[end] - frequencies_hz[end - 1]) / 2
    amplitude = hann_burst_spectrum(frequencies_hz ./ carrier_hz; cycles)
    weights = amplitude .* sqrt.(quadrature)
    weights ./= norm(weights)
    weights, quadrature
end

function power_whiten(work, frequency_hz; negative_tolerance=1.0e-8, floor_fraction=1.0e-12)
    size(work, 1) == size(work, 2) || throw(DimensionMismatch("work matrix must be square"))
    omega = 2pi * frequency_hz
    power = Matrix(Hermitian(
        (-omega / 2) .* ((work - work') ./ (2im)),
    ))
    decomposition = eigen(Hermitian(power))
    largest = maximum(decomposition.values)
    largest > 0 || error("accepted-power metric has no positive eigenvalue")
    minimum(decomposition.values) >= -negative_tolerance * largest || error(
        "accepted-power metric has a material negative eigenvalue",
    )
    floor_value = max(floor_fraction * largest, eps(largest))
    clipped = max.(decomposition.values, floor_value)
    inverse_root = decomposition.vectors * Diagonal(inv.(sqrt.(clipped))) *
        decomposition.vectors'
    whitening_residual = norm(inverse_root' * power * inverse_root -
        Matrix{ComplexF64}(I, size(power))) / sqrt(size(power, 1))
    (
        power,
        inverse_root,
        eigenvalues=decomposition.values,
        regularized_count=count(<(floor_value), decomposition.values),
        whitening_residual,
    )
end

function load_role(role)
    raw_transfer = Array{ComplexF64}(undef, 15, 15, length(FREQUENCIES_HZ))
    projected_transfer = similar(raw_transfer)
    raw_work = similar(raw_transfer)
    projected_work = similar(raw_transfer)
    whitened_transfer = similar(raw_transfer)
    frequency_rows = NamedTuple[]
    coordinates_mm = nothing
    for (frequency_index, frequency_hz) in enumerate(FREQUENCIES_HZ)
        path = result_path(role, frequency_hz)
        isfile(path) || error("missing broadband matrix: $path")
        data = JLD2.load(path)
        data["format_version"] == 1 || error("unsupported matrix format")
        H = ComplexF64.(data["output_ux_m"])
        W = ComplexF64.(data["port_work_compliance"])
        size(H) == (15, 15) || throw(DimensionMismatch("expected 15 x 15 transfer"))
        size(W) == (15, 15) || throw(DimensionMismatch("expected 15 x 15 work matrix"))
        y_mm = Float64.(data["output_y_mm"])
        isnothing(coordinates_mm) ? (coordinates_mm = y_mm) :
            y_mm == coordinates_mm || error("output coordinates changed across frequency")
        Hs = reflection_project(H)
        Ws = reflection_project(W)
        whitened = power_whiten(Ws, frequency_hz)
        Hw = Hs * whitened.inverse_root
        raw_transfer[:, :, frequency_index] .= H
        projected_transfer[:, :, frequency_index] .= Hs
        raw_work[:, :, frequency_index] .= W
        projected_work[:, :, frequency_index] .= Ws
        whitened_transfer[:, :, frequency_index] .= Hw

        even, odd = parity_bases(15)
        even_transfer = even' * Hw * even
        odd_transfer = odd' * Hw * odd
        parity_cross = hypot(norm(even' * Hw * odd), norm(odd' * Hw * even)) /
            max(norm(Hw), eps(Float64))
        full_summary = singular_summary(Hw, DYNAMIC_RANGES_DB)
        even_summary = singular_summary(even_transfer, DYNAMIC_RANGES_DB)
        odd_summary = singular_summary(odd_transfer, DYNAMIC_RANGES_DB)
        projected_reflection = norm(
            reverse(reverse(Hs; dims=1); dims=2) - Hs,
        ) / max(norm(Hs), eps(Float64))
        push!(frequency_rows, (
            role, frequency_hz,
            raw_reflection_residual=data["reflection_residual"],
            projected_reflection_residual=projected_reflection,
            reciprocity_residual=data["reciprocity_residual"],
            minimum_diagonal_accepted_work=minimum(-imag.(diag(Ws))),
            power_minimum_eigenvalue=minimum(whitened.eigenvalues),
            power_maximum_eigenvalue=maximum(whitened.eigenvalues),
            power_condition_number=maximum(whitened.eigenvalues) /
                max(minimum(whitened.eigenvalues), eps(Float64)),
            power_regularized_modes=whitened.regularized_count,
            power_whitening_residual=whitened.whitening_residual,
            parity_cross_residual=parity_cross,
            stable_rank=full_summary.stable_rank,
            modes_20db=full_summary.observable_modes[20.0],
            even_stable_rank=even_summary.stable_rank,
            even_modes_20db=even_summary.observable_modes[20.0],
            odd_stable_rank=odd_summary.stable_rank,
            odd_modes_20db=odd_summary.observable_modes[20.0],
            raw_reflection_gate_passed=data["reflection_residual"] <= 0.02,
            core_numerical_gate_passed=projected_reflection <= 1.0e-12 &&
                data["reciprocity_residual"] <= 1.0e-8 &&
                minimum(-imag.(diag(Ws))) > 0 &&
                minimum(whitened.eigenvalues) >=
                    -1.0e-8 * maximum(whitened.eigenvalues),
        ))
    end
    (
        role,
        coordinates_mm,
        raw_transfer,
        projected_transfer,
        raw_work,
        projected_work,
        whitened_transfer,
        frequency_rows,
    )
end

function refinement_rows()
    rows = NamedTuple[]
    frequency_hz = last(FREQUENCIES_HZ)
    for role in ROLES
        base_path = result_path(role, frequency_hz)
        refined_path = joinpath(REFINEMENT_ROOT, "$(role)_matrix.jld2")
        isfile(base_path) || error("missing base worst-case matrix: $base_path")
        isfile(refined_path) || error("missing refined worst-case matrix: $refined_path")
        base = JLD2.load(base_path)
        refined = JLD2.load(refined_path)
        base_H = reflection_project(ComplexF64.(base["output_ux_m"]))
        refined_H = reflection_project(ComplexF64.(refined["output_ux_m"]))
        base_W = reflection_project(ComplexF64.(base["port_work_compliance"]))
        refined_W = reflection_project(ComplexF64.(refined["port_work_compliance"]))
        base_power = power_whiten(base_W, frequency_hz)
        refined_power = power_whiten(refined_W, frequency_hz)
        base_whitened = base_H * base_power.inverse_root
        refined_whitened = refined_H * refined_power.inverse_root
        base_kernel = base_whitened * base_whitened'
        refined_kernel = refined_whitened * refined_whitened'
        base_kernel ./= real(tr(base_kernel))
        refined_kernel ./= real(tr(refined_kernel))
        kernel_error = norm(base_kernel - refined_kernel) /
            max(norm(refined_kernel), eps(Float64))
        base_singular = svdvals(base_whitened)
        refined_singular = svdvals(refined_whitened)
        base_singular ./= first(base_singular)
        refined_singular ./= first(refined_singular)
        singular_spectrum_error = norm(base_singular - refined_singular) /
            max(norm(refined_singular), eps(Float64))
        power_positive = minimum(refined_power.eigenvalues) >=
            -1.0e-8 * maximum(refined_power.eigenvalues)
        push!(rows, (
            role, frequency_hz,
            base_raw_reflection_residual=base["reflection_residual"],
            refined_raw_reflection_residual=refined["reflection_residual"],
            refined_reciprocity_residual=refined["reciprocity_residual"],
            refined_power_minimum_eigenvalue=minimum(refined_power.eigenvalues),
            refined_power_maximum_eigenvalue=maximum(refined_power.eigenvalues),
            normalized_kernel_relative_error=kernel_error,
            normalized_singular_spectrum_relative_error=singular_spectrum_error,
            passed=refined["reflection_residual"] <= 0.02 &&
                refined["reciprocity_residual"] <= 1.0e-8 &&
                power_positive && kernel_error <= 0.10 &&
                singular_spectrum_error <= 0.10,
        ))
    end
    rows
end

function band_operator(transfer, weights, output_basis, input_basis)
    blocks = [
        weights[index] .* (output_basis' * transfer[:, :, index] * input_basis)
        for index in eachindex(weights)
    ]
    hcat(blocks...)
end

function analyze_role(data, weights)
    even, odd = parity_bases(15)
    identity_basis = Matrix{Float64}(I, 15, 15)
    operators = Dict(
        :full => band_operator(data.whitened_transfer, weights, identity_basis, identity_basis),
        :even => band_operator(data.whitened_transfer, weights, even, even),
        :odd => band_operator(data.whitened_transfer, weights, odd, odd),
    )
    summaries = Dict(
        sector => singular_summary(operator, DYNAMIC_RANGES_DB)
        for (sector, operator) in operators
    )
    kernel = time_reversal_kernel(operators[:full])
    spacing_lambda = (data.coordinates_mm[2] - data.coordinates_mm[1]) /
        PRESSURE_WAVELENGTH_MM
    focuses = [
        focus_metrics(kernel, focal_bin; bin_spacing=spacing_lambda)
        for focal_bin in eachindex(data.coordinates_mm)
    ]
    (; data, operators, summaries, kernel, focuses)
end

function summary_rows(results)
    rows = NamedTuple[]
    for role in ROLES, sector in (:full, :even, :odd)
        summary = results[role].summaries[sector]
        focus = results[role].focuses
        push!(rows, (
            role, sector,
            output_dimension=size(results[role].operators[sector], 1),
            operator_columns=size(results[role].operators[sector], 2),
            stable_rank=summary.stable_rank,
            entropy_rank=summary.entropy_rank,
            modes_20db=summary.observable_modes[20.0],
            modes_30db=summary.observable_modes[30.0],
            modes_40db=summary.observable_modes[40.0],
            modes_60db=summary.observable_modes[60.0],
            median_fwhm_lambda=sector == :full ?
                median(getproperty.(focus, :fwhm)) : NaN,
            median_focus_contrast_db=sector == :full ?
                median(getproperty.(focus, :focus_contrast_db)) : NaN,
            median_peak_outside_fwhm_db=sector == :full ?
                median(getproperty.(focus, :peak_outside_fwhm_db)) : NaN,
            all_targets_are_global_maxima=sector == :full ?
                all(getproperty.(focus, :target_is_global_maximum)) : false,
        ))
    end
    rows
end

function write_csv(path, rows)
    isempty(rows) && return
    columns = propertynames(first(rows))
    csv_value(value) = value isa AbstractFloat && !isfinite(value) ? "NaN" : string(value)
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((csv_value(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function plot_spectra(results)
    panels = Any[]
    for sector in (:full, :even, :odd)
        panel = plot(
            xlabel="mode index", ylabel="relative singular value / dB",
            title=String(sector), ylim=(-65, 2), legend=:bottomleft,
        )
        for (role, label, color) in (
            (:sinusoidal, "sinusoidal", :darkorange),
            (:straight, "straight", :steelblue),
        )
            values = results[role].summaries[sector].singular_value_db
            plot!(panel, eachindex(values), values;
                marker=:circle, linewidth=2.2, label, color)
        end
        hline!(panel, [-20.0, -40.0]; color=:gray50, linestyle=:dash,
            label="dynamic range")
        push!(panels, panel)
    end
    savefig(plot(
        panels...; layout=(1, 3), size=(1500, 470),
        left_margin=5Plots.mm, bottom_margin=5Plots.mm,
    ), joinpath(OUTPUT_ROOT, "broadband_singular_spectra.png"))
end

function plot_kernels(results, coordinates_mm)
    panels = Any[]
    for (role, title) in ((:sinusoidal, "sinusoidal"), (:straight, "straight"))
        push!(panels, heatmap(
            coordinates_mm, coordinates_mm, coherence_kernel(results[role].kernel);
            xlabel="target y / mm", ylabel="observation y / mm", title,
            aspect_ratio=:equal, color=:viridis, clims=(0, 1), colorbar_title="|C|",
        ))
    end
    savefig(plot(
        panels...; layout=(1, 2), size=(1160, 510),
        left_margin=6Plots.mm, bottom_margin=5Plots.mm,
    ), joinpath(OUTPUT_ROOT, "broadband_time_reversal_coherence.png"))
end

function plot_focus_profiles(results, coordinates_mm)
    target_bins = (8, 12)
    panels = Any[]
    for focal_bin in target_bins
        panel = plot(
            xlabel="observation y / mm", ylabel="normalized intensity / dB",
            title="target y = $(coordinates_mm[focal_bin]) mm",
            ylim=(-55, 2), legend=:bottomleft,
        )
        for (role, label, color) in (
            (:sinusoidal, "sinusoidal", :darkorange),
            (:straight, "straight", :steelblue),
        )
            profile = normalized_focus_intensity(results[role].kernel, focal_bin)
            plot!(panel, coordinates_mm, 10 .* log10.(max.(profile, 10.0^(-5.5)));
                marker=:circle, linewidth=2.2, label, color)
        end
        vline!(panel, [coordinates_mm[focal_bin]];
            color=:black, linestyle=:dot, label="target")
        push!(panels, panel)
    end
    savefig(plot(
        panels...; layout=(1, 2), size=(1160, 480),
        left_margin=6Plots.mm, bottom_margin=5Plots.mm,
    ), joinpath(OUTPUT_ROOT, "broadband_focus_profiles.png"))
end

function plot_frequency_modes(role_data)
    figure = plot(
        xlabel="frequency / kHz", ylabel="stable rank",
        title="Power-normalized spatial rank by frequency",
        legend=:best, size=(920, 560),
    )
    for (role, label, color) in (
        (:sinusoidal, "sinusoidal", :darkorange),
        (:straight, "straight", :steelblue),
    )
        rows = role_data[role].frequency_rows
        plot!(figure, getproperty.(rows, :frequency_hz) ./ 1.0e3,
            getproperty.(rows, :stable_rank);
            marker=:circle, linewidth=2.2, label, color)
    end
    savefig(figure, joinpath(OUTPUT_ROOT, "frequency_stable_rank.png"))
end

function write_readme(rows, comparison, role_data, refinement)
    table = Dict((row.role, row.sector) => row for row in rows)
    lens = table[(:sinusoidal, :full)]
    straight = table[(:straight, :full)]
    core_gates = all(
        row.core_numerical_gate_passed
        for role in ROLES for row in role_data[role].frequency_rows
    )
    raw_failures = count(
        !row.raw_reflection_gate_passed
        for role in ROLES for row in role_data[role].frequency_rows
    )
    refinement_passed = all(getproperty.(refinement, :passed))
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# Full broadband response matrix: sinusoidal versus straight Al waveguides")
        println(io)
        println(io, "## Result")
        println(io)
        println(io, "All projected-symmetry, reciprocity, and accepted-power gates passed: `$(core_gates)`. The base grid has `$(raw_failures)` raw reflection failures above 2%; the preregistered upper-edge refinement control passed: `$(refinement_passed)`. Each frequency uses one FEM factorization and 15 source backsolves. Input channels are whitened by the full accepted-power quadratic form before the five-cycle band is assembled.")
        println(io)
        println(io, "| geometry | sector | stable rank | entropy rank | modes @20 dB | @30 dB | @40 dB | median FWHM / lambdaP | focus contrast / dB |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for role in ROLES, sector in (:full, :even, :odd)
            row = table[(role, sector)]
            println(io, "| $(role) | $(sector) | $(round(row.stable_rank; digits=3)) | $(round(row.entropy_rank; digits=3)) | $(row.modes_20db) | $(row.modes_30db) | $(row.modes_40db) | $(sector == :full ? round(row.median_fwhm_lambda; digits=3) : "-") | $(sector == :full ? round(row.median_focus_contrast_db; digits=2) : "-") |")
        end
        println(io)
        println(io, "The sinusoidal/full stable-rank ratio is `$(round(comparison.stable_rank_ratio; digits=3))`; its 30 dB mode-count change is `$(comparison.modes_30db_change)`, median FWHM ratio is `$(round(comparison.fwhm_ratio; digits=3))`, and focus-contrast change is `$(round(comparison.focus_contrast_change_db; digits=2)) dB`.")
        println(io)
        if comparison.stable_rank_ratio > 1
            println(io, "After exact accepted-power normalization and broadband combination, the sinusoidal geometry adds effective spatial diversity relative to the straight control.")
        else
            println(io, "Even after exact accepted-power normalization and broadband combination, the sinusoidal geometry does not add effective spatial diversity relative to the straight control.")
        end
        println(io)
        println(io, "## Boundaries")
        println(io)
        println(io, "This is a deterministic full-wave channel and time-reversal audit at a plane five millimetres after the strips. It includes signed positions and both parity sectors, but no stochastic electronics noise, clutter, defect ensemble, or new Pd estimate. The result diagnoses the waveguide operator; it does not reopen geometry optimization.")
    end
end

function analyze_broadband_matrix()
    mkpath(OUTPUT_ROOT)
    weights, quadrature_hz = trapezoid_pulse_weights(FREQUENCIES_HZ)
    role_data = Dict(role => load_role(role) for role in ROLES)
    all(row.core_numerical_gate_passed for role in ROLES for row in role_data[role].frequency_rows) ||
        error("at least one projected/core numerical gate failed")
    refinement = refinement_rows()
    all(getproperty.(refinement, :passed)) || error("worst-case refinement control failed")
    role_data[:sinusoidal].coordinates_mm == role_data[:straight].coordinates_mm ||
        error("geometry output coordinates differ")
    results = Dict(
        role => analyze_role(role_data[role], weights) for role in ROLES
    )
    rows = summary_rows(results)
    table = Dict((row.role, row.sector) => row for row in rows)
    lens = table[(:sinusoidal, :full)]
    straight = table[(:straight, :full)]
    comparison = (
        stable_rank_ratio=lens.stable_rank / straight.stable_rank,
        entropy_rank_ratio=lens.entropy_rank / straight.entropy_rank,
        modes_20db_change=lens.modes_20db - straight.modes_20db,
        modes_30db_change=lens.modes_30db - straight.modes_30db,
        modes_40db_change=lens.modes_40db - straight.modes_40db,
        fwhm_ratio=lens.median_fwhm_lambda / straight.median_fwhm_lambda,
        focus_contrast_change_db=
            lens.median_focus_contrast_db - straight.median_focus_contrast_db,
        peak_outside_fwhm_change_db=
            lens.median_peak_outside_fwhm_db - straight.median_peak_outside_fwhm_db,
    )
    frequency_rows = reduce(vcat, [role_data[role].frequency_rows for role in ROLES])
    pulse_rows = [
        (
            frequency_hz=FREQUENCIES_HZ[index],
            trapezoid_width_hz=quadrature_hz[index],
            pulse_amplitude_weight=weights[index],
        )
        for index in eachindex(weights)
    ]
    focus_rows = NamedTuple[]
    for role in ROLES, (focal_bin, metrics) in enumerate(results[role].focuses)
        push!(focus_rows, (
            role, focal_bin,
            focal_y_mm=role_data[role].coordinates_mm[focal_bin],
            peak_bin=metrics.peak_bin,
            target_is_global_maximum=metrics.target_is_global_maximum,
            fwhm_bins=metrics.fwhm_bins,
            fwhm_lambda=metrics.fwhm,
            focus_contrast_db=metrics.focus_contrast_db,
            peak_outside_fwhm_db=metrics.peak_outside_fwhm_db,
        ))
    end
    singular_rows = NamedTuple[]
    for role in ROLES, sector in (:full, :even, :odd)
        summary = results[role].summaries[sector]
        for mode_index in eachindex(summary.singular_values)
            push!(singular_rows, (
                role, sector, mode_index,
                singular_value=summary.singular_values[mode_index],
                relative_db=summary.singular_value_db[mode_index],
            ))
        end
    end

    write_csv(joinpath(OUTPUT_ROOT, "architecture_summary.csv"), rows)
    write_csv(joinpath(OUTPUT_ROOT, "comparison.csv"), [comparison])
    write_csv(joinpath(OUTPUT_ROOT, "frequency_metrics.csv"), frequency_rows)
    write_csv(joinpath(OUTPUT_ROOT, "refinement_control.csv"), refinement)
    write_csv(joinpath(OUTPUT_ROOT, "pulse_weights.csv"), pulse_rows)
    write_csv(joinpath(OUTPUT_ROOT, "focus_metrics.csv"), focus_rows)
    write_csv(joinpath(OUTPUT_ROOT, "singular_values.csv"), singular_rows)
    plot_spectra(results)
    plot_kernels(results, role_data[:sinusoidal].coordinates_mm)
    plot_focus_profiles(results, role_data[:sinusoidal].coordinates_mm)
    plot_frequency_modes(role_data)
    write_readme(rows, comparison, role_data, refinement)
    jldsave(
        joinpath(OUTPUT_ROOT, "broadband_response_matrix_audit.jld2");
        frequencies_hz=FREQUENCIES_HZ, weights, quadrature_hz,
        coordinates_mm=role_data[:sinusoidal].coordinates_mm,
        role_data, results, rows, comparison, refinement,
    )
    println("[+] Full stable-rank ratio: $(comparison.stable_rank_ratio)")
    println("[+] 30 dB mode-count change: $(comparison.modes_30db_change)")
    println("[+] FWHM ratio: $(comparison.fwhm_ratio)")
    println("[+] Focus-contrast change: $(comparison.focus_contrast_change_db) dB")
    rows
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    analyze_broadband_matrix()
end

end
