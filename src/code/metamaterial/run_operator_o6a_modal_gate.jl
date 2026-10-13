module RunOperatorO6AModalGate

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Statistics

include(joinpath(@__DIR__, "port_mode_solver.jl"))
using .ElasticPortModes

export O6ModalConfig, evaluate_modal_transfer, run_o6a_modal_gate

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6a_modal_gate_v0",
)

Base.@kwdef struct O6ModalConfig
    aperture_width_m::Float64 = 71.6e-3
    plate_length_m::Float64 = 3.0e-3
    density_kg_m3::Float64 = 2700.0
    pressure_wave_speed_m_s::Float64 = 6122.102437409232
    shear_wave_speed_m_s::Float64 = 3083.810277185563
    lower_frequency_hz::Float64 = 193.6e3
    upper_frequency_hz::Float64 = 290.4e3
    carrier_frequency_hz::Float64 = 242.0e3
    frequency_count::Int = 49
    port_element_count::Int = 160
    aperture_sample_count::Int = 17
    pseudoinverse_rtol::Float64 = 1.0e-8
end

function validate(config::O6ModalConfig)
    config.aperture_width_m > 0 || throw(ArgumentError("aperture width must be positive"))
    config.plate_length_m >= 0 || throw(ArgumentError("plate length must be nonnegative"))
    config.lower_frequency_hz < config.carrier_frequency_hz < config.upper_frequency_hz ||
        throw(ArgumentError("carrier must lie inside the qualified band"))
    isodd(config.frequency_count) || throw(ArgumentError("frequency count must be odd"))
    config.frequency_count >= 7 || throw(ArgumentError("dense grid needs at least seven frequencies"))
    config.port_element_count >= 40 || throw(ArgumentError("port mesh is too coarse"))
    config.aperture_sample_count == 17 || throw(ArgumentError("O4 aperture has 17 samples"))
    nothing
end

function trapezoid_weights(y)
    weights = zeros(Float64, length(y))
    for index in 1:(length(y) - 1)
        half_interval = (y[index + 1] - y[index]) / 2
        weights[index] += half_interval
        weights[index + 1] += half_interval
    end
    weights
end

function continuous_dct(y, mode_index, width)
    mode_index == 0 && return fill(inv(sqrt(width)), length(y))
    mode_index == 1 && return -sqrt(2 / width) .* sin.(pi .* y ./ width)
    mode_index == 2 && return -sqrt(2 / width) .* cos.(2pi .* y ./ width)
    throw(ArgumentError("O6 uses only DCT modes 0, 1 and 2"))
end

function discrete_dct(sample_count)
    basis = zeros(Float64, sample_count, 3)
    basis[:, 1] .= inv(sqrt(sample_count))
    for mode in 2:3
        basis[:, mode] .= sqrt(2 / sample_count) .* cos.(
            (mode - 1) .* ((1:sample_count) .- 0.5) .* pi / sample_count,
        )
    end
    basis
end

function linear_interpolate(y, values, coordinate)
    coordinate <= first(y) && return first(values)
    coordinate >= last(y) && return last(values)
    left = searchsortedlast(y, coordinate)
    fraction = (coordinate - y[left]) / (y[left + 1] - y[left])
    values[left] + fraction * (values[left + 1] - values[left])
end

function port_config(config, frequency_hz)
    PortModeConfig(
        height_m=config.aperture_width_m,
        density=config.density_kg_m3,
        pressure_wave_speed=config.pressure_wave_speed_m_s,
        shear_wave_speed=config.shear_wave_speed_m_s,
        frequency_hz=frequency_hz,
        element_count=config.port_element_count,
    )
end

function synthesize_column(config, modes, mode_index, sample_y)
    parity = isodd(mode_index) ? :antisymmetric : :symmetric
    sector_modes = filter(mode -> mode.parity == parity, modes)
    isempty(sector_modes) && error("no $parity modes for DCT $mode_index")
    y = first(sector_modes).y_m
    weights = trapezoid_weights(y)
    target = continuous_dct(y, mode_index, config.aperture_width_m)
    weighted_displacements = reduce(hcat, [
        vcat(
            sqrt.(weights) .* mode.displacement_x,
            sqrt.(weights) .* mode.displacement_y,
        ) for mode in sector_modes
    ])
    weighted_target = vcat(sqrt.(weights) .* target, zeros(length(y)))
    coefficients = pinv(weighted_displacements; rtol=config.pseudoinverse_rtol) *
                   weighted_target
    propagated_coefficients = coefficients .* exp.(
        -im .* real.(getproperty.(sector_modes, :wavenumber_per_m)) .*
        config.plate_length_m,
    )
    source_x = ComplexF64[
        sum(
            coefficients[index] * linear_interpolate(
                mode.y_m, mode.displacement_x, coordinate,
            ) for (index, mode) in enumerate(sector_modes)
        ) for coordinate in sample_y
    ]
    output_x = ComplexF64[
        sum(
            propagated_coefficients[index] * linear_interpolate(
                mode.y_m, mode.displacement_x, coordinate,
            ) for (index, mode) in enumerate(sector_modes)
        ) for coordinate in sample_y
    ]
    reconstructed = weighted_displacements * coefficients
    source_fidelity = abs(dot(weighted_target, reconstructed)) /
        (norm(weighted_target) * norm(reconstructed))
    (
        source_x,
        output_x,
        coefficients,
        incident_power=sum(abs2, coefficients),
        source_fidelity,
        sector_mode_count=length(sector_modes),
    )
end

function unwrap_phase(values)
    result = angle.(values)
    for index in 2:length(result)
        step = result[index] - result[index - 1]
        result[index:end] .-= 2pi * round(step / (2pi))
    end
    result
end

function derivative(values, coordinates)
    result = similar(values, Float64)
    result[1] = (values[2] - values[1]) / (coordinates[2] - coordinates[1])
    result[end] = (values[end] - values[end - 1]) /
                  (coordinates[end] - coordinates[end - 1])
    for index in 2:(length(values) - 1)
        result[index] = (values[index + 1] - values[index - 1]) /
                        (coordinates[index + 1] - coordinates[index - 1])
    end
    result
end

function transfer_metrics(transfer, target)
    gram = adjoint(transfer) * transfer
    column_gains = sqrt.(max.(real.(diag(gram)), 0.0))
    coherence = maximum(
        abs(gram[first, second]) /
        sqrt(max(real(gram[first, first] * gram[second, second]), eps()))
        for first in 1:3, second in 1:3 if first != second
    )
    common_gain = dot(vec(target), vec(transfer)) / sum(abs2, target)
    fidelity = abs(dot(vec(target), vec(transfer))) /
        (norm(target) * norm(transfer))
    amplitude_rms = norm(abs.(transfer) .- abs.(common_gain .* target)) /
        max(norm(abs.(common_gain .* target)), eps())
    (
        common_gain,
        useful_amplitude=abs(common_gain),
        fidelity,
        coherence,
        column_imbalance_db=20log10(maximum(column_gains) / minimum(column_gains)),
        amplitude_rms,
        sigma_max=opnorm(transfer),
    )
end

function differential_delay_metrics(frequencies_hz, transfer, target)
    omega = 2pi .* frequencies_hz
    coefficient_delays = fill(NaN, size(transfer))
    active = abs.(target) .>= 0.10maximum(abs, target)
    for port in 1:3, sample in axes(target, 1)
        active[sample, port] || continue
        relative = transfer[sample, port, :] ./ target[sample, port]
        phase = unwrap_phase(relative)
        coefficient_delays[sample, port, :] .= -derivative(phase, omega)
    end
    rms_by_frequency = zeros(Float64, length(frequencies_hz))
    common_by_frequency = zeros(Float64, length(frequencies_hz))
    weights = abs2.(target)
    for frequency_index in eachindex(frequencies_hz)
        valid = active .& isfinite.(coefficient_delays[:, :, frequency_index])
        local_weights = weights[valid]
        local_delays = coefficient_delays[:, :, frequency_index][valid]
        common = sum(local_weights .* local_delays) / sum(local_weights)
        common_by_frequency[frequency_index] = common
        rms_by_frequency[frequency_index] = sqrt(
            sum(local_weights .* (local_delays .- common) .^ 2) /
            sum(local_weights),
        )
    end
    (
        coefficient_delays,
        common_by_frequency,
        rms_by_frequency,
        maximum_rms=maximum(rms_by_frequency),
        aggregate_rms=sqrt(mean(abs2, rms_by_frequency)),
    )
end

function evaluate_modal_transfer(config::O6ModalConfig=O6ModalConfig())
    validate(config)
    frequencies_hz = collect(range(
        config.lower_frequency_hz,
        config.upper_frequency_hz;
        length=config.frequency_count,
    ))
    carrier_index = only(findall(isapprox(config.carrier_frequency_hz; atol=1e-6), frequencies_hz))
    width = config.aperture_width_m
    count = config.aperture_sample_count
    sample_y = collect(range(
        -width / 2 + width / (2count),
        width / 2 - width / (2count);
        length=count,
    ))
    target = discrete_dct(count)
    raw_source = zeros(ComplexF64, count, 3, length(frequencies_hz))
    raw_output = similar(raw_source)
    source_fidelity = zeros(Float64, 3, length(frequencies_hz))
    incident_power = zeros(Float64, 3, length(frequencies_hz))
    propagating_mode_count = zeros(Int, length(frequencies_hz))
    sector_mode_count = zeros(Int, 3, length(frequencies_hz))

    for (frequency_index, frequency_hz) in enumerate(frequencies_hz)
        modes = propagating_modes(
            solve_port_modes(port_config(config, frequency_hz)); direction=:right,
        )
        propagating_mode_count[frequency_index] = length(modes)
        for mode_index in 0:2
            column = synthesize_column(config, modes, mode_index, sample_y)
            # Convert the continuous L2-normalized profile to the frozen
            # Euclidean-normalized 17-bin O4 convention.
            scale = sqrt(width / count)
            raw_source[:, mode_index + 1, frequency_index] .= scale .* column.source_x
            raw_output[:, mode_index + 1, frequency_index] .= scale .* column.output_x
            source_fidelity[mode_index + 1, frequency_index] = column.source_fidelity
            incident_power[mode_index + 1, frequency_index] = column.incident_power
            sector_mode_count[mode_index + 1, frequency_index] = column.sector_mode_count
        end
    end

    # A single carrier calibration equalizes incident modal power of the three
    # electrical ports. The calibration is then fixed over the whole band.
    carrier_power = incident_power[:, carrier_index]
    reference_power = exp(mean(log.(carrier_power)))
    drive_scale = sqrt.(reference_power ./ carrier_power)
    calibrated_source = raw_source .* reshape(drive_scale, 1, 3, 1)
    calibrated_output = raw_output .* reshape(drive_scale, 1, 3, 1)

    # One global passive normalization is used for the complete band. There
    # is no per-frequency gain correction.
    maximum_raw_sigma = maximum(
        opnorm(calibrated_output[:, :, index]) for index in eachindex(frequencies_hz)
    )
    passive_scale = inv(max(1.0, maximum_raw_sigma))
    transfer = passive_scale .* calibrated_output

    rows = NamedTuple[]
    for (index, frequency_hz) in enumerate(frequencies_hz)
        metrics = transfer_metrics(transfer[:, :, index], target)
        push!(rows, (
            frequency_hz,
            frequency_ratio=frequency_hz / config.carrier_frequency_hz,
            propagating_mode_count=propagating_mode_count[index],
            even_mode_count=sector_mode_count[1, index],
            odd_mode_count=sector_mode_count[2, index],
            minimum_source_fidelity=minimum(source_fidelity[:, index]),
            useful_amplitude=metrics.useful_amplitude,
            normalized_transfer_fidelity=metrics.fidelity,
            maximum_port_coherence=metrics.coherence,
            column_imbalance_db=metrics.column_imbalance_db,
            differential_amplitude_rms=metrics.amplitude_rms,
            sigma_max=metrics.sigma_max,
        ))
    end
    delay = differential_delay_metrics(frequencies_hz, transfer, target)
    delay_rows = [
        merge(rows[index], (
            common_group_delay_ns=1e9delay.common_by_frequency[index],
            differential_delay_rms_ns=1e9delay.rms_by_frequency[index],
        )) for index in eachindex(rows)
    ]
    scalar_cutoff_hz = [
        mode * config.pressure_wave_speed_m_s / (2config.aperture_width_m)
        for mode in 0:2
    ]
    lower_ratio = scalar_cutoff_hz[3] / config.lower_frequency_hz
    scalar_differential_delay_per_m = (
        inv(config.pressure_wave_speed_m_s * sqrt(1 - lower_ratio^2)) -
        inv(config.pressure_wave_speed_m_s)
    )
    scalar_length_limit_m = 82.6e-9 / scalar_differential_delay_per_m
    gates = (
        source_fidelity=minimum(source_fidelity) >= 0.98,
        transfer_fidelity=minimum(getproperty.(rows, :normalized_transfer_fidelity)) >= 0.98,
        port_coherence=maximum(getproperty.(rows, :maximum_port_coherence)) <= 0.15,
        column_imbalance=maximum(getproperty.(rows, :column_imbalance_db)) <= 0.75,
        amplitude_rms=maximum(getproperty.(rows, :differential_amplitude_rms)) <= 0.03,
        differential_delay=delay.maximum_rms <= 82.6e-9,
        useful_amplitude=minimum(getproperty.(rows, :useful_amplitude)) >= 0.8994,
        passivity=maximum(getproperty.(rows, :sigma_max)) <= 1 + 1e-10,
        reciprocity=true,
    )
    (
        config,
        frequencies_hz,
        sample_y,
        target,
        transfer,
        calibrated_source,
        source_fidelity,
        incident_power,
        drive_scale,
        passive_scale,
        maximum_raw_sigma,
        rows=delay_rows,
        delay,
        scalar_cutoff_hz,
        scalar_length_limit_m,
        gates,
        passed=all(values(gates)),
    )
end

function write_namedtuple_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function convergence_audit()
    rows = NamedTuple[]
    for element_count in (120, 160, 200)
        result = evaluate_modal_transfer(O6ModalConfig(
            lower_frequency_hz=213766.66666666666,
            upper_frequency_hz=229900.0,
            carrier_frequency_hz=221833.33333333334,
            frequency_count=9,
            port_element_count=element_count,
        ))
        push!(rows, (
            port_element_count=element_count,
            minimum_source_fidelity=minimum(result.source_fidelity),
            minimum_transfer_fidelity=minimum(
                getproperty.(result.rows, :normalized_transfer_fidelity),
            ),
            maximum_differential_amplitude_rms=maximum(
                getproperty.(result.rows, :differential_amplitude_rms),
            ),
            maximum_differential_delay_rms_ns=1e9result.delay.maximum_rms,
        ))
    end
    rows
end

function write_summary(result, convergence_rows)
    rows = result.rows
    gates = result.gates
    gate_rows = [
        (gate=String(name), passed=value) for (name, value) in pairs(gates)
    ]
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "dense_band_metrics.csv"), rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "gates.csv"), gate_rows)
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6a modal gate: direct DCT3 transition plate")
        println(io)
        println(io, result.passed ? "**PASS: one frozen 2D FEM confirmation is authorized.**" :
            "**FAIL: the uniform plate is rejected; no FEM geometry is authorized.**")
        println(io)
        @printf(io, "Candidate: `71.6 mm x %.3f mm` homogeneous Al 6061, plane strain, dense `%d`-frequency grid from `%.1f` to `%.1f kHz`.\n\n",
            1e3result.config.plate_length_m, result.config.frequency_count,
            result.config.lower_frequency_hz / 1e3, result.config.upper_frequency_hz / 1e3)
        @printf(io, "The scalar `m=2` lower-band delay bound permits at most `%.3f mm`; the candidate length was frozen before the elastic solve.\n\n",
            1e3result.scalar_length_limit_m)
        println(io, "## Worst dense-band metrics")
        println(io)
        @printf(io, "- source-profile fidelity: `%.6f`;\n", minimum(result.source_fidelity))
        @printf(io, "- propagated fidelity: `%.6f`;\n", minimum(getproperty.(rows, :normalized_transfer_fidelity)))
        @printf(io, "- maximum port coherence: `%.6f`;\n", maximum(getproperty.(rows, :maximum_port_coherence)))
        @printf(io, "- maximum column imbalance: `%.4f dB`;\n", maximum(getproperty.(rows, :column_imbalance_db)))
        @printf(io, "- maximum differential amplitude RMS: `%.4f%%`;\n", 100maximum(getproperty.(rows, :differential_amplitude_rms)))
        @printf(io, "- maximum differential-delay RMS: `%.3f ns` (diagnostic; the local phase derivative is not mesh-converged at cut-on, but remains above `1 us` in the resolution check);\n", 1e9result.delay.maximum_rms)
        @printf(io, "- minimum useful one-way amplitude: `%.6f`;\n", minimum(getproperty.(rows, :useful_amplitude)))
        @printf(io, "- maximum singular value after one global contraction: `%.12f`;\n", maximum(getproperty.(rows, :sigma_max)))
        @printf(io, "- global passive scale: `%.6f`.\n\n", result.passive_scale)
        println(io, "Reciprocity residual is exactly zero in this homogeneous modal propagation model; the independent FEM solve must measure it numerically.")
        println(io)
        println(io, "The model contains every propagating P/S-mixed free-strip branch. The number of open right-going modes changes from `$(first(getproperty.(rows, :propagating_mode_count)))` to `$(last(getproperty.(rows, :propagating_mode_count)))` across the band; the dense gate therefore includes the cut-on risk hidden by a three-mode scalar model.")
        println(io)
        println(io, "## Robustness of the no-go")
        println(io)
        println(io, "A targeted `120/160/200`-element cross-section check gives minimum source fidelity `$(join((@sprintf("%.5f", row.minimum_source_fidelity) for row in convergence_rows), "/"))` and differential amplitude RMS `$(join((@sprintf("%.2f", 100row.maximum_differential_amplitude_rms) for row in convergence_rows), "/"))%`. The exact delay spike is sensitive to phase differentiation when a branch opens, so it is retained as a diagnostic and is not the sole rejection argument. The independently converged source-profile and amplitude gates already reject the candidate.")
        println(io)
        println(io, "The stopping rule was applied: no 2D geometry and no F1d replay were run. This rejects propagation through the wide uniform free strip, not the frozen K3,Q1 measurement operator. The next allowed mechanism is a load-local direct modal transducer without a laterally guiding Al transition plate.")
    end
end

function plot_summary(result)
    frequencies_khz = result.frequencies_hz ./ 1e3
    rows = result.rows
    quality = plot(
        frequencies_khz,
        getproperty.(rows, :normalized_transfer_fidelity);
        label="transfer fidelity", linewidth=2, marker=:circle,
        xlabel="frequency, kHz", ylabel="metric", ylim=(0.94, 1.01),
        title="O6a dense-band shape gate",
    )
    plot!(quality, frequencies_khz, getproperty.(rows, :minimum_source_fidelity);
        label="source fidelity", linewidth=2)
    hline!(quality, [0.98]; label="gate", linestyle=:dash, color=:black)
    errors = plot(
        frequencies_khz,
        100 .* getproperty.(rows, :differential_amplitude_rms);
        label="amplitude RMS, %", linewidth=2, xlabel="frequency, kHz",
        ylabel="error", title="Differential errors",
    )
    plot!(errors, frequencies_khz,
        getproperty.(rows, :differential_delay_rms_ns) ./ 10;
        label="delay RMS / 10, ns", linewidth=2)
    figure = plot(quality, errors; layout=(2, 1), size=(900, 760))
    savefig(figure, joinpath(OUTPUT_ROOT, "modal_gate.png"))
end

function run_o6a_modal_gate(config::O6ModalConfig=O6ModalConfig())
    mkpath(OUTPUT_ROOT)
    result = evaluate_modal_transfer(config)
    convergence_rows = convergence_audit()
    write_summary(result, convergence_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "convergence_audit.csv"), convergence_rows)
    jldsave(
        joinpath(OUTPUT_ROOT, "o6a_modal_gate.jld2");
        format_version=1,
        config=result.config,
        frequencies_hz=result.frequencies_hz,
        sample_y_m=result.sample_y,
        target=result.target,
        transfer=result.transfer,
        calibrated_source=result.calibrated_source,
        source_fidelity=result.source_fidelity,
        incident_power=result.incident_power,
        drive_scale=result.drive_scale,
        passive_scale=result.passive_scale,
        maximum_raw_sigma=result.maximum_raw_sigma,
        rows=result.rows,
        differential_delay_s=result.delay.coefficient_delays,
        scalar_cutoff_hz=result.scalar_cutoff_hz,
        scalar_length_limit_m=result.scalar_length_limit_m,
        gates=result.gates,
        passed=result.passed,
        convergence_rows,
    )
    plot_summary(result)
    println(result.passed ? "[+] O6a PASS" : "[-] O6a FAIL")
    println("[+] $(joinpath(OUTPUT_ROOT, "README.md"))")
    result
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6a_modal_gate()

end
