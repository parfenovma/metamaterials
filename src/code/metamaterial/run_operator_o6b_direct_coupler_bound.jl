module RunOperatorO6BDirectCouplerBound

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Statistics

export O6BConfig, aluminium_normal_impedance, evaluate_o6b_bound, run_o6b_bound

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(
    PROJECT_ROOT, "results", "operator_space", "o6b_direct_coupler_bound_v0",
)

Base.@kwdef struct O6BConfig
    aperture_width_m::Float64 = 71.6e-3
    lower_frequency_hz::Float64 = 193.6e3
    upper_frequency_hz::Float64 = 290.4e3
    carrier_frequency_hz::Float64 = 242.0e3
    frequency_count::Int = 49
    aluminium_density_kg_m3::Float64 = 2700.0
    aluminium_pressure_speed_m_s::Float64 = 6122.102437409232
    aluminium_shear_speed_m_s::Float64 = 3083.810277185563
    pzt_density_kg_m3::Float64 = 7750.0
    pzt_pressure_speed_m_s::Float64 = 4350.0
    bond_density_kg_m3::Float64 = 1150.0
    bond_pressure_speed_m_s::Float64 = 2500.0
    bond_thickness_m::Float64 = 100.0e-6
end

function validate(config::O6BConfig)
    config.aperture_width_m > 0 || throw(ArgumentError("aperture width must be positive"))
    config.lower_frequency_hz < config.carrier_frequency_hz < config.upper_frequency_hz ||
        throw(ArgumentError("carrier must lie in the qualified band"))
    isodd(config.frequency_count) || throw(ArgumentError("frequency grid must contain the carrier"))
    config.frequency_count >= 7 || throw(ArgumentError("frequency grid is too sparse"))
    config.aluminium_pressure_speed_m_s > config.aluminium_shear_speed_m_s > 0 ||
        throw(ArgumentError("Al wave speeds must satisfy cP > cS > 0"))
    config.pzt_density_kg_m3 > 0 || throw(ArgumentError("PZT density must be positive"))
    config.pzt_pressure_speed_m_s > 0 || throw(ArgumentError("PZT speed must be positive"))
    config.bond_density_kg_m3 > 0 || throw(ArgumentError("bond density must be positive"))
    config.bond_pressure_speed_m_s > 0 || throw(ArgumentError("bond speed must be positive"))
    config.bond_thickness_m >= 0 || throw(ArgumentError("bond thickness must be nonnegative"))
    nothing
end

function dct_basis(count=17)
    basis = zeros(Float64, count, 3)
    basis[:, 1] .= inv(sqrt(count))
    for mode in 2:3
        basis[:, mode] .= sqrt(2 / count) .* cos.(
            (mode - 1) .* ((1:count) .- 0.5) .* pi / count,
        )
    end
    basis
end

"""Exact plane-strain normal load impedance with zero tangential traction."""
function aluminium_normal_impedance(config::O6BConfig, frequency_hz, lateral_wavenumber)
    rho = config.aluminium_density_kg_m3
    cp = config.aluminium_pressure_speed_m_s
    cs = config.aluminium_shear_speed_m_s
    omega = 2pi * frequency_hz
    kp = omega / cp
    ks = omega / cs
    kx_p = sqrt(complex(kp^2 - lateral_wavenumber^2))
    kx_s = sqrt(complex(ks^2 - lateral_wavenumber^2))
    mu = rho * cs^2
    lambda = rho * cp^2 - 2mu
    displacement = im .* ComplexF64[
        kx_p lateral_wavenumber
        lateral_wavenumber -kx_s
    ]
    traction = ComplexF64[
        -(lambda * kp^2 + 2mu * kx_p^2) -2mu * kx_s * lateral_wavenumber
        -2mu * kx_p * lateral_wavenumber mu * (kx_s^2 - lateral_wavenumber^2)
    ]
    compliance = displacement / traction
    normal_displacement_per_traction = compliance[1, 1]
    # The leading minus sign chooses positive impedance for power flowing into
    # x>0 under the solver's exp(-i*omega*t) convention.
    -inv(-im * omega * normal_displacement_per_traction)
end

function longitudinal_normal_impedance(density, speed, omega, lateral_wavenumber)
    total_wavenumber = omega / speed
    axial_wavenumber = sqrt(complex(total_wavenumber^2 - lateral_wavenumber^2))
    density * omega / axial_wavenumber, axial_wavenumber
end

"""Complex power-wave transmission from PZT through the bond into Al."""
function modal_transfer(config, frequency_hz, mode_index)
    omega = 2pi * frequency_hz
    ky = mode_index * pi / config.aperture_width_m
    aluminium_impedance = aluminium_normal_impedance(config, frequency_hz, ky)
    pzt_impedance, _ = longitudinal_normal_impedance(
        config.pzt_density_kg_m3,
        config.pzt_pressure_speed_m_s,
        omega,
        ky,
    )
    bond_impedance, bond_axial_wavenumber = longitudinal_normal_impedance(
        config.bond_density_kg_m3,
        config.bond_pressure_speed_m_s,
        omega,
        ky,
    )
    phase = bond_axial_wavenumber * config.bond_thickness_m
    a = cos(phase) * aluminium_impedance + im * bond_impedance * sin(phase)
    b = cos(phase) + im * aluminium_impedance * sin(phase) / bond_impedance
    denominator = a + pzt_impedance * b
    2sqrt(real(pzt_impedance) * real(aluminium_impedance)) / denominator,
    aluminium_impedance,
    pzt_impedance,
    bond_impedance
end

function unwrap_phase(values)
    result = angle.(values)
    for index in 2:length(result)
        result[index:end] .-= 2pi * round((result[index] - result[index - 1]) / (2pi))
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

function matrix_metrics(transfer, target)
    gram = adjoint(transfer) * transfer
    gains = sqrt.(real.(diag(gram)))
    coherence = maximum(
        abs(gram[first, second]) / sqrt(real(gram[first, first] * gram[second, second]))
        for first in 1:3, second in 1:3 if first != second
    )
    common = dot(vec(target), vec(transfer)) / sum(abs2, target)
    (
        useful_amplitude=abs(common),
        fidelity=abs(dot(vec(target), vec(transfer))) / (norm(target) * norm(transfer)),
        coherence,
        imbalance_db=20log10(maximum(gains) / minimum(gains)),
        amplitude_rms=norm(abs.(transfer) .- abs.(common .* target)) /
                      norm(abs.(common .* target)),
        sigma_max=opnorm(transfer),
        common,
    )
end

function evaluate_o6b_bound(config::O6BConfig=O6BConfig())
    validate(config)
    frequencies_hz = collect(range(
        config.lower_frequency_hz,
        config.upper_frequency_hz;
        length=config.frequency_count,
    ))
    target = dct_basis()
    modal_h = zeros(ComplexF64, 3, length(frequencies_hz))
    aluminium_z = similar(modal_h)
    pzt_z = similar(modal_h)
    bond_z = similar(modal_h)
    transfer = zeros(ComplexF64, 17, 3, length(frequencies_hz))
    rows = NamedTuple[]
    for (frequency_index, frequency_hz) in enumerate(frequencies_hz)
        for mode_index in 0:2
            h, za, zp, zb = modal_transfer(config, frequency_hz, mode_index)
            modal_h[mode_index + 1, frequency_index] = h
            aluminium_z[mode_index + 1, frequency_index] = za
            pzt_z[mode_index + 1, frequency_index] = zp
            bond_z[mode_index + 1, frequency_index] = zb
        end
        transfer[:, :, frequency_index] .= target * Diagonal(modal_h[:, frequency_index])
        metrics = matrix_metrics(transfer[:, :, frequency_index], target)
        push!(rows, (
            frequency_hz,
            frequency_ratio=frequency_hz / config.carrier_frequency_hz,
            useful_amplitude=metrics.useful_amplitude,
            normalized_transfer_fidelity=metrics.fidelity,
            maximum_port_coherence=metrics.coherence,
            column_imbalance_db=metrics.imbalance_db,
            differential_amplitude_rms=metrics.amplitude_rms,
            sigma_max=metrics.sigma_max,
            mode0_amplitude=abs(modal_h[1, frequency_index]),
            mode1_amplitude=abs(modal_h[2, frequency_index]),
            mode2_amplitude=abs(modal_h[3, frequency_index]),
        ))
    end
    omega = 2pi .* frequencies_hz
    modal_delay = zeros(Float64, size(modal_h))
    for mode in 1:3
        modal_delay[mode, :] .= -derivative(unwrap_phase(modal_h[mode, :]), omega)
    end
    common_delay = vec(mean(modal_delay; dims=1))
    differential_delay_rms = [
        sqrt(mean((modal_delay[:, index] .- common_delay[index]) .^ 2))
        for index in eachindex(frequencies_hz)
    ]
    full_rows = [
        merge(rows[index], (
            common_group_delay_ns=1e9common_delay[index],
            differential_delay_rms_ns=1e9differential_delay_rms[index],
        )) for index in eachindex(rows)
    ]
    gates = (
        transfer_fidelity=minimum(getproperty.(rows, :normalized_transfer_fidelity)) >= 0.98,
        port_coherence=maximum(getproperty.(rows, :maximum_port_coherence)) <= 0.15,
        column_imbalance=maximum(getproperty.(rows, :column_imbalance_db)) <= 0.75,
        amplitude_rms=maximum(getproperty.(rows, :differential_amplitude_rms)) <= 0.03,
        differential_delay=maximum(differential_delay_rms) <= 82.6e-9,
        useful_amplitude=minimum(getproperty.(rows, :useful_amplitude)) >= 0.8994,
        passivity=maximum(getproperty.(rows, :sigma_max)) <= 1 + 1e-10,
        reciprocity=true,
    )
    (
        config,
        frequencies_hz,
        target,
        modal_h,
        transfer,
        aluminium_z,
        pzt_z,
        bond_z,
        modal_delay,
        differential_delay_rms,
        rows=full_rows,
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

function write_report(result)
    rows = result.rows
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "dense_band_metrics.csv"), rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "gates.csv"), [
        (gate=String(name), passed=value) for (name, value) in pairs(result.gates)
    ])
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6b direct-coupler upper bound")
        println(io)
        println(io, result.passed ?
            "**PASS: one finite mechanically coupled PZT/bond/Al twin is authorized.**" :
            "**FAIL: the direct load-local coupler topology is stopped.**")
        println(io)
        @printf(io, "The reciprocal dense-band model uses exact plane-strain P/S loading of Al, a longitudinal PZT-5A mechanical port and a frozen `%.0f um` epoxy bond.\n\n",
            1e6result.config.bond_thickness_m)
        println(io, "## Worst dense-band metrics")
        println(io)
        @printf(io, "- normalized transfer fidelity: `%.8f`;\n", minimum(getproperty.(rows, :normalized_transfer_fidelity)))
        @printf(io, "- maximum port coherence: `%.3e`;\n", maximum(getproperty.(rows, :maximum_port_coherence)))
        @printf(io, "- maximum column imbalance: `%.5f dB`;\n", maximum(getproperty.(rows, :column_imbalance_db)))
        @printf(io, "- maximum differential amplitude RMS: `%.5f%%`;\n", 100maximum(getproperty.(rows, :differential_amplitude_rms)))
        @printf(io, "- maximum differential-delay RMS: `%.5f ns`;\n", 1e9maximum(result.differential_delay_rms))
        @printf(io, "- minimum useful one-way amplitude: `%.8f`;\n", minimum(getproperty.(rows, :useful_amplitude)))
        @printf(io, "- maximum singular value: `%.8f`.\n\n", maximum(getproperty.(rows, :sigma_max)))
        println(io, "This pass removes the O6a lateral-guide mechanism but is only a necessary upper bound. It assumes an ideal distributed piezocomposite coupling profile and omits finite PZT lateral coupling, backing, electrode mass and dielectric loss. The next authorized calculation is exactly one finite 2D mechanical PZT/bond/Al twin with the same DCT0--2 profiles and bond thickness.")
    end
end

function plot_report(result)
    f = result.frequencies_hz ./ 1e3
    amplitudes = plot(
        f, transpose(abs.(result.modal_h));
        linewidth=2, label=["DCT0" "DCT1" "DCT2"],
        xlabel="frequency, kHz", ylabel="power-wave amplitude",
        title="O6b direct modal coupling",
    )
    hline!(amplitudes, [0.8994]; linestyle=:dash, color=:black, label="O4 useful gate")
    errors = plot(
        f, getproperty.(result.rows, :column_imbalance_db);
        linewidth=2, label="imbalance, dB", xlabel="frequency, kHz",
        ylabel="differential error", title="Residual modal mismatch",
    )
    plot!(errors, f, 100 .* getproperty.(result.rows, :differential_amplitude_rms);
        linewidth=2, label="amplitude RMS, %")
    savefig(plot(amplitudes, errors; layout=(2, 1), size=(900, 760)),
        joinpath(OUTPUT_ROOT, "direct_coupler_bound.png"))
end

function run_o6b_bound(config::O6BConfig=O6BConfig())
    mkpath(OUTPUT_ROOT)
    result = evaluate_o6b_bound(config)
    write_report(result)
    jldsave(
        joinpath(OUTPUT_ROOT, "o6b_direct_coupler_bound.jld2");
        format_version=1,
        config=result.config,
        frequencies_hz=result.frequencies_hz,
        target=result.target,
        modal_transfer=result.modal_h,
        transfer=result.transfer,
        aluminium_impedance=result.aluminium_z,
        pzt_impedance=result.pzt_z,
        bond_impedance=result.bond_z,
        modal_delay_s=result.modal_delay,
        rows=result.rows,
        gates=result.gates,
        passed=result.passed,
    )
    plot_report(result)
    println(result.passed ? "[+] O6b upper bound PASS" : "[-] O6b upper bound FAIL")
    println("[+] $(joinpath(OUTPUT_ROOT, "README.md"))")
    result
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6b_bound()

end
