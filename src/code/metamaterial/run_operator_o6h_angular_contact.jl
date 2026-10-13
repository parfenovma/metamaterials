module RunOperatorO6HAngularContact

ENV["GKSwstype"] = "100"

using LinearAlgebra
using Plots
using Printf

include(joinpath(@__DIR__, "run_operator_o6h_layered_contact.jl"))
using .RunOperatorO6HLayeredContact

export elastic_layer_scattering,
       direct_interface_scattering,
       run_o6h_angular_contact

const C = RunOperatorO6HLayeredContact
const M = C.M
const OUTPUT_ROOT = joinpath(
    M.PROJECT_ROOT, "results", "operator_space", "o6h_angular_contact_v0",
)
const ANGLES_DEG = (0.0, 15.0, 30.0, 45.0, 55.0)
const POWER_RATIO_GATE = 0.80
const S_FRACTION_CHANGE_GATE = 0.10
const ENERGY_BALANCE_TOLERANCE = 1e-8

function solid_material(density, pressure_speed, shear_speed)
    shear_modulus = density * shear_speed^2
    lame_first = density * pressure_speed^2 - 2shear_modulus
    (
        density=Float64(density),
        pressure_speed=Float64(pressure_speed),
        shear_speed=Float64(shear_speed),
        lame_first,
        shear_modulus,
    )
end

const PROBE = solid_material(
    C.PROBE_IMPEDANCE_RAYL / C.AL_PRESSURE_SPEED_M_S,
    C.AL_PRESSURE_SPEED_M_S,
    3083.8103,
)
const SAMPLE = solid_material(
    C.AL_DENSITY_KG_M3,
    C.AL_PRESSURE_SPEED_M_S,
    3083.8103,
)

function vertical_wavenumber(total_wavenumber, tangential_wavenumber, direction)
    root = sqrt(complex(total_wavenumber^2 - tangential_wavenumber^2))
    imag(root) < 0 && (root = -root)
    direction * root
end

function solid_mode(material, omega, tangential_wavenumber, kind, direction)
    speed = kind == :P ? material.pressure_speed : material.shear_speed
    total_wavenumber = omega / speed
    normal_wavenumber = vertical_wavenumber(
        total_wavenumber, tangential_wavenumber, direction,
    )
    wavevector = ComplexF64[normal_wavenumber, tangential_wavenumber]
    velocity = kind == :P ?
        wavevector / total_wavenumber :
        ComplexF64[-tangential_wavenumber, normal_wavenumber] / total_wavenumber
    divergence = dot(wavevector, velocity)
    stress_xx = -(
        material.lame_first * divergence +
        2material.shear_modulus * normal_wavenumber * velocity[1]
    ) / omega
    stress_xy = -material.shear_modulus * (
        normal_wavenumber * velocity[2] + tangential_wavenumber * velocity[1]
    ) / omega
    normal_flux = -0.5real(
        stress_xx * conj(velocity[1]) + stress_xy * conj(velocity[2]),
    )
    (; velocity, stress_xx, stress_xy, normal_flux, normal_wavenumber)
end

function elastic_layer_scattering(
    frequency_hz,
    incidence_angle_deg,
    thickness_m,
    fluid_impedance,
    fluid_speed;
    input_material=PROBE,
    output_material=SAMPLE,
)
    omega = 2pi * frequency_hz
    incidence_angle = deg2rad(incidence_angle_deg)
    tangential_wavenumber = omega / input_material.pressure_speed *
        sin(incidence_angle)
    incident = solid_mode(input_material, omega, tangential_wavenumber, :P, 1)
    reflected_p = solid_mode(input_material, omega, tangential_wavenumber, :P, -1)
    reflected_s = solid_mode(input_material, omega, tangential_wavenumber, :S, -1)
    transmitted_p = solid_mode(output_material, omega, tangential_wavenumber, :P, 1)
    transmitted_s = solid_mode(output_material, omega, tangential_wavenumber, :S, 1)
    fluid_density = fluid_impedance / fluid_speed
    fluid_total_wavenumber = omega / fluid_speed
    fluid_normal_wavenumber = vertical_wavenumber(
        fluid_total_wavenumber, tangential_wavenumber, 1,
    )
    fluid_admittance = fluid_normal_wavenumber / (omega * fluid_density)
    forward_phase = exp(im * fluid_normal_wavenumber * thickness_m)
    backward_phase = exp(-im * fluid_normal_wavenumber * thickness_m)
    matrix = zeros(ComplexF64, 6, 6)
    right = zeros(ComplexF64, 6)
    # Unknowns: reflected P, reflected SV, fluid +, fluid -, transmitted P, SV.
    matrix[1, :] .= [
        reflected_p.velocity[1], reflected_s.velocity[1],
        -fluid_admittance, fluid_admittance, 0, 0,
    ]
    right[1] = -incident.velocity[1]
    matrix[2, :] .= [
        reflected_p.stress_xx, reflected_s.stress_xx, 1, 1, 0, 0,
    ]
    right[2] = -incident.stress_xx
    matrix[3, :] .= [
        reflected_p.stress_xy, reflected_s.stress_xy, 0, 0, 0, 0,
    ]
    right[3] = -incident.stress_xy
    matrix[4, :] .= [
        0, 0,
        fluid_admittance * forward_phase,
        -fluid_admittance * backward_phase,
        -transmitted_p.velocity[1], -transmitted_s.velocity[1],
    ]
    matrix[5, :] .= [
        0, 0, forward_phase, backward_phase,
        transmitted_p.stress_xx, transmitted_s.stress_xx,
    ]
    matrix[6, :] .= [
        0, 0, 0, 0, transmitted_p.stress_xy, transmitted_s.stress_xy,
    ]
    coefficients = matrix \ right
    reflected_power = -(
        abs2(coefficients[1]) * reflected_p.normal_flux +
        abs2(coefficients[2]) * reflected_s.normal_flux
    ) / incident.normal_flux
    transmitted_p_power = abs2(coefficients[5]) *
        transmitted_p.normal_flux / incident.normal_flux
    transmitted_s_power = abs2(coefficients[6]) *
        transmitted_s.normal_flux / incident.normal_flux
    transmitted_power = transmitted_p_power + transmitted_s_power
    energy_balance = reflected_power + transmitted_power
    (
        reflected_p=coefficients[1],
        reflected_s=coefficients[2],
        transmitted_p=coefficients[5],
        transmitted_s=coefficients[6],
        reflected_power,
        transmitted_p_power,
        transmitted_s_power,
        transmitted_power,
        transmitted_s_fraction=transmitted_power > eps() ?
            transmitted_s_power / transmitted_power : 0.0,
        energy_balance,
        condition_number=cond(matrix),
    )
end

function direct_interface_scattering(
    frequency_hz,
    incidence_angle_deg;
    input_material=PROBE,
    output_material=SAMPLE,
)
    omega = 2pi * frequency_hz
    tangential_wavenumber = omega / input_material.pressure_speed *
        sin(deg2rad(incidence_angle_deg))
    incident = solid_mode(input_material, omega, tangential_wavenumber, :P, 1)
    reflected_p = solid_mode(input_material, omega, tangential_wavenumber, :P, -1)
    reflected_s = solid_mode(input_material, omega, tangential_wavenumber, :S, -1)
    transmitted_p = solid_mode(output_material, omega, tangential_wavenumber, :P, 1)
    transmitted_s = solid_mode(output_material, omega, tangential_wavenumber, :S, 1)
    matrix = [
        reflected_p.velocity[1] reflected_s.velocity[1] -transmitted_p.velocity[1] -transmitted_s.velocity[1]
        reflected_p.velocity[2] reflected_s.velocity[2] -transmitted_p.velocity[2] -transmitted_s.velocity[2]
        reflected_p.stress_xx reflected_s.stress_xx -transmitted_p.stress_xx -transmitted_s.stress_xx
        reflected_p.stress_xy reflected_s.stress_xy -transmitted_p.stress_xy -transmitted_s.stress_xy
    ]
    right = -ComplexF64[
        incident.velocity[1], incident.velocity[2],
        incident.stress_xx, incident.stress_xy,
    ]
    coefficients = matrix \ right
    reflected_power = -(
        abs2(coefficients[1]) * reflected_p.normal_flux +
        abs2(coefficients[2]) * reflected_s.normal_flux
    ) / incident.normal_flux
    transmitted_p_power = abs2(coefficients[3]) *
        transmitted_p.normal_flux / incident.normal_flux
    transmitted_s_power = abs2(coefficients[4]) *
        transmitted_s.normal_flux / incident.normal_flux
    transmitted_power = transmitted_p_power + transmitted_s_power
    (
        reflected_p=coefficients[1],
        reflected_s=coefficients[2],
        transmitted_p=coefficients[3],
        transmitted_s=coefficients[4],
        reflected_power,
        transmitted_p_power,
        transmitted_s_power,
        transmitted_power,
        transmitted_s_fraction=transmitted_power > eps() ?
            transmitted_s_power / transmitted_power : 0.0,
        energy_balance=reflected_power + transmitted_power,
    )
end

function write_rows(path, rows)
    M.write_rows(path, rows)
end

function run_o6h_angular_contact()
    mkpath(OUTPUT_ROOT)
    config = M.SurrogateConfig()
    scenarios = (
        (
            id=:minus_sigma,
            thickness=C.NOMINAL_THICKNESS_M - C.THICKNESS_RMS_M,
            impedance=C.COUPLANT_IMPEDANCE_RAYL * exp(-C.IMPEDANCE_LOG_RMS),
            speed=C.COUPLANT_SPEED_M_S * exp(-C.SPEED_LOG_RMS),
        ),
        (
            id=:nominal,
            thickness=C.NOMINAL_THICKNESS_M,
            impedance=C.COUPLANT_IMPEDANCE_RAYL,
            speed=C.COUPLANT_SPEED_M_S,
        ),
        (
            id=:plus_sigma,
            thickness=C.NOMINAL_THICKNESS_M + C.THICKNESS_RMS_M,
            impedance=C.COUPLANT_IMPEDANCE_RAYL * exp(C.IMPEDANCE_LOG_RMS),
            speed=C.COUPLANT_SPEED_M_S * exp(C.SPEED_LOG_RMS),
        ),
    )
    rows = NamedTuple[]
    maximum_energy_error = 0.0
    maximum_normal_error = 0.0
    for scenario in scenarios, ratio in config.frequency_ratios, angle in ANGLES_DEG
        frequency = C.CARRIER_FREQUENCY_HZ * ratio
        layer = elastic_layer_scattering(
            frequency,
            angle,
            scenario.thickness,
            scenario.impedance,
            scenario.speed,
        )
        direct = direct_interface_scattering(frequency, angle)
        energy_error = abs(layer.energy_balance - 1)
        direct_energy_error = abs(direct.energy_balance - 1)
        maximum_energy_error = max(
            maximum_energy_error, energy_error, direct_energy_error,
        )
        if scenario.id == :nominal && angle == 0.0
            expected = abs2(C.layer_power_transmission(
                frequency,
                scenario.thickness,
                scenario.impedance,
                scenario.speed,
            ))
            maximum_normal_error = max(
                maximum_normal_error, abs(layer.transmitted_power - expected),
            )
        end
        push!(rows, (
            scenario=scenario.id,
            frequency_ratio=ratio,
            incidence_angle_deg=angle,
            thickness_um=1e6scenario.thickness,
            fluid_impedance_mrayl=scenario.impedance / 1e6,
            fluid_speed_m_s=scenario.speed,
            layer_reflected_power=layer.reflected_power,
            layer_transmitted_p_power=layer.transmitted_p_power,
            layer_transmitted_s_power=layer.transmitted_s_power,
            layer_transmitted_power=layer.transmitted_power,
            layer_s_fraction=layer.transmitted_s_fraction,
            direct_transmitted_power=direct.transmitted_power,
            direct_s_fraction=direct.transmitted_s_fraction,
            layer_to_direct_power=layer.transmitted_power /
                max(direct.transmitted_power, eps(Float64)),
            s_fraction_change=abs(
                layer.transmitted_s_fraction - direct.transmitted_s_fraction,
            ),
            energy_balance_error=energy_error,
            condition_number=layer.condition_number,
        ))
    end
    maximum_energy_error <= ENERGY_BALANCE_TOLERANCE ||
        error("2D layer energy balance failed: $maximum_energy_error")
    maximum_normal_error <= 1e-8 ||
        error("2D normal incidence does not reproduce 1D: $maximum_normal_error")
    summaries = NamedTuple[]
    for scenario in scenarios
        qualification = filter(row ->
            row.scenario == scenario.id && row.incidence_angle_deg <= 45.0, rows)
        diagnostic = filter(row -> row.scenario == scenario.id, rows)
        minimum_power = minimum(getproperty.(qualification, :layer_to_direct_power))
        maximum_s_change = maximum(getproperty.(qualification, :s_fraction_change))
        gates = (
            transmitted_power=minimum_power >= POWER_RATIO_GATE,
            mode_conversion=maximum_s_change <= S_FRACTION_CHANGE_GATE,
        )
        push!(summaries, (
            scenario=scenario.id,
            minimum_layer_to_direct_power_0_45=minimum_power,
            maximum_s_fraction_change_0_45=maximum_s_change,
            minimum_layer_to_direct_power_0_55=
                minimum(getproperty.(diagnostic, :layer_to_direct_power)),
            maximum_s_fraction_change_0_55=
                maximum(getproperty.(diagnostic, :s_fraction_change)),
            maximum_condition_number=maximum(getproperty.(diagnostic, :condition_number)),
            power_gate=gates.transmitted_power,
            mode_conversion_gate=gates.mode_conversion,
            joint_pass=all(values(gates)),
        ))
    end
    branch_passed = all(getproperty.(summaries, :joint_pass))
    write_rows(joinpath(OUTPUT_ROOT, "angular_contact_rows.csv"), rows)
    write_rows(joinpath(OUTPUT_ROOT, "angular_contact_summary.csv"), summaries)
    write_rows(joinpath(OUTPUT_ROOT, "positive_controls.csv"), [(
        maximum_energy_balance_error=maximum_energy_error,
        maximum_normal_incidence_power_error=maximum_normal_error,
        energy_balance_tolerance=ENERGY_BALANCE_TOLERANCE,
        passed=maximum_energy_error <= ENERGY_BALANCE_TOLERANCE &&
            maximum_normal_error <= 1e-8,
    )])
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6h-C exact 2D angular contact spot check")
        println(io)
        println(io, branch_passed ?
            "**PASS for nominal and paired one-sigma contact anchors through 45 deg.**" :
            "**FAIL: the normal-incidence P12 contact qualification does not survive 2D angular P/S physics.**")
        println(io)
        @printf(io,
            "The exact solid--fluid--solid boundary solve conserves lossless power to `%.3e` and reproduces the 1D normal-incidence transmitted power to `%.3e`.\n\n",
            maximum_energy_error, maximum_normal_error,
        )
        println(io, "| Scenario | Min layer/direct power, 0--45 deg | Max Delta S fraction | 55-deg min power | Pass |")
        println(io, "|---|---:|---:|---:|---|")
        for row in summaries
            @printf(io, "| `%s` | `%.4f` | `%.4f` | `%.4f` | `%s` |\n",
                String(row.scenario),
                row.minimum_layer_to_direct_power_0_45,
                row.maximum_s_fraction_change_0_45,
                row.minimum_layer_to_direct_power_0_55,
                row.joint_pass,
            )
        end
        println(io)
        println(io, "This is an exact 2D plane-wave spot check over the operator band and ROI-relevant angles. It still omits lateral thickness texture, squeeze-out edges and a finite contact footprint; a failure is therefore a physical blocker, while a pass is only permission for a finite-footprint prototype model.")
    end
    power_plot = plot(
        xlabel="incidence angle, deg",
        ylabel="minimum layer/direct transmitted power",
        title="O6h-C angular contact transmission",
        legend=:bottomleft,
    )
    for (scenario, color) in zip(scenarios, (:steelblue, :black, :darkorange))
        values = [minimum(row.layer_to_direct_power for row in rows
            if row.scenario == scenario.id && row.incidence_angle_deg == angle)
            for angle in ANGLES_DEG]
        plot!(power_plot, collect(ANGLES_DEG), values;
            marker=:circle, linewidth=2, color, label=String(scenario.id))
    end
    hline!(power_plot, [POWER_RATIO_GATE];
        linestyle=:dash, color=:red, label="0.80 gate")
    savefig(power_plot, joinpath(OUTPUT_ROOT, "angular_contact_power.png"))
    println(branch_passed ? "[+] O6h-C angular contact PASS" : "[-] O6h-C angular contact FAIL")
    println("[+] wrote O6h-C to $OUTPUT_ROOT")
    (; branch_passed, summaries)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6h_angular_contact()

end
