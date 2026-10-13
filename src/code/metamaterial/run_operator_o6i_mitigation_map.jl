module RunOperatorO6IMitigationMap

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf

include(joinpath(@__DIR__, "run_operator_o6h_angular_contact.jl"))
using .RunOperatorO6HAngularContact

export solid_layer_scattering,
       run_o6i_mitigation_map

const A = RunOperatorO6HAngularContact
const C = A.C
const M = A.M
const OUTPUT_ROOT = joinpath(
    M.PROJECT_ROOT, "results", "operator_space", "o6i_mitigation_map_v0",
)
const IMPEDANCE_RATIOS = (0.10, 0.25, 0.50, 0.75, 1.00)
const PRESSURE_SPEED_RATIOS = (0.25, 0.50, 0.75, 1.00)
const SHEAR_SPEED_RATIOS = (0.10, 0.25, 0.50)
const THICKNESS_WAVELENGTHS = (0.0025, 0.005, 0.01, 0.02)
const MINIMUM_DENSITY = 800.0
const MAXIMUM_DENSITY = 8000.0
const ROBUST_POWER_GATE = 0.85
const ROBUST_S_CHANGE_GATE = 0.075
const INTRINSIC_Q = 35.0
const USEFUL_AMPLITUDE_GATE = 0.8994

function solid_layer_scattering(
    frequency_hz,
    incidence_angle_deg,
    thickness_m,
    layer_material;
    input_material=A.PROBE,
    output_material=A.SAMPLE,
)
    omega = 2pi * frequency_hz
    tangential_wavenumber = omega / input_material.pressure_speed *
        sin(deg2rad(incidence_angle_deg))
    incident = A.solid_mode(input_material, omega, tangential_wavenumber, :P, 1)
    reflected_p = A.solid_mode(input_material, omega, tangential_wavenumber, :P, -1)
    reflected_s = A.solid_mode(input_material, omega, tangential_wavenumber, :S, -1)
    layer_p_forward = A.solid_mode(layer_material, omega, tangential_wavenumber, :P, 1)
    layer_s_forward = A.solid_mode(layer_material, omega, tangential_wavenumber, :S, 1)
    layer_p_backward = A.solid_mode(layer_material, omega, tangential_wavenumber, :P, -1)
    layer_s_backward = A.solid_mode(layer_material, omega, tangential_wavenumber, :S, -1)
    transmitted_p = A.solid_mode(output_material, omega, tangential_wavenumber, :P, 1)
    transmitted_s = A.solid_mode(output_material, omega, tangential_wavenumber, :S, 1)
    layer_modes = (
        layer_p_forward, layer_s_forward, layer_p_backward, layer_s_backward,
    )
    phases = [exp(im * mode.normal_wavenumber * thickness_m)
        for mode in layer_modes]
    matrix = zeros(ComplexF64, 8, 8)
    right = zeros(ComplexF64, 8)
    # Unknowns: reflected P/S, layer P+/S+/P-/S-, transmitted P/S.
    for (row, field) in enumerate((:vx, :vy, :sxx, :sxy))
        reflected_values = field == :vx ?
            (reflected_p.velocity[1], reflected_s.velocity[1]) :
            field == :vy ?
            (reflected_p.velocity[2], reflected_s.velocity[2]) :
            field == :sxx ?
            (reflected_p.stress_xx, reflected_s.stress_xx) :
            (reflected_p.stress_xy, reflected_s.stress_xy)
        layer_values = field == :vx ? getindex.(getproperty.(layer_modes, :velocity), 1) :
            field == :vy ? getindex.(getproperty.(layer_modes, :velocity), 2) :
            field == :sxx ? getproperty.(layer_modes, :stress_xx) :
            getproperty.(layer_modes, :stress_xy)
        incident_value = field == :vx ? incident.velocity[1] :
            field == :vy ? incident.velocity[2] :
            field == :sxx ? incident.stress_xx : incident.stress_xy
        matrix[row, 1:2] .= reflected_values
        matrix[row, 3:6] .= -collect(layer_values)
        right[row] = -incident_value
    end
    for (local_row, field) in enumerate((:vx, :vy, :sxx, :sxy))
        row = local_row + 4
        layer_values = field == :vx ? getindex.(getproperty.(layer_modes, :velocity), 1) :
            field == :vy ? getindex.(getproperty.(layer_modes, :velocity), 2) :
            field == :sxx ? getproperty.(layer_modes, :stress_xx) :
            getproperty.(layer_modes, :stress_xy)
        transmitted_values = field == :vx ?
            (transmitted_p.velocity[1], transmitted_s.velocity[1]) :
            field == :vy ?
            (transmitted_p.velocity[2], transmitted_s.velocity[2]) :
            field == :sxx ?
            (transmitted_p.stress_xx, transmitted_s.stress_xx) :
            (transmitted_p.stress_xy, transmitted_s.stress_xy)
        matrix[row, 3:6] .= collect(layer_values) .* phases
        matrix[row, 7:8] .= -collect(transmitted_values)
    end
    coefficients = matrix \ right
    reflected_power = -(
        abs2(coefficients[1]) * reflected_p.normal_flux +
        abs2(coefficients[2]) * reflected_s.normal_flux
    ) / incident.normal_flux
    transmitted_p_power = abs2(coefficients[7]) *
        transmitted_p.normal_flux / incident.normal_flux
    transmitted_s_power = abs2(coefficients[8]) *
        transmitted_s.normal_flux / incident.normal_flux
    transmitted_power = transmitted_p_power + transmitted_s_power
    (
        reflected_p=coefficients[1],
        reflected_s=coefficients[2],
        transmitted_p=coefficients[7],
        transmitted_s=coefficients[8],
        reflected_power,
        transmitted_p_power,
        transmitted_s_power,
        transmitted_power,
        transmitted_s_fraction=transmitted_power > eps() ?
            transmitted_s_power / transmitted_power : 0.0,
        energy_balance=reflected_power + transmitted_power,
        condition_number=cond(matrix),
    )
end

function candidate_summary(
    config,
    impedance_ratio,
    pressure_speed_ratio,
    shear_speed_ratio,
    thickness_wavelengths,
)
    reference_impedance = sqrt(
        C.PROBE_IMPEDANCE_RAYL * C.AL_IMPEDANCE_RAYL,
    )
    layer_impedance = impedance_ratio * reference_impedance
    pressure_speed = pressure_speed_ratio * C.AL_PRESSURE_SPEED_M_S
    shear_speed = shear_speed_ratio * pressure_speed
    density = layer_impedance / pressure_speed
    layer = A.solid_material(density, pressure_speed, shear_speed)
    thickness_m = thickness_wavelengths * pressure_speed / C.CARRIER_FREQUENCY_HZ
    rows = NamedTuple[]
    for ratio in config.frequency_ratios, angle in A.ANGLES_DEG
        frequency = C.CARRIER_FREQUENCY_HZ * ratio
        result = solid_layer_scattering(
            frequency, angle, thickness_m, layer,
        )
        direct = A.direct_interface_scattering(frequency, angle)
        push!(rows, (
            frequency_ratio=ratio,
            incidence_angle_deg=angle,
            layer_to_direct_power=result.transmitted_power /
                max(direct.transmitted_power, eps(Float64)),
            s_fraction_change=abs(
                result.transmitted_s_fraction - direct.transmitted_s_fraction,
            ),
            energy_balance_error=abs(result.energy_balance - 1),
            condition_number=result.condition_number,
        ))
    end
    qualification = filter(row -> row.incidence_angle_deg <= 45.0, rows)
    minimum_power = minimum(getproperty.(qualification, :layer_to_direct_power))
    maximum_s_change = maximum(getproperty.(qualification, :s_fraction_change))
    maximum_energy_error = maximum(getproperty.(rows, :energy_balance_error))
    (
        impedance_ratio,
        pressure_speed_ratio,
        shear_speed_ratio,
        thickness_wavelengths,
        layer_impedance_mrayl=layer_impedance / 1e6,
        pressure_speed_m_s=pressure_speed,
        shear_speed_m_s=shear_speed,
        density_kg_m3=density,
        thickness_um=1e6thickness_m,
        minimum_layer_to_direct_power_0_45=minimum_power,
        maximum_s_fraction_change_0_45=maximum_s_change,
        minimum_layer_to_direct_power_0_55=
            minimum(getproperty.(rows, :layer_to_direct_power)),
        maximum_s_fraction_change_0_55=
            maximum(getproperty.(rows, :s_fraction_change)),
        maximum_energy_balance_error=maximum_energy_error,
        maximum_condition_number=maximum(getproperty.(rows, :condition_number)),
        density_valid=MINIMUM_DENSITY <= density <= MAXIMUM_DENSITY,
        power_gate=minimum_power >= A.POWER_RATIO_GATE,
        mode_conversion_gate=maximum_s_change <= A.S_FRACTION_CHANGE_GATE,
        joint_pass=minimum_power >= A.POWER_RATIO_GATE &&
            maximum_s_change <= A.S_FRACTION_CHANGE_GATE &&
            maximum_energy_error <= A.ENERGY_BALANCE_TOLERANCE,
        robust_interior=minimum_power >= ROBUST_POWER_GATE &&
            maximum_s_change <= ROBUST_S_CHANGE_GATE &&
            impedance_ratio in (0.25, 0.50, 0.75) &&
            MINIMUM_DENSITY <= density <= MAXIMUM_DENSITY,
    )
end

function choose_handoff(rows)
    candidates = filter(row -> row.robust_interior, rows)
    isempty(candidates) && return nothing
    density_center = sqrt(MINIMUM_DENSITY * MAXIMUM_DENSITY)
    sort!(candidates; by=row -> (
        -row.thickness_um,
        abs(log(row.density_kg_m3 / density_center)),
        row.impedance_ratio,
        row.pressure_speed_ratio,
        row.shear_speed_ratio,
        row.thickness_wavelengths,
    ))
    first(candidates)
end

function write_rows(path, rows)
    M.write_rows(path, rows)
end

function run_o6i_mitigation_map()
    mkpath(OUTPUT_ROOT)
    o6h = JLD2.load(joinpath(
        M.PROJECT_ROOT, "results", "operator_space", "o6h_bvd_loading_v0",
        "o6h_bvd_loading.jld2",
    ))
    required_q = Float64(o6h["quality_ceiling"])
    damping_efficiency = required_q / INTRINSIC_Q
    damping_amplitude_ceiling = sqrt(damping_efficiency)
    damping_passed = damping_amplitude_ceiling >= USEFUL_AMPLITUDE_GATE
    damping_row = (
        intrinsic_q=INTRINSIC_Q,
        required_effective_q=required_q,
        resistance_multiplier=INTRINSIC_Q / required_q,
        maximum_useful_power_fraction=damping_efficiency,
        maximum_useful_amplitude=damping_amplitude_ceiling,
        useful_amplitude_gate=USEFUL_AMPLITUDE_GATE,
        passed=damping_passed,
    )
    println(damping_passed ?
        "[+] O6i-A resistive damping PASS" :
        "[-] O6i-A resistive damping no-go")

    config = M.SurrogateConfig()
    rows = NamedTuple[]
    for impedance_ratio in IMPEDANCE_RATIOS,
        pressure_speed_ratio in PRESSURE_SPEED_RATIOS,
        shear_speed_ratio in SHEAR_SPEED_RATIOS,
        thickness_wavelengths in THICKNESS_WAVELENGTHS
        candidate = candidate_summary(
            config,
            impedance_ratio,
            pressure_speed_ratio,
            shear_speed_ratio,
            thickness_wavelengths,
        )
        candidate.density_valid && push!(rows, candidate)
    end
    maximum(getproperty.(rows, :maximum_energy_balance_error)) <=
        A.ENERGY_BALANCE_TOLERANCE || error("solid-layer energy balance failed")
    passing = filter(row -> row.joint_pass, rows)
    handoff = choose_handoff(rows)
    isnothing(handoff) && error("no robust interior solid-contact handoff exists")
    write_rows(joinpath(OUTPUT_ROOT, "q8_resistive_damping.csv"), [damping_row])
    write_rows(joinpath(OUTPUT_ROOT, "solid_contact_map.csv"), rows)
    write_rows(joinpath(OUTPUT_ROOT, "solid_contact_passing.csv"), passing)
    write_rows(joinpath(OUTPUT_ROOT, "solid_contact_handoff.csv"), [handoff])
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6i mitigation-class map")
        println(io)
        println(io, damping_passed ?
            "**q8 resistive broadening PASS.**" :
            "**q8 resistive broadening STOP; lossless/multiresonant matching remains required.**")
        @printf(io,
            " Reducing intrinsic Q `%.1f` to `%.5f` by added series resistance requires `%.2fx` total resistance and leaves at most `%.3f%%` useful power (`%.4f` amplitude), below the frozen `%.4f` amplitude gate.\n\n",
            INTRINSIC_Q, required_q, INTRINSIC_Q / required_q,
            100damping_efficiency, damping_amplitude_ceiling,
            USEFUL_AMPLITUDE_GATE,
        )
        println(io, "**A shear-transmitting solid-contact class exists.**")
        @printf(io,
            " `%d/%d` density-valid grid candidates pass the O6h-C angular gates; `%d` satisfy the stricter interior rule. The frozen handoff is:\n\n",
            length(passing), length(rows), count(row -> row.robust_interior, rows),
        )
        println(io, "| ZP ratio | cP ratio | cS/cP | d/lambdaP | rho | thickness | min power | max Delta S |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|---:|")
        @printf(io,
            "| `%.2f` | `%.2f` | `%.2f` | `%.4f` | `%.1f kg/m3` | `%.2f um` | `%.4f` | `%.4f` |\n\n",
            handoff.impedance_ratio,
            handoff.pressure_speed_ratio,
            handoff.shear_speed_ratio,
            handoff.thickness_wavelengths,
            handoff.density_kg_m3,
            handoff.thickness_um,
            handoff.minimum_layer_to_direct_power_0_45,
            handoff.maximum_s_fraction_change_0_45,
        )
        println(io, "This is an elastic-property handoff, not a named commercial film. O6j must apply attachment-property/thickness variations to this frozen interior point and replay P12 detector Pfa before any prototype recommendation.")
    end
    contact_plot = scatter(
        getproperty.(rows, :minimum_layer_to_direct_power_0_45),
        getproperty.(rows, :maximum_s_fraction_change_0_45);
        marker_z=getproperty.(rows, :thickness_um),
        color=:viridis,
        xlabel="minimum layer/direct power, 0--45 deg",
        ylabel="maximum Delta S fraction",
        title="O6i-B solid-contact property map",
        label="density-valid candidates",
    )
    vline!(contact_plot, [A.POWER_RATIO_GATE];
        linestyle=:dash, color=:black, label="power gate")
    hline!(contact_plot, [A.S_FRACTION_CHANGE_GATE];
        linestyle=:dash, color=:red, label="mode gate")
    scatter!(contact_plot,
        [handoff.minimum_layer_to_direct_power_0_45],
        [handoff.maximum_s_fraction_change_0_45];
        marker=:star5, markersize=10, color=:gold, label="frozen handoff")
    savefig(contact_plot, joinpath(OUTPUT_ROOT, "solid_contact_map.png"))
    println("[+] O6i-B $(length(passing))/$(length(rows)) solid contacts pass")
    println("[+] wrote O6i to $OUTPUT_ROOT")
    (; damping_row, rows, passing, handoff)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6i_mitigation_map()

end
