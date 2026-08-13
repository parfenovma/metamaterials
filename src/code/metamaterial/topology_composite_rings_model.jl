module TopologyCompositeRingsModel

using FFTW
using LinearAlgebra
using Statistics

include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .ImpulseRiskAnalysis

export CompositeProperty,
       CompositeRingsConfig,
       LossCase,
       BondCase,
       RingZone,
       measured_dp190_w_epoxy,
       interpolate_composite,
       design_ring_zones,
       required_delay_s,
       layer_transfer,
       zone_transfer,
       aperture_transfer,
       evaluate_case,
       pulse_spectrum

"""Measured or interpolated normal-incidence longitudinal property set."""
Base.@kwdef struct CompositeProperty
    label::String
    tungsten_volume_fraction::Float64
    density_kg_m3::Float64
    pressure_wave_speed_m_s::Float64
    impedance_mrayl::Float64
end

"""
Autoclave DP190/W anchors from Powell et al. (2026), DOI
10.1177/07316844261436540, Table 5. Pure DP190 is the same paper's reference.
Only interpolation inside the measured 9.92--20.22 vol% interval is used for
the focusing core; the 31.52 vol% point is used as the common face-sheet
candidate.
"""
function measured_dp190_w_epoxy()
    CompositeProperty[
        CompositeProperty(
            label="DP190 neat, atmosphere-cured reference",
            tungsten_volume_fraction=0.0,
            density_kg_m3=1110.0,
            pressure_wave_speed_m_s=2382.0,
            impedance_mrayl=2.61,
        ),
        CompositeProperty(
            label="DP190 + 9.92 vol% W + 1 vol% glass beads, 4 bar autoclave",
            tungsten_volume_fraction=0.0992,
            density_kg_m3=3000.0,
            pressure_wave_speed_m_s=1982.0,
            impedance_mrayl=5.99,
        ),
        CompositeProperty(
            label="DP190 + 20.22 vol% W + 1 vol% glass beads, 4 bar autoclave",
            tungsten_volume_fraction=0.2022,
            density_kg_m3=4950.0,
            pressure_wave_speed_m_s=1698.0,
            impedance_mrayl=8.41,
        ),
        CompositeProperty(
            label="DP190 + 31.52 vol% W + 1 vol% glass beads, 4 bar autoclave",
            tungsten_volume_fraction=0.3152,
            density_kg_m3=5940.0,
            pressure_wave_speed_m_s=1688.0,
            impedance_mrayl=10.02,
        ),
    ]
end

function validate(property::CompositeProperty)
    0 <= property.tungsten_volume_fraction < 1 ||
        throw(ArgumentError("invalid tungsten volume fraction"))
    property.density_kg_m3 > 0 || throw(ArgumentError("density must be positive"))
    property.pressure_wave_speed_m_s > 0 || throw(ArgumentError("speed must be positive"))
    property.impedance_mrayl > 0 || throw(ArgumentError("impedance must be positive"))
    predicted = property.density_kg_m3 * property.pressure_wave_speed_m_s / 1.0e6
    isapprox(predicted, property.impedance_mrayl; rtol=0.015) ||
        throw(ArgumentError("density, speed and impedance are inconsistent"))
    property
end

"""Piecewise-linear interpolation; extrapolation is deliberately forbidden."""
function interpolate_composite(volume_fraction::Real; anchors=measured_dp190_w_epoxy())
    fraction = Float64(volume_fraction)
    sorted = sort(collect(anchors); by=x -> x.tungsten_volume_fraction)
    first(sorted).tungsten_volume_fraction <= fraction <= last(sorted).tungsten_volume_fraction ||
        throw(ArgumentError("requested composition is outside measured anchors"))
    right = searchsortedfirst(getproperty.(sorted, :tungsten_volume_fraction), fraction)
    right == 1 && return validate(sorted[1])
    right > length(sorted) && return validate(sorted[end])
    left = right - 1
    a, b = sorted[left], sorted[right]
    weight = (fraction - a.tungsten_volume_fraction) /
             (b.tungsten_volume_fraction - a.tungsten_volume_fraction)
    result = CompositeProperty(
        label="linear interpolation $(round(100fraction, digits=3)) vol% W",
        tungsten_volume_fraction=fraction,
        density_kg_m3=(1 - weight) * a.density_kg_m3 + weight * b.density_kg_m3,
        pressure_wave_speed_m_s=(1 - weight) * a.pressure_wave_speed_m_s +
                                weight * b.pressure_wave_speed_m_s,
        impedance_mrayl=(1 - weight) * a.impedance_mrayl + weight * b.impedance_mrayl,
    )
    # Interpolating all three measured observables independently introduces a
    # small closure error. Preserve measured density and speed, and use their
    # product in the wave model.
    CompositeProperty(
        label=result.label,
        tungsten_volume_fraction=result.tungsten_volume_fraction,
        density_kg_m3=result.density_kg_m3,
        pressure_wave_speed_m_s=result.pressure_wave_speed_m_s,
        impedance_mrayl=result.density_kg_m3 * result.pressure_wave_speed_m_s / 1.0e6,
    )
end

Base.@kwdef struct CompositeRingsConfig
    center_frequency_hz::Float64 = 242.0e3
    pulse_cycles::Float64 = 5.0
    aperture_diameter_mm::Float64 = 71.6
    focus_depth_mm::Float64 = 60.0
    aluminium_speed_m_s::Float64 = 6122.102437409232
    aluminium_density_kg_m3::Float64 = 2700.0
    core_thickness_mm::Float64 = 18.0
    ring_count::Int = 7
    outer_tungsten_volume_fraction::Float64 = 0.0992
    face_tungsten_volume_fraction::Float64 = 0.2022
    face_thickness_mm::Float64 = 1.65
    aperture_samples::Int = 401
end

Base.@kwdef struct LossCase
    label::String
    core_db_per_cm_at_f0::Float64
    frequency_power::Float64 = 1.0
    face_multiplier::Float64 = 1.35
end

Base.@kwdef struct BondCase
    label::String
    thickness_mm::Float64
    property::CompositeProperty = first(measured_dp190_w_epoxy())
    loss_db_per_cm_at_f0::Float64 = 0.5
end

struct RingZone
    index::Int
    inner_radius_mm::Float64
    outer_radius_mm::Float64
    representative_radius_mm::Float64
    target_delay_s::Float64
    target_speed_m_s::Float64
    property::CompositeProperty
end

function validate(config::CompositeRingsConfig)
    config.center_frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    config.aperture_diameter_mm > 0 || throw(ArgumentError("aperture must be positive"))
    config.focus_depth_mm > 0 || throw(ArgumentError("focus must be positive"))
    config.core_thickness_mm > 0 || throw(ArgumentError("core thickness must be positive"))
    config.ring_count >= 3 || throw(ArgumentError("at least three rings are required"))
    config.aperture_samples >= 101 || throw(ArgumentError("aperture sampling is too coarse"))
    config.core_thickness_mm + 2config.face_thickness_mm <= 55.0 ||
        throw(ArgumentError("device exceeds the benchmark depth"))
    config
end

"""Extra central delay relative to radius `reference_radius_mm`."""
function required_delay_s(radius_mm::Real, reference_radius_mm::Real, config::CompositeRingsConfig)
    path = hypot(config.focus_depth_mm, Float64(radius_mm))
    reference_path = hypot(config.focus_depth_mm, Float64(reference_radius_mm))
    (reference_path - path) * 1.0e-3 / config.aluminium_speed_m_s
end

"""Area centroid of a 3-D annulus, used to choose one recipe per physical ring."""
function annulus_centroid_radius_mm(inner::Real, outer::Real)
    a, b = Float64(inner), Float64(outer)
    iszero(a) && return 2b / 3
    (2 / 3) * (b^3 - a^3) / (b^2 - a^2)
end

function fraction_for_speed(speed_m_s::Real; anchors=measured_dp190_w_epoxy())
    # The manufacturable core is intentionally constrained to the measured,
    # non-porous 9.92--20.22 vol% autoclave interval.
    low_fraction = 0.0992
    high_fraction = 0.2022
    fast = interpolate_composite(low_fraction; anchors)
    slow = interpolate_composite(high_fraction; anchors)
    slow.pressure_wave_speed_m_s <= speed_m_s <= fast.pressure_wave_speed_m_s ||
        throw(ArgumentError("target speed lies outside the measured autoclave interval"))
    low_fraction + (fast.pressure_wave_speed_m_s - Float64(speed_m_s)) /
                   (fast.pressure_wave_speed_m_s - slow.pressure_wave_speed_m_s) *
                   (high_fraction - low_fraction)
end

function design_ring_zones(config::CompositeRingsConfig=CompositeRingsConfig())
    validate(config)
    radius = config.aperture_diameter_mm / 2
    boundaries = collect(range(0.0, radius; length=config.ring_count + 1))
    representatives = [
        annulus_centroid_radius_mm(boundaries[index], boundaries[index + 1])
        for index in 1:config.ring_count
    ]
    reference_radius = last(representatives)
    outer = interpolate_composite(config.outer_tungsten_volume_fraction)
    base_transit_s = config.core_thickness_mm * 1.0e-3 / outer.pressure_wave_speed_m_s
    [begin
        delay = required_delay_s(representatives[index], reference_radius, config)
        target_speed = config.core_thickness_mm * 1.0e-3 / (base_transit_s + delay)
        fraction = fraction_for_speed(target_speed)
        RingZone(
            index,
            boundaries[index],
            boundaries[index + 1],
            representatives[index],
            delay,
            target_speed,
            interpolate_composite(fraction),
        )
    end for index in 1:config.ring_count]
end

function attenuation_np_per_m(
    db_per_cm_at_f0::Real,
    frequency_hz::Real,
    center_frequency_hz::Real,
    frequency_power::Real,
)
    frequency_hz <= 0 && return 0.0
    db_per_m = 100Float64(db_per_cm_at_f0) *
               (Float64(frequency_hz) / Float64(center_frequency_hz))^Float64(frequency_power)
    db_per_m * log(10) / 20
end


function layer_matrix(
    property::CompositeProperty,
    thickness_mm::Real,
    frequency_hz::Real;
    loss_db_per_cm_at_f0::Real=0.0,
    center_frequency_hz::Real=242.0e3,
    frequency_power::Real=1.0,
)
    validate(property)
    thickness_m = Float64(thickness_mm) * 1.0e-3
    omega = 2pi * Float64(frequency_hz)
    alpha = attenuation_np_per_m(
        loss_db_per_cm_at_f0,
        frequency_hz,
        center_frequency_hz,
        frequency_power,
    )
    k = omega / property.pressure_wave_speed_m_s - im * alpha
    z = property.density_kg_m3 * property.pressure_wave_speed_m_s
    phase = k * thickness_m
    ComplexF64[
        cos(phase) im * z * sin(phase)
        im * sin(phase) / z cos(phase)
    ]
end

"""Pressure transmission/reflection for a layer stack between real half-spaces."""
function layer_transfer(matrix, input_impedance_rayl::Real, output_impedance_rayl::Real)
    z0, zl = Float64(input_impedance_rayl), Float64(output_impedance_rayl)
    a, b, c, d = matrix[1, 1], matrix[1, 2], matrix[2, 1], matrix[2, 2]
    denominator = a + b / zl + z0 * c + z0 * d / zl
    transmission = 2 / denominator
    reflection = (a + b / zl) * transmission - 1
    power_transmission = abs2(transmission) * z0 / zl
    (; transmission, reflection, power_transmission)
end

function zone_transfer(
    zone_property::CompositeProperty,
    frequency_hz::Real,
    config::CompositeRingsConfig,
    loss::LossCase,
    bond::BondCase;
    face_property=interpolate_composite(config.face_tungsten_volume_fraction),
)
    z_al = config.aluminium_density_kg_m3 * config.aluminium_speed_m_s
    matrix = Matrix{ComplexF64}(I, 2, 2)
    if bond.thickness_mm > 0
        matrix *= layer_matrix(
            bond.property,
            bond.thickness_mm,
            frequency_hz;
            loss_db_per_cm_at_f0=bond.loss_db_per_cm_at_f0,
            center_frequency_hz=config.center_frequency_hz,
            frequency_power=loss.frequency_power,
        )
    end
    matrix *= layer_matrix(
        face_property,
        config.face_thickness_mm,
        frequency_hz;
        loss_db_per_cm_at_f0=loss.face_multiplier * loss.core_db_per_cm_at_f0,
        center_frequency_hz=config.center_frequency_hz,
        frequency_power=loss.frequency_power,
    )
    matrix *= layer_matrix(
        zone_property,
        config.core_thickness_mm,
        frequency_hz;
        loss_db_per_cm_at_f0=loss.core_db_per_cm_at_f0 *
            (1 + 0.5 * (zone_property.tungsten_volume_fraction - 0.0992) / (0.2022 - 0.0992)),
        center_frequency_hz=config.center_frequency_hz,
        frequency_power=loss.frequency_power,
    )
    matrix *= layer_matrix(
        face_property,
        config.face_thickness_mm,
        frequency_hz;
        loss_db_per_cm_at_f0=loss.face_multiplier * loss.core_db_per_cm_at_f0,
        center_frequency_hz=config.center_frequency_hz,
        frequency_power=loss.frequency_power,
    )
    if bond.thickness_mm > 0
        matrix *= layer_matrix(
            bond.property,
            bond.thickness_mm,
            frequency_hz;
            loss_db_per_cm_at_f0=bond.loss_db_per_cm_at_f0,
            center_frequency_hz=config.center_frequency_hz,
            frequency_power=loss.frequency_power,
        )
    end
    layer_transfer(matrix, z_al, z_al)
end

function zone_index(radius_mm::Real, zones)
    value = abs(Float64(radius_mm))
    for zone in zones
        value <= zone.outer_radius_mm + 1.0e-12 && return zone.index
    end
    last(zones).index
end

function straight_property(zones)
    # Preserve total tungsten inventory in the actual axisymmetric part.
    total_area = last(zones).outer_radius_mm^2
    fraction = sum(zones) do zone
        weight = (zone.outer_radius_mm^2 - zone.inner_radius_mm^2) / total_area
        weight * zone.property.tungsten_volume_fraction
    end
    interpolate_composite(fraction)
end

"""
2-D plane-strain-style scalar Huygens surrogate. It is deliberately named a
surrogate: it omits elastic mode conversion and finite-defect scattering.
"""
function aperture_transfer(
    frequency_hz::Real,
    config::CompositeRingsConfig,
    zones,
    loss::LossCase,
    bond::BondCase;
    role::Symbol=:candidate,
    observation_x_mm::Real=config.focus_depth_mm,
    observation_y_mm::Real=0.0,
)
    role in (:candidate, :straight, :bare) || throw(ArgumentError("unknown role"))
    frequency_hz <= 0 && return (field=0.0 + 0.0im, aperture_reflection=0.0 + 0.0im,
                                 weighted_power_transmission=0.0)
    radius = config.aperture_diameter_mm / 2
    y_mm = collect(range(-radius, radius; length=config.aperture_samples))
    dy_m = (y_mm[2] - y_mm[1]) * 1.0e-3
    wavelength_m = config.aluminium_speed_m_s / Float64(frequency_hz)
    omega = 2pi * Float64(frequency_hz)
    straight = straight_property(zones)
    field = 0.0 + 0.0im
    reflection = 0.0 + 0.0im
    power = 0.0
    for y in y_mm
        local_transfer = if role == :bare
            (transmission=1.0 + 0im, reflection=0.0 + 0im, power_transmission=1.0)
        elseif role == :straight
            zone_transfer(straight, frequency_hz, config, loss, bond)
        else
            zone_transfer(zones[zone_index(y, zones)].property, frequency_hz, config, loss, bond)
        end
        distance_m = hypot(Float64(observation_x_mm), Float64(observation_y_mm) - y) * 1.0e-3
        green = dy_m / sqrt(wavelength_m * distance_m) *
                cis(-omega * distance_m / config.aluminium_speed_m_s)
        field += local_transfer.transmission * green
        reflection += local_transfer.reflection / length(y_mm)
        power += local_transfer.power_transmission / length(y_mm)
    end
    (; field, aperture_reflection=reflection, weighted_power_transmission=power)
end

function active_indices(pulse; lower_hz=80.0e3, upper_hz=420.0e3, relative_floor=1.0e-5)
    findall(
        (pulse.frequency_hz .>= lower_hz) .&
        (pulse.frequency_hz .<= upper_hz) .&
        (pulse.relative_amplitude .>= relative_floor),
    )
end

function waveform(pulse, transfer)
    irfft(pulse.spectrum .* transfer, length(pulse.signal))
end

function group_delay_s(frequency_hz, response)
    raw_phases = angle.(response)
    phases = Float64[first(raw_phases)]
    for phase in Iterators.drop(raw_phases, 1)
        push!(phases, last(phases) + mod(phase - last(phases) + pi, 2pi) - pi)
    end
    omega = 2pi .* frequency_hz
    result = similar(phases)
    result[1] = -(phases[2] - phases[1]) / (omega[2] - omega[1])
    for index in 2:(length(result) - 1)
        result[index] = -(phases[index + 1] - phases[index - 1]) /
                        (omega[index + 1] - omega[index - 1])
    end
    result[end] = -(phases[end] - phases[end - 1]) / (omega[end] - omega[end - 1])
    result
end

function normalized_peak(signal)
    maximum(abs, signal)
end

function evaluate_case(
    config::CompositeRingsConfig,
    loss::LossCase,
    bond::BondCase;
    fft_length::Int=65536,
)
    zones = design_ring_zones(config)
    pulse = pulse_spectrum(PulseConfig(
        center_frequency_hz=config.center_frequency_hz,
        cycles=config.pulse_cycles,
        samples_per_period=80,
        fft_length=fft_length,
    ))
    candidate = zeros(ComplexF64, length(pulse.spectrum))
    straight = zeros(ComplexF64, length(candidate))
    bare = zeros(ComplexF64, length(candidate))
    candidate_reflection = zeros(ComplexF64, length(candidate))
    straight_reflection = zeros(ComplexF64, length(candidate))
    candidate_power = fill(NaN, length(candidate))
    straight_power = fill(NaN, length(candidate))
    active = active_indices(pulse)
    for index in active
        frequency = pulse.frequency_hz[index]
        c = aperture_transfer(frequency, config, zones, loss, bond; role=:candidate)
        s = aperture_transfer(frequency, config, zones, loss, bond; role=:straight)
        b = aperture_transfer(frequency, config, zones, loss, bond; role=:bare)
        candidate[index] = c.field
        straight[index] = s.field
        bare[index] = b.field
        candidate_reflection[index] = c.aperture_reflection
        straight_reflection[index] = s.aperture_reflection
        candidate_power[index] = c.weighted_power_transmission
        straight_power[index] = s.weighted_power_transmission
    end
    candidate_wave = waveform(pulse, candidate)
    straight_wave = waveform(pulse, straight)
    bare_wave = waveform(pulse, bare)
    ideal = bare_wave # conservative pulse-fidelity reference
    candidate_vs_bare = ImpulseRiskAnalysis.pulse_metrics(
        candidate_wave, bare_wave, ideal, pulse.time_s[2] - pulse.time_s[1],
    )
    candidate_vs_straight = ImpulseRiskAnalysis.pulse_metrics(
        candidate_wave, straight_wave, ideal, pulse.time_s[2] - pulse.time_s[1],
    )
    weights = abs2.(pulse.spectrum)
    finite = active
    weighted_candidate_power = sum(weights[finite] .* candidate_power[finite]) / sum(weights[finite])
    weighted_straight_power = sum(weights[finite] .* straight_power[finite]) / sum(weights[finite])

    # Strictly labelled point-scatterer proxy. H is normalized by the bare
    # carrier response before squaring, so it cannot be confused with a FEM
    # void result or an absolute CNR/POD prediction.
    carrier_index = argmin(abs.(pulse.frequency_hz .- config.center_frequency_hz))
    normalization = abs(bare[carrier_index])
    scatter_weight = zeros(Float64, length(candidate))
    scatter_weight[active] .= (pulse.frequency_hz[active] ./ config.center_frequency_hz).^2
    candidate_echo_transfer = scatter_weight .* (candidate ./ normalization).^2
    straight_echo_transfer = scatter_weight .* (straight ./ normalization).^2
    bare_echo_transfer = scatter_weight .* (bare ./ normalization).^2
    candidate_echo = waveform(pulse, candidate_echo_transfer)
    straight_echo = waveform(pulse, straight_echo_transfer)
    bare_echo = waveform(pulse, bare_echo_transfer)
    candidate_ringing = waveform(pulse, candidate_reflection)
    straight_ringing = waveform(pulse, straight_reflection)
    echo_peak_gain_bare_proxy = normalized_peak(candidate_echo) / normalized_peak(bare_echo)
    echo_energy_gain_bare_proxy = sum(abs2, candidate_echo) / sum(abs2, bare_echo)

    band = findall((pulse.relative_amplitude .>= 0.1) .& (pulse.frequency_hz .> 0))
    relative_zone_response = [
        ComplexF64[
            zone_transfer(zone.property, pulse.frequency_hz[index], config, loss, bond).transmission
            for index in band
        ] for zone in zones
    ]
    delays = [group_delay_s(pulse.frequency_hz[band], response) for response in relative_zone_response]
    delay_spread_s = [maximum(map(values -> values[index], delays)) -
                      minimum(map(values -> values[index], delays)) for index in eachindex(band)]

    (; config, loss, bond, zones, pulse, active,
       candidate_transfer=candidate,
       straight_transfer=straight,
       bare_transfer=bare,
       candidate_reflection,
       straight_reflection,
       candidate_wave,
       straight_wave,
       bare_wave,
       candidate_ringing,
       straight_ringing,
       candidate_echo,
       straight_echo,
       bare_echo,
       candidate_vs_bare,
       candidate_vs_straight,
       weighted_candidate_power,
       weighted_straight_power,
       echo_peak_gain_bare_proxy,
       echo_energy_gain_bare_proxy,
       band_frequency_hz=pulse.frequency_hz[band],
       zone_group_delay_s=delays,
       delay_spread_s)
end

end
