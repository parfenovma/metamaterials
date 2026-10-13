module TopologyFresnelRingsModel

using FFTW
using LinearAlgebra

include(joinpath(@__DIR__, "ndt_benchmark.jl"))
using .NDTBenchmark

export AcousticMaterial,
       FresnelRingConfig,
       RingLevel,
       FresnelDesign,
       aluminium,
       peek,
       polyurethane,
       tungsten_epoxy,
       design_fresnel_rings,
       required_delay_s,
       phase_span_rad,
       halfwave_zone_radii_mm,
       layer_transmission,
       level_transfer,
       fresnel_pulse_metrics,
       radial_transfer,
       scalar_aperture_transfer,
       lateral_profile,
       axial_peak_mm,
       weighted_mean

Base.@kwdef struct AcousticMaterial
    name::Symbol
    density_kg_m3::Float64
    speed_m_s::Float64
    attenuation_db_cm_at_fc::Float64 = 0.0
end

impedance(material::AcousticMaterial) = material.density_kg_m3 * material.speed_m_s

aluminium() = AcousticMaterial(
    name=:al6061,
    density_kg_m3=2700.0,
    speed_m_s=6122.102437409232,
)

peek(; attenuation_db_cm_at_fc=0.8) = AcousticMaterial(
    name=:peek,
    density_kg_m3=1300.0,
    speed_m_s=2530.0,
    attenuation_db_cm_at_fc=attenuation_db_cm_at_fc,
)

polyurethane(; attenuation_db_cm_at_fc=2.0) = AcousticMaterial(
    name=:cast_polyurethane,
    density_kg_m3=1100.0,
    speed_m_s=1500.0,
    attenuation_db_cm_at_fc=attenuation_db_cm_at_fc,
)

tungsten_epoxy(; attenuation_db_cm_at_fc=0.5) = AcousticMaterial(
    name=:tungsten_epoxy_match,
    density_kg_m3=3000.0,
    speed_m_s=2500.0,
    attenuation_db_cm_at_fc=attenuation_db_cm_at_fc,
)

Base.@kwdef struct FresnelRingConfig
    stand::NDTStandConfig = NDTStandConfig()
    level_count::Int = 3
    core_depth_mm::Float64 = 8.0
    matching_layer_depth_mm::Float64 = 2500.0 / (4 * 242.0e3) * 1.0e3
    minimum_feature_mm::Float64 = 0.8
end

struct RingLevel
    id::Int
    inner_radius_mm::Float64
    outer_radius_mm::Float64
    target_delay_s::Float64
    wrapped_delay_s::Float64
    polyurethane_depth_mm::Float64
    peek_depth_mm::Float64
end

struct FresnelDesign
    config::FresnelRingConfig
    levels::Vector{RingLevel}
    maximum_delay_s::Float64
    carrier_period_s::Float64
    carrier_phase_span_rad::Float64
    full_wrap_count::Int
    first_halfwave_radius_mm::Float64
    matching_layer_depth_mm::Float64
    total_device_depth_mm::Float64
end

aperture_radius_mm(config::FresnelRingConfig) = config.stand.aperture_height_mm / 2

function path_excess_mm(radius_mm::Real, focus_mm::Real)
    hypot(Float64(focus_mm), Float64(radius_mm)) - Float64(focus_mm)
end

"""Delay needed at radius r relative to the aperture edge for equal arrival at focus."""
function required_delay_s(radius_mm::Real, config::FresnelRingConfig=FresnelRingConfig())
    radius = Float64(radius_mm)
    edge = aperture_radius_mm(config)
    0 <= radius <= edge || throw(ArgumentError("radius must lie in aperture"))
    focus = config.stand.nominal_focus_depth_mm
    path_delta_mm = path_excess_mm(edge, focus) - path_excess_mm(radius, focus)
    path_delta_mm * 1.0e-3 / config.stand.pressure_wave_speed_m_s
end

phase_span_rad(config::FresnelRingConfig=FresnelRingConfig()) =
    2pi * config.stand.frequency_hz * required_delay_s(0.0, config)

function halfwave_zone_radii_mm(config::FresnelRingConfig=FresnelRingConfig())
    wavelength = pressure_wavelength_mm(config.stand)
    focus = config.stand.nominal_focus_depth_mm
    edge = aperture_radius_mm(config)
    radii = Float64[]
    order = 1
    while true
        path = order * wavelength / 2
        radius = sqrt((focus + path)^2 - focus^2)
        radius > edge && break
        push!(radii, radius)
        order += 1
    end
    radii
end

function radius_for_delay_s(delay_s::Real, config::FresnelRingConfig)
    edge = aperture_radius_mm(config)
    focus = config.stand.nominal_focus_depth_mm
    edge_excess = path_excess_mm(edge, focus)
    excess = edge_excess - Float64(delay_s) * config.stand.pressure_wave_speed_m_s * 1.0e3
    excess = clamp(excess, 0.0, edge_excess)
    sqrt((focus + excess)^2 - focus^2)
end

function design_fresnel_rings(config::FresnelRingConfig=FresnelRingConfig())
    2 <= config.level_count <= 8 || throw(ArgumentError("level_count must be in 2:8"))
    config.core_depth_mm > 0 || throw(ArgumentError("core depth must be positive"))
    max_delay = required_delay_s(0.0, config)
    period = inv(config.stand.frequency_hz)
    phase_span = 2pi * max_delay / period
    delays = collect(range(max_delay, 0.0; length=config.level_count))
    boundaries = Float64[0.0]
    for index in 1:(config.level_count - 1)
        midpoint = (delays[index] + delays[index + 1]) / 2
        push!(boundaries, radius_for_delay_s(midpoint, config))
    end
    push!(boundaries, aperture_radius_mm(config))

    slow = polyurethane(attenuation_db_cm_at_fc=0.0)
    fast = peek(attenuation_db_cm_at_fc=0.0)
    delay_per_m = inv(slow.speed_m_s) - inv(fast.speed_m_s)
    levels = RingLevel[]
    for index in eachindex(delays)
        wrapped = mod(delays[index], period)
        slow_depth_mm = wrapped / delay_per_m * 1.0e3
        slow_depth_mm <= config.core_depth_mm + 1.0e-10 ||
            throw(ArgumentError("core is too shallow for requested delay"))
        push!(levels, RingLevel(
            index,
            boundaries[index],
            boundaries[index + 1],
            delays[index],
            wrapped,
            slow_depth_mm,
            config.core_depth_mm - slow_depth_mm,
        ))
    end

    first_path = pressure_wavelength_mm(config.stand) / 2
    focus = config.stand.nominal_focus_depth_mm
    first_radius = sqrt((focus + first_path)^2 - focus^2)
    matching_depth = config.matching_layer_depth_mm
    total_depth = config.core_depth_mm + 2matching_depth
    total_depth <= config.stand.maximum_device_depth_mm ||
        throw(ArgumentError("device exceeds stand depth budget"))
    FresnelDesign(
        config,
        levels,
        max_delay,
        period,
        phase_span,
        floor(Int, phase_span / 2pi),
        first_radius,
        matching_depth,
        total_depth,
    )
end

function material_with_loss(material::AcousticMaterial, loss_scale::Real)
    AcousticMaterial(
        name=material.name,
        density_kg_m3=material.density_kg_m3,
        speed_m_s=material.speed_m_s,
        attenuation_db_cm_at_fc=material.attenuation_db_cm_at_fc * Float64(loss_scale),
    )
end

function layer_matrix(material::AcousticMaterial, depth_mm::Real, frequency_hz::Real, fc_hz::Real)
    depth_m = Float64(depth_mm) * 1.0e-3
    depth_m == 0 && return Matrix{ComplexF64}(I, 2, 2)
    frequency = Float64(frequency_hz)
    frequency == 0 && return Matrix{ComplexF64}(I, 2, 2)
    alpha_db_m = material.attenuation_db_cm_at_fc * 100 * frequency / fc_hz
    alpha_np_m = alpha_db_m * log(10) / 20
    wavenumber = 2pi * frequency / material.speed_m_s - im * alpha_np_m
    phase = wavenumber * depth_m
    z = impedance(material)
    ComplexF64[
        cos(phase) im * z * sin(phase)
        im * sin(phase) / z cos(phase)
    ]
end

"""Pressure transmission/reflection for a normal-incidence layer stack in Al."""
function layer_transmission(
    layers,
    frequency_hz::Real;
    incident::AcousticMaterial=aluminium(),
    transmitted::AcousticMaterial=incident,
    fc_hz::Real=242.0e3,
)
    matrix = Matrix{ComplexF64}(I, 2, 2)
    for (material, depth_mm) in layers
        matrix *= layer_matrix(material, depth_mm, frequency_hz, fc_hz)
    end
    zl = impedance(incident)
    zr = impedance(transmitted)
    denominator = matrix[1, 1] + matrix[1, 2] / zr + zl * matrix[2, 1] +
                  zl * matrix[2, 2] / zr
    transmission = 2 / denominator
    reflection = (matrix[1, 1] + matrix[1, 2] / zr - zl * matrix[2, 1] -
                  zl * matrix[2, 2] / zr) / denominator
    power_transmission = abs2(transmission) * zr / zl
    (transmission=transmission, reflection=reflection,
     power_transmission=power_transmission, power_reflection=abs2(reflection))
end

function level_layers(
    polyurethane_depth_mm::Real,
    design::FresnelDesign;
    loss_scale::Real=1.0,
)
    fast = material_with_loss(peek(), loss_scale)
    slow = material_with_loss(polyurethane(), loss_scale)
    matching = material_with_loss(tungsten_epoxy(), loss_scale)
    slow_depth = Float64(polyurethane_depth_mm)
    fast_depth = design.config.core_depth_mm - slow_depth
    [
        (matching, design.matching_layer_depth_mm),
        (fast, fast_depth),
        (slow, slow_depth),
        (matching, design.matching_layer_depth_mm),
    ]
end

function level_transfer(
    level::RingLevel,
    design::FresnelDesign,
    frequency_hz::Real;
    loss_scale::Real=1.0,
)
    layer_transmission(
        level_layers(level.polyurethane_depth_mm, design; loss_scale),
        frequency_hz;
        fc_hz=design.config.stand.frequency_hz,
    )
end

function continuous_transfer(
    radius_mm::Real,
    design::FresnelDesign,
    frequency_hz::Real;
    loss_scale::Real=1.0,
)
    slow = polyurethane(attenuation_db_cm_at_fc=0.0)
    fast = peek(attenuation_db_cm_at_fc=0.0)
    delay_per_m = inv(slow.speed_m_s) - inv(fast.speed_m_s)
    slow_depth_mm = required_delay_s(radius_mm, design.config) / delay_per_m * 1.0e3
    layer_transmission(
        level_layers(slow_depth_mm, design; loss_scale),
        frequency_hz;
        fc_hz=design.config.stand.frequency_hz,
    )
end

function level_at_radius(radius_mm::Real, design::FresnelDesign)
    radius = Float64(radius_mm)
    for level in design.levels
        if level.inner_radius_mm <= radius <= level.outer_radius_mm + 1.0e-10
            return level
        end
    end
    throw(ArgumentError("radius outside aperture"))
end

function weighted_mean(values::AbstractVector, weights::AbstractVector)
    length(values) == length(weights) || throw(DimensionMismatch())
    sum(values .* weights) / sum(weights)
end

function hann_burst(time_s::AbstractVector, frequency_hz::Real, cycles::Real)
    duration = Float64(cycles) / Float64(frequency_hz)
    [
        0 <= time <= duration ?
            0.5 * (1 - cos(2pi * time / duration)) * sin(2pi * frequency_hz * time) : 0.0
        for time in time_s
    ]
end

function energy_width_samples(signal::AbstractVector, fraction::Real=0.9)
    energy = abs2.(signal)
    total = sum(energy)
    total == 0 && return 0
    cumulative = cumsum(energy) / total
    lower = searchsortedfirst(cumulative, (1 - fraction) / 2)
    upper = searchsortedfirst(cumulative, 1 - (1 - fraction) / 2)
    max(1, upper - lower + 1)
end

function fresnel_pulse_metrics(
    transfer::Function;
    frequency_hz::Real=242.0e3,
    cycles::Real=5.0,
    dt_s::Real=0.1e-6,
    sample_count::Int=8192,
)
    time = collect(0:(sample_count - 1)) .* Float64(dt_s)
    input = hann_burst(time, frequency_hz, cycles)
    spectrum = rfft(input)
    frequencies = collect(0:length(spectrum)-1) ./ (sample_count * Float64(dt_s))
    response = ComplexF64[transfer(f) for f in frequencies]
    output = irfft(spectrum .* response, sample_count)
    broadening = energy_width_samples(output) / energy_width_samples(input)
    correlation = maximum(abs.(ifft(conj.(fft(input)) .* fft(output)))) /
                  (norm(input) * norm(output))

    peak_index = argmax(abs.(output))
    main_halfwidth = ceil(Int, cycles / frequency_hz / dt_s / 2)
    post_start = min(sample_count + 1, peak_index + main_halfwidth + ceil(Int, 1 / frequency_hz / dt_s))
    postcursor = post_start <= sample_count ? maximum(abs, @view output[post_start:end]) : 0.0
    main_peak = maximum(abs, output)
    postcursor_ratio = main_peak > 0 ? postcursor / main_peak : Inf
    (
        broadening_ratio=broadening,
        correlation=correlation,
        postcursor_ratio=postcursor_ratio,
        input=input,
        output=output,
        time_s=time,
    )
end

function radial_transfer(
    radius_mm::Real,
    design::FresnelDesign,
    frequency_hz::Real,
    topology::Symbol;
    loss_scale::Real=1.0,
)
    if topology == :bare
        return 1.0 + 0im
    elseif topology == :ideal_unwrapped_delay
        return exp(-2pi * im * Float64(frequency_hz) *
                   required_delay_s(radius_mm, design.config))
    elseif topology == :fresnel_quantized
        return level_transfer(level_at_radius(radius_mm, design), design, frequency_hz; loss_scale).transmission
    elseif topology == :unwrapped_continuous
        return continuous_transfer(radius_mm, design, frequency_hz; loss_scale).transmission
    elseif topology == :straight_footprint
        edge = aperture_radius_mm(design.config)
        # Area-weighted mean filler depth preserves cassette material inventory.
        slow_depth = sum(
            level.polyurethane_depth_mm *
            (level.outer_radius_mm^2 - level.inner_radius_mm^2) / edge^2
            for level in design.levels
        )
        return layer_transmission(
            level_layers(slow_depth, design; loss_scale),
            frequency_hz;
            fc_hz=design.config.stand.frequency_hz,
        ).transmission
    end
    throw(ArgumentError("unknown topology: $topology"))
end

"""Axisymmetric scalar bulk-P aperture proxy, never labelled elastic FEM."""
function scalar_aperture_transfer(
    design::FresnelDesign,
    frequency_hz::Real,
    topology::Symbol;
    axial_mm::Real=design.config.stand.nominal_focus_depth_mm,
    lateral_mm::Real=0.0,
    radial_samples::Int=160,
    azimuth_samples::Int=(lateral_mm == 0 ? 1 : 96),
    loss_scale::Real=1.0,
)
    radius_max = aperture_radius_mm(design.config)
    dr = radius_max / radial_samples
    dtheta = 2pi / azimuth_samples
    z = Float64(axial_mm)
    y = Float64(lateral_mm)
    k = 2pi * Float64(frequency_hz) / design.config.stand.pressure_wave_speed_m_s * 1.0e-3
    integral = 0.0 + 0im
    for radial_index in 1:radial_samples
        radius = (radial_index - 0.5) * dr
        local_transfer = radial_transfer(radius, design, frequency_hz, topology; loss_scale)
        angular = 0.0 + 0im
        for theta_index in 1:azimuth_samples
            theta = (theta_index - 0.5) * dtheta
            distance = sqrt(z^2 + radius^2 + y^2 - 2radius * y * cos(theta))
            # Rayleigh--Sommerfeld normal derivative kernel (constant prefactors
            # cancel in equal-input comparisons at a fixed frequency).
            angular += exp(-im * k * distance) *
                       (1 + im * k * distance) * z / distance^3
        end
        angular /= azimuth_samples
        integral += 2radius * dr / radius_max^2 * local_transfer * angular
    end
    integral
end

function lateral_profile(
    design::FresnelDesign,
    topology::Symbol;
    frequency_hz::Real=design.config.stand.frequency_hz,
    maximum_mm::Real=60.0,
    step_mm::Real=0.5,
    loss_scale::Real=1.0,
)
    positions = collect(0.0:Float64(step_mm):Float64(maximum_mm))
    amplitude = abs.([
        scalar_aperture_transfer(
            design,
            frequency_hz,
            topology;
            lateral_mm=position,
            loss_scale,
        )
        for position in positions
    ])
    (position_mm=positions, amplitude=amplitude)
end

function axial_peak_mm(
    design::FresnelDesign,
    topology::Symbol,
    frequency_hz::Real;
    range_mm=35.0:0.5:90.0,
    loss_scale::Real=1.0,
)
    topology == :bare && return NaN
    positions = collect(range_mm)
    amplitudes = [
        abs(scalar_aperture_transfer(design, frequency_hz, topology; axial_mm=z, loss_scale)) /
        abs(scalar_aperture_transfer(design, frequency_hz, :bare; axial_mm=z))
        for z in positions
    ]
    positions[argmax(amplitudes)]
end

end
