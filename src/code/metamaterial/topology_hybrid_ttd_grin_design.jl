module TopologyHybridTTDGRINDesign

export ElasticConstituent,
       HybridTTDGRINConfig,
       aluminium_6061,
       copper_c110,
       channel_centers_mm,
       geometric_delay_s,
       maximum_delay_s,
       guided_material,
       guide_width_mm,
       slow_insert_profile,
       smooth_arc_length_mm,
       smooth_amplitude_mm,
       minimum_inner_radius_mm,
       design_summary,
       validate

Base.@kwdef struct ElasticConstituent
    name::Symbol
    density_kg_m3::Float64
    young_modulus_pa::Float64
    poisson_ratio::Float64
end

Base.@kwdef struct HybridTTDGRINConfig
    frequency_hz::Float64 = 242.0e3
    pressure_wave_speed_m_s::Float64 = 6122.102437409232
    calibrated_aluminium_guide_speed_m_s::Float64 = 4793.0
    calibrated_path_delay_speed_m_s::Float64 = 2930.2208499011267
    slow_delay_calibration_factor::Float64 = 0.755190146001182
    aperture_height_mm::Float64 = 71.6
    channel_count::Int = 15
    guide_width_mm::Float64 = 3.2
    device_depth_mm::Float64 = 55.0
    lead_length_mm::Float64 = 6.0
    focal_depth_mm::Float64 = 60.0
    path_delay_fraction::Float64 = 0.60
    transition_length_mm::Float64 = 2.0
    minimum_feature_mm::Float64 = 0.8
    minimum_inner_radius_mm::Float64 = 3.93
    profile_samples::Int = 8001
end

aluminium_6061() = ElasticConstituent(
    name=:al6061,
    density_kg_m3=2700.0,
    young_modulus_pa=68.3e9,
    poisson_ratio=0.33,
)

copper_c110() = ElasticConstituent(
    name=:copper_c110,
    density_kg_m3=8960.0,
    young_modulus_pa=110.0e9,
    poisson_ratio=0.34,
)

function validate(config::HybridTTDGRINConfig)
    config.frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    isodd(config.channel_count) || throw(ArgumentError("channel count must be odd"))
    config.channel_count >= 3 || throw(ArgumentError("at least three channels are required"))
    config.aperture_height_mm > config.guide_width_mm > 0 ||
        throw(ArgumentError("guide must fit inside the aperture"))
    2config.lead_length_mm < config.device_depth_mm ||
        throw(ArgumentError("leads leave no active route"))
    0 < config.path_delay_fraction < 1 ||
        throw(ArgumentError("hybrid split must lie strictly between zero and one"))
    config.calibrated_path_delay_speed_m_s > 0 ||
        throw(ArgumentError("path-delay calibration speed must be positive"))
    0 < config.slow_delay_calibration_factor <= 1 ||
        throw(ArgumentError("slow-delay calibration factor must lie in (0,1]"))
    config.transition_length_mm >= config.minimum_feature_mm ||
        throw(ArgumentError("graded transition is below the minimum feature"))
    config.profile_samples >= 1001 || throw(ArgumentError("route integration is undersampled"))
    config
end

function channel_centers_mm(config::HybridTTDGRINConfig=HybridTTDGRINConfig())
    validate(config)
    pitch = (config.aperture_height_mm - config.guide_width_mm) /
            (config.channel_count - 1)
    half = (config.channel_count - 1) ÷ 2
    collect((-half):half) .* pitch
end

function geometric_delay_s(y_mm::Real, config::HybridTTDGRINConfig=HybridTTDGRINConfig())
    centers = channel_centers_mm(config)
    edge_distance_mm = hypot(config.focal_depth_mm, maximum(abs, centers))
    local_distance_mm = hypot(config.focal_depth_mm, Float64(y_mm))
    (edge_distance_mm - local_distance_mm) * 1e-3 / config.pressure_wave_speed_m_s
end

maximum_delay_s(config::HybridTTDGRINConfig=HybridTTDGRINConfig()) =
    geometric_delay_s(0.0, config)

"""
Voigt Al--Cu mixture used as an explicit screening assumption, not as a
validated constitutive law for a particular sintered coupon. The guided speed
is scaled from the measured project Al-strip speed by sqrt((E/rho)/(E_Al/rho_Al)).
"""
function guided_material(copper_fraction::Real, config::HybridTTDGRINConfig=HybridTTDGRINConfig())
    validate(config)
    fraction = Float64(copper_fraction)
    0 <= fraction <= 1 || throw(ArgumentError("copper fraction must lie in [0,1]"))
    al = aluminium_6061()
    cu = copper_c110()
    density = (1 - fraction) * al.density_kg_m3 + fraction * cu.density_kg_m3
    young = (1 - fraction) * al.young_modulus_pa + fraction * cu.young_modulus_pa
    poisson = (1 - fraction) * al.poisson_ratio + fraction * cu.poisson_ratio
    speed = config.calibrated_aluminium_guide_speed_m_s *
            sqrt((young / density) / (al.young_modulus_pa / al.density_kg_m3))
    (
        copper_fraction=fraction,
        density_kg_m3=density,
        young_modulus_pa=young,
        poisson_ratio=poisson,
        guided_speed_m_s=speed,
        line_impedance_per_area=density * speed,
    )
end

"""Area-equivalent 2D guide width that keeps rho*c*A equal to the Al lead."""
function guide_width_mm(copper_fraction::Real, config::HybridTTDGRINConfig=HybridTTDGRINConfig())
    material = guided_material(copper_fraction, config)
    al = guided_material(0.0, config)
    config.guide_width_mm * al.line_impedance_per_area / material.line_impedance_per_area
end

function slow_slowness_excess_s_m(copper_fraction, config)
    1 / guided_material(copper_fraction, config).guided_speed_m_s -
    1 / config.calibrated_aluminium_guide_speed_m_s
end

function trapezoid_average(function_value; samples=4001)
    fractions = range(0.0, 1.0; length=samples)
    values = function_value.(fractions)
    (sum(values) - (first(values) + last(values)) / 2) / (samples - 1)
end

"""
Return the symmetric linear-gradient / pure-Cu / linear-gradient insert.
Lengths are centreline lengths; the two graded ramps are always retained.
"""
function slow_insert_profile(
    delay_s::Real=maximum_delay_s(),
    config::HybridTTDGRINConfig=HybridTTDGRINConfig();
    transition_steps::Integer=5,
)
    validate(config)
    transition_steps >= 2 || throw(ArgumentError("at least two transition steps are required"))
    target_slow_delay = (1 - config.path_delay_fraction) * Float64(delay_s) /
                        config.slow_delay_calibration_factor
    target_slow_delay >= 0 || throw(ArgumentError("delay must be non-negative"))
    if iszero(target_slow_delay)
        return (
            target_delay_s=0.0,
            transition_length_mm=0.0,
            maximum_copper_fraction=0.0,
            plateau_length_mm=0.0,
            total_length_mm=0.0,
            minimum_width_mm=config.guide_width_mm,
            segments=NamedTuple[],
        )
    end
    ramp_average = trapezoid_average(
        fraction -> slow_slowness_excess_s_m(fraction, config),
    )
    full_ramp_delay = 2config.transition_length_mm * 1e-3 * ramp_average
    copper_excess = slow_slowness_excess_s_m(1.0, config)
    maximum_fraction = 1.0
    plateau_length_mm = 0.0
    if target_slow_delay < full_ramp_delay
        lower, upper = 0.0, 1.0
        for _ in 1:70
            middle = (lower + upper) / 2
            average = trapezoid_average(
                coordinate -> slow_slowness_excess_s_m(middle * coordinate, config),
            )
            trial = 2config.transition_length_mm * 1e-3 * average
            trial < target_slow_delay ? (lower = middle) : (upper = middle)
        end
        maximum_fraction = (lower + upper) / 2
    else
        plateau_length_mm = (target_slow_delay - full_ramp_delay) /
                            copper_excess * 1e3
    end
    fractions = collect(range(0.0, maximum_fraction; length=transition_steps + 1))
    step_length_mm = config.transition_length_mm / transition_steps
    segments = NamedTuple[]
    for fraction in fractions[2:end]
        push!(segments, (
            copper_fraction=fraction,
            length_mm=step_length_mm,
            width_mm=guide_width_mm(fraction, config),
        ))
    end
    plateau_length_mm > 1e-10 && push!(segments, (
        copper_fraction=maximum_fraction,
        length_mm=plateau_length_mm,
        width_mm=guide_width_mm(maximum_fraction, config),
    ))
    for fraction in reverse(fractions[1:(end - 1)])
        push!(segments, (
            copper_fraction=fraction,
            length_mm=step_length_mm,
            width_mm=guide_width_mm(fraction, config),
        ))
    end
    (
        target_delay_s=target_slow_delay,
        transition_length_mm=config.transition_length_mm,
        maximum_copper_fraction=maximum_fraction,
        plateau_length_mm,
        total_length_mm=2config.transition_length_mm + plateau_length_mm,
        minimum_width_mm=minimum(getproperty.(segments, :width_mm)),
        segments,
    )
end

smooth_centerline_mm(x_mm, axial_length_mm, amplitude_mm) =
    Float64(amplitude_mm) * sinpi(clamp(Float64(x_mm) / axial_length_mm, 0.0, 1.0))^4

smooth_slope(x_mm, axial_length_mm, amplitude_mm) = begin
    t = clamp(Float64(x_mm) / axial_length_mm, 0.0, 1.0)
    4Float64(amplitude_mm) * pi / axial_length_mm * sinpi(t)^3 * cospi(t)
end

smooth_second(x_mm, axial_length_mm, amplitude_mm) = begin
    t = clamp(Float64(x_mm) / axial_length_mm, 0.0, 1.0)
    4Float64(amplitude_mm) * pi^2 / axial_length_mm^2 *
    (3sinpi(t)^2 * cospi(t)^2 - sinpi(t)^4)
end

function smooth_arc_length_mm(axial_length_mm::Real, amplitude_mm::Real; samples::Integer=8001)
    x = range(0.0, Float64(axial_length_mm); length=samples)
    values = hypot.(1.0, smooth_slope.(x, axial_length_mm, amplitude_mm))
    step = Float64(axial_length_mm) / (samples - 1)
    step * (sum(values) - (first(values) + last(values)) / 2)
end

function smooth_amplitude_mm(
    axial_length_mm::Real,
    extra_path_mm::Real;
    samples::Integer=8001,
)
    extra_path_mm >= 0 || throw(ArgumentError("extra path must be non-negative"))
    iszero(extra_path_mm) && return 0.0
    target = Float64(axial_length_mm + extra_path_mm)
    lower, upper = 0.0, max(Float64(axial_length_mm), 1.0)
    while smooth_arc_length_mm(axial_length_mm, upper; samples) < target
        upper *= 2
    end
    for _ in 1:70
        middle = (lower + upper) / 2
        if smooth_arc_length_mm(axial_length_mm, middle; samples) < target
            lower = middle
        else
            upper = middle
        end
    end
    (lower + upper) / 2
end

function minimum_inner_radius_mm(
    axial_length_mm::Real,
    amplitude_mm::Real,
    width_mm::Real;
    samples::Integer=20001,
)
    iszero(amplitude_mm) && return Inf
    x = range(0.0, Float64(axial_length_mm); length=samples)
    slope = smooth_slope.(x, axial_length_mm, amplitude_mm)
    second = smooth_second.(x, axial_length_mm, amplitude_mm)
    curvature = maximum(abs.(second) ./ (1 .+ slope .^ 2) .^ 1.5)
    inv(curvature) - Float64(width_mm) / 2
end

function route_for_extra_path(extra_path_mm, config)
    active_length = config.device_depth_mm - 2config.lead_length_mm
    amplitude = smooth_amplitude_mm(
        active_length,
        extra_path_mm;
        samples=config.profile_samples,
    )
    (
        active_length_mm=active_length,
        extra_path_mm=Float64(extra_path_mm),
        amplitude_mm=amplitude,
        inner_radius_mm=minimum_inner_radius_mm(
            active_length,
            amplitude,
            config.guide_width_mm,
        ),
    )
end

function minimum_active_length_mm(extra_path_mm, config)
    lower = max(Float64(extra_path_mm), config.guide_width_mm)
    upper = config.device_depth_mm * 2
    margin(length_mm) = begin
        amplitude = smooth_amplitude_mm(length_mm, extra_path_mm; samples=config.profile_samples)
        minimum_inner_radius_mm(length_mm, amplitude, config.guide_width_mm) -
        config.minimum_inner_radius_mm
    end
    while margin(upper) < 0
        upper *= 1.5
    end
    for _ in 1:60
        middle = (lower + upper) / 2
        margin(middle) < 0 ? (lower = middle) : (upper = middle)
    end
    (lower + upper) / 2
end

function design_summary(config::HybridTTDGRINConfig=HybridTTDGRINConfig())
    validate(config)
    delay = maximum_delay_s(config)
    pure_extra = delay * config.calibrated_path_delay_speed_m_s * 1e3
    hybrid_extra = config.path_delay_fraction * pure_extra
    pure_route = route_for_extra_path(pure_extra, config)
    hybrid_route = route_for_extra_path(hybrid_extra, config)
    pure_minimum_active = minimum_active_length_mm(pure_extra, config)
    pure_compliant_amplitude = smooth_amplitude_mm(
        pure_minimum_active,
        pure_extra;
        samples=config.profile_samples,
    )
    insert = slow_insert_profile(delay, config)
    (
        frequency_hz=config.frequency_hz,
        aperture_height_mm=config.aperture_height_mm,
        channel_count=config.channel_count,
        pitch_mm=channel_centers_mm(config)[2] - channel_centers_mm(config)[1],
        maximum_channel_delay_us=delay * 1e6,
        pure_ttd_extra_path_mm=pure_extra,
        hybrid_extra_path_mm=hybrid_extra,
        extra_path_reduction_fraction=1 - hybrid_extra / pure_extra,
        slow_wave_delay_us=(1 - config.path_delay_fraction) * delay * 1e6,
        slow_insert_length_mm=insert.total_length_mm,
        copper_plateau_length_mm=insert.plateau_length_mm,
        minimum_insert_width_mm=insert.minimum_width_mm,
        hybrid_amplitude_mm=hybrid_route.amplitude_mm,
        hybrid_inner_radius_mm=hybrid_route.inner_radius_mm,
        pure_same_depth_amplitude_mm=pure_route.amplitude_mm,
        pure_same_depth_inner_radius_mm=pure_route.inner_radius_mm,
        pure_curvature_compliant_device_depth_mm=
            pure_minimum_active + 2config.lead_length_mm,
        pure_curvature_compliant_amplitude_mm=pure_compliant_amplitude,
        out_of_plane_excursion_reduction_fraction=
            1 - hybrid_route.amplitude_mm / pure_compliant_amplitude,
        hybrid_manufacturable=
            hybrid_route.inner_radius_mm >= config.minimum_inner_radius_mm &&
            insert.minimum_width_mm >= config.minimum_feature_mm &&
            insert.total_length_mm < hybrid_route.active_length_mm,
    )
end

end
