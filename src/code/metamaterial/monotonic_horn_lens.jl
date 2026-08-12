module MonotonicHornLens

export MonotonicHornLensConfig,
       aperture_centers_mm,
       geometric_delays_s,
       normalized_delays,
       smooth_centerline_mm,
       smooth_slope,
       smooth_second_derivative,
       smooth_arc_length_mm,
       smooth_amplitude_mm,
       minimum_inner_radius_mm,
       minimum_common_axial_length_mm,
       synthesize_monotonic_geometry

Base.@kwdef struct MonotonicHornLensConfig
    element_count::Int = 15
    pitch_mm::Float64 = 4.8
    throat_width_mm::Float64 = 1.6
    focal_distance_mm::Float64 = 35.0
    receiver_speed_m_s::Float64 = 2340.0
    minimum_inner_radius_mm::Float64 = 3.93
    arc_samples::Int = 4001
    curvature_samples::Int = 20001
end

function validate(config::MonotonicHornLensConfig)
    isodd(config.element_count) || throw(ArgumentError("element count must be odd"))
    config.element_count >= 3 || throw(ArgumentError("at least three channels are required"))
    config.pitch_mm > config.throat_width_mm > 0 ||
        throw(ArgumentError("the throat must fit inside one pitch"))
    config.focal_distance_mm > 0 || throw(ArgumentError("focus distance must be positive"))
    config.receiver_speed_m_s > 0 || throw(ArgumentError("receiver speed must be positive"))
    config.minimum_inner_radius_mm > 0 ||
        throw(ArgumentError("minimum inner radius must be positive"))
    config.arc_samples >= 101 || throw(ArgumentError("arc integration is undersampled"))
    config.curvature_samples >= 1001 || throw(ArgumentError("curvature is undersampled"))
    nothing
end

function aperture_centers_mm(config::MonotonicHornLensConfig=MonotonicHornLensConfig())
    validate(config)
    half = (config.element_count - 1) ÷ 2
    collect((-half):half) .* config.pitch_mm
end

"""Unwrapped non-negative delay that equalizes propagation to the requested focus."""
function geometric_delays_s(config::MonotonicHornLensConfig=MonotonicHornLensConfig())
    centers_m = aperture_centers_mm(config) .* 1e-3
    focus_m = config.focal_distance_mm * 1e-3
    distance_m = hypot.(focus_m, centers_m)
    (maximum(distance_m) .- distance_m) ./ config.receiver_speed_m_s
end

function normalized_delays(config::MonotonicHornLensConfig=MonotonicHornLensConfig())
    delays = geometric_delays_s(config)
    delays ./ maximum(delays)
end

function smooth_centerline_mm(x_mm::Real, axial_length_mm::Real, amplitude_mm::Real)
    t = clamp(Float64(x_mm) / axial_length_mm, 0.0, 1.0)
    Float64(amplitude_mm) * sinpi(t)^4
end

function smooth_slope(x_mm::Real, axial_length_mm::Real, amplitude_mm::Real)
    t = clamp(Float64(x_mm) / axial_length_mm, 0.0, 1.0)
    4Float64(amplitude_mm) * pi / axial_length_mm * sinpi(t)^3 * cospi(t)
end

function smooth_second_derivative(x_mm::Real, axial_length_mm::Real, amplitude_mm::Real)
    t = clamp(Float64(x_mm) / axial_length_mm, 0.0, 1.0)
    4Float64(amplitude_mm) * pi^2 / axial_length_mm^2 *
    (3sinpi(t)^2 * cospi(t)^2 - sinpi(t)^4)
end

function smooth_arc_length_mm(
    axial_length_mm::Real,
    amplitude_mm::Real;
    samples::Integer=4001,
)
    axial_length_mm > 0 || throw(ArgumentError("axial length must be positive"))
    amplitude_mm >= 0 || throw(ArgumentError("amplitude must be non-negative"))
    samples >= 101 || throw(ArgumentError("arc integration is undersampled"))
    x_mm = range(0.0, Float64(axial_length_mm); length=samples)
    integrand = hypot.(1.0, smooth_slope.(x_mm, axial_length_mm, amplitude_mm))
    step_mm = Float64(axial_length_mm) / (samples - 1)
    step_mm * (sum(integrand) - (first(integrand) + last(integrand)) / 2)
end

function smooth_amplitude_mm(
    axial_length_mm::Real,
    extra_path_mm::Real;
    samples::Integer=4001,
)
    extra_path_mm >= 0 || throw(ArgumentError("extra path must be non-negative"))
    iszero(extra_path_mm) && return 0.0
    target_mm = Float64(axial_length_mm + extra_path_mm)
    lower_mm, upper_mm = 0.0, max(Float64(axial_length_mm), 1.0)
    while smooth_arc_length_mm(axial_length_mm, upper_mm; samples) < target_mm
        upper_mm *= 2
    end
    for _ in 1:70
        middle_mm = (lower_mm + upper_mm) / 2
        if smooth_arc_length_mm(axial_length_mm, middle_mm; samples) < target_mm
            lower_mm = middle_mm
        else
            upper_mm = middle_mm
        end
    end
    (lower_mm + upper_mm) / 2
end

function minimum_inner_radius_mm(
    axial_length_mm::Real,
    amplitude_mm::Real,
    throat_width_mm::Real;
    samples::Integer=20001,
)
    iszero(amplitude_mm) && return Inf
    x_mm = range(0.0, Float64(axial_length_mm); length=samples)
    slope = smooth_slope.(x_mm, axial_length_mm, amplitude_mm)
    second = smooth_second_derivative.(x_mm, axial_length_mm, amplitude_mm)
    maximum_curvature_per_mm = maximum(abs.(second) ./ (1 .+ slope .^ 2) .^ 1.5)
    inv(maximum_curvature_per_mm) - Float64(throat_width_mm) / 2
end

function minimum_common_axial_length_mm(
    maximum_extra_path_mm::Real;
    config::MonotonicHornLensConfig=MonotonicHornLensConfig(),
)
    validate(config)
    maximum_extra_path_mm >= 0 || throw(ArgumentError("extra path must be non-negative"))
    iszero(maximum_extra_path_mm) && return 0.0
    function radius_margin(axial_length_mm)
        amplitude_mm = smooth_amplitude_mm(
            axial_length_mm,
            maximum_extra_path_mm;
            samples=config.arc_samples,
        )
        minimum_inner_radius_mm(
            axial_length_mm,
            amplitude_mm,
            config.throat_width_mm;
            samples=config.curvature_samples,
        ) - config.minimum_inner_radius_mm
    end
    lower_mm = max(Float64(maximum_extra_path_mm), config.throat_width_mm)
    upper_mm = max(2lower_mm, 20.0)
    while radius_margin(upper_mm) < 0
        upper_mm *= 1.5
    end
    for _ in 1:60
        middle_mm = (lower_mm + upper_mm) / 2
        radius_margin(middle_mm) < 0 ? (lower_mm = middle_mm) : (upper_mm = middle_mm)
    end
    (lower_mm + upper_mm) / 2
end

function synthesize_monotonic_geometry(
    maximum_extra_path_mm::Real;
    config::MonotonicHornLensConfig=MonotonicHornLensConfig(),
    axial_length_mm::Union{Nothing, Real}=nothing,
)
    validate(config)
    fractions = normalized_delays(config)
    extra_path_mm = Float64(maximum_extra_path_mm) .* fractions
    common_axial_mm = isnothing(axial_length_mm) ?
        minimum_common_axial_length_mm(maximum_extra_path_mm; config) :
        Float64(axial_length_mm)
    amplitude_mm = [
        smooth_amplitude_mm(common_axial_mm, extra; samples=config.arc_samples)
        for extra in extra_path_mm
    ]
    inner_radius_mm = [
        minimum_inner_radius_mm(
            common_axial_mm,
            amplitude,
            config.throat_width_mm;
            samples=config.curvature_samples,
        )
        for amplitude in amplitude_mm
    ]
    minimum(inner_radius_mm) + 1e-8 >= config.minimum_inner_radius_mm ||
        error("synthesized guide violates the curvature constraint")
    (
        centers_mm=aperture_centers_mm(config),
        delay_s=geometric_delays_s(config),
        delay_fraction=fractions,
        extra_path_mm,
        axial_length_mm=common_axial_mm,
        path_length_mm=common_axial_mm .+ extra_path_mm,
        amplitude_mm,
        inner_radius_mm,
    )
end

end
