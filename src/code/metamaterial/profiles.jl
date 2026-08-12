module MetamaterialProfiles

export WallProfile,
       LegacySinusoidalProfile,
       SinusoidalProfile,
       ExponentialProfile,
       PowerProfile,
       RoundedNotchProfile,
       CoupledConstrictionProfile,
       GeometryConfig,
       amplitude_mm,
       indentation,
       generate_boundary_points,
       minimum_gap_mm,
       profile_slug,
       validate_geometry

abstract type WallProfile end

"""Profile reproducing the wall formula used by the original project."""
struct LegacySinusoidalProfile <: WallProfile
    amplitude_mm::Float64
    correction_percent::Float64

    function LegacySinusoidalProfile(amplitude_mm::Real; correction_percent::Real=4.0)
        amplitude_mm >= 0 || throw(ArgumentError("amplitude must be non-negative"))
        correction_percent >= 0 || throw(ArgumentError("correction_percent must be non-negative"))
        new(Float64(amplitude_mm), Float64(correction_percent))
    end
end

"""Smooth periodic raised-cosine profile with a prescribed number of periods."""
struct SinusoidalProfile <: WallProfile
    amplitude_mm::Float64
    periods::Int

    function SinusoidalProfile(amplitude_mm::Real; periods::Integer=2)
        amplitude_mm >= 0 || throw(ArgumentError("amplitude must be non-negative"))
        periods > 0 || throw(ArgumentError("periods must be positive"))
        new(Float64(amplitude_mm), Int(periods))
    end
end

"""
Smooth periodic exponential deformation of a raised cosine.

`sharpness -> 0` recovers `SinusoidalProfile`; positive values localize the
indentation near its maxima.
"""
struct ExponentialProfile <: WallProfile
    amplitude_mm::Float64
    periods::Int
    sharpness::Float64

    function ExponentialProfile(amplitude_mm::Real; periods::Integer=2, sharpness::Real=1.0)
        amplitude_mm >= 0 || throw(ArgumentError("amplitude must be non-negative"))
        periods > 0 || throw(ArgumentError("periods must be positive"))
        isfinite(sharpness) || throw(ArgumentError("sharpness must be finite"))
        new(Float64(amplitude_mm), Int(periods), Float64(sharpness))
    end
end

"""Smooth raised-cosine power family. Powers below one are rejected to retain C1 joins."""
struct PowerProfile <: WallProfile
    amplitude_mm::Float64
    periods::Int
    power::Float64

    function PowerProfile(amplitude_mm::Real; periods::Integer=2, power::Real=2.0)
        amplitude_mm >= 0 || throw(ArgumentError("amplitude must be non-negative"))
        periods > 0 || throw(ArgumentError("periods must be positive"))
        power >= 1 || throw(ArgumentError("power must be at least one for smooth periodic joins"))
        new(Float64(amplitude_mm), Int(periods), Float64(power))
    end
end

"""
Alternating one-sided U-notches with vertical mouths and semicircular tips.

Unlike the smooth periodic families, only one wall is indented at each notch.
The alternating bottom/top sequence creates the sharp zig-zag channel sketched
by the user while `notch_width_mm/2` supplies a finite inner-tip radius.
"""
struct RoundedNotchProfile <: WallProfile
    amplitude_mm::Float64
    notch_count::Int
    notch_width_mm::Float64
    end_margin_mm::Float64
    start_from_bottom::Bool

    function RoundedNotchProfile(
        amplitude_mm::Real;
        notch_count::Integer=4,
        notch_width_mm::Real=1.4,
        end_margin_mm::Real=1.2,
        start_from_bottom::Bool=true,
    )
        amplitude_mm >= 0 || throw(ArgumentError("amplitude must be non-negative"))
        notch_count >= 2 || throw(ArgumentError("at least two alternating notches are required"))
        iseven(notch_count) || throw(ArgumentError("notch count must be even for equal wall loading"))
        notch_width_mm > 0 || throw(ArgumentError("notch width must be positive"))
        end_margin_mm > 0 || throw(ArgumentError("notch end margin must be positive"))
        iszero(amplitude_mm) || amplitude_mm >= notch_width_mm / 2 ||
            throw(ArgumentError("notch depth must be at least the rounded-tip radius"))
        new(
            Float64(amplitude_mm),
            Int(notch_count),
            Float64(notch_width_mm),
            Float64(end_margin_mm),
            start_from_bottom,
        )
    end
end

const StaggeredPeriodicProfile = Union{
    SinusoidalProfile,
    ExponentialProfile,
    PowerProfile,
}

"""
One or more localized, independently positioned symmetric constrictions.

Each entry describes the full gap between the upper and lower walls at the
centre of a raised-cosine constriction. Centres and widths are absolute
millimetre coordinates, which makes their spacing independent of total length.
"""
struct CoupledConstrictionProfile <: WallProfile
    minimum_gaps_mm::Vector{Float64}
    widths_mm::Vector{Float64}
    centers_mm::Vector{Float64}
    nominal_height_mm::Float64

    function CoupledConstrictionProfile(
        minimum_gaps_mm::AbstractVector{<:Real},
        widths_mm::AbstractVector{<:Real},
        centers_mm::AbstractVector{<:Real};
        height_mm::Real=7.0,
    )
        count = length(minimum_gaps_mm)
        count > 0 || throw(ArgumentError("at least one constriction is required"))
        length(widths_mm) == count || throw(DimensionMismatch("one width is required per constriction"))
        length(centers_mm) == count || throw(DimensionMismatch("one centre is required per constriction"))
        height = Float64(height_mm)
        height > 0 || throw(ArgumentError("height_mm must be positive"))
        gaps = Float64.(minimum_gaps_mm)
        widths = Float64.(widths_mm)
        centers = Float64.(centers_mm)
        all(0 .< gaps .< height) || throw(ArgumentError("gaps must lie between zero and height_mm"))
        all(widths .> 0) || throw(ArgumentError("widths must be positive"))
        all(isfinite, centers) || throw(ArgumentError("centres must be finite"))
        new(gaps, widths, centers, height)
    end
end

Base.@kwdef struct GeometryConfig
    length_mm::Float64 = 17.0
    height_mm::Float64 = 7.0
    end_margin_mm::Float64 = 0.4
    samples::Int = 200
end

raised_cosine(phase::Real) = 0.5 * (1.0 - cos(2.0 * pi * phase))

amplitude_mm(profile::WallProfile) = profile.amplitude_mm
amplitude_mm(profile::CoupledConstrictionProfile) =
    maximum((profile.nominal_height_mm .- profile.minimum_gaps_mm) ./ 2)

function indentation(profile::SinusoidalProfile, ξ::Real)
    profile.amplitude_mm * raised_cosine(profile.periods * ξ)
end

function indentation(profile::ExponentialProfile, ξ::Real)
    s = raised_cosine(profile.periods * ξ)
    κ = profile.sharpness
    factor = abs(κ) < 1.0e-8 ? s : expm1(κ * s) / expm1(κ)
    profile.amplitude_mm * factor
end

function indentation(profile::PowerProfile, ξ::Real)
    profile.amplitude_mm * raised_cosine(profile.periods * ξ)^profile.power
end

function notch_centers_mm(profile::RoundedNotchProfile, config::GeometryConfig=GeometryConfig())
    active_span = config.length_mm - 2profile.end_margin_mm
    [
        profile.end_margin_mm + active_span * (index - 0.5) / profile.notch_count
        for index in 1:profile.notch_count
    ]
end

function rounded_notch_depth_at_x(
    profile::RoundedNotchProfile,
    config::GeometryConfig,
    x_mm::Real,
    from_bottom::Bool,
)
    iszero(profile.amplitude_mm) && return 0.0
    radius = profile.notch_width_mm / 2
    stem = profile.amplitude_mm - radius
    local_depth = 0.0
    for (index, center) in enumerate(notch_centers_mm(profile, config))
        belongs_to_wall = isodd(index) == profile.start_from_bottom
        belongs_to_wall == from_bottom || continue
        offset = abs(Float64(x_mm) - center)
        if offset <= radius
            local_depth = max(
                local_depth,
                stem + sqrt(max(0.0, radius^2 - offset^2)),
            )
        end
    end
    local_depth
end

function localized_depth_mm(profile::CoupledConstrictionProfile, x_mm::Real)
    maximum(zip(profile.minimum_gaps_mm, profile.widths_mm, profile.centers_mm); init=0.0) do (gap, width, center)
        offset = abs(Float64(x_mm) - center)
        offset >= width / 2 && return 0.0
        depth = (profile.nominal_height_mm - gap) / 2
        depth * 0.5 * (1 + cos(2pi * offset / width))
    end
end

function indentation(profile::LegacySinusoidalProfile, ξ::Real)
    angle_deg = 4.0 * 180.0 * ξ
    correction_mm = profile.amplitude_mm * profile.correction_percent / 100.0
    if ξ <= 0.25 || ξ >= 0.75
        return 0.5 * profile.amplitude_mm * (1.0 - cosd(angle_deg))
    end
    (0.5 * profile.amplitude_mm - 0.5 * correction_mm) *
        (1.0 - cosd(angle_deg)) + correction_mm
end

function legacy_indentation(profile::LegacySinusoidalProfile, i::Integer, samples::Integer)
    angle_deg = 4.0 * 180.0 * i / (samples + 1.0)
    correction_mm = profile.amplitude_mm * profile.correction_percent / 100.0
    if i <= samples / 4 || i >= 3 * samples / 4
        return 0.5 * profile.amplitude_mm * (1.0 - cosd(angle_deg))
    end
    (0.5 * profile.amplitude_mm - 0.5 * correction_mm) *
        (1.0 - cosd(angle_deg)) + correction_mm
end

function active_span_mm(profile::LegacySinusoidalProfile, config::GeometryConfig)
    4.0 * (config.length_mm - 2.0 * config.end_margin_mm) / 5.0
end

function active_span_mm(profile::StaggeredPeriodicProfile, config::GeometryConfig)
    available_span = config.length_mm - 2.0 * config.end_margin_mm
    available_span * 2.0 * profile.periods / (2.0 * profile.periods + 1.0)
end

function half_period_shift_mm(profile::StaggeredPeriodicProfile, config::GeometryConfig)
    active_span_mm(profile, config) / (2.0 * profile.periods)
end

function depth_at_x(profile::WallProfile, config::GeometryConfig, x_mm::Real)
    active_span = active_span_mm(profile, config)
    x_start = config.end_margin_mm
    x_end = x_start + active_span
    (x_mm < x_start || x_mm > x_end) && return 0.0

    ξ_position = (x_mm - x_start) / active_span
    if profile isa LegacySinusoidalProfile
        # Continuous counterpart of the original i/N position and i/(N+1)
        # phase conventions.
        i = ξ_position * config.samples
        angle_deg = 4.0 * 180.0 * i / (config.samples + 1.0)
        correction_mm = profile.amplitude_mm * profile.correction_percent / 100.0
        if i <= config.samples / 4 || i >= 3 * config.samples / 4
            return 0.5 * profile.amplitude_mm * (1.0 - cosd(angle_deg))
        end
        return (0.5 * profile.amplitude_mm - 0.5 * correction_mm) *
               (1.0 - cosd(angle_deg)) + correction_mm
    end
    indentation(profile, ξ_position)
end


depth_at_x(profile::RoundedNotchProfile, config::GeometryConfig, x_mm::Real) =
    rounded_notch_depth_at_x(profile, config, x_mm, true)


function top_depth_at_x(
    profile::LegacySinusoidalProfile,
    config::GeometryConfig,
    x_mm::Real,
)
    depth_at_x(profile, config, config.length_mm - x_mm)
end


top_depth_at_x(
    profile::RoundedNotchProfile,
    config::GeometryConfig,
    x_mm::Real,
) = rounded_notch_depth_at_x(profile, config, x_mm, false)


function top_depth_at_x(
    profile::StaggeredPeriodicProfile,
    config::GeometryConfig,
    x_mm::Real,
)
    amplitude = amplitude_mm(profile)
    iszero(amplitude) && return 0.0

    bottom_start = config.end_margin_mm
    active_span = active_span_mm(profile, config)
    shift = half_period_shift_mm(profile, config)
    bottom_end = bottom_start + active_span
    top_start = bottom_start + shift
    top_end = bottom_end + shift
    (x_mm < top_start || x_mm > top_end) && return 0.0

    # Over the common span, y_top = y_bottom + (H - A), so the channel has
    # exactly constant thickness H - A for every periodic profile family.
    # After the bottom profile ends, replay its first half-period to taper the
    # top wall smoothly back to the straight right port.
    source_x = x_mm <= bottom_end ? x_mm : bottom_start + (x_mm - bottom_end)
    clamp(amplitude - depth_at_x(profile, config, source_x), 0.0, amplitude)
end


depth_at_x(profile::CoupledConstrictionProfile, config::GeometryConfig, x_mm::Real) =
    localized_depth_mm(profile, x_mm)

function minimum_gap_mm(profile::WallProfile, config::GeometryConfig=GeometryConfig())
    sample_count = max(4 * config.samples, 256)
    active_span = active_span_mm(profile, config)
    profile_nodes = [
        config.end_margin_mm + active_span * i / config.samples
        for i in 0:config.samples
    ]
    top_nodes = if profile isa StaggeredPeriodicProfile
        shift = half_period_shift_mm(profile, config)
        profile_nodes .+ shift
    else
        config.length_mm .- profile_nodes
    end
    candidates = vcat(
        collect(range(0.0, config.length_mm; length=sample_count + 1)),
        profile_nodes,
        top_nodes,
    )
    minimum(
        config.height_mm -
        depth_at_x(profile, config, x_mm) -
        top_depth_at_x(profile, config, x_mm)
        for x_mm in candidates
    )
end


function minimum_gap_mm(profile::CoupledConstrictionProfile, config::GeometryConfig=GeometryConfig())
    height_offset = config.height_mm - profile.nominal_height_mm
    minimum(profile.minimum_gaps_mm) + height_offset
end


minimum_gap_mm(profile::RoundedNotchProfile, config::GeometryConfig=GeometryConfig()) =
    config.height_mm - profile.amplitude_mm

function validate_geometry(profile::WallProfile, config::GeometryConfig)
    config.length_mm > 0 || throw(ArgumentError("length_mm must be positive"))
    config.height_mm > 0 || throw(ArgumentError("height_mm must be positive"))
    config.end_margin_mm >= 0 || throw(ArgumentError("end_margin_mm must be non-negative"))
    2.0 * config.end_margin_mm < config.length_mm ||
        throw(ArgumentError("end margins leave no active profile span"))
    config.samples >= 8 || throw(ArgumentError("at least eight samples are required"))

    minimum_gap_mm(profile, config) > 1.0e-9 ||
        throw(ArgumentError("opposing indentations intersect or close the waveguide"))
    nothing
end


function validate_geometry(profile::CoupledConstrictionProfile, config::GeometryConfig)
    config.length_mm > 0 || throw(ArgumentError("length_mm must be positive"))
    config.height_mm ≈ profile.nominal_height_mm ||
        throw(ArgumentError("profile and geometry heights must agree"))
    config.samples >= 8 || throw(ArgumentError("at least eight samples are required"))
    for (width, center) in zip(profile.widths_mm, profile.centers_mm)
        center - width / 2 > 0 || throw(ArgumentError("a constriction reaches the left port"))
        center + width / 2 < config.length_mm || throw(ArgumentError("a constriction reaches the right port"))
    end
    ordering = sortperm(profile.centers_mm)
    for (left_index, right_index) in zip(ordering[1:end-1], ordering[2:end])
        left_edge = profile.centers_mm[left_index] + profile.widths_mm[left_index] / 2
        right_edge = profile.centers_mm[right_index] - profile.widths_mm[right_index] / 2
        left_edge <= right_edge || throw(ArgumentError("localized constrictions overlap"))
    end
    minimum_gap_mm(profile, config) > 1.0e-9 ||
        throw(ArgumentError("opposing indentations intersect or close the waveguide"))
    nothing
end


function validate_geometry(profile::RoundedNotchProfile, config::GeometryConfig)
    config.length_mm > 0 || throw(ArgumentError("length_mm must be positive"))
    config.height_mm > 0 || throw(ArgumentError("height_mm must be positive"))
    config.samples >= 8 || throw(ArgumentError("at least eight samples are required"))
    2profile.end_margin_mm < config.length_mm ||
        throw(ArgumentError("rounded-notch end margins leave no active span"))
    profile.amplitude_mm < config.height_mm ||
        throw(ArgumentError("rounded notches close the channel"))
    spacing = (config.length_mm - 2profile.end_margin_mm) / profile.notch_count
    profile.notch_width_mm < spacing ||
        throw(ArgumentError("neighbouring rounded-notch mouths overlap"))
    minimum_gap_mm(profile, config) > 1.0e-9 ||
        throw(ArgumentError("rounded notches close the waveguide"))
    nothing
end

"""
Return bottom and top boundary points in millimetres.

The legacy profile deliberately preserves the active-span convention of the
original mesher. New periodic profiles use staggered walls: the top starts half
a period after the bottom and is its vertical translate over the common span.
This retains a constant channel thickness of `height_mm - amplitude_mm`.
"""
function generate_boundary_points(profile::WallProfile, config::GeometryConfig=GeometryConfig())
    validate_geometry(profile, config)

    active_span = active_span_mm(profile, config)

    bottom = Tuple{Float64, Float64}[(0.0, 0.0)]
    for i in 0:config.samples
        ξ_position = i / config.samples
        x = config.end_margin_mm + active_span * ξ_position

        # The original code used N in the x coordinate and N+1 in the phase.
        # Keep that convention only for legacy datasets.
        depth = if profile isa LegacySinusoidalProfile
            legacy_indentation(profile, i, config.samples)
        else
            indentation(profile, ξ_position)
        end
        push!(bottom, (x, depth))
    end
    push!(bottom, (config.length_mm, 0.0))

    top = if profile isa LegacySinusoidalProfile
        [(config.length_mm - x, config.height_mm - y) for (x, y) in bottom]
    else
        shift = half_period_shift_mm(profile, config)
        top_start = config.end_margin_mm + shift
        top_x = [top_start + active_span * i / config.samples for i in 0:config.samples]
        bottom_end = config.end_margin_mm + active_span
        push!(top_x, bottom_end)
        sort!(unique!(top_x))
        top_increasing = Tuple{Float64, Float64}[(0.0, config.height_mm)]
        append!(top_increasing, [
            (x, config.height_mm - top_depth_at_x(profile, config, x))
            for x in top_x
        ])
        push!(top_increasing, (config.length_mm, config.height_mm))
        reverse(top_increasing)
    end
    bottom, top
end


function generate_boundary_points(
    profile::CoupledConstrictionProfile,
    config::GeometryConfig=GeometryConfig(),
)
    validate_geometry(profile, config)
    x_values = collect(range(0.0, config.length_mm; length=config.samples + 1))
    # Include exact edges and centres so the requested gap and width do not
    # depend on the background sampling grid.
    append!(x_values, profile.centers_mm)
    append!(x_values, profile.centers_mm .- profile.widths_mm ./ 2)
    append!(x_values, profile.centers_mm .+ profile.widths_mm ./ 2)
    sort!(unique!(x_values))
    bottom = [(x, localized_depth_mm(profile, x)) for x in x_values]
    top = reverse([(x, config.height_mm - localized_depth_mm(profile, x)) for x in x_values])
    bottom, top
end


function generate_boundary_points(
    profile::RoundedNotchProfile,
    config::GeometryConfig=GeometryConfig(),
)
    validate_geometry(profile, config)
    depth = profile.amplitude_mm
    iszero(depth) && return (
        [(0.0, 0.0), (config.length_mm, 0.0)],
        [(config.length_mm, config.height_mm), (0.0, config.height_mm)],
    )
    radius = profile.notch_width_mm / 2
    stem = depth - radius
    centers = notch_centers_mm(profile, config)
    bottom_centers = [
        center for (index, center) in enumerate(centers)
        if isodd(index) == profile.start_from_bottom
    ]
    top_centers = [
        center for (index, center) in enumerate(centers)
        if isodd(index) != profile.start_from_bottom
    ]
    arc_samples = max(12, config.samples ÷ (4profile.notch_count))

    bottom = Tuple{Float64, Float64}[(0.0, 0.0)]
    for center in bottom_centers
        left, right = center - radius, center + radius
        push!(bottom, (left, 0.0), (left, stem))
        for theta in range(pi, 0.0; length=arc_samples + 1)
            push!(bottom, (
                center + radius * cos(theta),
                stem + radius * sin(theta),
            ))
        end
        push!(bottom, (right, 0.0))
    end
    push!(bottom, (config.length_mm, 0.0))

    top = Tuple{Float64, Float64}[(config.length_mm, config.height_mm)]
    for center in reverse(top_centers)
        right, left = center + radius, center - radius
        push!(top, (right, config.height_mm), (right, config.height_mm - stem))
        for theta in range(0.0, -pi; length=arc_samples + 1)
            push!(top, (
                center + radius * cos(theta),
                config.height_mm - stem + radius * sin(theta),
            ))
        end
        push!(top, (left, config.height_mm))
    end
    push!(top, (0.0, config.height_mm))
    bottom, top
end

function profile_slug(profile::LegacySinusoidalProfile)
    "A_$(profile.amplitude_mm)"
end

function profile_slug(profile::SinusoidalProfile)
    base = "sin_A_$(profile.amplitude_mm)_N_$(profile.periods)"
    iszero(profile.amplitude_mm) ? base : "$(base)_STG"
end

function profile_slug(profile::ExponentialProfile)
    base = "exp_A_$(profile.amplitude_mm)_N_$(profile.periods)_K_$(profile.sharpness)"
    iszero(profile.amplitude_mm) ? base : "$(base)_STG"
end


function profile_slug(profile::PowerProfile)
    base = "pow_A_$(profile.amplitude_mm)_N_$(profile.periods)_P_$(profile.power)"
    iszero(profile.amplitude_mm) ? base : "$(base)_STG"
end


function profile_slug(profile::RoundedNotchProfile)
    side = profile.start_from_bottom ? "B" : "T"
    "notch_A_$(profile.amplitude_mm)_N_$(profile.notch_count)_W_$(profile.notch_width_mm)_$(side)_STG"
end


function profile_slug(profile::CoupledConstrictionProfile)
    label(values) = join(replace.(string.(values), "." => "p"), "-")
    "coupled_G_$(label(profile.minimum_gaps_mm))_W_$(label(profile.widths_mm))_C_$(label(profile.centers_mm))"
end

end
