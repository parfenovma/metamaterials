module SinusoidalMaterialLens

using FFTW
using LinearAlgebra: I

export ElasticMaterial,
       AbstractLensCell,
       SinusoidalCell,
       RoundedNotchCell,
       MaterialLensConfig,
       LensDesignCandidate,
       RoundedLensDesignCandidate,
       photopolymer,
       aluminium_6061,
       geometric_matching_material,
       quarter_wave_thickness_mm,
       cell_gap_mm,
       cell_transfer,
       focus_spectrum,
       pulse_response,
       optimize_configuration,
       optimize_rounded_configuration,
       design_rows,
       rounded_design_rows,
       notch_centers_mm

"""Isotropic material used by the reduced longitudinal-wave screening model."""
Base.@kwdef struct ElasticMaterial
    name::String
    density_kg_m3::Float64
    pressure_wave_speed_m_s::Float64
    shear_wave_speed_m_s::Float64
    rayleigh_alpha_s_inv::Float64 = 0.0
    rayleigh_beta_s::Float64 = 0.0
end

"""
Nominal project photopolymer. These are the values used by the existing FEM
datasets, not a claim that every printable resin has the same properties.
"""
photopolymer() = ElasticMaterial(
    name="project_photopolymer",
    density_kg_m3=1210.0,
    pressure_wave_speed_m_s=2340.0,
    shear_wave_speed_m_s=1170.0,
    rayleigh_alpha_s_inv=79560.0,
    rayleigh_beta_s=2.5e-9,
)

"""
Nominal 6061 aluminium. The wave speeds correspond to E=68.3 GPa,
rho=2700 kg/m^3 and nu=0.33.
"""
aluminium_6061() = ElasticMaterial(
    name="aluminium_6061",
    density_kg_m3=2700.0,
    pressure_wave_speed_m_s=6122.102437409232,
    shear_wave_speed_m_s=3083.810277185563,
)

abstract type AbstractLensCell end

"""Single manufacturable sinusoidal strip with flat, equal-height end sections."""
Base.@kwdef struct SinusoidalCell <: AbstractLensCell
    length_mm::Float64
    height_mm::Float64 = 7.0
    minimum_gap_mm::Float64 = 7.0
    periods::Int = 1
    end_margin_mm::Float64 = 0.4
end

"""
Alternating one-sided U-notches with sharp mouths and semicircular inner ends.

`minimum_gap_mm` is the ligament left between a notch tip and the opposite
face. Mirrored aperture elements swap `start_from_bottom` to preserve global
symmetry.
"""
Base.@kwdef struct RoundedNotchCell <: AbstractLensCell
    length_mm::Float64
    height_mm::Float64 = 7.0
    minimum_gap_mm::Float64 = 4.0
    notch_count::Int = 4
    notch_width_mm::Float64 = 1.4
    end_margin_mm::Float64 = 1.2
    start_from_bottom::Bool = true
end

cell_gap_mm(cell::SinusoidalCell) = cell.minimum_gap_mm
cell_amplitude_mm(cell::SinusoidalCell) = (cell.height_mm - cell.minimum_gap_mm) / 2
cell_gap_mm(cell::RoundedNotchCell) = cell.minimum_gap_mm
notch_depth_mm(cell::RoundedNotchCell) = cell.height_mm - cell.minimum_gap_mm
notch_tip_radius_mm(cell::RoundedNotchCell) = cell.notch_width_mm / 2

function notch_centers_mm(cell::RoundedNotchCell)
    active = cell.length_mm - 2cell.end_margin_mm
    [
        cell.end_margin_mm + active * (index - 0.5) / cell.notch_count
        for index in 1:cell.notch_count
    ]
end

Base.@kwdef struct MaterialLensConfig
    center_frequency_hz::Float64 = 242.0e3
    pulse_cycles::Float64 = 4.0
    element_count::Int = 15
    element_height_mm::Float64 = 7.0
    slot_width_mm::Float64 = 1.2
    focal_distance_mm::Float64 = 60.0
    quadrature_points::Int = 7
    transfer_segments::Int = 120
    samples_per_period::Int = 32
    record_duration_s::Float64 = 120.0e-6
end

struct LensDesignCandidate
    configuration::Symbol
    lens_material::ElasticMaterial
    output_material::ElasticMaterial
    matching_material::Union{Nothing, ElasticMaterial}
    matching_thickness_mm::Float64
    cells::Vector{SinusoidalCell}
    harmonic_focus_amplitude::Float64
    pulse_peak_amplitude::Float64
    pulse_peak_time_s::Float64
    uniform_pulse_peak_amplitude::Float64
end

struct RoundedLensDesignCandidate
    configuration::Symbol
    lens_material::ElasticMaterial
    output_material::ElasticMaterial
    matching_material::Union{Nothing, ElasticMaterial}
    matching_thickness_mm::Float64
    cells::Vector{RoundedNotchCell}
    harmonic_focus_amplitude::Float64
    pulse_peak_amplitude::Float64
    pulse_peak_time_s::Float64
    uniform_pulse_peak_amplitude::Float64
end

function validate(material::ElasticMaterial)
    material.density_kg_m3 > 0 || throw(ArgumentError("density must be positive"))
    material.pressure_wave_speed_m_s > 0 || throw(ArgumentError("P-wave speed must be positive"))
    material.shear_wave_speed_m_s > 0 || throw(ArgumentError("S-wave speed must be positive"))
    material.pressure_wave_speed_m_s > sqrt(2.0) * material.shear_wave_speed_m_s ||
        throw(ArgumentError("material gives a non-positive first Lame constant"))
    nothing
end

function validate(cell::SinusoidalCell)
    cell.length_mm > 0 || throw(ArgumentError("cell length must be positive"))
    cell.height_mm > 0 || throw(ArgumentError("cell height must be positive"))
    0 < cell.minimum_gap_mm <= cell.height_mm ||
        throw(ArgumentError("minimum gap must lie in (0,height]"))
    cell.periods > 0 || throw(ArgumentError("period count must be positive"))
    2cell.end_margin_mm < cell.length_mm || throw(ArgumentError("end margins close the cell"))
    nothing
end

function validate(cell::RoundedNotchCell)
    cell.length_mm > 0 || throw(ArgumentError("cell length must be positive"))
    cell.height_mm > 0 || throw(ArgumentError("cell height must be positive"))
    0 < cell.minimum_gap_mm <= cell.height_mm ||
        throw(ArgumentError("minimum gap must lie in (0,height]"))
    cell.notch_count > 0 || throw(ArgumentError("notch count must be positive"))
    cell.notch_width_mm > 0 || throw(ArgumentError("notch width must be positive"))
    2cell.end_margin_mm < cell.length_mm || throw(ArgumentError("end margins close the cell"))
    spacing = (cell.length_mm - 2cell.end_margin_mm) / cell.notch_count
    cell.notch_width_mm < spacing || throw(ArgumentError("neighbouring notch mouths overlap"))
    depth = notch_depth_mm(cell)
    if depth > 0
        depth >= notch_tip_radius_mm(cell) ||
            throw(ArgumentError("notch depth is smaller than the rounded tip radius"))
    end
    nothing
end

function validate(config::MaterialLensConfig)
    isodd(config.element_count) || throw(ArgumentError("element count must be odd"))
    config.element_count >= 3 || throw(ArgumentError("at least three elements are required"))
    config.center_frequency_hz > 0 || throw(ArgumentError("center frequency must be positive"))
    config.focal_distance_mm > 0 || throw(ArgumentError("focal distance must be positive"))
    config.transfer_segments >= 8 || throw(ArgumentError("at least eight transfer segments are required"))
    config.samples_per_period >= 8 || throw(ArgumentError("at least eight time samples per period are required"))
    nothing
end

"""
Create the normal-incidence geometric-mean impedance layer. Density and both
wave speeds are geometric means, so the layer matches P and S impedances in
this isotropic approximation while retaining a physical Poisson ratio.
"""
function geometric_matching_material(left::ElasticMaterial, right::ElasticMaterial)
    validate(left)
    validate(right)
    ElasticMaterial(
        name="geometric_impedance_match",
        density_kg_m3=sqrt(left.density_kg_m3 * right.density_kg_m3),
        pressure_wave_speed_m_s=sqrt(
            left.pressure_wave_speed_m_s * right.pressure_wave_speed_m_s,
        ),
        shear_wave_speed_m_s=sqrt(
            left.shear_wave_speed_m_s * right.shear_wave_speed_m_s,
        ),
    )
end

quarter_wave_thickness_mm(material::ElasticMaterial, frequency_hz::Real) =
    1.0e3 * material.pressure_wave_speed_m_s / (4Float64(frequency_hz))

function gap_at(cell::SinusoidalCell, x_mm::Real)
    validate(cell)
    active_length = cell.length_mm - 2cell.end_margin_mm
    if x_mm <= cell.end_margin_mm || x_mm >= cell.length_mm - cell.end_margin_mm
        return cell.height_mm
    end
    xi = (Float64(x_mm) - cell.end_margin_mm) / active_length
    amplitude = cell_amplitude_mm(cell)
    cell.height_mm - 2amplitude * 0.5 * (1 - cos(2pi * cell.periods * xi))
end

function gap_at(cell::RoundedNotchCell, x_mm::Real)
    validate(cell)
    depth = notch_depth_mm(cell)
    iszero(depth) && return cell.height_mm
    radius = notch_tip_radius_mm(cell)
    local_depth = 0.0
    for center in notch_centers_mm(cell)
        offset = abs(Float64(x_mm) - center)
        if offset <= radius
            local_depth = max(local_depth, depth - radius + sqrt(max(0.0, radius^2 - offset^2)))
        end
    end
    cell.height_mm - local_depth
end

function complex_wavenumber_per_m(material::ElasticMaterial, frequency_hz::Real)
    omega = 2pi * Float64(frequency_hz)
    iszero(omega) && return 0.0 + 0.0im
    loss_factor = material.rayleigh_alpha_s_inv / omega + material.rayleigh_beta_s * omega
    (omega / material.pressure_wave_speed_m_s) * (1 - 0.5im * loss_factor)
end

function segment_matrix(
    material::ElasticMaterial,
    frequency_hz::Real,
    length_m::Real,
    area_m2::Real,
)
    area_m2 > 0 || throw(ArgumentError("segment area must be positive"))
    k = complex_wavenumber_per_m(material, frequency_hz)
    impedance = material.density_kg_m3 * material.pressure_wave_speed_m_s / area_m2
    phase = k * Float64(length_m)
    ComplexF64[
        cos(phase) im * impedance * sin(phase)
        im * sin(phase) / impedance cos(phase)
    ]
end

function cell_matrix(
    cell::AbstractLensCell,
    material::ElasticMaterial,
    frequency_hz::Real;
    segments::Integer=120,
)
    validate(cell)
    validate(material)
    segments >= 8 || throw(ArgumentError("at least eight segments are required"))
    dx_mm = cell.length_mm / segments
    matrix = Matrix{ComplexF64}(I, 2, 2)
    for index in 1:segments
        x_mm = (index - 0.5) * dx_mm
        area_m2 = gap_at(cell, x_mm) * 1.0e-3 # unit out-of-plane thickness
        matrix *= segment_matrix(material, frequency_hz, dx_mm * 1.0e-3, area_m2)
    end
    matrix
end

function layer_matrix(
    material::ElasticMaterial,
    frequency_hz::Real,
    thickness_mm::Real,
    height_mm::Real,
)
    thickness_mm >= 0 || throw(ArgumentError("layer thickness must be non-negative"))
    iszero(thickness_mm) && return Matrix{ComplexF64}(I, 2, 2)
    segment_matrix(
        material,
        frequency_hz,
        Float64(thickness_mm) * 1.0e-3,
        Float64(height_mm) * 1.0e-3,
    )
end

"""
Return stress and particle-velocity transmission for unit incident stress.

The transfer matrices use the state `[stress; volume_velocity]`. The velocity
coefficient is the quantity propagated by the reduced focusing model.
"""
function cell_transfer(
    cell::AbstractLensCell,
    frequency_hz::Real,
    lens_material::ElasticMaterial,
    output_material::ElasticMaterial;
    matching_material::Union{Nothing, ElasticMaterial}=nothing,
    matching_thickness_mm::Real=0.0,
    segments::Integer=120,
)
    validate(lens_material)
    validate(output_material)
    matrix = cell_matrix(cell, lens_material, frequency_hz; segments)
    if !isnothing(matching_material) && matching_thickness_mm > 0
        validate(matching_material)
        matrix *= layer_matrix(
            matching_material,
            frequency_hz,
            matching_thickness_mm,
            cell.height_mm,
        )
    end

    input_area_m2 = cell.height_mm * 1.0e-3
    output_area_m2 = input_area_m2
    input_impedance = lens_material.density_kg_m3 *
                      lens_material.pressure_wave_speed_m_s / input_area_m2
    output_impedance = output_material.density_kg_m3 *
                       output_material.pressure_wave_speed_m_s / output_area_m2
    denominator = matrix[1, 1] * output_impedance + matrix[1, 2] +
                  input_impedance * (matrix[2, 1] * output_impedance + matrix[2, 2])
    stress = 2output_impedance / denominator
    velocity = 2input_impedance / denominator
    (stress=ComplexF64(stress), velocity=ComplexF64(velocity))
end

function lens_centers_mm(config::MaterialLensConfig)
    pitch = config.element_height_mm + config.slot_width_mm
    half = (config.element_count - 1) ÷ 2
    collect((-half):half) .* pitch
end

function aperture_kernel(
    center_y_mm::Real,
    frequency_hz::Real,
    output_material::ElasticMaterial,
    config::MaterialLensConfig,
)
    iszero(frequency_hz) && return 0.0 + 0.0im
    dy = config.element_height_mm / config.quadrature_points
    y0 = center_y_mm - config.element_height_mm / 2 + dy / 2
    wavelength_mm = 1.0e3 * output_material.pressure_wave_speed_m_s / frequency_hz
    k_per_mm = 2pi / wavelength_mm
    dy * sum(1:config.quadrature_points) do index
        source_y = y0 + (index - 1) * dy
        distance = hypot(config.focal_distance_mm, source_y)
        exp(-im * (k_per_mm * distance - pi / 4)) / sqrt(wavelength_mm * distance)
    end
end

function focus_spectrum(
    cells::AbstractVector{<:AbstractLensCell},
    frequencies_hz::AbstractVector{<:Real},
    lens_material::ElasticMaterial,
    output_material::ElasticMaterial,
    config::MaterialLensConfig;
    matching_material::Union{Nothing, ElasticMaterial}=nothing,
    matching_thickness_mm::Real=0.0,
)
    length(cells) == config.element_count || throw(DimensionMismatch("one cell per element is required"))
    centers = lens_centers_mm(config)
    unique_cells = unique(cells)
    transfers = Dict(
        cell => ComplexF64[
            cell_transfer(
                cell,
                frequency_hz,
                lens_material,
                output_material;
                matching_material,
                matching_thickness_mm,
                segments=config.transfer_segments,
            ).stress
            for frequency_hz in frequencies_hz
        ]
        for cell in unique_cells
    )
    ComplexF64[
        sum(zip(cells, centers)) do (cell, center)
            # Unit incident stress is common to all three experiments.  The
            # transmitted stress coefficient therefore also gives output
            # particle velocity relative to the all-aluminium reference
            # velocity p_inc/Z_Al.
            transfers[cell][frequency_index] *
            aperture_kernel(center, frequency_hz, output_material, config)
        end
        for (frequency_index, frequency_hz) in enumerate(frequencies_hz)
    ]
end

function pulse_time_grid(config::MaterialLensConfig)
    dt = 1 / (config.center_frequency_hz * config.samples_per_period)
    sample_count = nextpow(2, Int(ceil(config.record_duration_s / dt)))
    collect(0:(sample_count - 1)) .* dt
end

function incident_pulse(time_s::AbstractVector, config::MaterialLensConfig)
    duration = config.pulse_cycles / config.center_frequency_hz
    [
        t < duration ?
        0.5 * (1 - cos(2pi * t / duration)) * sin(2pi * config.center_frequency_hz * t) :
        0.0
        for t in time_s
    ]
end

function pulse_response(
    cells::AbstractVector{<:AbstractLensCell},
    lens_material::ElasticMaterial,
    output_material::ElasticMaterial,
    config::MaterialLensConfig;
    matching_material::Union{Nothing, ElasticMaterial}=nothing,
    matching_thickness_mm::Real=0.0,
)
    validate(config)
    time_s = pulse_time_grid(config)
    incident = incident_pulse(time_s, config)
    incident_spectrum = rfft(incident)
    frequencies_hz = collect(0:(length(incident_spectrum) - 1)) ./
                     (length(time_s) * (time_s[2] - time_s[1]))
    active = findall(abs.(incident_spectrum) .>= 1.0e-4 * maximum(abs, incident_spectrum))
    transfer = zeros(ComplexF64, length(frequencies_hz))
    transfer[active] .= focus_spectrum(
        cells,
        frequencies_hz[active],
        lens_material,
        output_material,
        config;
        matching_material,
        matching_thickness_mm,
    )
    focused = irfft(incident_spectrum .* transfer, length(time_s))
    peak_index = argmax(abs.(focused))
    (
        time_s,
        incident,
        focused,
        frequencies_hz,
        transfer,
        peak_amplitude=abs(focused[peak_index]),
        peak_time_s=time_s[peak_index],
    )
end

function symmetric_cells(group_cells::AbstractVector{<:AbstractLensCell}, config::MaterialLensConfig)
    group_count = (config.element_count + 1) ÷ 2
    length(group_cells) == group_count || throw(DimensionMismatch("wrong symmetric group count"))
    [group_cells[abs(index - group_count) + 1] for index in 1:config.element_count]
end

function mirrored(cell::RoundedNotchCell)
    RoundedNotchCell(
        length_mm=cell.length_mm,
        height_mm=cell.height_mm,
        minimum_gap_mm=cell.minimum_gap_mm,
        notch_count=cell.notch_count,
        notch_width_mm=cell.notch_width_mm,
        end_margin_mm=cell.end_margin_mm,
        start_from_bottom=!cell.start_from_bottom,
    )
end

function symmetric_rounded_cells(
    group_cells::AbstractVector{RoundedNotchCell},
    config::MaterialLensConfig,
)
    half_count = (config.element_count - 1) ÷ 2
    length(group_cells) == half_count + 1 || throw(DimensionMismatch("wrong symmetric group count"))
    RoundedNotchCell[
        let cell=group_cells[abs(offset) + 1]
            offset < 0 ? mirrored(cell) : cell
        end
        for offset in (-half_count):half_count
    ]
end

function target_frequency_shortlist(
    cell_options,
    lens_material,
    output_material,
    config;
    matching_material=nothing,
    matching_thickness_mm=0.0,
    keep::Integer=24,
    fixed_center_index::Union{Nothing, Integer}=nothing,
)
    centers = lens_centers_mm(config)
    absolute_centers = sort(unique(abs.(centers)))
    frequency = config.center_frequency_hz
    contributions = [
        cell_transfer(
            cell,
            frequency,
            lens_material,
            output_material;
            matching_material,
            matching_thickness_mm,
            segments=config.transfer_segments,
        ).stress * (
            iszero(center) ?
            aperture_kernel(0.0, frequency, output_material, config) :
            aperture_kernel(center, frequency, output_material, config) +
            aperture_kernel(-center, frequency, output_material, config)
        )
        for center in absolute_centers, cell in cell_options
    ]

    # For a trial direction phi, maximizing the projection of the total field
    # is separable by symmetric group. Scanning phi therefore avoids an
    # exponential Cartesian product for the 15-element (8-group) aperture.
    designs = Dict{Tuple{Vararg{Int}}, Float64}()
    for phi in range(0.0, 2pi; length=1441)[1:end-1]
        rotation = cis(-phi)
        indices = Tuple(
            !isnothing(fixed_center_index) && group == 1 ?
            Int(fixed_center_index) :
            argmax(real.(view(contributions, group, :) .* rotation))
            for group in axes(contributions, 1)
        )
        amplitude = abs(sum(contributions[group, indices[group]] for group in eachindex(indices)))
        designs[indices] = max(get(designs, indices, -Inf), amplitude)
    end
    ranked = sort(collect(designs); by=last, rev=true)
    [
        (amplitude, collect(indices))
        for (indices, amplitude) in Iterators.take(ranked, min(keep, length(ranked)))
    ]
end

function optimize_rounded_configuration(
    configuration::Symbol;
    config::MaterialLensConfig=MaterialLensConfig(transfer_segments=180),
    lengths_mm::AbstractVector{<:Real}=[20.0, 24.0, 28.0, 32.0, 36.0],
    gaps_mm::AbstractVector{<:Real}=[2.0, 3.0, 4.0, 5.0],
    notch_counts::AbstractVector{<:Integer}=[3, 4, 5],
    notch_widths_mm::AbstractVector{<:Real}=[1.2, 1.6, 2.0],
    matching_thicknesses_mm::AbstractVector{<:Real}=[1.0, 1.4, 1.8, 2.2],
    shortlist_count::Integer=20,
)
    validate(config)
    lens_material, output_material, matching_material = configuration_materials(configuration)
    thicknesses = configuration == :polymer_matched ? Float64.(matching_thicknesses_mm) : [0.0]
    best = nothing

    for length_mm in lengths_mm, thickness_mm in thicknesses
        cell_options = RoundedNotchCell[
            RoundedNotchCell(
                length_mm=Float64(length_mm),
                height_mm=config.element_height_mm,
                minimum_gap_mm=config.element_height_mm,
                notch_count=1,
                notch_width_mm=1.2,
            ),
        ]
        for notch_count in notch_counts, width_mm in notch_widths_mm, gap_mm in gaps_mm
            cell = RoundedNotchCell(
                length_mm=Float64(length_mm),
                height_mm=config.element_height_mm,
                minimum_gap_mm=Float64(gap_mm),
                notch_count=Int(notch_count),
                notch_width_mm=Float64(width_mm),
            )
            try
                validate(cell)
                push!(cell_options, cell)
            catch error
                error isa ArgumentError || rethrow()
            end
        end
        shortlist = target_frequency_shortlist(
            cell_options,
            lens_material,
            output_material,
            config;
            matching_material,
            matching_thickness_mm=thickness_mm,
            keep=shortlist_count,
            fixed_center_index=1,
        )
        for (harmonic_amplitude, indices) in shortlist
            cells = symmetric_rounded_cells(cell_options[indices], config)
            response = pulse_response(
                cells,
                lens_material,
                output_material,
                config;
                matching_material,
                matching_thickness_mm=thickness_mm,
            )
            if isnothing(best) || response.peak_amplitude > best.pulse_peak_amplitude
                best = RoundedLensDesignCandidate(
                    configuration,
                    lens_material,
                    output_material,
                    matching_material,
                    Float64(thickness_mm),
                    cells,
                    harmonic_amplitude,
                    response.peak_amplitude,
                    response.peak_time_s,
                    NaN,
                )
            end
        end
    end
    isnothing(best) && error("empty rounded-notch design search")
    uniform = RoundedNotchCell(
        length_mm=first(best.cells).length_mm,
        height_mm=config.element_height_mm,
        minimum_gap_mm=config.element_height_mm,
        notch_count=1,
        notch_width_mm=1.2,
    )
    uniform_response = pulse_response(
        fill(uniform, config.element_count),
        lens_material,
        output_material,
        config;
        matching_material,
        matching_thickness_mm=best.matching_thickness_mm,
    )
    RoundedLensDesignCandidate(
        best.configuration,
        best.lens_material,
        best.output_material,
        best.matching_material,
        best.matching_thickness_mm,
        best.cells,
        best.harmonic_focus_amplitude,
        best.pulse_peak_amplitude,
        best.pulse_peak_time_s,
        uniform_response.peak_amplitude,
    )
end

function configuration_materials(configuration::Symbol)
    polymer = photopolymer()
    aluminium = aluminium_6061()
    if configuration == :polymer_direct
        return polymer, aluminium, nothing
    elseif configuration == :polymer_matched
        return polymer, aluminium, geometric_matching_material(polymer, aluminium)
    elseif configuration == :aluminium
        return aluminium, aluminium, nothing
    end
    throw(ArgumentError("unknown configuration: $configuration"))
end

function rounded_design_rows(candidate::RoundedLensDesignCandidate; slot_width_mm::Real=1.2)
    count = length(candidate.cells)
    pitch = first(candidate.cells).height_mm + Float64(slot_width_mm)
    centers = collect(-((count - 1) ÷ 2):((count - 1) ÷ 2)) .* pitch
    [
        (
            element_index=index,
            center_y_mm=centers[index],
            length_mm=cell.length_mm,
            minimum_gap_mm=cell.minimum_gap_mm,
            notch_count=cell.notch_count,
            notch_width_mm=cell.notch_width_mm,
            start_from_bottom=cell.start_from_bottom,
        )
        for (index, cell) in enumerate(candidate.cells)
    ]
end

"""
Optimize only peak amplitude at the prescribed point. The search is exact at
the carrier frequency; the best carrier candidates are re-ranked by the actual
finite-cycle pulse peak, so a narrow harmonic optimum cannot win automatically.
"""
function optimize_configuration(
    configuration::Symbol;
    config::MaterialLensConfig=MaterialLensConfig(),
    lengths_mm::AbstractVector{<:Real}=[12.0, 17.0, 22.0, 28.0, 34.0],
    gaps_mm::AbstractVector{<:Real}=[1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0],
    periods::AbstractVector{<:Integer}=[1, 2],
    matching_thicknesses_mm::AbstractVector{<:Real}=Float64[],
    shortlist_count::Integer=24,
)
    validate(config)
    lens_material, output_material, matching_material = configuration_materials(configuration)
    thicknesses = if configuration == :polymer_matched
        isempty(matching_thicknesses_mm) ? collect(2.0:0.5:6.0) : Float64.(matching_thicknesses_mm)
    else
        [0.0]
    end

    best = nothing
    for length_mm in lengths_mm, thickness_mm in thicknesses
        cell_options = SinusoidalCell[]
        for period_count in periods, gap_mm in gaps_mm
            push!(cell_options, SinusoidalCell(
                length_mm=Float64(length_mm),
                height_mm=config.element_height_mm,
                minimum_gap_mm=Float64(gap_mm),
                periods=Int(period_count),
            ))
        end
        shortlist = target_frequency_shortlist(
            cell_options,
            lens_material,
            output_material,
            config;
            matching_material,
            matching_thickness_mm=thickness_mm,
            keep=shortlist_count,
        )
        for (harmonic_amplitude, indices) in shortlist
            group_cells = cell_options[indices]
            cells = symmetric_cells(group_cells, config)
            response = pulse_response(
                cells,
                lens_material,
                output_material,
                config;
                matching_material,
                matching_thickness_mm=thickness_mm,
            )
            if isnothing(best) || response.peak_amplitude > best.pulse_peak_amplitude
                best = LensDesignCandidate(
                    configuration,
                    lens_material,
                    output_material,
                    matching_material,
                    Float64(thickness_mm),
                    cells,
                    harmonic_amplitude,
                    response.peak_amplitude,
                    response.peak_time_s,
                    NaN,
                )
            end
        end
    end
    isnothing(best) && error("empty design search")
    uniform_cell = SinusoidalCell(
        length_mm=first(best.cells).length_mm,
        height_mm=config.element_height_mm,
        minimum_gap_mm=config.element_height_mm,
        periods=1,
    )
    uniform_response = pulse_response(
        fill(uniform_cell, config.element_count),
        lens_material,
        output_material,
        config;
        matching_material,
        matching_thickness_mm=best.matching_thickness_mm,
    )
    LensDesignCandidate(
        best.configuration,
        best.lens_material,
        best.output_material,
        best.matching_material,
        best.matching_thickness_mm,
        best.cells,
        best.harmonic_focus_amplitude,
        best.pulse_peak_amplitude,
        best.pulse_peak_time_s,
        uniform_response.peak_amplitude,
    )
end

function design_rows(candidate::LensDesignCandidate; slot_width_mm::Real=1.2)
    count = length(candidate.cells)
    pitch = first(candidate.cells).height_mm + Float64(slot_width_mm)
    centers = collect(-((count - 1) ÷ 2):((count - 1) ÷ 2)) .* pitch
    [
        (
            element_index=index,
            center_y_mm=centers[index],
            length_mm=cell.length_mm,
            minimum_gap_mm=cell.minimum_gap_mm,
            amplitude_mm=cell_amplitude_mm(cell),
            periods=cell.periods,
        )
        for (index, cell) in enumerate(candidate.cells)
    ]
end

end
