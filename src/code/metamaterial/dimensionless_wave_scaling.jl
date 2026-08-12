module DimensionlessWaveScaling

export WaveScale,
       wavelength_mm,
       pressure_wavelength_mm,
       shear_wavelength_mm,
       normalized_frequency,
       length_over_pressure_wavelength,
       length_over_shear_wavelength,
       physical_pressure_length_mm,
       physical_shear_length_mm,
       mesh_resolution,
       delay_cycles,
       delay_phase_rad,
       impedance_ratio

"""Reference scales used internally by wave and mesh models."""
struct WaveScale
    reference_frequency_hz::Float64
    pressure_wave_speed_m_s::Float64
    shear_wave_speed_m_s::Float64
end

function WaveScale(material; reference_frequency_hz::Real)
    reference_frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    material.pressure_wave_speed_m_s > 0 ||
        throw(ArgumentError("P-wave speed must be positive"))
    material.shear_wave_speed_m_s > 0 ||
        throw(ArgumentError("S-wave speed must be positive"))
    WaveScale(
        Float64(reference_frequency_hz),
        Float64(material.pressure_wave_speed_m_s),
        Float64(material.shear_wave_speed_m_s),
    )
end

wavelength_mm(wave_speed_m_s::Real, frequency_hz::Real) =
    1.0e3 * Float64(wave_speed_m_s) / Float64(frequency_hz)

pressure_wavelength_mm(scale::WaveScale; frequency_hz=scale.reference_frequency_hz) =
    wavelength_mm(scale.pressure_wave_speed_m_s, frequency_hz)

shear_wavelength_mm(scale::WaveScale; frequency_hz=scale.reference_frequency_hz) =
    wavelength_mm(scale.shear_wave_speed_m_s, frequency_hz)

normalized_frequency(frequency_hz::Real, scale::WaveScale) =
    Float64(frequency_hz) / scale.reference_frequency_hz

length_over_pressure_wavelength(
    length_mm::Real,
    scale::WaveScale;
    frequency_hz=scale.reference_frequency_hz,
) = Float64(length_mm) / pressure_wavelength_mm(scale; frequency_hz)

length_over_shear_wavelength(
    length_mm::Real,
    scale::WaveScale;
    frequency_hz=scale.reference_frequency_hz,
) = Float64(length_mm) / shear_wavelength_mm(scale; frequency_hz)

physical_pressure_length_mm(
    normalized_length::Real,
    scale::WaveScale;
    frequency_hz=scale.reference_frequency_hz,
) = Float64(normalized_length) * pressure_wavelength_mm(scale; frequency_hz)

physical_shear_length_mm(
    normalized_length::Real,
    scale::WaveScale;
    frequency_hz=scale.reference_frequency_hz,
) = Float64(normalized_length) * shear_wavelength_mm(scale; frequency_hz)

function mesh_resolution(
    mesh_size_mm::Real,
    scale::WaveScale;
    frequency_hz=scale.reference_frequency_hz,
)
    h_over_lambda_s = length_over_shear_wavelength(
        mesh_size_mm,
        scale;
        frequency_hz,
    )
    (;
        h_over_lambda_s,
        k_s_h=2pi * h_over_lambda_s,
        shear_elements_per_wavelength=inv(h_over_lambda_s),
    )
end

delay_cycles(delay_s::Real, frequency_hz::Real) = Float64(delay_s) * Float64(frequency_hz)
delay_phase_rad(delay_s::Real, frequency_hz::Real) =
    2pi * delay_cycles(delay_s, frequency_hz)

impedance_ratio(
    left_density_kg_m3::Real,
    left_wave_speed_m_s::Real,
    right_density_kg_m3::Real,
    right_wave_speed_m_s::Real,
) = (Float64(right_density_kg_m3) * Float64(right_wave_speed_m_s)) /
    (Float64(left_density_kg_m3) * Float64(left_wave_speed_m_s))

end

