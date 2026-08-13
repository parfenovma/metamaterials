module NDTBenchmark

export NDTStandConfig,
       NDTDefectCase,
       validate,
       required_defects,
       extended_defects,
       all_defects,
       burst_duration_s,
       round_trip_time_s,
       nominal_dead_zone_depth_mm,
       pressure_wavelength_mm,
       shear_wavelength_mm,
       metric_columns,
       topology_result_directory

Base.@kwdef struct NDTStandConfig
    frequency_hz::Float64 = 242.0e3
    pulse_cycles::Float64 = 5.0
    density_kg_m3::Float64 = 2700.0
    pressure_wave_speed_m_s::Float64 = 6122.102437409232
    shear_wave_speed_m_s::Float64 = 3083.810277185563
    aperture_height_mm::Float64 = 71.6
    maximum_device_depth_mm::Float64 = 55.0
    specimen_length_mm::Float64 = 125.0
    specimen_height_mm::Float64 = 71.6
    nominal_focus_depth_mm::Float64 = 60.0
    mesh_size_mm::Float64 = 1.5
    defect_mesh_size_mm::Float64 = 0.35
    defect_transition_distance_mm::Float64 = 7.0
end

Base.@kwdef struct NDTDefectCase
    id::Symbol
    kind::Symbol = :circular_void
    center_x_mm::Float64
    center_y_mm::Float64 = 0.0
    radius_mm::Float64 = 0.0
    major_radius_mm::Float64 = 0.0
    minor_radius_mm::Float64 = 0.0
    angle_deg::Float64 = 0.0
end

function validate(config::NDTStandConfig)
    config.frequency_hz > 0 || throw(ArgumentError("frequency must be positive"))
    config.pulse_cycles > 0 || throw(ArgumentError("pulse cycle count must be positive"))
    config.density_kg_m3 > 0 || throw(ArgumentError("density must be positive"))
    config.pressure_wave_speed_m_s > 0 || throw(ArgumentError("P-wave speed must be positive"))
    config.shear_wave_speed_m_s > 0 || throw(ArgumentError("S-wave speed must be positive"))
    config.aperture_height_mm > 0 || throw(ArgumentError("aperture must be positive"))
    config.specimen_length_mm > 0 || throw(ArgumentError("specimen length must be positive"))
    config.specimen_height_mm >= config.aperture_height_mm ||
        throw(ArgumentError("specimen must contain the full aperture"))
    0 < config.nominal_focus_depth_mm < config.specimen_length_mm ||
        throw(ArgumentError("focus must lie inside the specimen"))
    config.mesh_size_mm > 0 || throw(ArgumentError("mesh size must be positive"))
    config.defect_mesh_size_mm > 0 || throw(ArgumentError("defect mesh size must be positive"))
    config.defect_transition_distance_mm > 0 ||
        throw(ArgumentError("defect transition distance must be positive"))
    config
end

function defect_extents_mm(defect::NDTDefectCase)
    defect.kind == :none && return (0.0, 0.0)
    defect.kind == :circular_void && return (defect.radius_mm, defect.radius_mm)
    defect.kind == :slit && return (
        max(defect.major_radius_mm, defect.minor_radius_mm),
        max(defect.major_radius_mm, defect.minor_radius_mm),
    )
    throw(ArgumentError("unsupported defect kind: $(defect.kind)"))
end

function validate(defect::NDTDefectCase, config::NDTStandConfig=NDTStandConfig())
    validate(config)
    defect.kind in (:none, :circular_void, :slit) ||
        throw(ArgumentError("unsupported defect kind: $(defect.kind)"))
    defect.kind == :none && return defect
    if defect.kind == :circular_void
        defect.radius_mm > 0 || throw(ArgumentError("void radius must be positive"))
    else
        defect.major_radius_mm > defect.minor_radius_mm > 0 ||
            throw(ArgumentError("slit semi-axes must satisfy major > minor > 0"))
    end
    extent_x, extent_y = defect_extents_mm(defect)
    defect.center_x_mm - extent_x > 0 || throw(ArgumentError("defect intersects input face"))
    defect.center_x_mm + extent_x < config.specimen_length_mm ||
        throw(ArgumentError("defect intersects downstream boundary"))
    abs(defect.center_y_mm) + extent_y < config.specimen_height_mm / 2 ||
        throw(ArgumentError("defect intersects a lateral boundary"))
    defect
end

function required_defects(config::NDTStandConfig=NDTStandConfig())
    validate(config)
    NDTDefectCase[
        NDTDefectCase(id=:none, kind=:none, center_x_mm=config.nominal_focus_depth_mm),
        NDTDefectCase(
            id=:void_r2_focus,
            center_x_mm=config.nominal_focus_depth_mm,
            radius_mm=2.0,
        ),
    ]
end

function extended_defects(config::NDTStandConfig=NDTStandConfig())
    validate(config)
    focus = config.nominal_focus_depth_mm
    NDTDefectCase[
        NDTDefectCase(id=:void_r1_focus, center_x_mm=focus, radius_mm=1.0),
        NDTDefectCase(id=:void_r4_focus, center_x_mm=focus, radius_mm=4.0),
        NDTDefectCase(id=:void_r2_near, center_x_mm=35.0, radius_mm=2.0),
        NDTDefectCase(id=:void_r2_deep, center_x_mm=85.0, radius_mm=2.0),
        NDTDefectCase(id=:void_r2_offset, center_x_mm=focus, center_y_mm=10.0, radius_mm=2.0),
        [
            NDTDefectCase(
                id=Symbol("slit_$(angle)_deg"),
                kind=:slit,
                center_x_mm=focus,
                major_radius_mm=4.0,
                minor_radius_mm=0.20,
                angle_deg=Float64(angle),
            )
            for angle in (0, 15, 30)
        ]...,
    ]
end

all_defects(config::NDTStandConfig=NDTStandConfig()) =
    vcat(required_defects(config), extended_defects(config))

burst_duration_s(config::NDTStandConfig=NDTStandConfig()) =
    config.pulse_cycles / config.frequency_hz

round_trip_time_s(depth_mm::Real, config::NDTStandConfig=NDTStandConfig()) =
    2 * Float64(depth_mm) * 1.0e-3 / config.pressure_wave_speed_m_s

nominal_dead_zone_depth_mm(config::NDTStandConfig=NDTStandConfig()) =
    0.5 * config.pressure_wave_speed_m_s * burst_duration_s(config) * 1.0e3

pressure_wavelength_mm(config::NDTStandConfig=NDTStandConfig()) =
    config.pressure_wave_speed_m_s / config.frequency_hz * 1.0e3

shear_wavelength_mm(config::NDTStandConfig=NDTStandConfig()) =
    config.shear_wave_speed_m_s / config.frequency_hz * 1.0e3

metric_columns() = (
    :topology,
    :case_id,
    :reference,
    :frequency_hz,
    :input_energy_j_per_m,
    :echo_mf_gain_bare,
    :echo_energy_gain_bare,
    :echo_to_ringing_db,
    :equal_noise_snr_proxy_db,
    :dead_zone_depth_mm,
    :focus_peak_gain_bare,
    :focus_energy_gain_bare,
    :fwhm_mm,
    :dof_minus3db_mm,
    :sidelobe_ratio,
    :pulse_broadening_ratio,
    :pulse_correlation,
    :postcursor_ratio,
    :weighted_transmission,
    :weighted_reflection,
    :weighted_mode_conversion,
    :minimum_feature_mm,
    :verdict,
)

function topology_result_directory(name::AbstractString)
    occursin(r"^[a-z0-9_]+$", name) ||
        throw(ArgumentError("topology name must use lowercase ASCII, digits and underscores"))
    joinpath("results", "topology_benchmark", name)
end

end
