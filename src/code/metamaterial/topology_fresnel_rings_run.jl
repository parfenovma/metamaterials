module TopologyFresnelRingsRun

using Printf
using Statistics

include(joinpath(@__DIR__, "topology_fresnel_rings_model.jl"))
using .TopologyFresnelRingsModel

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RESULT_DIRECTORY = joinpath(PROJECT_ROOT, "results", "topology_benchmark", "fresnel_rings")
const RAW_DIRECTORY = joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "fresnel_rings")

function burst_weights(frequencies_hz; fc_hz=242.0e3, cycles=5.0)
    duration = cycles / fc_hz
    time = collect(range(0.0, duration; length=1201))
    dt = time[2] - time[1]
    signal = @. 0.5 * (1 - cos(2pi * time / duration)) * sin(2pi * fc_hz * time)
    weights = Float64[]
    for frequency in frequencies_hz
        spectrum = sum(signal .* exp.(-2pi * im * frequency .* time)) * dt
        push!(weights, abs2(spectrum))
    end
    weights ./ sum(weights)
end

function unwrap_phase(values)
    output = copy(Float64.(values))
    for index in 2:length(output)
        jump = output[index] - output[index - 1]
        while jump > pi
            output[index:end] .-= 2pi
            jump -= 2pi
        end
        while jump < -pi
            output[index:end] .+= 2pi
            jump += 2pi
        end
    end
    output
end

function weighted_group_delay_s(transmission, frequencies_hz, weights)
    phase = unwrap_phase(angle.(transmission))
    omega = 2pi .* frequencies_hz
    delay = similar(phase)
    delay[1] = -(phase[2] - phase[1]) / (omega[2] - omega[1])
    delay[end] = -(phase[end] - phase[end - 1]) / (omega[end] - omega[end - 1])
    for index in 2:(length(delay) - 1)
        delay[index] = -(phase[index + 1] - phase[index - 1]) /
                       (omega[index + 1] - omega[index - 1])
    end
    mean_delay = weighted_mean(delay, weights)
    ripple = sqrt(weighted_mean((delay .- mean_delay) .^ 2, weights))
    (mean_s=mean_delay, ripple_s=ripple)
end

function level_metrics(level, design, frequencies_hz, weights; loss_scale=1.0)
    transfer = [level_transfer(level, design, f; loss_scale) for f in frequencies_hz]
    complex_t = getproperty.(transfer, :transmission)
    group = weighted_group_delay_s(complex_t, frequencies_hz, weights)
    pulse = fresnel_pulse_metrics(f -> level_transfer(level, design, f; loss_scale).transmission)
    (
        weighted_transmission=weighted_mean(getproperty.(transfer, :power_transmission), weights),
        weighted_reflection=weighted_mean(getproperty.(transfer, :power_reflection), weights),
        group_delay_us=group.mean_s * 1e6,
        group_delay_ripple_us=group.ripple_s * 1e6,
        pulse_broadening_ratio=pulse.broadening_ratio,
        pulse_correlation=pulse.correlation,
        postcursor_ratio=pulse.postcursor_ratio,
    )
end

function screen_designs(frequencies_hz, weights)
    rows = NamedTuple[]
    for level_count in 2:4, core_depth in 6.5:0.5:18.0, matching_depth in 1.5:0.125:3.5
        config = FresnelRingConfig(
            level_count=level_count,
            core_depth_mm=core_depth,
            matching_layer_depth_mm=matching_depth,
        )
        design = try
            design_fresnel_rings(config)
        catch error
            error isa ArgumentError || rethrow()
            continue
        end
        metrics = [
            level_metrics(level, design, frequencies_hz, weights; loss_scale=0.0)
            for level in design.levels
        ]
        push!(rows, (
            level_count,
            core_depth_mm=core_depth,
            matching_depth_mm=matching_depth,
            worst_weighted_transmission=minimum(getproperty.(metrics, :weighted_transmission)),
            mean_weighted_transmission=mean(getproperty.(metrics, :weighted_transmission)),
            worst_pulse_correlation=minimum(getproperty.(metrics, :pulse_correlation)),
            worst_broadening=maximum(getproperty.(metrics, :pulse_broadening_ratio)),
            worst_postcursor=maximum(getproperty.(metrics, :postcursor_ratio)),
        ))
    end
    rows
end

function write_namedtuple_csv(path, rows)
    isempty(rows) && error("cannot write empty CSV")
    keys_order = collect(keys(first(rows)))
    open(path, "w") do io
        println(io, join(string.(keys_order), ','))
        for row in rows
            println(io, join((getproperty(row, key) for key in keys_order), ','))
        end
    end
    path
end

function profile_metrics(profile)
    amplitude = profile.amplitude
    positions = profile.position_mm
    peak = amplitude[1]
    half_index = findfirst(value -> value <= peak / 2, amplitude)
    fwhm = isnothing(half_index) ? NaN : 2positions[half_index]
    minimum_index = nothing
    for index in 2:(length(amplitude) - 1)
        if amplitude[index] <= amplitude[index - 1] && amplitude[index] < amplitude[index + 1]
            minimum_index = index
            break
        end
    end
    sidelobe = isnothing(minimum_index) || minimum_index == length(amplitude) ? NaN :
               maximum(@view amplitude[(minimum_index + 1):end]) / peak
    (fwhm_mm=fwhm, sidelobe_ratio=sidelobe)
end

function aperture_power_metrics(design, topology, frequencies_hz, weights; loss_scale=1.0)
    radius_max = design.config.stand.aperture_height_mm / 2
    radial_samples = 160
    dr = radius_max / radial_samples
    transmission = Float64[]
    for frequency in frequencies_hz
        local_power = 0.0
        for index in 1:radial_samples
            radius = (index - 0.5) * dr
            amplitude = radial_transfer(radius, design, frequency, topology; loss_scale)
            local_power += 2radius * dr / radius_max^2 * abs2(amplitude)
        end
        push!(transmission, local_power)
    end
    (weighted_transmission=weighted_mean(transmission, weights),)
end

function phase_coherence_at_focus(design, topology, frequency_hz; loss_scale=1.0)
    radius_max = design.config.stand.aperture_height_mm / 2
    radial_samples = 240
    dr = radius_max / radial_samples
    focus = design.config.stand.nominal_focus_depth_mm
    k = 2pi * frequency_hz / design.config.stand.pressure_wave_speed_m_s * 1e-3
    coherent = 0.0 + 0im
    incoherent = 0.0
    for index in 1:radial_samples
        radius = (index - 0.5) * dr
        distance = hypot(focus, radius)
        transfer = radial_transfer(radius, design, frequency_hz, topology; loss_scale)
        weight = 2radius * dr / radius_max^2 *
                 abs((1 + im * k * distance) * focus / distance^3)
        contribution = weight * transfer * exp(-im * k * distance)
        coherent += contribution
        incoherent += abs(contribution)
    end
    abs(coherent) / incoherent
end

function aperture_metrics(design, topology; loss_scale=1.0)
    fc = design.config.stand.frequency_hz
    response = scalar_aperture_transfer(design, fc, topology; loss_scale)
    bare = scalar_aperture_transfer(design, fc, :bare)
    pulse = fresnel_pulse_metrics(
        f -> scalar_aperture_transfer(
            design,
            f,
            topology;
            radial_samples=80,
            loss_scale,
        ),
    )
    profile = lateral_profile(design, topology; loss_scale, maximum_mm=55.0, step_mm=0.5)
    spatial = profile_metrics(profile)
    (
        topology=String(topology),
        focus_peak_gain_bare=abs(response) / abs(bare),
        focus_energy_gain_bare=abs2(response) / abs2(bare),
        fwhm_mm=spatial.fwhm_mm,
        sidelobe_ratio=spatial.sidelobe_ratio,
        pulse_broadening_ratio=pulse.broadening_ratio,
        pulse_correlation=pulse.correlation,
        postcursor_ratio=pulse.postcursor_ratio,
    )
end

function common_metric_row(design, aperture, power; verdict="stop")
    minimum_axial_feature = minimum(
        vcat(
            [design.matching_layer_depth_mm],
            getproperty.(design.levels, :polyurethane_depth_mm),
            getproperty.(design.levels, :peek_depth_mm),
        )[vcat(
            [design.matching_layer_depth_mm],
            getproperty.(design.levels, :polyurethane_depth_mm),
            getproperty.(design.levels, :peek_depth_mm),
        ) .> 1.0e-9]
    )
    minimum_radial_feature = minimum(
        level.outer_radius_mm - level.inner_radius_mm for level in design.levels
    )
    (
        topology="fresnel_rings",
        case_id="none",
        reference=aperture.topology,
        frequency_hz=design.config.stand.frequency_hz,
        input_energy_j_per_m=NaN,
        echo_mf_gain_bare=NaN,
        echo_energy_gain_bare=NaN,
        echo_to_ringing_db=NaN,
        equal_noise_snr_proxy_db=NaN,
        dead_zone_depth_mm=NaN,
        focus_peak_gain_bare=aperture.focus_peak_gain_bare,
        focus_energy_gain_bare=aperture.focus_energy_gain_bare,
        fwhm_mm=aperture.fwhm_mm,
        dof_minus3db_mm=NaN,
        sidelobe_ratio=aperture.sidelobe_ratio,
        pulse_broadening_ratio=aperture.pulse_broadening_ratio,
        pulse_correlation=aperture.pulse_correlation,
        postcursor_ratio=aperture.postcursor_ratio,
        weighted_transmission=power.weighted_transmission,
        weighted_reflection=NaN,
        weighted_mode_conversion=NaN,
        minimum_feature_mm=min(minimum_axial_feature, minimum_radial_feature),
        verdict,
    )
end

function write_g0(path, design)
    open(path, "w") do io
        println(io, "quantity,value,unit")
        println(io, "aperture_radius,$(design.config.stand.aperture_height_mm / 2),mm")
        println(io, "focus_depth,$(design.config.stand.nominal_focus_depth_mm),mm")
        wavelength_mm = design.config.stand.pressure_wave_speed_m_s /
                        design.config.stand.frequency_hz * 1.0e3
        println(io, "wavelength,$wavelength_mm,mm")
        println(io, "edge_center_path_spread,$(design.maximum_delay_s * design.config.stand.pressure_wave_speed_m_s * 1e3),mm")
        println(io, "maximum_required_delay,$(design.maximum_delay_s * 1e6),us")
        println(io, "carrier_period,$(design.carrier_period_s * 1e6),us")
        println(io, "carrier_phase_span,$(design.carrier_phase_span_rad),rad")
        println(io, "full_2pi_wrap_count,$(design.full_wrap_count),count")
        println(io, "first_halfwave_radius,$(design.first_halfwave_radius_mm),mm")
        println(io, "total_device_depth,$(design.total_device_depth_mm),mm")
        for level in design.levels
            println(io, "ring_$(level.id)_inner_radius,$(level.inner_radius_mm),mm")
            println(io, "ring_$(level.id)_outer_radius,$(level.outer_radius_mm),mm")
            println(io, "ring_$(level.id)_polyurethane_depth,$(level.polyurethane_depth_mm),mm")
        end
    end
end

function main()
    mkpath(RESULT_DIRECTORY)
    mkpath(RAW_DIRECTORY)
    frequencies = collect(range(162.2e3, 321.7e3; step=500.0))
    weights = burst_weights(frequencies)

    println("[G0/G1 screen] 2--4 levels, core and matching depth")
    screen = screen_designs(frequencies, weights)
    write_namedtuple_csv(joinpath(RAW_DIRECTORY, "manufacturing_screen.csv"), screen)
    best_by_count = map(2:4) do count
        subset = filter(row -> row.level_count == count, screen)
        subset[argmax(getproperty.(subset, :worst_weighted_transmission))]
    end
    write_namedtuple_csv(joinpath(RESULT_DIRECTORY, "screen_best_by_level_count.csv"), best_by_count)

    # Three levels retain a genuine intermediate delay while the lossless
    # screen passes T>=0.8; the four-level cassette fails that gate even before
    # material attenuation, while two levels leave a larger phase error.
    best = only(filter(row -> row.level_count == 3, best_by_count))
    config = FresnelRingConfig(
        level_count=best.level_count,
        core_depth_mm=best.core_depth_mm,
        matching_layer_depth_mm=best.matching_depth_mm,
    )
    design = design_fresnel_rings(config)
    write_g0(joinpath(RESULT_DIRECTORY, "g0_design.csv"), design)

    println("[G1] broadband local transfer under attenuation bracket")
    g1_rows = NamedTuple[]
    nominal_metrics = NamedTuple[]
    for (bracket, loss_scale) in (("lossless", 0.0), ("optimistic", 0.25),
                                  ("nominal", 1.0), ("pessimistic", 3.0))
        for level in design.levels
            metric = level_metrics(level, design, frequencies, weights; loss_scale)
            push!(g1_rows, merge((bracket=bracket, level=level.id), metric))
            bracket == "nominal" && push!(nominal_metrics, metric)
        end
    end
    write_namedtuple_csv(joinpath(RESULT_DIRECTORY, "g1_level_metrics.csv"), g1_rows)

    println("[G2] equal-input scalar aperture controls")
    apertures = [
        aperture_metrics(design, topology; loss_scale=1.0)
        for topology in (:bare, :straight_footprint, :unwrapped_continuous, :fresnel_quantized)
    ]
    write_namedtuple_csv(joinpath(RESULT_DIRECTORY, "g2_aperture_metrics.csv"), apertures)

    chromatic_rows = NamedTuple[]
    for topology in (:bare, :unwrapped_continuous, :fresnel_quantized),
        frequency in (162.2e3, 193.8e3, 242.0e3, 290.1e3, 321.7e3)
        response = scalar_aperture_transfer(design, frequency, topology; loss_scale=1.0)
        bare_response = scalar_aperture_transfer(design, frequency, :bare)
        spatial = profile_metrics(lateral_profile(
            design,
            topology;
            frequency_hz=frequency,
            loss_scale=1.0,
            maximum_mm=55.0,
            step_mm=1.0,
        ))
        push!(chromatic_rows, (
            topology=String(topology),
            frequency_hz=frequency,
            nominal_focus_gain_bare=abs(response) / abs(bare_response),
            phase_coherence=phase_coherence_at_focus(design, topology, frequency; loss_scale=1.0),
            fwhm_mm=spatial.fwhm_mm,
            sidelobe_ratio=spatial.sidelobe_ratio,
        ))
    end
    write_namedtuple_csv(joinpath(RESULT_DIRECTORY, "g2_chromatic_focus.csv"), chromatic_rows)

    powers = Dict(
        aperture.topology => aperture_power_metrics(
            design,
            Symbol(aperture.topology),
            frequencies,
            weights;
            loss_scale=1.0,
        )
        for aperture in apertures
    )
    common = [
        common_metric_row(design, aperture, powers[aperture.topology]; verdict="stop")
        for aperture in apertures
    ]
    write_namedtuple_csv(joinpath(RESULT_DIRECTORY, "common_metrics.csv"), common)

    println("[done] primary three-level design")
    @printf("  core %.3f mm, matching skins %.3f mm each\n",
            design.config.core_depth_mm, design.matching_layer_depth_mm)
    @printf("  phase span %.4f rad, wraps %d, first pi-zone %.3f mm\n",
            design.carrier_phase_span_rad, design.full_wrap_count,
            design.first_halfwave_radius_mm)
    for row in g1_rows
        row.bracket == "nominal" || continue
        @printf("  ring %d: Tbar %.3f, Bt %.3f, rho %.3f, post %.3f\n",
                row.level, row.weighted_transmission, row.pulse_broadening_ratio,
                row.pulse_correlation, row.postcursor_ratio)
    end
    for row in apertures
        @printf("  %-20s focus gain %.3f, FWHM %.2f mm, sidelobe %.3f\n",
                row.topology, row.focus_peak_gain_bare, row.fwhm_mm,
                row.sidelobe_ratio)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
