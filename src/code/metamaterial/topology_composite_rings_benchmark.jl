module TopologyCompositeRingsBenchmark

ENV["GKSwstype"] = "100"

using FFTW
using Plots
using Statistics

include(joinpath(@__DIR__, "topology_composite_rings_model.jl"))
using .TopologyCompositeRingsModel
include(joinpath(@__DIR__, "ndt_benchmark.jl"))
using .NDTBenchmark

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_COMPOSITE_RINGS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "composite_rings"),
)

const CONFIG = CompositeRingsConfig(
    center_frequency_hz=242.0e3,
    pulse_cycles=5.0,
    aperture_diameter_mm=71.6,
    focus_depth_mm=60.0,
    core_thickness_mm=18.0,
    ring_count=7,
    outer_tungsten_volume_fraction=0.0992,
    face_tungsten_volume_fraction=0.2022,
    face_thickness_mm=1.65,
    aperture_samples=401,
)

const LOSS_CASES = LossCase[
    LossCase(label="lossless_bound", core_db_per_cm_at_f0=0.0, frequency_power=1.0, face_multiplier=1.25),
    LossCase(label="low_0p3_n1", core_db_per_cm_at_f0=0.3, frequency_power=1.0, face_multiplier=1.25),
    LossCase(label="optimistic_0p53_n2", core_db_per_cm_at_f0=0.53, frequency_power=2.0, face_multiplier=1.25),
    LossCase(label="threshold_0p65_n1", core_db_per_cm_at_f0=0.65, frequency_power=1.0, face_multiplier=1.25),
    LossCase(label="nominal_1p0_n1", core_db_per_cm_at_f0=1.0, frequency_power=1.0, face_multiplier=1.25),
    LossCase(label="linear_from_1mhz_2p2_n1", core_db_per_cm_at_f0=2.2, frequency_power=1.0, face_multiplier=1.25),
    LossCase(label="high_3p0_n1", core_db_per_cm_at_f0=3.0, frequency_power=1.0, face_multiplier=1.25),
]

const BOND_CASES = BondCase[
    BondCase(label="perfect_cocure", thickness_mm=0.0),
    BondCase(label="controlled_0p10mm", thickness_mm=0.10, loss_db_per_cm_at_f0=0.5),
    BondCase(label="thick_0p25mm", thickness_mm=0.25, loss_db_per_cm_at_f0=0.5),
]

function slug(value)
    replace(String(value), r"[^A-Za-z0-9_]+" => "_")
end

function write_design(zones)
    path = joinpath(OUTPUT_ROOT, "ring_design.csv")
    open(path, "w") do io
        println(io, "ring,inner_radius_mm,outer_radius_mm,representative_radius_mm,target_delay_us,target_speed_m_s,tungsten_volume_fraction,density_kg_m3,impedance_mrayl")
        for zone in zones
            p = zone.property
            println(io, join((
                zone.index, zone.inner_radius_mm, zone.outer_radius_mm,
                zone.representative_radius_mm, zone.target_delay_s * 1.0e6,
                zone.target_speed_m_s, p.tungsten_volume_fraction,
                p.density_kg_m3, p.impedance_mrayl,
            ), ','))
        end
    end
    path
end

function finite_weighted_reflection(result, role)
    reflection = role == :candidate ? result.candidate_reflection : result.straight_reflection
    weights = abs2.(result.pulse.spectrum[result.active])
    sum(weights .* abs2.(reflection[result.active])) / sum(weights)
end

peak(signal) = maximum(abs, signal)

function echo_window_ratio_db(result, role)
    echo = role == :candidate ? result.candidate_echo : result.straight_echo
    ringing = role == :candidate ? result.candidate_ringing : result.straight_ringing
    centre = argmax(abs.(echo))
    half_window = round(Int, CONFIG.pulse_cycles / CONFIG.center_frequency_hz /
                             (result.pulse.time_s[2] - result.pulse.time_s[1]))
    left = max(firstindex(echo), centre - half_window)
    right = min(lastindex(echo), centre + half_window)
    ringing_peak = maximum(abs, ringing[left:right])
    20log10(maximum(abs, echo[left:right]) / max(ringing_peak, eps(Float64)))
end

function write_sweep(rows)
    path = joinpath(OUTPUT_ROOT, "loss_bond_sweep.csv")
    open(path, "w") do io
        println(io, "loss_case,core_loss_db_cm_f0,frequency_power,bond_case,bond_thickness_mm,focus_peak_gain_bare,focus_peak_gain_straight,weighted_transmission,weighted_reflection,point_scatterer_echo_peak_gain_bare_proxy,point_scatterer_echo_energy_gain_bare_proxy,pulse_broadening_ratio,pulse_correlation,postcursor_ratio,echo_to_ringing_proxy_db,echo_to_ringing_delta_straight_db")
        for row in rows
            println(io, join((
                row.loss.label, row.loss.core_db_per_cm_at_f0, row.loss.frequency_power,
                row.bond.label, row.bond.thickness_mm,
                row.result.candidate_vs_bare.gain_peak,
                row.result.candidate_vs_straight.gain_peak,
                row.result.weighted_candidate_power,
                finite_weighted_reflection(row.result, :candidate),
                row.result.echo_peak_gain_bare_proxy,
                row.result.echo_energy_gain_bare_proxy,
                row.result.candidate_vs_bare.broadening_ratio,
                row.result.candidate_vs_bare.pulse_correlation,
                row.result.candidate_vs_bare.postcursor_ratio,
                row.echo_to_ringing_candidate_db,
                row.echo_to_ringing_candidate_db - row.echo_to_ringing_straight_db,
            ), ','))
        end
    end
    path
end

function carrier_profile(zones, loss, bond, role; direction=:transverse)
    coordinates = direction == :transverse ? collect(-40.0:0.25:40.0) : collect(25.0:0.25:95.0)
    amplitude = Float64[]
    for coordinate in coordinates
        response = if direction == :transverse
            aperture_transfer(
                CONFIG.center_frequency_hz, CONFIG, zones, loss, bond;
                role, observation_x_mm=CONFIG.focus_depth_mm, observation_y_mm=coordinate,
            )
        else
            aperture_transfer(
                CONFIG.center_frequency_hz, CONFIG, zones, loss, bond;
                role, observation_x_mm=coordinate, observation_y_mm=0.0,
            )
        end
        push!(amplitude, abs(response.field))
    end
    (; coordinates, amplitude)
end

function central_width(profile; threshold=inv(sqrt(2.0)))
    peak_index = argmax(profile.amplitude)
    level = threshold * profile.amplitude[peak_index]
    left = peak_index
    right = peak_index
    while left > 1 && profile.amplitude[left - 1] >= level
        left -= 1
    end
    while right < length(profile.amplitude) && profile.amplitude[right + 1] >= level
        right += 1
    end
    # A point just outside the -3 dB lobe is not a sidelobe. Exclude through
    # the first local minimum on each side of the peak.
    left_minimum = left
    while left_minimum > 2
        profile.amplitude[left_minimum - 1] > profile.amplitude[left_minimum] && break
        left_minimum -= 1
    end
    right_minimum = right
    while right_minimum < length(profile.amplitude) - 1
        profile.amplitude[right_minimum + 1] > profile.amplitude[right_minimum] && break
        right_minimum += 1
    end
    outside = vcat(
        profile.amplitude[1:max(0, left_minimum - 1)],
        profile.amplitude[min(length(profile.amplitude) + 1, right_minimum + 1):end],
    )
    sidelobe = isempty(outside) ? NaN : maximum(outside) / profile.amplitude[peak_index]
    (
        width_mm=profile.coordinates[right] - profile.coordinates[left],
        peak_coordinate_mm=profile.coordinates[peak_index],
        peak_amplitude=profile.amplitude[peak_index],
        sidelobe_ratio=sidelobe,
    )
end

function write_profiles(zones, loss, bond)
    transverse = Dict(role => carrier_profile(zones, loss, bond, role; direction=:transverse)
                      for role in (:candidate, :straight, :bare))
    axial = Dict(role => carrier_profile(zones, loss, bond, role; direction=:axial)
                 for role in (:candidate, :straight, :bare))
    transverse_path = joinpath(OUTPUT_ROOT, "carrier_transverse_profile.csv")
    open(transverse_path, "w") do io
        println(io, "y_mm,candidate,straight,bare")
        for index in eachindex(transverse[:candidate].coordinates)
            println(io, join((
                transverse[:candidate].coordinates[index],
                transverse[:candidate].amplitude[index],
                transverse[:straight].amplitude[index],
                transverse[:bare].amplitude[index],
            ), ','))
        end
    end
    axial_path = joinpath(OUTPUT_ROOT, "carrier_axial_profile.csv")
    open(axial_path, "w") do io
        println(io, "x_mm,candidate,straight,bare")
        for index in eachindex(axial[:candidate].coordinates)
            println(io, join((
                axial[:candidate].coordinates[index],
                axial[:candidate].amplitude[index],
                axial[:straight].amplitude[index],
                axial[:bare].amplitude[index],
            ), ','))
        end
    end
    metrics = Dict(role => (
        transverse=central_width(transverse[role]),
        axial=central_width(axial[role]),
    ) for role in keys(transverse))
    p1 = plot(
        transverse[:candidate].coordinates, transverse[:candidate].amplitude ./ maximum(transverse[:bare].amplitude);
        label="rings", lw=2.5, xlabel="y, mm", ylabel="|H| / max(|H_bare|)",
        title="242 kHz transverse profile",
    )
    plot!(p1, transverse[:straight].coordinates, transverse[:straight].amplitude ./ maximum(transverse[:bare].amplitude);
          label="straight", lw=2, ls=:dash)
    plot!(p1, transverse[:bare].coordinates, transverse[:bare].amplitude ./ maximum(transverse[:bare].amplitude);
          label="bare", lw=2, ls=:dot)
    p2 = plot(
        axial[:candidate].coordinates, axial[:candidate].amplitude ./ maximum(axial[:bare].amplitude);
        label="rings", lw=2.5, xlabel="x from interface, mm", ylabel="|H| / max(|H_bare|)",
        title="242 kHz axial profile",
    )
    plot!(p2, axial[:straight].coordinates, axial[:straight].amplitude ./ maximum(axial[:bare].amplitude);
          label="straight", lw=2, ls=:dash)
    plot!(p2, axial[:bare].coordinates, axial[:bare].amplitude ./ maximum(axial[:bare].amplitude);
          label="bare", lw=2, ls=:dot)
    plot_path = joinpath(OUTPUT_ROOT, "carrier_profiles.png")
    savefig(plot(p1, p2; layout=(1, 2), size=(1250, 470)), plot_path)
    (; transverse, axial, metrics, transverse_path, axial_path, plot_path)
end

function write_common_metrics(primary, profiles)
    path = joinpath(OUTPUT_ROOT, "common_metrics.csv")
    columns = metric_columns()
    open(path, "w") do io
        println(io, join(columns, ','))
        for role in (:candidate, :straight, :bare), case_id in (:none, :void_r2_focus)
            unavailable_defect = case_id == :void_r2_focus
            focus_gain = role == :candidate ? primary.candidate_vs_bare.gain_peak :
                         role == :straight ? peak(primary.straight_wave) / peak(primary.bare_wave) : 1.0
            weighted_t = role == :candidate ? primary.weighted_candidate_power :
                         role == :straight ? primary.weighted_straight_power : 1.0
            weighted_r = role == :bare ? 0.0 : finite_weighted_reflection(primary, role)
            profile = profiles.metrics[role]
            broadening = role == :candidate ? primary.candidate_vs_bare.broadening_ratio : NaN
            correlation = role == :candidate ? primary.candidate_vs_bare.pulse_correlation : NaN
            postcursor = role == :candidate ? primary.candidate_vs_bare.postcursor_ratio : NaN
            values = Dict{Symbol, Any}(
                :topology => role == :candidate ? "composite_rings" : String(role),
                :case_id => String(case_id),
                :reference => role == :candidate ? "bare_and_straight" : "control",
                :frequency_hz => CONFIG.center_frequency_hz,
                :input_energy_j_per_m => NaN,
                :echo_mf_gain_bare => NaN,
                :echo_energy_gain_bare => NaN,
                :echo_to_ringing_db => NaN,
                :equal_noise_snr_proxy_db => NaN,
                :dead_zone_depth_mm => NaN,
                :focus_peak_gain_bare => unavailable_defect ? NaN : focus_gain,
                :focus_energy_gain_bare => unavailable_defect ? NaN : focus_gain^2,
                :fwhm_mm => unavailable_defect ? NaN : profile.transverse.width_mm,
                # The scalar field has no isolated axial maximum at 60 mm;
                # reporting its global -3 dB width as DOF would be misleading.
                :dof_minus3db_mm => NaN,
                :sidelobe_ratio => unavailable_defect ? NaN : profile.transverse.sidelobe_ratio,
                :pulse_broadening_ratio => unavailable_defect ? NaN : broadening,
                :pulse_correlation => unavailable_defect ? NaN : correlation,
                :postcursor_ratio => unavailable_defect ? NaN : postcursor,
                :weighted_transmission => unavailable_defect ? NaN : weighted_t,
                :weighted_reflection => unavailable_defect ? NaN : weighted_r,
                :weighted_mode_conversion => NaN,
                :minimum_feature_mm => 0.10,
                :verdict => "conditional",
            )
            println(io, join((values[column] for column in columns), ','))
        end
    end
    path
end

function write_local_transfer(zones, loss, bond)
    frequencies = Float64[
        162_200, 178_160, 194_120, 210_080, 226_040, 242_000,
        257_940, 273_880, 289_820, 305_760, 321_700,
    ]
    responses = [
        [zone_transfer(zone.property, frequency, CONFIG, loss, bond)
         for frequency in frequencies]
        for zone in zones
    ]
    phases = Vector{Vector{Float64}}()
    for response in responses
        raw = angle.(getproperty.(response, :transmission))
        phase = Float64[first(raw)]
        for value in Iterators.drop(raw, 1)
            push!(phase, last(phase) + mod(value - last(phase) + pi, 2pi) - pi)
        end
        push!(phases, phase)
    end
    delays = [begin
        omega = 2pi .* frequencies
        delay = similar(phase)
        delay[1] = -(phase[2] - phase[1]) / (omega[2] - omega[1])
        for index in 2:(length(phase) - 1)
            delay[index] = -(phase[index + 1] - phase[index - 1]) /
                           (omega[index + 1] - omega[index - 1])
        end
        delay[end] = -(phase[end] - phase[end - 1]) / (omega[end] - omega[end - 1])
        delay
    end for phase in phases]
    path = joinpath(OUTPUT_ROOT, "local_transfer.csv")
    open(path, "w") do io
        println(io, "ring,frequency_hz,s21_real,s21_imag,s21_abs,s11_real,s11_imag,s11_abs,power_transmission,unwrapped_phase_rad,group_delay_us")
        for ring in eachindex(zones), index in eachindex(frequencies)
            response = responses[ring][index]
            println(io, join((
                ring, frequencies[index],
                real(response.transmission), imag(response.transmission), abs(response.transmission),
                real(response.reflection), imag(response.reflection), abs(response.reflection),
                response.power_transmission, phases[ring][index], delays[ring][index] * 1.0e6,
            ), ','))
        end
    end
    path
end

function make_summary_plot(rows)
    controlled = filter(row -> row.bond.label == "controlled_0p10mm", rows)
    losses = getproperty.(getproperty.(controlled, :loss), :core_db_per_cm_at_f0)
    focus_gain = getproperty.(getproperty.(getproperty.(controlled, :result), :candidate_vs_bare), :gain_peak)
    echo_gain = getproperty.(getproperty.(controlled, :result), :echo_peak_gain_bare_proxy)
    transmission = getproperty.(getproperty.(controlled, :result), :weighted_candidate_power)
    p = plot(losses, focus_gain; marker=:circle, lw=2.5, label="focus peak / bare",
             xlabel="core attenuation at 242 kHz, dB/cm", ylabel="gain or power fraction",
             title="Loss gate, 0.10 mm bond")
    plot!(p, losses, echo_gain; marker=:square, lw=2.5, label="point-scatterer echo proxy / bare")
    plot!(p, losses, transmission; marker=:diamond, lw=2.5, label="weighted one-way power T")
    hline!(p, [1.0]; color=:black, ls=:dash, label="break-even")
    path = joinpath(OUTPUT_ROOT, "loss_gate.png")
    savefig(p, path)
    path
end

function main()
    mkpath(OUTPUT_ROOT)
    zones = design_ring_zones(CONFIG)
    design_path = write_design(zones)
    rows = NamedTuple[]
    for loss in LOSS_CASES, bond in BOND_CASES
        println("[i] $(loss.label), $(bond.label)")
        result = evaluate_case(CONFIG, loss, bond; fft_length=32768)
        push!(rows, (
            loss,
            bond,
            result,
            echo_to_ringing_candidate_db=echo_window_ratio_db(result, :candidate),
            echo_to_ringing_straight_db=echo_window_ratio_db(result, :straight),
        ))
    end
    sweep_path = write_sweep(rows)
    primary_row = only(filter(
        row -> row.loss.label == "linear_from_1mhz_2p2_n1" &&
               row.bond.label == "controlled_0p10mm",
        rows,
    ))
    local_path = write_local_transfer(zones, primary_row.loss, primary_row.bond)
    profiles = write_profiles(zones, primary_row.loss, primary_row.bond)
    common_path = write_common_metrics(primary_row.result, profiles)
    plot_path = make_summary_plot(rows)
    println("[+] design: $design_path")
    println("[+] sweep: $sweep_path")
    println("[+] common metrics: $common_path")
    println("[+] local transfer: $local_path")
    println("[+] loss plot: $plot_path")
    println("[+] profile plot: $(profiles.plot_path)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
