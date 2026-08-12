module WorkingBandBoundaryStudy

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PILOT_ROOT = joinpath(PROJECT_ROOT, "tmp", "pilot_220khz")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "working_band_boundary_study")
const BAND_SCAN_FREQUENCIES_HZ = collect(120.0e3:10.0e3:160.0e3)
const FINAL_TIME_S = 120.0e-6
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

if !isdefined(@__MODULE__, :MetamaterialProfiles)
    include(joinpath(@__DIR__, "profiles.jl"))
end
using .MetamaterialProfiles

if !isnothing(REQUESTED_STAGE) &&
   (startswith(REQUESTED_STAGE, "band-one-") ||
    startswith(REQUESTED_STAGE, "reflect-one-"))
    include(joinpath(@__DIR__, "step2_solver.jl"))
elseif REQUESTED_STAGE == "reflect-solve"
    using JLD2
elseif REQUESTED_STAGE in ("band-analyze", "analyze")
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
    include(joinpath(@__DIR__, "spectral_analysis.jl"))
    using .SpectralAnalysis
    include(joinpath(@__DIR__, "plot_wavelet_analysis.jl"))
end

const TOPOLOGIES = [
    (
        id="rectangular",
        label="rectangular",
        profile=SinusoidalProfile(0.0; periods=2),
    ),
    (
        id="sinusoidal",
        label="sinusoidal",
        profile=SinusoidalProfile(2.5; periods=2),
    ),
    (
        id="exponential_k0p5",
        label="soft exponential κ=0.5",
        profile=ExponentialProfile(2.5; periods=2, sharpness=0.5),
    ),
    (
        id="exponential_k0p5_n1",
        label="stretched exponential κ=0.5, N=1",
        profile=ExponentialProfile(2.5; periods=1, sharpness=0.5),
    ),
    (
        id="exponential_k1",
        label="exponential κ=1",
        profile=ExponentialProfile(2.5; periods=2, sharpness=1.0),
    ),
    (
        id="exponential_k3",
        label="exponential κ=3",
        profile=ExponentialProfile(2.5; periods=2, sharpness=3.0),
    ),
    (
        id="power_p2",
        label="power p=2",
        profile=PowerProfile(2.5; periods=2, power=2.0),
    ),
    (
        id="power_p4",
        label="power p=4",
        profile=PowerProfile(2.5; periods=2, power=4.0),
    ),
    (
        id="rounded_notch",
        label="rounded U-notch",
        profile=RoundedNotchProfile(
            2.5;
            notch_count=4,
            notch_width_mm=1.4,
            end_margin_mm=1.2,
        ),
    ),
    (
        id="legacy",
        label="legacy",
        profile=LegacySinusoidalProfile(2.5),
    ),
]

frequency_khz(frequency_hz::Real) = Float64(frequency_hz) / 1.0e3
frequency_label(frequency_hz::Real) = string(frequency_khz(frequency_hz))

function absorbing_signal_path(topology, frequency_hz::Real)
    joinpath(
        PILOT_ROOT,
        "signals",
        "data_$(profile_slug(topology.profile))_F_$(frequency_label(frequency_hz)).jld2",
    )
end

function reflecting_signal_path(topology, frequency_hz::Real)
    joinpath(
        OUTPUT_ROOT,
        "signals",
        "data_$(profile_slug(topology.profile))_F_$(frequency_label(frequency_hz))_BC_free_reflecting.jld2",
    )
end

function study_simulation(frequency_hz::Real, boundary_condition::Symbol)
    SimulationConfig(
        frequencies_hz=[Float64(frequency_hz)],
        final_time_s=FINAL_TIME_S,
        samples_per_period=30,
        pulse_cycles=4.0,
        save_vtk=false,
        skip_existing=true,
        right_boundary_condition=boundary_condition,
        model_dir=joinpath(PILOT_ROOT, "models"),
        vtk_dir=joinpath(OUTPUT_ROOT, "vtk"),
        signal_dir=boundary_condition == :absorbing ?
                   joinpath(PILOT_ROOT, "signals") :
                   joinpath(OUTPUT_ROOT, "signals"),
    )
end

function run_absorbing_one(topology_index::Integer, frequency_index::Integer)
    topology = TOPOLOGIES[topology_index]
    frequency_hz = BAND_SCAN_FREQUENCIES_HZ[frequency_index]
    run_acoustic_simulation(
        topology.profile,
        frequency_hz;
        simulation=study_simulation(frequency_hz, :absorbing),
    )
end

function selection_path()
    joinpath(OUTPUT_ROOT, "working_band_selection.jld2")
end

function selected_carrier_hz()
    isfile(selection_path()) || error(
        "working-band selection is missing; run --stage=band-analyze first",
    )
    Float64(JLD2.load(selection_path())["selected_carrier_hz"])
end

function run_reflecting_one(topology_index::Integer)
    topology = TOPOLOGIES[topology_index]
    frequency_hz = selected_carrier_hz()
    run_acoustic_simulation(
        topology.profile,
        frequency_hz;
        simulation=study_simulation(frequency_hz, :free_reflecting),
    )
end

function process_count(job_count::Integer)
    requested = parse(Int, get(ENV, "METAMATERIALS_BOUNDARY_JOBS", "5"))
    requested > 0 || error("METAMATERIALS_BOUNDARY_JOBS must be positive")
    min(requested, job_count, Sys.CPU_THREADS)
end

function child_command(stage::AbstractString)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_parallel_stages(stages, description)
    isempty(stages) && return println("[=] No pending $description jobs")
    jobs = process_count(length(stages))
    println("=== $description: $(length(stages)) jobs, $jobs processes ===")
    semaphore = Base.Semaphore(jobs)
    @sync for stage in stages
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command(stage))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function run_band_solver()
    pending = [
        "band-one-$(topology_index)-$(frequency_index)"
        for topology_index in eachindex(TOPOLOGIES)
        for frequency_index in eachindex(BAND_SCAN_FREQUENCIES_HZ)
        if !isfile(absorbing_signal_path(
            TOPOLOGIES[topology_index],
            BAND_SCAN_FREQUENCIES_HZ[frequency_index],
        ))
    ]
    run_parallel_stages(pending, "absorbing working-band scan")
end

function run_reflecting_solver()
    frequency_hz = selected_carrier_hz()
    pending = [
        "reflect-one-$(topology_index)"
        for topology_index in eachindex(TOPOLOGIES)
        if !isfile(reflecting_signal_path(TOPOLOGIES[topology_index], frequency_hz))
    ]
    run_parallel_stages(pending, "free-boundary topology scan")
end

function load_boundary_signal(path::AbstractString)
    JLD2.jldopen(path, "r") do file
        (
            time_s=Float64.(file["time_s"]),
            drive_mpa=Float64.(file["source_drive_mpa"]),
            right_pressure_mpa=.-Float64.(file["right_normal_traction_mpa"]),
            right_velocity_m_per_s=Float64.(file["right_normal_velocity_m_per_s"]),
        )
    end
end

function reliable_calibrated_transfer(sample, reference)
    spectrum_config = SpectrumConfig(input_floor_relative=1.0e-2)
    sample_transfer = analyze_transfer(
        sample.time_s,
        sample.drive_mpa,
        sample.right_pressure_mpa;
        config=spectrum_config,
    )
    reference_transfer = analyze_transfer(
        reference.time_s,
        reference.drive_mpa,
        reference.right_pressure_mpa;
        config=spectrum_config,
    )
    calibrated = relative_transfer(sample_transfer, reference_transfer)
    reference_floor = 1.0e-2 * maximum(reference_transfer.amplitude[reference_transfer.valid])
    reliable = calibrated.valid .&
               reference_transfer.valid .&
               (reference_transfer.amplitude .>= reference_floor)
    calibrated, reliable
end

function carrier_point(sample, reference, requested_frequency_hz::Real)
    calibrated, reliable = reliable_calibrated_transfer(sample, reference)
    index = argmin(abs.(calibrated.frequency_hz .- requested_frequency_hz))
    (
        sampled_frequency_hz=calibrated.frequency_hz[index],
        amplitude=calibrated.amplitude[index],
        magnitude_squared=calibrated.magnitude_squared[index],
        phase_rad=angle(calibrated.transfer[index]),
        group_delay_s=calibrated.group_delay_s[index],
        reliable=reliable[index],
    )
end

function unwrap_sampled_phases(rows)
    result = copy(rows)
    for topology_id in unique(getproperty.(rows, :topology))
        indices = sort(
            findall(row -> row.topology == topology_id, rows);
            by=index -> rows[index].requested_frequency_hz,
        )
        previous_raw = rows[first(indices)].relative_phase_rad
        unwrapped = previous_raw
        result[first(indices)] = merge(rows[first(indices)], (relative_phase_rad=unwrapped,))
        for index in Iterators.drop(indices, 1)
            raw = rows[index].relative_phase_rad
            unwrapped += mod(raw - previous_raw + pi, 2pi) - pi
            result[index] = merge(rows[index], (relative_phase_rad=unwrapped,))
            previous_raw = raw
        end
    end
    result
end

function band_scan_rows()
    rows = NamedTuple[]
    reference_topology = first(TOPOLOGIES)
    for frequency_hz in BAND_SCAN_FREQUENCIES_HZ
        reference = load_boundary_signal(absorbing_signal_path(reference_topology, frequency_hz))
        for topology in TOPOLOGIES
            sample = load_boundary_signal(absorbing_signal_path(topology, frequency_hz))
            point = carrier_point(sample, reference, frequency_hz)
            push!(rows, (
                topology=topology.id,
                label=topology.label,
                requested_frequency_hz=frequency_hz,
                sampled_frequency_hz=point.sampled_frequency_hz,
                amplitude=point.amplitude,
                magnitude_squared=point.magnitude_squared,
                relative_phase_rad=point.phase_rad,
                excess_group_delay_s=point.group_delay_s,
                excess_envelope_delay_s=envelope_delay(
                    sample.time_s,
                    reference.right_pressure_mpa,
                    sample.right_pressure_mpa,
                ),
                reliable=point.reliable,
            ))
        end
    end
    unwrap_sampled_phases(rows)
end

function contiguous_band_indices(values::AbstractVector, peak_index::Integer, threshold::Real)
    left = peak_index
    right = peak_index
    while left > firstindex(values) && values[left - 1] >= threshold
        left -= 1
    end
    while right < lastindex(values) && values[right + 1] >= threshold
        right += 1
    end
    left:right
end

function threshold_crossing(x1, y1, x2, y2, threshold)
    y1 == y2 && return (x1 + x2) / 2
    x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
end

function interpolated_band_edges(frequencies, values, peak_index, threshold)
    sampled_band = contiguous_band_indices(values, peak_index, threshold)
    left = first(sampled_band)
    right = last(sampled_band)
    low = left == firstindex(values) ? frequencies[left] : threshold_crossing(
        frequencies[left - 1],
        values[left - 1],
        frequencies[left],
        values[left],
        threshold,
    )
    high = right == lastindex(values) ? frequencies[right] : threshold_crossing(
        frequencies[right],
        values[right],
        frequencies[right + 1],
        values[right + 1],
        threshold,
    )
    (
        sampled_low_hz=frequencies[left],
        sampled_high_hz=frequencies[right],
        estimated_low_hz=low,
        estimated_high_hz=high,
    )
end

function topology_band_summary(rows, topology)
    current = sort(
        filter(row -> row.topology == topology.id && row.reliable, rows);
        by=row -> row.requested_frequency_hz,
    )
    length(current) == length(BAND_SCAN_FREQUENCIES_HZ) ||
        error("incomplete reliable band scan for $(topology.id)")
    amplitudes = getproperty.(current, :amplitude)
    peak_index = argmax(amplitudes)
    threshold = amplitudes[peak_index] / sqrt(2.0)
    frequencies = getproperty.(current, :requested_frequency_hz)
    band = interpolated_band_edges(frequencies, amplitudes, peak_index, threshold)
    (
        topology=topology.id,
        peak_frequency_hz=current[peak_index].requested_frequency_hz,
        peak_amplitude=amplitudes[peak_index],
        minus3db_threshold=threshold,
        sampled_band_low_hz=band.sampled_low_hz,
        sampled_band_high_hz=band.sampled_high_hz,
        estimated_band_low_hz=band.estimated_low_hz,
        estimated_band_high_hz=band.estimated_high_hz,
        group_delay_at_peak_s=current[peak_index].excess_group_delay_s,
        envelope_delay_at_peak_s=current[peak_index].excess_envelope_delay_s,
    )
end

function common_band_summary(rows)
    shaped = Iterators.drop(TOPOLOGIES, 1)
    scores = Float64[]
    minimum_amplitudes = Float64[]
    for frequency_hz in BAND_SCAN_FREQUENCIES_HZ
        amplitudes = [
            only(filter(rows) do row
                row.topology == topology.id &&
                row.requested_frequency_hz == frequency_hz &&
                row.reliable
            end).amplitude
            for topology in shaped
        ]
        push!(scores, exp(sum(log.(max.(amplitudes, eps(Float64)))) / length(amplitudes)))
        push!(minimum_amplitudes, minimum(amplitudes))
    end
    peak_index = argmax(scores)
    threshold = scores[peak_index] / sqrt(2.0)
    band = interpolated_band_edges(
        BAND_SCAN_FREQUENCIES_HZ,
        scores,
        peak_index,
        threshold,
    )
    (
        selected_carrier_hz=BAND_SCAN_FREQUENCIES_HZ[peak_index],
        geometric_mean_amplitude=scores[peak_index],
        minimum_topology_amplitude=minimum_amplitudes[peak_index],
        minus3db_threshold=threshold,
        sampled_band_low_hz=band.sampled_low_hz,
        sampled_band_high_hz=band.sampled_high_hz,
        estimated_band_low_hz=band.estimated_low_hz,
        estimated_band_high_hz=band.estimated_high_hz,
        scores,
        minimum_amplitudes,
    )
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
end

function save_band_scan_plot(rows, common)
    colors = (:gray40, :royalblue, :deepskyblue3, :navy, :darkorange, :firebrick, :seagreen, :purple, :brown3, :goldenrod)
    figure = plot(
        xlabel="Carrier frequency, kHz",
        ylabel="Calibrated amplitude |Hcal|",
        title="Working-band scan for staggered wall topologies",
        xticks=frequency_khz.(BAND_SCAN_FREQUENCIES_HZ),
        gridalpha=0.25,
        legend=:outertopright,
        size=(1150, 650),
        margin=5Plots.mm,
    )
    for (index, topology) in enumerate(TOPOLOGIES)
        current = sort(filter(row -> row.topology == topology.id, rows); by=row -> row.requested_frequency_hz)
        plot!(
            figure,
            frequency_khz.(getproperty.(current, :requested_frequency_hz)),
            getproperty.(current, :amplitude);
            marker=:circle,
            linewidth=2.2,
            color=colors[index],
            label=topology.label,
        )
    end
    vline!(
        figure,
        [frequency_khz(common.selected_carrier_hz)];
        color=:black,
        linestyle=:dash,
        linewidth=1.5,
        label="selected carrier",
    )
    vspan!(
        figure,
        frequency_khz.([common.estimated_band_low_hz, common.estimated_band_high_hz]);
        color=:gray,
        alpha=0.10,
        label="estimated common -3 dB band",
    )
    path = joinpath(OUTPUT_ROOT, "working_band_scan.png")
    savefig(figure, path)
    println("[+] $path")
    path
end

function save_band_delay_plot(rows)
    colors = (:gray40, :royalblue, :deepskyblue3, :navy, :darkorange, :firebrick, :seagreen, :purple, :brown3, :goldenrod)
    phase_panel = plot(
        xlabel="Carrier frequency, kHz",
        ylabel="relative phase, rad",
        title="Phase relative to rectangular reference",
        gridalpha=0.25,
        legend=false,
    )
    group_panel = plot(
        xlabel="Carrier frequency, kHz",
        ylabel="excess group delay, μs",
        title="Local phase-slope delay (shown only for |Hcal| >= 0.5)",
        gridalpha=0.25,
        legend=false,
    )
    envelope_panel = plot(
        xlabel="Carrier frequency, kHz",
        ylabel="envelope delay, μs",
        title="Output-envelope correlation delay",
        gridalpha=0.25,
        legend=:outertopright,
    )
    for (index, topology) in enumerate(TOPOLOGIES)
        current = sort(filter(row -> row.topology == topology.id, rows); by=row -> row.requested_frequency_hz)
        frequencies = getproperty.(current, :requested_frequency_hz) ./ 1e3
        trusted = getproperty.(current, :reliable) .& (getproperty.(current, :amplitude) .>= 0.5)
        group_delay_us = getproperty.(current, :excess_group_delay_s) .* 1e6
        envelope_delay_us = getproperty.(current, :excess_envelope_delay_s) .* 1e6
        group_delay_us[.!trusted] .= NaN
        envelope_delay_us[.!trusted] .= NaN
        plot!(phase_panel, frequencies, getproperty.(current, :relative_phase_rad);
              color=colors[index], marker=:circle, linewidth=2, label=false)
        plot!(group_panel, frequencies, group_delay_us;
              color=colors[index], marker=:circle, linewidth=2, label=false)
        plot!(envelope_panel, frequencies, envelope_delay_us;
              color=colors[index], marker=:circle, linewidth=2, label=topology.label)
    end
    figure = plot(
        phase_panel,
        group_panel,
        envelope_panel;
        layout=(3, 1),
        size=(1200, 1150),
        margin=5Plots.mm,
        plot_title="Working-band phase and delay diagnostics",
    )
    path = joinpath(OUTPUT_ROOT, "working_band_phase_delay.png")
    savefig(figure, path)
    println("[+] $path")
    path
end

function run_band_analysis()
    mkpath(OUTPUT_ROOT)
    rows = band_scan_rows()
    summaries = [topology_band_summary(rows, topology) for topology in TOPOLOGIES]
    common = common_band_summary(rows)
    common_row = (
        topology="common_geometric_mean",
        peak_frequency_hz=common.selected_carrier_hz,
        peak_amplitude=common.geometric_mean_amplitude,
        minus3db_threshold=common.minus3db_threshold,
        sampled_band_low_hz=common.sampled_band_low_hz,
        sampled_band_high_hz=common.sampled_band_high_hz,
        estimated_band_low_hz=common.estimated_band_low_hz,
        estimated_band_high_hz=common.estimated_band_high_hz,
        group_delay_at_peak_s=NaN,
        envelope_delay_at_peak_s=NaN,
    )
    scan_csv = joinpath(OUTPUT_ROOT, "working_band_scan.csv")
    summary_csv = joinpath(OUTPUT_ROOT, "working_band_summary.csv")
    write_csv(scan_csv, rows)
    write_csv(summary_csv, vcat(summaries, [common_row]))
    JLD2.jldsave(
        selection_path();
        format_version=1,
        scan_frequencies_hz=BAND_SCAN_FREQUENCIES_HZ,
        selected_carrier_hz=common.selected_carrier_hz,
        common_sampled_band_low_hz=common.sampled_band_low_hz,
        common_sampled_band_high_hz=common.sampled_band_high_hz,
        common_estimated_band_low_hz=common.estimated_band_low_hz,
        common_estimated_band_high_hz=common.estimated_band_high_hz,
        common_geometric_mean_amplitude=common.geometric_mean_amplitude,
        common_minimum_topology_amplitude=common.minimum_topology_amplitude,
        rows,
        summaries,
    )
    save_band_scan_plot(rows, common)
    save_band_delay_plot(rows)
    println("[+] $scan_csv")
    println("[+] $summary_csv")
    println(
        "[+] selected carrier = $(frequency_khz(common.selected_carrier_hz)) kHz, ",
        "estimated common -3 dB band = ",
        "$(round(frequency_khz(common.estimated_band_low_hz); digits=1))--",
        "$(round(frequency_khz(common.estimated_band_high_hz); digits=1)) kHz",
    )
    common
end

function coi_peak_ridge_s(result, power)
    [
        begin
            masked_power = copy(view(power, frequency_index, :))
            masked_power[.!view(result.cone_of_influence_valid, frequency_index, :)] .= 0.0
            result.time_s[argmax(masked_power)]
        end
        for frequency_index in eachindex(result.frequency_hz)
    ]
end

function scalogram_panel(result, power_db, ridge_us; title, carrier_khz, colorbar=true)
    time_us = result.time_s .* 1.0e6
    frequency_khz_values = result.frequency_hz ./ 1.0e3
    scales_s = morlet_scale_s.(result.frequency_hz)
    coi_margin_s = sqrt(2.0) .* scales_s
    left_coi_us = (first(result.time_s) .+ coi_margin_s) .* 1.0e6
    right_coi_us = (last(result.time_s) .- coi_margin_s) .* 1.0e6
    panel = heatmap(
        frequency_khz_values,
        time_us,
        permutedims(power_db);
        xlabel="Frequency, kHz",
        ylabel="Time, μs",
        title,
        xlims=(first(frequency_khz_values), last(frequency_khz_values)),
        ylims=(first(time_us), last(time_us)),
        xticks=[50.0, 200.0, 400.0, 600.0, 800.0, 1000.0],
        color=:viridis,
        clims=(-50.0, 0.0),
        colorbar,
        colorbar_title=colorbar ? "dB" : "",
        framestyle=:box,
        grid=false,
        legend=:topright,
        left_margin=8Plots.mm,
        right_margin=5Plots.mm,
        bottom_margin=6Plots.mm,
    )
    plot!(panel, frequency_khz_values, ridge_us; color=:white, linewidth=1.8, label="peak ridge")
    plot!(panel, frequency_khz_values, left_coi_us; color=:white, linestyle=:dash, label="COI")
    plot!(panel, frequency_khz_values, right_coi_us; color=:white, linestyle=:dash, label="")
    vline!(panel, [carrier_khz]; color=:white, linestyle=:dot, label="carrier")
    panel
end

function wavelet_pair(topology, frequency_hz)
    absorbing_path = absorbing_signal_path(topology, frequency_hz)
    reflecting_path = reflecting_signal_path(topology, frequency_hz)
    absorbing = load_boundary_signal(absorbing_path)
    reflecting = load_boundary_signal(reflecting_path)
    absorbing.time_s == reflecting.time_s || error("time-grid mismatch for $(topology.id)")
    config = recommended_morlet_config(absorbing.time_s)
    absorbing_cwt = continuous_wavelet_transform(
        absorbing.time_s,
        absorbing.right_velocity_m_per_s;
        config,
    )
    reflecting_cwt = continuous_wavelet_transform(
        reflecting.time_s,
        reflecting.right_velocity_m_per_s;
        config,
    )
    absorbing_power = wavelet_power(absorbing_cwt)
    reflecting_power = wavelet_power(reflecting_cwt)
    normalization = max(maximum(absorbing_power), maximum(reflecting_power))
    absorbing_db = 10.0 .* log10.(absorbing_power ./ normalization .+ eps(Float64))
    reflecting_db = 10.0 .* log10.(reflecting_power ./ normalization .+ eps(Float64))
    absorbing_ridge_us = energy_masked_ridge_us(
        absorbing_cwt,
        coi_peak_ridge_s(absorbing_cwt, absorbing_power),
        1.0e-2,
    )
    reflecting_ridge_us = energy_masked_ridge_us(
        reflecting_cwt,
        coi_peak_ridge_s(reflecting_cwt, reflecting_power),
        1.0e-2,
    )
    dt_s = absorbing.time_s[2] - absorbing.time_s[1]
    late_mask = absorbing.time_s .>= 60.0e-6
    absorbing_energy = sum(abs2, absorbing.right_velocity_m_per_s) * dt_s
    reflecting_energy = sum(abs2, reflecting.right_velocity_m_per_s) * dt_s
    absorbing_late = sum(abs2, absorbing.right_velocity_m_per_s[late_mask]) /
                     sum(abs2, absorbing.right_velocity_m_per_s)
    reflecting_late = sum(abs2, reflecting.right_velocity_m_per_s[late_mask]) /
                      sum(abs2, reflecting.right_velocity_m_per_s)
    (;
        topology,
        absorbing_path,
        reflecting_path,
        absorbing,
        reflecting,
        config,
        absorbing_cwt,
        reflecting_cwt,
        absorbing_db,
        reflecting_db,
        absorbing_ridge_us,
        reflecting_ridge_us,
        absorbing_energy,
        reflecting_energy,
        absorbing_late,
        reflecting_late,
    )
end

function plot_all_boundary_scalograms(frequency_hz::Real)
    output_dir = joinpath(OUTPUT_ROOT, "scalograms")
    mkpath(output_dir)
    pairs = wavelet_pair.(TOPOLOGIES, Ref(frequency_hz))
    metric_rows = NamedTuple[]
    gallery_panels = Any[]
    carrier_khz = frequency_khz(frequency_hz)
    for pair in pairs
        topology = pair.topology
        absorbing_title = "$(topology.label): absorbing boundary"
        reflecting_title = "$(topology.label): free reflecting boundary"
        absorbing_panel = scalogram_panel(
            pair.absorbing_cwt,
            pair.absorbing_db,
            pair.absorbing_ridge_us;
            title=absorbing_title,
            carrier_khz,
        )
        reflecting_panel = scalogram_panel(
            pair.reflecting_cwt,
            pair.reflecting_db,
            pair.reflecting_ridge_us;
            title=reflecting_title,
            carrier_khz,
        )
        figure = plot(absorbing_panel, reflecting_panel; layout=(2, 1), size=(1200, 1050))
        figure_path = joinpath(
            output_dir,
            "scalograms_$(topology.id)_absorbing_vs_reflecting.png",
        )
        data_path = joinpath(output_dir, "wavelets_$(topology.id).jld2")
        savefig(figure, figure_path)
        JLD2.jldsave(
            data_path;
            carrier_frequency_hz=frequency_hz,
            topology=topology.id,
            observable="right_normal_velocity_m_per_s",
            absorbing_source=abspath(pair.absorbing_path),
            reflecting_source=abspath(pair.reflecting_path),
            morlet_config=pair.config,
            absorbing_cwt=pair.absorbing_cwt,
            reflecting_cwt=pair.reflecting_cwt,
        )
        push!(metric_rows, (
            topology=topology.id,
            carrier_frequency_hz=frequency_hz,
            absorbing_peak_velocity=maximum(abs, pair.absorbing.right_velocity_m_per_s),
            reflecting_peak_velocity=maximum(abs, pair.reflecting.right_velocity_m_per_s),
            reflecting_over_absorbing_peak=
                maximum(abs, pair.reflecting.right_velocity_m_per_s) /
                maximum(abs, pair.absorbing.right_velocity_m_per_s),
            absorbing_integrated_velocity_squared=pair.absorbing_energy,
            reflecting_integrated_velocity_squared=pair.reflecting_energy,
            reflecting_over_absorbing_integrated=pair.reflecting_energy / pair.absorbing_energy,
            absorbing_late_fraction_after_60us=pair.absorbing_late,
            reflecting_late_fraction_after_60us=pair.reflecting_late,
        ))
        push!(gallery_panels, scalogram_panel(
            pair.absorbing_cwt,
            pair.absorbing_db,
            pair.absorbing_ridge_us;
            title=absorbing_title,
            carrier_khz,
            colorbar=false,
        ))
        push!(gallery_panels, scalogram_panel(
            pair.reflecting_cwt,
            pair.reflecting_db,
            pair.reflecting_ridge_us;
            title=reflecting_title,
            carrier_khz,
            colorbar=false,
        ))
        println("[+] $figure_path")
    end
    gallery = plot(
        gallery_panels...;
        layout=(length(TOPOLOGIES), 2),
        size=(1800, 420 * length(TOPOLOGIES)),
        plot_title="Absorbing vs reflecting right boundary at $(carrier_khz) kHz",
    )
    gallery_path = joinpath(output_dir, "scalograms_all_topologies.png")
    metrics_path = joinpath(OUTPUT_ROOT, "boundary_metrics_all_topologies.csv")
    savefig(gallery, gallery_path)
    write_csv(metrics_path, metric_rows)
    println("[+] $gallery_path")
    println("[+] $metrics_path")
    gallery_path
end

function plot_transmission_at_selected_carrier(frequency_hz::Real)
    reference_topology = first(TOPOLOGIES)
    reference = load_boundary_signal(absorbing_signal_path(reference_topology, frequency_hz))
    colors = (:gray40, :royalblue, :deepskyblue3, :navy, :darkorange, :firebrick, :seagreen, :purple, :brown3, :goldenrod)
    figure = plot(
        xlabel="Frequency, kHz",
        ylabel="Calibrated amplitude |Hcal|",
        title="Transmission near the selected working band",
        xlims=(50.0, 500.0),
        gridalpha=0.25,
        legend=:outertopright,
        size=(1200, 650),
        margin=5Plots.mm,
    )
    rows = NamedTuple[]
    for (index, topology) in enumerate(TOPOLOGIES)
        sample = load_boundary_signal(absorbing_signal_path(topology, frequency_hz))
        calibrated, reliable = reliable_calibrated_transfer(sample, reference)
        amplitude = copy(calibrated.amplitude)
        amplitude[.!reliable] .= NaN
        plot!(
            figure,
            calibrated.frequency_hz ./ 1.0e3,
            amplitude;
            linewidth=2.2,
            color=colors[index],
            label=topology.label,
        )
        append!(rows, [
            (
                topology=topology.id,
                frequency_hz=calibrated.frequency_hz[i],
                amplitude=calibrated.amplitude[i],
                reliable=reliable[i],
            )
            for i in eachindex(calibrated.frequency_hz)
        ])
    end
    vline!(figure, [frequency_khz(frequency_hz)]; color=:black, linestyle=:dash, label="carrier")
    figure_path = joinpath(OUTPUT_ROOT, "transmission_selected_working_band.png")
    csv_path = joinpath(OUTPUT_ROOT, "transmission_selected_working_band.csv")
    savefig(figure, figure_path)
    write_csv(csv_path, rows)
    println("[+] $figure_path")
    println("[+] $csv_path")
    figure_path
end

function run_analysis()
    frequency_hz = selected_carrier_hz()
    plot_transmission_at_selected_carrier(frequency_hz)
    plot_all_boundary_scalograms(frequency_hz)
end

function run_child_stage(stage)
    if stage == "band-solve"
        run_band_solver()
    elseif startswith(stage, "band-one-")
        parts = split(stage, '-')
        length(parts) == 4 || error("invalid band stage: $stage")
        run_absorbing_one(parse(Int, parts[3]), parse(Int, parts[4]))
    elseif stage == "band-analyze"
        run_band_analysis()
    elseif stage == "reflect-solve"
        run_reflecting_solver()
    elseif startswith(stage, "reflect-one-")
        parts = split(stage, '-')
        length(parts) == 3 || error("invalid reflecting stage: $stage")
        run_reflecting_one(parse(Int, parts[3]))
    elseif stage == "analyze"
        run_analysis()
    else
        error("unknown stage: $stage")
    end
end

function run_pipeline()
    for stage in ("band-solve", "band-analyze", "reflect-solve", "analyze")
        println("\n=== Stage: $stage ===")
        run(child_command(stage))
    end
end

function main(args=ARGS)
    if isempty(args)
        run_pipeline()
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error(
            "usage: julia run_working_band_boundary_study.jl ",
            "[--stage=band-solve|band-analyze|reflect-solve|analyze]",
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
