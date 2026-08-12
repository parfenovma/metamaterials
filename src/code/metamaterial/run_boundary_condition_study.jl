const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PILOT_ROOT = joinpath(PROJECT_ROOT, "tmp", "pilot_220khz")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "boundary_condition_study_220khz")
const CARRIER_FREQUENCY_HZ = 220.0e3

const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

if REQUESTED_STAGE == "solve"
    include(joinpath(@__DIR__, "step2_solver.jl"))
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
    include(joinpath(@__DIR__, "spectral_analysis.jl"))
    using .SpectralAnalysis
    include(joinpath(@__DIR__, "plot_wavelet_analysis.jl"))
elseif REQUESTED_STAGE == "meshes"
    include(joinpath(@__DIR__, "plot_reference_meshes.jl"))
end

const TRANSMISSION_CASES = [
    (id="rectangular", label="rectangular", slug="sin_A_0.0_N_2"),
    (id="sinusoidal", label="sinusoidal", slug="sin_A_2.5_N_2_STG"),
    (id="exponential_k1", label="exponential κ=1", slug="exp_A_2.5_N_2_K_1.0_STG"),
    (id="exponential_k3", label="exponential κ=3", slug="exp_A_2.5_N_2_K_3.0_STG"),
    (id="power_p2", label="power p=2", slug="pow_A_2.5_N_2_P_2.0_STG"),
    (id="power_p4", label="power p=4", slug="pow_A_2.5_N_2_P_4.0_STG"),
    (id="legacy", label="legacy", slug="A_2.5"),
]

function absorbing_signal_path(slug)
    joinpath(PILOT_ROOT, "signals", "data_$(slug)_F_220.0.jld2")
end

function reflecting_signal_path()
    joinpath(
        OUTPUT_ROOT,
        "signals",
        "data_sin_A_2.5_N_2_STG_F_220.0_BC_free_reflecting.jld2",
    )
end

function run_reflecting_signal()
    simulation = SimulationConfig(
        frequencies_hz=[CARRIER_FREQUENCY_HZ],
        final_time_s=120.0e-6,
        samples_per_period=30,
        pulse_cycles=4.0,
        save_vtk=false,
        skip_existing=true,
        right_boundary_condition=:free_reflecting,
        model_dir=joinpath(PILOT_ROOT, "models"),
        vtk_dir=joinpath(OUTPUT_ROOT, "vtk"),
        signal_dir=joinpath(OUTPUT_ROOT, "signals"),
    )
    run_acoustic_simulation(
        SinusoidalProfile(2.5; periods=2),
        CARRIER_FREQUENCY_HZ;
        simulation,
    )
end

function load_boundary_signal(path)
    jldopen(path, "r") do file
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

function write_transmission_csv(path, rows)
    open(path, "w") do io
        println(io, "topology,frequency_hz,amplitude,magnitude_squared,reliable")
        for row in rows
            println(
                io,
                join((
                    row.topology,
                    row.frequency_hz,
                    row.amplitude,
                    row.magnitude_squared,
                    row.reliable,
                ), ','),
            )
        end
    end
end

function plot_transmission_coefficients()
    reference = load_boundary_signal(absorbing_signal_path("sin_A_0.0_N_2"))
    curves = NamedTuple[]
    rows = NamedTuple[]
    for case in TRANSMISSION_CASES
        sample = load_boundary_signal(absorbing_signal_path(case.slug))
        calibrated, reliable = reliable_calibrated_transfer(sample, reference)
        amplitude = copy(calibrated.amplitude)
        amplitude[.!reliable] .= NaN
        push!(curves, (
            id=case.id,
            label=case.label,
            frequency_khz=calibrated.frequency_hz ./ 1.0e3,
            amplitude,
        ))
        append!(rows, [
            (
                topology=case.id,
                frequency_hz=calibrated.frequency_hz[index],
                amplitude=calibrated.amplitude[index],
                magnitude_squared=calibrated.magnitude_squared[index],
                reliable=reliable[index],
            )
            for index in eachindex(calibrated.frequency_hz)
        ])
    end

    colors = (:gray40, :royalblue, :darkorange, :firebrick, :seagreen, :purple, :goldenrod)
    figure = plot(
        xlabel="Frequency, kHz",
        ylabel="Calibrated amplitude transmission |Hcal|",
        title="Transmission with absorbing right boundary: staggered walls",
        xlims=(50.0, 1000.0),
        xticks=[50.0, 200.0, 400.0, 600.0, 800.0, 1000.0],
        gridalpha=0.25,
        legend=:outertopright,
        size=(1250, 650),
        margin=5Plots.mm,
    )
    for (index, curve) in enumerate(curves)
        plot!(
            figure,
            curve.frequency_khz,
            curve.amplitude;
            linewidth=2.2,
            color=colors[index],
            label=curve.label,
        )
    end
    hline!(figure, [1.0]; color=:black, linestyle=:dash, linewidth=1.2, label="reference")
    mkpath(OUTPUT_ROOT)
    figure_path = joinpath(OUTPUT_ROOT, "transmission_absorbing_boundary.png")
    csv_path = joinpath(OUTPUT_ROOT, "transmission_absorbing_boundary.csv")
    savefig(figure, figure_path)
    write_transmission_csv(csv_path, rows)
    println("[+] $figure_path")
    println("[+] $csv_path")
    figure_path
end

function scalogram_panel(result, power_db, ridge_us; title, carrier_khz)
    time_us = result.time_s .* 1.0e6
    frequency_khz = result.frequency_hz ./ 1.0e3
    scales_s = morlet_scale_s.(result.frequency_hz)
    coi_margin_s = sqrt(2.0) .* scales_s
    left_coi_us = (first(result.time_s) .+ coi_margin_s) .* 1.0e6
    right_coi_us = (last(result.time_s) .- coi_margin_s) .* 1.0e6
    panel = heatmap(
        frequency_khz,
        time_us,
        permutedims(power_db);
        xlabel="Frequency, kHz",
        ylabel="Time, μs",
        title,
        xlims=(first(frequency_khz), last(frequency_khz)),
        ylims=(first(time_us), last(time_us)),
        xticks=[50.0, 200.0, 400.0, 600.0, 800.0, 1000.0],
        color=:viridis,
        clims=(-50.0, 0.0),
        colorbar_title="dB",
        framestyle=:box,
        grid=false,
        legend=:topright,
        left_margin=8Plots.mm,
        right_margin=5Plots.mm,
        bottom_margin=6Plots.mm,
    )
    plot!(panel, frequency_khz, ridge_us; color=:white, linewidth=1.8, label="peak ridge")
    plot!(panel, frequency_khz, left_coi_us; color=:white, linestyle=:dash, label="COI")
    plot!(panel, frequency_khz, right_coi_us; color=:white, linestyle=:dash, label="")
    vline!(panel, [carrier_khz]; color=:white, linestyle=:dot, label="carrier")
    panel
end

function coi_peak_ridge_s(result, power)
    [
        begin
            masked_power = copy(view(power, frequency_index, :))
            masked_power[.!view(
                result.cone_of_influence_valid,
                frequency_index,
                :,
            )] .= 0.0
            result.time_s[argmax(masked_power)]
        end
        for frequency_index in eachindex(result.frequency_hz)
    ]
end

function plot_boundary_scalograms()
    absorbing = load_boundary_signal(absorbing_signal_path("sin_A_2.5_N_2_STG"))
    reflecting = load_boundary_signal(reflecting_signal_path())
    absorbing.time_s == reflecting.time_s ||
        error("absorbing and reflecting signals use different time grids")
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
    dt_s = absorbing.time_s[2] - absorbing.time_s[1]
    absorbing_time_energy = sum(abs2, absorbing.right_velocity_m_per_s) * dt_s
    reflecting_time_energy = sum(abs2, reflecting.right_velocity_m_per_s) * dt_s
    late_mask = absorbing.time_s .>= 60.0e-6
    absorbing_late_fraction = sum(abs2, absorbing.right_velocity_m_per_s[late_mask]) /
                              sum(abs2, absorbing.right_velocity_m_per_s)
    reflecting_late_fraction = sum(abs2, reflecting.right_velocity_m_per_s[late_mask]) /
                               sum(abs2, reflecting.right_velocity_m_per_s)
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

    absorbing_panel = scalogram_panel(
        absorbing_cwt,
        absorbing_db,
        absorbing_ridge_us;
        title="Right-port velocity: absorbing boundary",
        carrier_khz=CARRIER_FREQUENCY_HZ / 1.0e3,
    )
    reflecting_panel = scalogram_panel(
        reflecting_cwt,
        reflecting_db,
        reflecting_ridge_us;
        title="Right-port velocity: free reflecting boundary",
        carrier_khz=CARRIER_FREQUENCY_HZ / 1.0e3,
    )
    figure = plot(
        absorbing_panel,
        reflecting_panel;
        layout=(2, 1),
        size=(1200, 1050),
    )
    mkpath(OUTPUT_ROOT)
    figure_path = joinpath(OUTPUT_ROOT, "scalograms_absorbing_vs_reflecting.png")
    data_path = joinpath(OUTPUT_ROOT, "wavelets_absorbing_vs_reflecting.jld2")
    metrics_path = joinpath(OUTPUT_ROOT, "boundary_condition_metrics.csv")
    savefig(figure, figure_path)
    open(metrics_path, "w") do io
        println(io, "boundary,peak_velocity_m_per_s,time_integrated_velocity_squared,late_energy_fraction_after_60us")
        println(io, join((
            "absorbing",
            maximum(abs, absorbing.right_velocity_m_per_s),
            absorbing_time_energy,
            absorbing_late_fraction,
        ), ','))
        println(io, join((
            "free_reflecting",
            maximum(abs, reflecting.right_velocity_m_per_s),
            reflecting_time_energy,
            reflecting_late_fraction,
        ), ','))
    end
    JLD2.jldsave(
        data_path;
        carrier_frequency_hz=CARRIER_FREQUENCY_HZ,
        observable="right_normal_velocity_m_per_s",
        morlet_config=config,
        absorbing_source=abspath(absorbing_signal_path("sin_A_2.5_N_2_STG")),
        reflecting_source=abspath(reflecting_signal_path()),
        absorbing_cwt,
        reflecting_cwt,
        absorbing_time_energy,
        reflecting_time_energy,
        absorbing_late_fraction,
        reflecting_late_fraction,
    )
    println("[+] $figure_path")
    println("[+] $data_path")
    println("[+] $metrics_path")
    figure_path
end

function run_analysis()
    plot_transmission_coefficients()
    plot_boundary_scalograms()
end

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("solve", "analyze", "meshes")
            run(child_command(stage))
        end
    elseif REQUESTED_STAGE == "solve"
        run_reflecting_signal()
    elseif REQUESTED_STAGE == "analyze"
        run_analysis()
    elseif REQUESTED_STAGE == "meshes"
        plot_reference_meshes()
    else
        error("usage: julia run_boundary_condition_study.jl [--stage=solve|analyze|meshes]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
