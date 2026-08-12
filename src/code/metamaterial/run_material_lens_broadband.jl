module RunMaterialLensBroadband

ENV["GKSwstype"] = "100"

using FFTW
using JLD2
using Plots
using Statistics

include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
using .ImpulseRiskAnalysis
include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
using .SinusoidalMaterialLens

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const FEM_DRIVER = joinpath(@__DIR__, "run_material_lens_fem.jl")
const DESIGN_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_BROADBAND_DESIGN",
    joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_5cycle_242khz"),
)
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_BROADBAND_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_material_lens_broadband_5cycle_242khz"),
)
const JOBS = parse(Int, get(ENV, "METAMATERIALS_MATERIAL_LENS_BROADBAND_JOBS", "2"))
const ROLES = (:selected, :uniform)
const FREQUENCIES_HZ = Float64[
    162_200, 178_160, 194_120, 202_100, 210_080, 218_060, 226_040,
    234_020, 242_000, 249_970, 257_940, 265_910, 273_880, 281_850,
    289_820, 305_760, 321_700,
]

frequency_tag(frequency_hz) = "f$(round(Int, frequency_hz))hz"
result_path(role, frequency_hz) = joinpath(
    OUTPUT_ROOT,
    "results",
    "aluminium_$(role)_$(frequency_tag(frequency_hz)).jld2",
)
mesh_path(role) = joinpath(OUTPUT_ROOT, "meshes", "aluminium_$(role).msh")
specs() = [(; role, frequency_hz) for role in ROLES for frequency_hz in FREQUENCIES_HZ]

function driver_command(action, role; frequency_hz=nothing)
    julia = Base.julia_cmd()
    command = `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $FEM_DRIVER $action aluminium:$role`
    environment = Pair{String, String}[
        "METAMATERIALS_MATERIAL_LENS_DESIGN" => DESIGN_ROOT,
        "METAMATERIALS_MATERIAL_LENS_FEM_OUTPUT" => OUTPUT_ROOT,
        "METAMATERIALS_MATERIAL_LENS_FEM_ORDER" => "1",
        "METAMATERIALS_MATERIAL_LENS_H_LENS_MM" => "0.42",
        "METAMATERIALS_MATERIAL_LENS_H_OUTPUT_MM" => "1.55",
        "METAMATERIALS_MATERIAL_LENS_SCAN" => "false",
    ]
    if !isnothing(frequency_hz)
        push!(environment, "METAMATERIALS_MATERIAL_LENS_FREQUENCY_HZ" => string(Float64(frequency_hz)))
    end
    addenv(command, environment...)
end

function build_meshes()
    for role in ROLES
        isfile(mesh_path(role)) && continue
        run(driver_command("--mesh", role))
    end
end

function solve_all()
    build_meshes()
    pending = filter(spec -> !isfile(result_path(spec.role, spec.frequency_hz)), specs())
    isempty(pending) && return println("[=] All broadband material-lens points exist")
    jobs = min(max(JOBS, 1), length(pending), Sys.CPU_THREADS)
    semaphore = Base.Semaphore(jobs)
    println("=== Al material-lens broadband: $(length(pending)) solves, $jobs processes ===")
    @sync for spec in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(driver_command("--solve", spec.role; frequency_hz=spec.frequency_hz))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function unwrap_phase(phases)
    result = Float64[first(phases)]
    for phase in Iterators.drop(phases, 1)
        push!(result, last(result) + mod(phase - last(result) + pi, 2pi) - pi)
    end
    result
end

function linear_interpolate(x, y, query)
    length(x) == length(y) || throw(DimensionMismatch("x and y differ"))
    first(x) <= query <= last(x) || throw(ArgumentError("query outside interpolation range"))
    right = searchsortedfirst(x, query)
    right == firstindex(x) && return y[right]
    right > lastindex(x) && return y[end]
    left = right - 1
    weight = (query - x[left]) / (x[right] - x[left])
    (1 - weight) * y[left] + weight * y[right]
end

function sampled_response(role)
    data = [load(result_path(role, frequency_hz)) for frequency_hz in FREQUENCIES_HZ]
    displacement = ComplexF64[item["focus_displacement"][1] for item in data]
    amplitude = abs.(displacement)
    phase = unwrap_phase(angle.(displacement))
    (; role, frequency_hz=copy(FREQUENCIES_HZ), displacement, amplitude, phase)
end

function interpolate_response(sample, frequency_hz)
    log_amplitude = log.(max.(sample.amplitude, eps(Float64)))
    exp(linear_interpolate(sample.frequency_hz, log_amplitude, frequency_hz)) *
    cis(linear_interpolate(sample.frequency_hz, sample.phase, frequency_hz))
end

function bandlimited_waveform(pulse, sample)
    active = (pulse.frequency_hz .>= first(FREQUENCIES_HZ)) .&
             (pulse.frequency_hz .<= last(FREQUENCIES_HZ))
    transfer = zeros(ComplexF64, length(pulse.frequency_hz))
    for index in findall(active)
        transfer[index] = interpolate_response(sample, pulse.frequency_hz[index])
    end
    waveform = irfft(pulse.spectrum .* transfer, length(pulse.signal))
    (; waveform, transfer, active)
end

function ideal_scalar_same_aperture(pulse)
    config = MaterialLensConfig(
        center_frequency_hz=242.0e3,
        pulse_cycles=5.0,
        element_count=15,
        element_height_mm=7.0,
        slot_width_mm=1.2,
        focal_distance_mm=60.0,
    )
    material = aluminium_6061()
    centers_mm = SinusoidalMaterialLens.lens_centers_mm(config)
    distances_mm = hypot.(config.focal_distance_mm, centers_mm)
    delays_s = (maximum(distances_mm) .- distances_mm) .* 1.0e-3 ./
               material.pressure_wave_speed_m_s
    active = (pulse.frequency_hz .>= first(FREQUENCIES_HZ)) .&
             (pulse.frequency_hz .<= last(FREQUENCIES_HZ))
    focused_transfer = zeros(ComplexF64, length(pulse.frequency_hz))
    uniform_transfer = zeros(ComplexF64, length(pulse.frequency_hz))
    for index in findall(active)
        frequency_hz = pulse.frequency_hz[index]
        kernels = ComplexF64[
            SinusoidalMaterialLens.aperture_kernel(center, frequency_hz, material, config)
            for center in centers_mm
        ]
        uniform_transfer[index] = sum(kernels)
        focused_transfer[index] = sum(
            kernels .* cis.(-2pi * frequency_hz .* delays_s),
        )
    end
    focused = irfft(pulse.spectrum .* focused_transfer, length(pulse.signal))
    uniform = irfft(pulse.spectrum .* uniform_transfer, length(pulse.signal))
    metrics = pulse_metrics(
        focused,
        uniform,
        uniform,
        pulse.time_s[2] - pulse.time_s[1],
    )
    (; focused, uniform, focused_transfer, uniform_transfer, delays_s, metrics)
end

function analytic_envelope(signal)
    count = length(signal)
    spectrum = fft(signal)
    multiplier = zeros(Float64, count)
    multiplier[1] = 1.0
    if iseven(count)
        multiplier[2:(count ÷ 2)] .= 2.0
        multiplier[count ÷ 2 + 1] = 1.0
    else
        multiplier[2:((count + 1) ÷ 2)] .= 2.0
    end
    abs.(ifft(spectrum .* multiplier))
end

function group_delay_s(frequency_hz, unwrapped_phase)
    omega = 2pi .* frequency_hz
    result = similar(unwrapped_phase)
    result[1] = -(unwrapped_phase[2] - unwrapped_phase[1]) / (omega[2] - omega[1])
    for index in 2:(length(result) - 1)
        result[index] = -(unwrapped_phase[index + 1] - unwrapped_phase[index - 1]) /
                        (omega[index + 1] - omega[index - 1])
    end
    result[end] = -(unwrapped_phase[end] - unwrapped_phase[end - 1]) /
                  (omega[end] - omega[end - 1])
    result
end

function write_frequency_csv(selected, uniform)
    ratio = selected.displacement ./ uniform.displacement
    ratio_phase = unwrap_phase(angle.(ratio))
    ratio_delay = group_delay_s(FREQUENCIES_HZ, ratio_phase)
    path = joinpath(OUTPUT_ROOT, "aluminium_broadband_focus_response.csv")
    open(path, "w") do io
        println(io, "frequency_hz,selected_focus_ux_real_m,selected_focus_ux_imag_m,selected_focus_ux_abs_m,uniform_focus_ux_real_m,uniform_focus_ux_imag_m,uniform_focus_ux_abs_m,relative_amplitude,relative_phase_rad,relative_group_delay_s")
        for index in eachindex(FREQUENCIES_HZ)
            println(io, join((
                FREQUENCIES_HZ[index],
                real(selected.displacement[index]), imag(selected.displacement[index]),
                selected.amplitude[index],
                real(uniform.displacement[index]), imag(uniform.displacement[index]),
                uniform.amplitude[index], abs(ratio[index]), ratio_phase[index],
                ratio_delay[index],
            ), ','))
        end
    end
    path
end

function make_plot(selected, uniform, pulse, selected_wave, uniform_wave, metrics, ideal_gain)
    response = plot(
        FREQUENCIES_HZ ./ 1.0e3, selected.amplitude .* 1.0e9;
        xlabel="frequency, kHz", ylabel="|uₓ(focus)|, nm",
        title="Full FEM focus response", marker=:circle, lw=2.5, label="lens",
    )
    plot!(response, FREQUENCIES_HZ ./ 1.0e3, uniform.amplitude .* 1.0e9;
          marker=:circle, lw=2.5, ls=:dash, label="straight")

    ratio = selected.displacement ./ uniform.displacement
    ratio_phase = unwrap_phase(angle.(ratio))
    relative = plot(
        FREQUENCIES_HZ ./ 1.0e3, abs.(ratio);
        xlabel="frequency, kHz", ylabel="|u_lens/u_straight|",
        title="Matched broadband gain", marker=:circle, lw=2.5,
        label="amplitude gain",
    )
    hline!(relative, [2.0]; ls=:dash, color=:black, label="target 2")
    phase_axis = twinx(relative)
    plot!(phase_axis, FREQUENCIES_HZ ./ 1.0e3, rad2deg.(ratio_phase);
          ylabel="unwrapped relative phase, deg", marker=:square, lw=2,
          color=:darkorange, label="phase")

    selected_envelope = analytic_envelope(selected_wave.waveform)
    uniform_envelope = analytic_envelope(uniform_wave.waveform)
    peak_index = argmax(selected_envelope)
    time_us = pulse.time_s .* 1.0e6
    left = max(firstindex(time_us), peak_index - 650)
    right = min(lastindex(time_us), peak_index + 900)
    waveform = plot(
        time_us[left:right], selected_wave.waveform[left:right] .* 1.0e9;
        xlabel="time, μs", ylabel="uₓ(focus), nm",
        title="Band-limited 5-cycle focus", lw=1.6, label="lens",
    )
    plot!(waveform, time_us[left:right], uniform_wave.waveform[left:right] .* 1.0e9;
          lw=1.6, ls=:dash, label="straight")
    plot!(waveform, time_us[left:right], selected_envelope[left:right] .* 1.0e9;
          lw=2.3, color=:navy, label="lens envelope")
    plot!(waveform, time_us[left:right], uniform_envelope[left:right] .* 1.0e9;
          lw=2.3, color=:red, ls=:dash, label="straight envelope")

    summary = bar(
        ["G_FEM", "G_ideal", "B_t", "ρ", "postcursor"],
        [metrics.gain_peak, ideal_gain, metrics.broadening_ratio,
         metrics.pulse_correlation, metrics.postcursor_ratio];
        ylabel="metric value", title="Impulse metrics", label=false,
        color=[:steelblue, :goldenrod, :darkorange, :seagreen, :mediumpurple],
    )
    hline!(summary, [1.0]; color=:gray45, ls=:dot, label=false)

    path = joinpath(OUTPUT_ROOT, "aluminium_broadband_impulse.png")
    savefig(plot(
        response, relative, waveform, summary;
        layout=(2, 2), size=(1450, 950), margin=5Plots.mm,
    ), path)
    path
end

function analyze()
    all(isfile(result_path(spec.role, spec.frequency_hz)) for spec in specs()) ||
        error("broadband FEM set is incomplete")
    mkpath(OUTPUT_ROOT)
    selected = sampled_response(:selected)
    uniform = sampled_response(:uniform)
    pulse = pulse_spectrum(PulseConfig(
        center_frequency_hz=242.0e3,
        cycles=5.0,
        samples_per_period=80,
        fft_length=65536,
    ))
    selected_wave = bandlimited_waveform(pulse, selected)
    uniform_wave = bandlimited_waveform(pulse, uniform)
    ideal = ideal_scalar_same_aperture(pulse)
    dt_s = pulse.time_s[2] - pulse.time_s[1]
    metrics = pulse_metrics(
        selected_wave.waveform,
        uniform_wave.waveform,
        uniform_wave.waveform,
        dt_s,
    )
    energy_fraction = spectral_energy_fraction(
        pulse;
        lower_hz=first(FREQUENCIES_HZ),
        upper_hz=last(FREQUENCIES_HZ),
    )
    carrier_index = argmin(abs.(FREQUENCIES_HZ .- 242.0e3))
    carrier_gain = selected.amplitude[carrier_index] / uniform.amplitude[carrier_index]
    passed = metrics.gain_peak >= 2.0 && metrics.broadening_ratio <= 1.25 &&
             metrics.pulse_correlation >= 0.90 && metrics.postcursor_ratio <= 0.10
    frequency_csv = write_frequency_csv(selected, uniform)
    summary_path = joinpath(OUTPUT_ROOT, "aluminium_broadband_impulse_summary.csv")
    open(summary_path, "w") do io
        println(io, "lower_frequency_hz,upper_frequency_hz,band_spectral_energy_fraction,selected_peak_m,uniform_peak_m,impulse_peak_gain,carrier_gain,scalar_ideal_same_aperture_peak_gain,broadening_ratio,pulse_correlation,postcursor_ratio,gate_passed")
        println(io, join((
            first(FREQUENCIES_HZ), last(FREQUENCIES_HZ), energy_fraction,
            metrics.focused_peak, metrics.uniform_peak, metrics.gain_peak,
            carrier_gain, ideal.metrics.gain_peak, metrics.broadening_ratio,
            metrics.pulse_correlation,
            metrics.postcursor_ratio, passed,
        ), ','))
    end
    plot_path = make_plot(
        selected,
        uniform,
        pulse,
        selected_wave,
        uniform_wave,
        metrics,
        ideal.metrics.gain_peak,
    )
    data_path = joinpath(OUTPUT_ROOT, "aluminium_broadband_impulse.jld2")
    jldsave(
        data_path;
        format_version=1,
        selected,
        uniform,
        pulse_time_s=pulse.time_s,
        pulse_input=pulse.signal,
        selected_waveform_m=selected_wave.waveform,
        uniform_waveform_m=uniform_wave.waveform,
        selected_transfer_m=selected_wave.transfer,
        uniform_transfer_m=uniform_wave.transfer,
        ideal_scalar_focused_waveform=ideal.focused,
        ideal_scalar_uniform_waveform=ideal.uniform,
        ideal_scalar_delay_s=ideal.delays_s,
        ideal_scalar_metrics=ideal.metrics,
        band_spectral_energy_fraction=energy_fraction,
        metrics,
        carrier_gain,
        gate_passed=passed,
    )
    println("[+] Band spectral energy fraction: $energy_fraction")
    println("[+] G_peak=$(metrics.gain_peak), carrier=$carrier_gain")
    println("[+] Ideal scalar same-aperture G_peak=$(ideal.metrics.gain_peak)")
    println("[+] Bt=$(metrics.broadening_ratio), rho=$(metrics.pulse_correlation), post=$(metrics.postcursor_ratio)")
    println("[+] Gate passed: $passed")
    println("[+] Frequency response: $frequency_csv")
    println("[+] Summary: $summary_path")
    println("[+] Plot: $plot_path")
    println("[+] Data: $data_path")
end

function main(args=ARGS)
    stage = isempty(args) ? "all" : only(args)
    if stage == "--mesh"
        build_meshes()
    elseif stage == "--solve"
        solve_all()
    elseif stage == "--analyze"
        analyze()
    elseif stage == "--all"
        solve_all()
        analyze()
    else
        error("usage: run_material_lens_broadband.jl --mesh|--solve|--analyze|--all")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
