module DistributedSlowWavePilot

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_SLOW_WAVE_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "distributed_slow_wave_pilot"),
)
const MAX_WORKERS = parse(Int, get(ENV, "METAMATERIALS_SLOW_WAVE_WORKERS", "3"))
const VARIANTS = (:solid, :mid, :max)
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "distributed_slow_wave_mesher.jl"))
using .DistributedSlowWaveMesher

is_solver_stage(stage) = any(
    prefix -> stage == prefix || startswith(stage, "$prefix-worker-"),
    ("lossless", "fine-lossless", "fine-lossy"),
) || stage == "validate"

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) &&
       (is_solver_stage(REQUESTED_STAGE) ||
        REQUESTED_STAGE in ("analyze", "analyze-fine", "analyze-lossy"))
    include(joinpath(@__DIR__, "modal_harmonic_solver.jl"))
    include(joinpath(@__DIR__, "modal_projection.jl"))
    using .ModalHarmonicElasticity
    using .ElasticModalProjection
    using JLD2
end

if REQUESTED_STAGE in ("analyze", "analyze-fine", "analyze-lossy")
    ENV["GKSwstype"] = "100"
    using Plots
    using FFTW
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    using .ImpulseRiskAnalysis
end

const CONFIG = DistributedSlowWaveConfig(
    lead_length_mm=parse(Float64, get(ENV, "METAMATERIALS_SLOW_WAVE_LEAD_MM", "12.0")),
    ligament_height_mm=parse(
        Float64,
        get(ENV, "METAMATERIALS_SLOW_WAVE_LIGAMENT_MM", "0.45"),
    ),
)
const MODEL_TOTAL_LENGTH_MM = parse(
    Float64,
    get(
        ENV,
        "METAMATERIALS_THREE_STATE_TOTAL_LENGTH_MM",
        string(total_length_mm(CONFIG)),
    ),
)
const FIGURE_TITLE = get(
    ENV,
    "METAMATERIALS_THREE_STATE_TITLE",
    "Distributed slow-wave cell — 3.2 mm symmetric port",
)

const LOSSLESS_FREQUENCIES_HZ = Float64[
    162.0e3,
    180.0e3,
    200.0e3,
    220.0e3,
    242.0e3,
    260.0e3,
    280.0e3,
    300.0e3,
    322.0e3,
]
const FINE_FREQUENCIES_HZ = collect(162.0e3:4.0e3:322.0e3)

mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")
model_path(variant) = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
point_path(stage, variant, frequency_hz) = joinpath(
    OUTPUT_ROOT,
    String(stage),
    "$(variant)_$(round(Int, frequency_hz))hz.jld2",
)

lossless_specs() = [
    (; variant, frequency_hz)
    for variant in VARIANTS
    for frequency_hz in LOSSLESS_FREQUENCIES_HZ
]

fine_specs() = [
    (; variant, frequency_hz)
    for variant in VARIANTS
    for frequency_hz in FINE_FREQUENCIES_HZ
]

stage_specs(stage) = stage == :lossless ? lossless_specs() : fine_specs()
stage_loss_scale(stage) = stage == :fine_lossy ? 1.0 : 0.0

function child_command(stage)
    # Use the resolved runtime rather than the juliaup launcher. Concurrent
    # child launches otherwise contend for juliaup's configuration lock.
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            build_distributed_slow_wave_mesh(mesh_path(variant); config=CONFIG, variant)
        end
    finally
        gmsh.finalize()
    end
end

function run_convert_stage()
    for variant in VARIANTS
        convert_mesh(mesh_path(variant); output_dir=dirname(model_path(variant)))
    end
end

function project_all(state, tag)
    right_modes = filter(
        mode -> mode.kind == :propagating && mode.direction == :right,
        state.modes,
    )
    [
        begin
            left_candidates = filter(
                mode -> mode.kind == :propagating &&
                        mode.direction == :left &&
                        mode.parity == right_mode.parity,
                state.modes,
            )
            left_mode = first(sort(left_candidates; by=mode -> abs(
                mode.wavenumber_per_m + right_mode.wavenumber_per_m,
            )))
            projection = project_boundary_mode_pair(
                state,
                tag,
                right_mode,
                left_mode;
                quadrature_degree=6,
            )
            (; right_mode, left_mode, projection...)
        end
        for right_mode in right_modes
    ]
end

function solve_scattering(variant, frequency_hz; loss_scale=0.0)
    harmonic_config = ModalHarmonicElasticity.HarmonicElasticity.HarmonicConfig(
        rayleigh_alpha=Float64(loss_scale) * 79560.0,
        rayleigh_beta=Float64(loss_scale) * 2.5e-9,
        element_order=2,
        quadrature_degree=4,
    )
    state = solve_symmetry_reduced_modal_state(
        model_path(variant),
        frequency_hz;
        config=harmonic_config,
        port_height_m=CONFIG.port_height_mm * 1e-3,
        port_element_count=80,
    )
    left = project_all(state, "Source")
    right = project_all(state, "Microphone")
    input_index = findfirst(pair -> pair.right_mode.parity == :symmetric, left)
    input_index === nothing && error("symmetric input mode was not found")
    input_left = left[input_index]
    input_right = right[input_index]
    incident = abs2(input_left.right_amplitude)
    reflection = abs2(input_left.left_amplitude) / incident
    transmission = abs2(input_right.right_amplitude) / incident
    right_incoming = abs2(input_right.left_amplitude) / incident
    conversion = sum(
        abs2(left_pair.left_amplitude) + abs2(right_pair.right_amplitude)
        for (index, (left_pair, right_pair)) in enumerate(zip(left, right))
        if index != input_index
    ) / incident
    transmission_amplitude = input_right.right_amplitude / input_left.right_amplitude
    length_m = MODEL_TOTAL_LENGTH_MM * 1e-3
    wavenumber = real(state.right_mode.wavenumber_per_m)
    deembedded = transmission_amplitude * exp(im * wavenumber * length_m)
    (
        variant,
        frequency_hz=Float64(frequency_hz),
        loss_scale=Float64(loss_scale),
        incident_power=incident,
        T00=transmission,
        R00=reflection,
        C0=conversion,
        right_incoming,
        transmission_amplitude=ComplexF64(transmission_amplitude),
        deembedded_amplitude=abs(deembedded),
        deembedded_phase_rad=angle(deembedded),
        propagating_mode_count=length(left),
        power_sum=transmission + reflection + conversion,
        variational_balance_residual=state.metrics.balance_residual_w_per_m / incident,
    )
end

function solve_and_save(stage, spec)
    row = solve_scattering(
        spec.variant,
        spec.frequency_hz;
        loss_scale=stage_loss_scale(stage),
    )
    path = point_path(stage, spec.variant, spec.frequency_hz)
    mkpath(dirname(path))
    JLD2.jldsave(path; row)
    println(
        "[+] $stage $(spec.variant) ", round(spec.frequency_hz / 1e3; digits=1),
        " kHz: T=", round(row.T00; digits=4),
        ", R=", round(row.R00; digits=4),
        ", C=", round(row.C0; digits=3),
        ", sum=", round(row.power_sum; digits=4),
    )
end

function run_worker_specs(stage, worker_index)
    specs = stage_specs(stage)
    for index in worker_index:MAX_WORKERS:length(specs)
        spec = specs[index]
        isfile(point_path(stage, spec.variant, spec.frequency_hz)) && continue
        solve_and_save(stage, spec)
    end
end

function run_parallel_specs(stage)
    specs = stage_specs(stage)
    any_missing = any(
        !isfile(point_path(stage, spec.variant, spec.frequency_hz)) for spec in specs
    )
    any_missing || return
    @sync for worker_index in 1:MAX_WORKERS
        stage_text = replace(String(stage), '_' => '-')
        @async run(child_command("$stage_text-worker-$worker_index"))
    end
end

function load_rows(stage)
    [
        JLD2.load(point_path(stage, spec.variant, spec.frequency_hz))["row"]
        for spec in stage_specs(stage)
    ]
end

function unwrap_phase(values)
    result = Float64[first(values)]
    for raw in Iterators.drop(values, 1)
        push!(result, last(result) + mod(raw - last(result) + pi, 2pi) - pi)
    end
    result
end

function phase_delay_rows(rows, variant)
    selected = sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
    references = Dict(
        row.frequency_hz => row for row in rows if row.variant == :solid
    )
    relative_transfer = [
        row.transmission_amplitude / references[row.frequency_hz].transmission_amplitude
        for row in selected
    ]
    phase = unwrap_phase(angle.(relative_transfer))
    omega = 2pi .* getproperty.(selected, :frequency_hz)
    delay = similar(phase)
    delay[1] = -(phase[2] - phase[1]) / (omega[2] - omega[1])
    for index in 2:(length(phase) - 1)
        delay[index] = -(phase[index + 1] - phase[index - 1]) /
                       (omega[index + 1] - omega[index - 1])
    end
    delay[end] = -(phase[end] - phase[end - 1]) / (omega[end] - omega[end - 1])
    [
        merge(
            row,
            (
                relative_amplitude=abs(relative_transfer[index]),
                relative_power=abs2(relative_transfer[index]),
                relative_phase_rad=phase[index],
                excess_group_delay_s=delay[index],
            ),
        )
        for (index, row) in enumerate(selected)
    ]
end

function linear_interpolate(x, values, query)
    query <= first(x) && return first(values)
    query >= last(x) && return last(values)
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function trapezoid_integral(x, values)
    sum((x[index + 1] - x[index]) * (values[index + 1] + values[index]) / 2
        for index in 1:(length(x) - 1))
end

function weighted_metrics(rows, pulse; calibrated=false)
    frequencies = getproperty.(rows, :frequency_hz)
    spectral_power = [
        abs2(linear_interpolate(pulse.frequency_hz, pulse.spectrum, frequency))
        for frequency in frequencies
    ]
    normalization = trapezoid_integral(frequencies, spectral_power)
    weighted(field) = trapezoid_integral(
        frequencies,
        spectral_power .* getproperty.(rows, field),
    ) / normalization
    (
        Tbar=weighted(calibrated ? :relative_power : :T00),
        Rbar=weighted(:R00),
        Cbar=weighted(:C0),
        relative_power_bar=weighted(:relative_power),
        carrier_delay_us=linear_interpolate(
            frequencies,
            getproperty.(rows, :excess_group_delay_s),
            pulse.config.center_frequency_hz,
        ) * 1e6,
        minimum_T=minimum(getproperty.(rows, calibrated ? :relative_power : :T00)),
        maximum_R=maximum(getproperty.(rows, :R00)),
        maximum_power_error=maximum(abs.(getproperty.(rows, :power_sum) .- 1)),
    )
end

function write_csv(rows, stem)
    path = joinpath(OUTPUT_ROOT, "$(stem).csv")
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
    println("[+] $path")
end

function save_figure(rows, stem)
    power_panel = plot(
        xlabel="frequency, kHz",
        ylabel="power fraction",
        title="Absolute lossless scattering",
        ylims=(0, 1.05),
        gridalpha=0.25,
    )
    phase_panel = plot(
        xlabel="frequency, kHz",
        ylabel="relative phase, rad",
        title="Phase relative to solid strip",
        gridalpha=0.25,
    )
    delay_panel = plot(
        xlabel="frequency, kHz",
        ylabel="excess group delay, μs",
        title="Coarse phase-slope diagnostic",
        gridalpha=0.25,
    )
    reflection_panel = plot(
        xlabel="frequency, kHz",
        ylabel="R00",
        title="Reflection",
        ylims=(0, 1.05),
        gridalpha=0.25,
    )
    styles = Dict(:solid => (:black, :circle), :mid => (:darkorange, :diamond), :max => (:royalblue, :square))
    for variant in VARIANTS
        current = sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
        color, marker = styles[variant]
        frequencies = getproperty.(current, :frequency_hz) ./ 1e3
        plot!(power_panel, frequencies, getproperty.(current, :T00);
              color, marker, linewidth=2, label=String(variant))
        plot!(reflection_panel, frequencies, getproperty.(current, :R00);
              color, marker, linewidth=2, label=String(variant))
        plot!(phase_panel, frequencies, getproperty.(current, :relative_phase_rad);
              color, marker, linewidth=2, label=String(variant))
        plot!(delay_panel, frequencies, getproperty.(current, :excess_group_delay_s) .* 1e6;
              color, marker, linewidth=2, label=String(variant))
    end
    for panel in (power_panel, reflection_panel, phase_panel, delay_panel)
        vline!(panel, [242.0]; color=:gray, linestyle=:dash, label=false)
    end
    figure = plot(
        power_panel,
        reflection_panel,
        phase_panel,
        delay_panel;
        layout=(2, 2),
        size=(1450, 900),
        margin=5Plots.mm,
        plot_title=FIGURE_TITLE,
    )
    path = joinpath(OUTPUT_ROOT, "$(stem).png")
    savefig(figure, path)
    println("[+] $path")
end

function reconstructed_transfer(rows, pulse; common_delay_s=100.0e-6)
    frequencies = getproperty.(rows, :frequency_hz)
    amplitude = getproperty.(rows, :relative_amplitude)
    phase = getproperty.(rows, :relative_phase_rad)
    ComplexF64[
        linear_interpolate(frequencies, amplitude, frequency_hz) *
        cis(linear_interpolate(frequencies, phase, frequency_hz)) *
        cis(-2pi * frequency_hz * common_delay_s)
        for frequency_hz in pulse.frequency_hz
    ]
end

function save_impulse_diagnostics(rows, pulse, stem)
    dt_s = pulse.time_s[2] - pulse.time_s[1]
    common_delay_s = 100.0e-6
    reference_transfer = cis.(-2pi .* pulse.frequency_hz .* common_delay_s)
    reference = irfft(pulse.spectrum .* reference_transfer, length(pulse.signal))
    metrics_rows = NamedTuple[]
    waveforms = Dict{Symbol, Vector{Float64}}(:solid => reference)
    for variant in (:mid, :max)
        current = sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
        transfer = reconstructed_transfer(current, pulse; common_delay_s)
        waveform = irfft(pulse.spectrum .* transfer, length(pulse.signal))
        waveforms[variant] = waveform
        metrics = pulse_metrics(waveform, reference, reference, dt_s)
        push!(metrics_rows, (variant=variant, metrics...))
    end

    csv_path = joinpath(OUTPUT_ROOT, "$(stem)_impulse_metrics.csv")
    columns = propertynames(first(metrics_rows))
    open(csv_path, "w") do io
        println(io, join(string.(columns), ','))
        for row in metrics_rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end

    window = (pulse.time_s .>= 65.0e-6) .& (pulse.time_s .<= 145.0e-6)
    figure = plot(
        pulse.time_s[window] .* 1e6,
        reference[window];
        linewidth=2,
        color=:black,
        label="solid reference",
        xlabel="time, μs",
        ylabel="normalized response",
        title="Provisional pulse reconstructed from fine S(ω)",
        gridalpha=0.25,
    )
    for (variant, color) in ((:mid, :darkorange), (:max, :royalblue))
        plot!(figure, pulse.time_s[window] .* 1e6, waveforms[variant][window];
              linewidth=2, color, label=String(variant))
    end
    png_path = joinpath(OUTPUT_ROOT, "$(stem)_impulse.png")
    savefig(figure, png_path)
    println("[+] $csv_path")
    println("[+] $png_path")
    for row in metrics_rows
        println(
            "[+] impulse $(row.variant): peak ratio=", round(row.gain_peak; digits=4),
            ", Bt=", round(row.broadening_ratio; digits=3),
            ", rho=", round(row.pulse_correlation; digits=4),
            ", postcursor=", round(row.postcursor_ratio; digits=4),
        )
    end
end

function run_analysis_stage(stage=:lossless)
    raw = load_rows(stage)
    rows = reduce(vcat, (phase_delay_rows(raw, variant) for variant in VARIANTS))
    stem = stage == :lossless ? "distributed_slow_wave_lossless" :
           "distributed_slow_wave_$(stage)"
    write_csv(rows, stem)
    save_figure(rows, stem)
    pulse = pulse_spectrum(PulseConfig(fft_length=65536))
    summary_path = joinpath(OUTPUT_ROOT, "$(stem)_summary.csv")
    open(summary_path, "w") do io
        println(io, "variant,Tbar,Rbar,Cbar,relative_power_bar,carrier_delay_us,minimum_T,maximum_R,maximum_power_error")
        for variant in VARIANTS
            current = sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
            metrics = weighted_metrics(current, pulse; calibrated=stage == :fine_lossy)
            println(io, join((variant, values(metrics)...), ','))
            println(
                "[+] $variant: Tbar=", round(metrics.Tbar; digits=4),
                ", Rbar=", round(metrics.Rbar; digits=4),
                ", delay@242=", round(metrics.carrier_delay_us; digits=3), " μs",
                ", min T=", round(metrics.minimum_T; digits=4),
                ", max |sum-1|=", round(metrics.maximum_power_error; sigdigits=3),
            )
        end
    end
    println("[+] $summary_path")
    stage in (:fine_lossless, :fine_lossy) && save_impulse_diagnostics(rows, pulse, stem)
end

function run_validation_stage()
    frequencies_hz = (242.0e3, 318.0e3)
    path = joinpath(OUTPUT_ROOT, "lead_validation.csv")
    open(path, "w") do io
        println(io, "lead_mm,frequency_hz,T00,R00,C0,relative_amplitude,relative_phase_rad,power_sum")
        for frequency_hz in frequencies_hz
            solid = solve_scattering(:solid, frequency_hz)
            trim = solve_scattering(:max, frequency_hz)
            relative = trim.transmission_amplitude / solid.transmission_amplitude
            println(io, join((
                CONFIG.lead_length_mm,
                frequency_hz,
                trim.T00,
                trim.R00,
                trim.C0,
                abs(relative),
                angle(relative),
                trim.power_sum,
            ), ','))
            println(
                "[+] lead=", CONFIG.lead_length_mm,
                " mm, f=", frequency_hz / 1e3,
                " kHz: T=", round(trim.T00; digits=5),
                ", R=", round(trim.R00; digits=5),
                ", phase=", round(angle(relative); digits=4),
                ", sum=", round(trim.power_sum; digits=5),
            )
        end
    end
    println("[+] $path")
end

function run_child_stage(stage)
    if stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage == "lossless"
        run_parallel_specs(:lossless)
    elseif stage == "fine-lossless"
        run_parallel_specs(:fine_lossless)
    elseif stage == "fine-lossy"
        run_parallel_specs(:fine_lossy)
    elseif any(prefix -> startswith(stage, "$prefix-worker-"),
               ("lossless", "fine-lossless", "fine-lossy"))
        prefix = first(filter(prefix -> startswith(stage, "$prefix-worker-"),
                              ("lossless", "fine-lossless", "fine-lossy")))
        worker_index = parse(Int, split(stage, "-worker-"; limit=2)[2])
        1 <= worker_index <= MAX_WORKERS || error("invalid worker index")
        run_worker_specs(Symbol(replace(prefix, '-' => '_')), worker_index)
    elseif stage == "analyze"
        run_analysis_stage(:lossless)
    elseif stage == "analyze-fine"
        run_analysis_stage(:fine_lossless)
    elseif stage == "analyze-lossy"
        run_analysis_stage(:fine_lossy)
    elseif stage == "validate"
        run_validation_stage()
    else
        error("unknown distributed slow-wave stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("mesh", "convert", "lossless", "analyze")
            println("\n=== Distributed slow-wave stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_distributed_slow_wave_pilot.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
