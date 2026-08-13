module TopologyHybridTTDGRINIsolated

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HYBRID_TTD_GRIN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "hybrid_ttd_grin", "isolated"),
)
const WORKERS = parse(Int, get(ENV, "METAMATERIALS_HYBRID_TTD_GRIN_WORKERS", "2"))
const VARIANTS = (:straight, :path_only, :hybrid, :pure_same_depth)
const FREQUENCIES_HZ = Float64[162.2e3, 180e3, 200e3, 220e3, 242e3, 260e3,
                               280e3, 300e3, 321.7e3]

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "analyze")

include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_design.jl"))
using .TopologyHybridTTDGRINDesign

if STAGE == "mesh"
    using Gmsh: gmsh
    include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_mesher.jl"))
    using .TopologyHybridTTDGRINMesher
elseif STAGE == "convert"
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif STAGE == "solve" || startswith(STAGE, "solve-worker-")
    using JLD2
    include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_harmonic.jl"))
    using .TopologyHybridTTDGRINHarmonic
elseif STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using FFTW
    using JLD2
    using Plots
    include(joinpath(@__DIR__, "topology_hybrid_ttd_grin_transfer.jl"))
    using .TopologyHybridTTDGRINTransfer
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    using .ImpulseRiskAnalysis
end

const DESIGN = HybridTTDGRINConfig()

mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")
model_path(variant) = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
point_path(variant, frequency_hz) = joinpath(
    OUTPUT_ROOT,
    "points",
    "$(variant)_$(round(Int, frequency_hz))hz.jld2",
)

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            build_hybrid_channel_mesh(
                mesh_path(variant);
                variant,
                config=HybridChannelMeshConfig(design=DESIGN),
            )
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

material_fractions(variant) = variant == :hybrid ?
    Float64[0.0, 0.2, 0.4, 0.6, 0.8, 1.0] : Float64[0.0]

specs() = [(; variant, frequency_hz) for variant in VARIANTS for frequency_hz in FREQUENCIES_HZ]

function solve_one(spec)
    path = point_path(spec.variant, spec.frequency_hz)
    isfile(path) && return
    row = solve_hybrid_channel_scattering(
        model_path(spec.variant),
        spec.frequency_hz,
        material_fractions(spec.variant);
        design=DESIGN,
    )
    row = merge(row, (variant=spec.variant,))
    mkpath(dirname(path))
    jldsave(path; row)
    println(
        "[+] $(spec.variant) $(spec.frequency_hz / 1e3) kHz: T=",
        round(row.T00; digits=4),
        ", R=", round(row.R00; digits=4),
        ", C=", round(row.C0; digits=4),
    )
end

function child_command(stage)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_worker(index)
    all_specs = specs()
    for spec_index in index:WORKERS:length(all_specs)
        solve_one(all_specs[spec_index])
    end
end

function run_solve_stage()
    @sync for worker in 1:WORKERS
        @async run(child_command("solve-worker-$worker"))
    end
end

function load_rows()
    [JLD2.load(point_path(spec.variant, spec.frequency_hz))["row"] for spec in specs()]
end

function unwrap_phase(values)
    result = Float64[first(values)]
    for value in Iterators.drop(values, 1)
        push!(result, last(result) + mod(value - last(result) + pi, 2pi) - pi)
    end
    result
end

function relative_rows(rows, variant)
    selected = sort(filter(row -> row.variant == variant, rows); by=row -> row.frequency_hz)
    reference = Dict(
        row.frequency_hz => row for row in rows if row.variant == :straight
    )
    transfer = [row.S21 / reference[row.frequency_hz].S21 for row in selected]
    phase = unwrap_phase(angle.(transfer))
    frequency = getproperty.(selected, :frequency_hz)
    delay = similar(phase)
    delay[1] = -(phase[2] - phase[1]) / (2pi * (frequency[2] - frequency[1]))
    for index in 2:(length(phase) - 1)
        delay[index] = -(phase[index + 1] - phase[index - 1]) /
                       (2pi * (frequency[index + 1] - frequency[index - 1]))
    end
    delay[end] = -(phase[end] - phase[end - 1]) /
                 (2pi * (frequency[end] - frequency[end - 1]))
    [merge(row, (
        relative_S21=transfer[index],
        relative_amplitude=abs(transfer[index]),
        relative_phase_rad=phase[index],
        excess_group_delay_s=delay[index],
    )) for (index, row) in enumerate(selected)]
end

function linear_interpolate(x, values, query)
    query <= first(x) && return first(values)
    query >= last(x) && return last(values)
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function transfer_from_rows(relative, pulse)
    frequencies = getproperty.(relative, :frequency_hz)
    amplitude = getproperty.(relative, :relative_amplitude)
    phase = getproperty.(relative, :relative_phase_rad)
    delay_low = relative[1].excess_group_delay_s
    delay_high = relative[end].excess_group_delay_s
    ComplexF64[
        if frequency < first(frequencies)
            first(amplitude) * cis(first(phase) - 2pi * (frequency - first(frequencies)) * delay_low)
        elseif frequency > last(frequencies)
            last(amplitude) * cis(last(phase) - 2pi * (frequency - last(frequencies)) * delay_high)
        else
            linear_interpolate(frequencies, amplitude, frequency) *
            cis(linear_interpolate(frequencies, phase, frequency))
        end
        for frequency in pulse.frequency_hz
    ]
end

function weighted(values, frequencies, pulse)
    weights = abs2.(ComplexF64[
        linear_interpolate(pulse.frequency_hz, pulse.spectrum, frequency)
        for frequency in frequencies
    ])
    sum(values .* weights) / sum(weights)
end

function pulse_result(relative, pulse)
    transfer = transfer_from_rows(relative, pulse)
    mean_delay = weighted(
        getproperty.(relative, :excess_group_delay_s),
        getproperty.(relative, :frequency_hz),
        pulse,
    )
    ideal_transfer = cis.(-2pi .* pulse.frequency_hz .* mean_delay)
    waveform = irfft(pulse.spectrum .* transfer, length(pulse.signal))
    ideal = irfft(pulse.spectrum .* ideal_transfer, length(pulse.signal))
    metrics = pulse_metrics(
        waveform,
        ideal,
        ideal,
        pulse.time_s[2] - pulse.time_s[1],
    )
    (; transfer, waveform, ideal, mean_delay, metrics)
end

function write_csv(path, rows)
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
end

function run_analysis_stage()
    rows = load_rows()
    relative = reduce(vcat, (relative_rows(rows, variant) for variant in VARIANTS))
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "isolated_frequency.csv"), relative)
    pulse = pulse_spectrum(PulseConfig(
        center_frequency_hz=DESIGN.frequency_hz,
        cycles=5.0,
        samples_per_period=80,
        fft_length=16384,
    ))
    summary_rows = NamedTuple[]
    pulse_results = Dict{Symbol, Any}()
    for variant in VARIANTS
        current = sort(filter(row -> row.variant == variant, relative); by=row -> row.frequency_hz)
        result = pulse_result(current, pulse)
        pulse_results[variant] = result
        frequencies = getproperty.(current, :frequency_hz)
        pulse_metrics_value = result.metrics
        push!(summary_rows, (
            variant,
            weighted_transmission=weighted(getproperty.(current, :T00), frequencies, pulse),
            weighted_reflection=weighted(getproperty.(current, :R00), frequencies, pulse),
            weighted_mode_conversion=weighted(getproperty.(current, :C0), frequencies, pulse),
            minimum_transmission=minimum(getproperty.(current, :T00)),
            maximum_reflection=maximum(getproperty.(current, :R00)),
            excess_group_delay_us=result.mean_delay * 1e6,
            delay_ripple_us=(maximum(getproperty.(current, :excess_group_delay_s)) -
                             minimum(getproperty.(current, :excess_group_delay_s))) * 1e6,
            pulse_peak_ratio=pulse_metrics_value.gain_peak,
            pulse_broadening_ratio=pulse_metrics_value.broadening_ratio,
            pulse_correlation=pulse_metrics_value.pulse_correlation,
            postcursor_ratio=pulse_metrics_value.postcursor_ratio,
            gate_pass=
                weighted(getproperty.(current, :T00), frequencies, pulse) >= 0.80 &&
                pulse_metrics_value.broadening_ratio <= 1.25 &&
                pulse_metrics_value.pulse_correlation >= 0.90 &&
                pulse_metrics_value.postcursor_ratio <= 0.10,
        ))
    end
    write_csv(joinpath(OUTPUT_ROOT, "isolated_summary.csv"), summary_rows)
    jldsave(joinpath(OUTPUT_ROOT, "isolated_analysis.jld2");
            rows, relative, summary_rows, pulse, pulse_results)

    power_panel = plot(
        xlabel="frequency, kHz", ylabel="power fraction",
        title="Isolated full-elastic scattering", ylims=(0, 1.08),
    )
    delay_panel = plot(
        xlabel="frequency, kHz", ylabel="excess group delay, μs",
        title="Delay relative to straight footprint",
    )
    waveform_panel = plot(
        xlabel="time, μs", ylabel="relative response",
        title="Five-cycle channel response",
    )
    colors = Dict(:straight => :black, :path_only => :darkorange,
                  :hybrid => :royalblue, :pure_same_depth => :firebrick)
    for variant in VARIANTS
        current = sort(filter(row -> row.variant == variant, relative); by=row -> row.frequency_hz)
        plot!(power_panel, getproperty.(current, :frequency_hz) ./ 1e3,
              getproperty.(current, :T00); color=colors[variant], marker=:circle,
              linewidth=2, label="T $(variant)")
        plot!(delay_panel, getproperty.(current, :frequency_hz) ./ 1e3,
              getproperty.(current, :excess_group_delay_s) .* 1e6;
              color=colors[variant], marker=:circle, linewidth=2, label=String(variant))
        response = pulse_results[variant].waveform
        window = pulse.time_s .<= 45e-6
        plot!(waveform_panel, pulse.time_s[window] .* 1e6,
              response[window]; color=colors[variant], linewidth=2, label=String(variant))
    end
    summary_panel = bar(
        String.(getproperty.(summary_rows, :variant)),
        getproperty.(summary_rows, :weighted_transmission);
        ylabel="pulse-weighted transmission", title="G1 transmission gate",
        label=false, color=[:gray50, :darkorange, :royalblue, :firebrick],
    )
    hline!(summary_panel, [0.8]; linestyle=:dash, color=:black, label="gate")
    savefig(plot(
        power_panel, delay_panel, waveform_panel, summary_panel;
        layout=(2, 2), size=(1500, 950), margin=5Plots.mm,
    ), joinpath(OUTPUT_ROOT, "isolated_gate.png"))
    for row in summary_rows
        println("[+] ", row)
    end
end

function main()
    STAGE == "mesh" && return run_mesh_stage()
    STAGE == "convert" && return run_convert_stage()
    STAGE == "solve" && return run_solve_stage()
    startswith(STAGE, "solve-worker-") &&
        return run_worker(parse(Int, split(STAGE, "solve-worker-"; limit=2)[2]))
    STAGE == "analyze" && return run_analysis_stage()
    error("unknown stage: $STAGE")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
