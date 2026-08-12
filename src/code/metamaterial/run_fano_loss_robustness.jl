module FanoLossRobustness

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SOURCE_ROOT = get(
    ENV,
    "METAMATERIALS_FANO_LOSS_SOURCE",
    joinpath(PROJECT_ROOT, "tmp", "fano_dark_mass_tuning"),
)
const SAMPLE_MODEL_NAME = get(
    ENV,
    "METAMATERIALS_FANO_LOSS_SAMPLE_MODEL",
    "model_bd_ld_3p06.json",
)
const SAMPLE_LABEL = get(ENV, "METAMATERIALS_FANO_LOSS_SAMPLE_LABEL", "ld=3.06 mm")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_FANO_LOSS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "fano_loss_robustness"),
)
const DAMPING_SCALES = parse.(
    Float64,
    split(get(ENV, "METAMATERIALS_FANO_DAMPING_SCALES", "0,0.25,0.5,1"), ','),
)
const FREQUENCIES_HZ = collect(
    parse(Float64, get(ENV, "METAMATERIALS_FANO_LOSS_START_HZ", "240500")):
    parse(Float64, get(ENV, "METAMATERIALS_FANO_LOSS_STEP_HZ", "125")):
    parse(Float64, get(ENV, "METAMATERIALS_FANO_LOSS_STOP_HZ", "243500")),
)
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

if !isnothing(REQUESTED_STAGE) && startswith(REQUESTED_STAGE, "solve-")
    include(joinpath(@__DIR__, "harmonic_solver.jl"))
    using .HarmonicElasticity
    using JLD2
elseif REQUESTED_STAGE == "analyze"
    ENV["GKSwstype"] = "100"
    using JLD2
    using Plots
end

scale_label(scale) = replace(string(Float64(scale)), "." => "p")
reference_model_path() = joinpath(SOURCE_ROOT, "models", "model_r0_reference.json")
sample_model_path() = joinpath(SOURCE_ROOT, "models", SAMPLE_MODEL_NAME)
result_path(kind, scale) = joinpath(
    OUTPUT_ROOT,
    "harmonic",
    "$(kind)_s$(scale_label(scale)).jld2",
)

function child_command(stage)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function solve_case(kind, scale)
    model = kind == "reference" ? reference_model_path() : sample_model_path()
    isfile(model) || error("run run_fano_dark_mass_tuning.jl first; missing $model")
    config = HarmonicConfig(
        rayleigh_alpha=Float64(scale) * 79560.0,
        rayleigh_beta=Float64(scale) * 2.5e-9,
    )
    points = solve_harmonic_sweep(model, FREQUENCIES_HZ; config)
    path = result_path(kind, scale)
    mkpath(dirname(path))
    JLD2.jldsave(
        path;
        frequency_hz=getproperty.(points, :frequency_hz),
        right_displacement=getproperty.(points, :right_displacement),
        source_power_w_per_m=getproperty.(points, :source_power_w_per_m),
        left_absorbed_power_w_per_m=getproperty.(points, :left_absorbed_power_w_per_m),
        right_absorbed_power_w_per_m=getproperty.(points, :right_absorbed_power_w_per_m),
        internal_dissipated_power_w_per_m=getproperty.(points, :internal_dissipated_power_w_per_m),
    )
    println("[+] $path")
end

function run_single_stage(specification)
    columns = split(specification, '-')
    length(columns) == 2 || error("expected kind-index")
    kind = columns[1]
    kind in ("reference", "sample") || error("unknown kind: $kind")
    index = parse(Int, columns[2])
    1 <= index <= length(DAMPING_SCALES) || error("invalid scale index")
    solve_case(kind, DAMPING_SCALES[index])
end

function run_solve_stage()
    stages = [
        "solve-$kind-$index"
        for index in eachindex(DAMPING_SCALES)
        for kind in ("reference", "sample")
    ]
    semaphore = Base.Semaphore(min(4, length(stages), Sys.CPU_THREADS))
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

function unwrap_phase(values)
    raw = angle.(values)
    result = Float64[first(raw)]
    for index in 2:length(raw)
        push!(result, last(result) + mod(raw[index] - raw[index - 1] + pi, 2pi) - pi)
    end
    result
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

function run_analysis_stage()
    rows = NamedTuple[]
    summary = NamedTuple[]
    amplitude_plot = plot(
        xlabel="frequency, kHz",
        ylabel="reference-calibrated |H|",
        title="Fano line for $(SAMPLE_LABEL) under Rayleigh damping",
        gridalpha=0.25,
        legend=:outertopright,
    )
    power_plot = plot(
        xlabel="frequency, kHz",
        ylabel="Pout / Pout,R0",
        gridalpha=0.25,
        legend=:outertopright,
    )
    colors = [:royalblue, :darkorange, :seagreen, :firebrick]
    for (scale_index, scale) in enumerate(DAMPING_SCALES)
        reference = JLD2.load(result_path("reference", scale))
        sample = JLD2.load(result_path("sample", scale))
        response = sample["right_displacement"] ./ reference["right_displacement"]
        amplitude = abs.(response)
        power_ratio = sample["right_absorbed_power_w_per_m"] ./
                      reference["right_absorbed_power_w_per_m"]
        phase = unwrap_phase(response)
        peak_index = argmax(amplitude)
        dip_index = argmin(amplitude)
        balance = sample["source_power_w_per_m"] .-
                  sample["left_absorbed_power_w_per_m"] .-
                  sample["right_absorbed_power_w_per_m"] .-
                  sample["internal_dissipated_power_w_per_m"]
        push!(summary, (
            damping_scale=scale,
            peak_frequency_hz=FREQUENCIES_HZ[peak_index],
            peak_amplitude=amplitude[peak_index],
            peak_power_ratio=power_ratio[peak_index],
            dip_frequency_hz=FREQUENCIES_HZ[dip_index],
            dip_amplitude=amplitude[dip_index],
            dip_power_ratio=power_ratio[dip_index],
            peak_to_dip_contrast=amplitude[peak_index] / amplitude[dip_index],
            max_relative_power_balance_error=maximum(abs.(balance) ./
                max.(abs.(sample["source_power_w_per_m"]), eps())),
        ))
        for frequency_index in eachindex(FREQUENCIES_HZ)
            push!(rows, (
                damping_scale=scale,
                frequency_hz=FREQUENCIES_HZ[frequency_index],
                amplitude=amplitude[frequency_index],
                phase_rad=phase[frequency_index],
                transmitted_power_ratio=power_ratio[frequency_index],
                internal_dissipated_power_w_per_m=
                    sample["internal_dissipated_power_w_per_m"][frequency_index],
                power_balance_error_w_per_m=balance[frequency_index],
            ))
        end
        label = "s=$(scale)"
        plot!(
            amplitude_plot,
            FREQUENCIES_HZ ./ 1e3,
            amplitude;
            color=colors[scale_index],
            linewidth=2.3,
            label,
        )
        plot!(
            power_plot,
            FREQUENCIES_HZ ./ 1e3,
            power_ratio;
            color=colors[scale_index],
            linewidth=2.3,
            label,
        )
    end
    mkpath(OUTPUT_ROOT)
    write_csv(joinpath(OUTPUT_ROOT, "fano_loss_response.csv"), rows)
    write_csv(joinpath(OUTPUT_ROOT, "fano_loss_summary.csv"), summary)
    figure = plot(
        amplitude_plot,
        power_plot;
        layout=(2, 1),
        size=(1150, 900),
        margin=5Plots.mm,
    )
    path = joinpath(OUTPUT_ROOT, "fano_loss_robustness.png")
    savefig(figure, path)
    println("[+] $path")
    foreach(println, summary)
end

function run_child_stage(stage)
    if stage == "solve"
        run_solve_stage()
    elseif startswith(stage, "solve-")
        run_single_stage(stage[7:end])
    elseif stage == "analyze"
        run_analysis_stage()
    else
        error("unknown loss-robustness stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("solve", "analyze")
            println("\n=== Fano loss stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_fano_loss_robustness.jl [--stage=...]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
