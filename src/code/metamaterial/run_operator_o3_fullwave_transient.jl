module RunOperatorO3FullWaveTransient

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf

include(joinpath(@__DIR__, "run_operator_o3_fullwave_carrier.jl"))
include(joinpath(@__DIR__, "ndt_fullwave_transient_solver.jl"))

using .RunOperatorO3FullWaveCarrier
using .NDTFullWaveTransientSolver

export project_time_response,
       transient_windows,
       analyze_transient,
       run_campaign

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_O3_F1_OUTPUT",
    joinpath(PROJECT_ROOT, "results", "operator_space", "fullwave_transient_v0"),
)
const MESH_ROOT = get(
    ENV,
    "METAMATERIALS_O3_F1_MESH_ROOT",
    joinpath(
        PROJECT_ROOT,
        "results",
        "operator_space",
        "fullwave_carrier_commonmesh_h1p0_v0",
        "meshes",
    ),
)
const CASES = (:none, :void_r2_focus)

mesh_path(case_id) = joinpath(MESH_ROOT, "$(case_id).msh")
response_path(case_id) = joinpath(OUTPUT_ROOT, "responses", "$(case_id)_transient.jld2")

function solver_config()
    NDTFullWaveTransientConfig(
        samples_per_cycle=parse(Int, get(ENV, "METAMATERIALS_O3_F1_SAMPLES_PER_CYCLE", "30")),
        final_time_s=parse(Float64, get(ENV, "METAMATERIALS_O3_F1_FINAL_TIME_US", "70.0")) * 1.0e-6,
    )
end

function solve_case(case_id; force=false)
    case_id in CASES || throw(ArgumentError("unknown F1 case: $case_id"))
    path = response_path(case_id)
    if isfile(path) && !force
        println("[=] Existing F1 response: $path")
        return path
    end
    isfile(mesh_path(case_id)) || error("F0 common mesh not found: $(mesh_path(case_id))")
    mkpath(dirname(path))
    config = solver_config()
    println("[+] F1 direct transient solve: $case_id")
    solved = solve_ndt_element_transient_matrix(mesh_path(case_id); config)
    jldsave(
        path;
        format_version=1,
        case_id,
        mesh_path=abspath(mesh_path(case_id)),
        time_s=solved.time_s,
        source_drive=solved.source_drive,
        normal_response_m2=solved.normal_response_m2,
        tangent_response_m2=solved.tangent_response_m2,
        dt_s=solved.dt_s,
        matrix_size=solved.matrix_size,
        mass_nonzeros=solved.mass_nonzeros,
        damping_nonzeros=solved.damping_nonzeros,
        stiffness_nonzeros=solved.stiffness_nonzeros,
        effective_nonzeros=solved.effective_nonzeros,
        element_lengths_m=solved.element_lengths_m,
        solver_config=(
            frequency_hz=config.frequency_hz,
            pulse_cycles=config.pulse_cycles,
            pressure_amplitude_pa=config.pressure_amplitude_pa,
            final_time_s=config.final_time_s,
            samples_per_cycle=config.samples_per_cycle,
            density_kg_m3=config.density_kg_m3,
            pressure_wave_speed_m_s=config.pressure_wave_speed_m_s,
            shear_wave_speed_m_s=config.shear_wave_speed_m_s,
            rayleigh_alpha_s_inv=config.rayleigh_alpha_s_inv,
            rayleigh_beta_s=config.rayleigh_beta_s,
            element_count=config.element_count,
            element_order=config.element_order,
            quadrature_degree=config.quadrature_degree,
            newmark_beta=config.newmark_beta,
            newmark_gamma=config.newmark_gamma,
        ),
    )
    println("[+] Saved F1 response: $path")
    path
end

function project_time_response(response, definition)
    time_count = size(response, 3)
    first_projection = RunOperatorO3FullWaveCarrier.project_operator(
        view(response, :, :, 1), definition,
    )
    projected = Matrix{ComplexF64}(undef, length(first_projection), time_count)
    projected[:, 1] .= first_projection
    for time_index in 2:time_count
        projected[:, time_index] .= RunOperatorO3FullWaveCarrier.project_operator(
            view(response, :, :, time_index), definition,
        )
    end
    projected
end

function transient_windows(time_s, config)
    dt = time_s[2] - time_s[1]
    guard = 2dt
    defect_x_m = 60.0e-3
    defect_radius_m = 2.0e-3
    specimen_length_m = 125.0e-3
    burst_duration = config.pulse_cycles / config.frequency_hz
    earliest_void = 2(defect_x_m - defect_radius_m) / config.pressure_wave_speed_m_s
    latest_void = 2(defect_x_m + defect_radius_m) / config.pressure_wave_speed_m_s + burst_duration
    backwall_arrival = 2specimen_length_m / config.pressure_wave_speed_m_s
    (
        precursor=findall(<(max(earliest_void - guard, 0.0)), time_s),
        void=findall(time -> earliest_void - guard <= time <= latest_void + guard, time_s),
        backwall=findall(time ->
            backwall_arrival - guard <= time <= min(last(time_s), backwall_arrival + burst_duration + guard),
            time_s,
        ),
        earliest_void_s=earliest_void,
        latest_void_s=latest_void,
        backwall_arrival_s=backwall_arrival,
        guard_s=guard,
    )
end

positive_floor(value::Real) = max(Float64(value), floatmin(Float64))
energy(signal, indices, dt) = dt * sum(abs2, view(signal, :, indices))
finite_signal(signal) = all(isfinite, real.(signal)) && all(isfinite, imag.(signal))

function relative_residual(first, second)
    norm(first - second) / positive_floor(norm(first))
end

mirror_response(response) = response[end:-1:1, end:-1:1, :]

function write_namedtuple_csv(path, rows)
    isempty(rows) && return path
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((getproperty(row, column) for column in columns), ','))
        end
    end
    path
end

function analyze_transient()
    free_data = JLD2.load(response_path(:none))
    void_data = JLD2.load(response_path(:void_r2_focus))
    time_s = Float64.(free_data["time_s"])
    time_s == Float64.(void_data["time_s"]) || error("paired time grids differ")
    config = solver_config()
    windows = transient_windows(time_s, config)
    dt = time_s[2] - time_s[1]
    free = Float64.(free_data["normal_response_m2"])
    void = Float64.(void_data["normal_response_m2"])
    free_tangent = Float64.(free_data["tangent_response_m2"])
    void_tangent = Float64.(void_data["tangent_response_m2"])
    paired = void - free
    paired_tangent = void_tangent - free_tangent

    reciprocity_residual = maximum((
        relative_residual(free, permutedims(free, (2, 1, 3))),
        relative_residual(void, permutedims(void, (2, 1, 3))),
        relative_residual(paired, permutedims(paired, (2, 1, 3))),
    ))
    mirror_residual = relative_residual(paired, mirror_response(paired))
    total_paired_energy = dt * sum(abs2, paired)
    precursor_energy = dt * sum(abs2, view(paired, :, :, windows.precursor))
    void_window_energy = dt * sum(abs2, view(paired, :, :, windows.void))
    backwall_free_energy = dt * sum(abs2, view(free, :, :, windows.backwall))
    tangential_over_normal = norm(paired_tangent) / positive_floor(norm(paired))
    precursor_fraction = precursor_energy / positive_floor(total_paired_energy)

    definitions = RunOperatorO3FullWaveCarrier.frozen_definitions()
    projected = Dict{Symbol,Any}()
    for definition in definitions
        free_features = project_time_response(free, definition)
        void_features = project_time_response(void, definition)
        projected[definition.id] = (
            free=free_features,
            void=void_features,
            paired=void_features - free_features,
        )
    end
    fmc_void_energy = energy(
        projected[:paut_fmc_oracle].paired, windows.void, dt,
    )
    architecture_rows = NamedTuple[]
    for definition in definitions
        signals = projected[definition.id]
        paired_void_energy = energy(signals.paired, windows.void, dt)
        free_void_energy = energy(signals.free, windows.void, dt)
        free_backwall_energy = energy(signals.free, windows.backwall, dt)
        paired_precursor_energy = energy(signals.paired, windows.precursor, dt)
        envelope = vec(sqrt.(sum(abs2, signals.paired; dims=1)))
        local_peak = windows.void[argmax(view(envelope, windows.void))]
        push!(architecture_rows, (
            architecture_id=definition.id,
            kind=definition.kind,
            electronic_channels=definition.electronic_channels,
            physical_states=definition.physical_states,
            transmit_events=definition.transmit_events,
            complex_observables=size(signals.paired, 1),
            paired_void_energy=paired_void_energy,
            free_void_window_energy=free_void_energy,
            paired_to_free_void_db=10log10(positive_floor(
                paired_void_energy / positive_floor(free_void_energy),
            )),
            paired_energy_vs_fmc=paired_void_energy / positive_floor(fmc_void_energy),
            paired_precursor_fraction=paired_precursor_energy /
                positive_floor(dt * sum(abs2, signals.paired)),
            free_backwall_energy=free_backwall_energy,
            paired_peak_time_us=1.0e6time_s[local_peak],
        ))
    end

    common_physics =
        free_data["solver_config"] == void_data["solver_config"] &&
        isapprox(
            free_data["element_lengths_m"], void_data["element_lengths_m"];
            rtol=1.0e-12, atol=1.0e-15,
        )
    all_finite = finite_signal(free) && finite_signal(void) &&
        finite_signal(free_tangent) && finite_signal(void_tangent)
    gate_rows = [
        (gate=:all_time_samples_finite, limit="true", value=string(all_finite), passed=all_finite),
        (gate=:maximum_reciprocity_residual, limit="1e-8", value=string(reciprocity_residual), passed=reciprocity_residual <= 1.0e-8),
        (gate=:paired_mirror_residual, limit="0.03", value=string(mirror_residual), passed=mirror_residual <= 0.03),
        (gate=:paired_precursor_fraction, limit="0.01", value=string(precursor_fraction), passed=precursor_fraction <= 0.01),
        (gate=:paired_void_window_energy, limit=">0", value=string(void_window_energy), passed=isfinite(void_window_energy) && void_window_energy > 0),
        (gate=:free_backwall_window_energy, limit=">0", value=string(backwall_free_energy), passed=isfinite(backwall_free_energy) && backwall_free_energy > 0),
        (gate=:paired_tangential_response_finite, limit="true", value=string(finite_signal(paired_tangent)), passed=finite_signal(paired_tangent)),
        (gate=:common_pulse_mesh_receive_chain, limit="true", value=string(common_physics), passed=common_physics),
    ]
    gate_passed = all(getproperty.(gate_rows, :passed))

    mkpath(OUTPUT_ROOT)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "f1a_gates.csv"), gate_rows)
    write_namedtuple_csv(joinpath(OUTPUT_ROOT, "architecture_pulse_metrics.csv"), architecture_rows)
    summary_path = joinpath(OUTPUT_ROOT, "fullwave_transient_summary.jld2")
    jldsave(
        summary_path;
        format_version=1,
        contract="docs/operator_o3_f1_contract.md",
        f1a_gate_passed=gate_passed,
        time_s,
        windows,
        reciprocity_residual,
        mirror_residual,
        precursor_fraction,
        tangential_over_normal,
        total_paired_energy,
        void_window_energy,
        backwall_free_energy,
        gate_rows,
        architecture_rows,
    )

    signal_plot = plot(
        xlabel="time, us", ylabel="||paired architecture signal(t)||",
        title="F1a direct transient paired responses", legend=:topright,
    )
    for definition in definitions
        signal = projected[definition.id].paired
        envelope = vec(sqrt.(sum(abs2, signal; dims=1)))
        plot!(signal_plot, 1.0e6 .* time_s, envelope; lw=2, label=string(definition.id))
    end
    vspan!(signal_plot, 1.0e6 .* [time_s[first(windows.void)], time_s[last(windows.void)]];
        alpha=0.08, color=:green, label="void window")
    vspan!(signal_plot, 1.0e6 .* [time_s[first(windows.backwall)], time_s[last(windows.backwall)]];
        alpha=0.06, color=:orange, label="backwall window")
    plot_path = joinpath(OUTPUT_ROOT, "paired_transient_architectures.png")
    savefig(signal_plot, plot_path)

    readme_path = joinpath(OUTPUT_ROOT, "README.md")
    open(readme_path, "w") do io
        println(io, "# O3 F1a: direct paired full-wave transient")
        println(io)
        println(io, "> F1a physics gate: **$(gate_passed ? "PASS" : "FAIL")**. F1b stochastic detection is still required.")
        println(io)
        println(io, @sprintf("Five-cycle 242 kHz pulse, %.3f us step, %d samples, direct multi-RHS Newmark.", 1.0e6dt, length(time_s)))
        println(io)
        println(io, "## Physics diagnostics")
        println(io)
        println(io, @sprintf("- reciprocity residual: `%.3e`", reciprocity_residual))
        println(io, @sprintf("- paired mirror residual: `%.3f%%`", 100mirror_residual))
        println(io, @sprintf("- precursor energy fraction: `%.3f%%`", 100precursor_fraction))
        println(io, @sprintf("- paired tangential/normal norm: `%.3f`", tangential_over_normal))
        println(io, @sprintf("- geometric void window: `%.2f--%.2f us`", 1.0e6time_s[first(windows.void)], 1.0e6time_s[last(windows.void)]))
        println(io, @sprintf("- backwall window: `%.2f--%.2f us`", 1.0e6time_s[first(windows.backwall)], 1.0e6time_s[last(windows.backwall)]))
        println(io)
        println(io, "## Frozen pulse projections")
        println(io)
        println(io, "| Architecture | E | Paired/free void, dB | Paired energy / FMC | Peak, us |")
        println(io, "|---|---:|---:|---:|---:|")
        for row in architecture_rows
            println(io, @sprintf("| `%s` | %d | %.2f | %.3f | %.2f |", row.architecture_id, row.transmit_events, row.paired_to_free_void_db, row.paired_energy_vs_fmc, row.paired_peak_time_us))
        end
        println(io)
        println(io, "These are deterministic raw-pulse transfer metrics, not Pd/Pfa. F1b must add independent nuisance/noise splits and a finite defect ROI ensemble.")
    end
    println("[+] F1a gate: $(gate_passed ? "PASS" : "FAIL")")
    println("[+] $readme_path")
    (; gate_passed, gate_rows, architecture_rows, summary_path, plot_path, readme_path)
end

function run_campaign(; force=false)
    mkpath(OUTPUT_ROOT)
    solve_case(:none; force)
    solve_case(:void_r2_focus; force)
    analyze_transient()
end

function main(args=ARGS)
    if "--analyze-only" in args
        analyze_transient()
    elseif any(startswith(argument, "--solve-one=") for argument in args)
        argument = only(filter(argument -> startswith(argument, "--solve-one="), args))
        solve_case(Symbol(split(argument, '='; limit=2)[2]); force="--force" in args)
    else
        run_campaign(; force="--force" in args)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

end

