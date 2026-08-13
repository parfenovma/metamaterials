module TopologyBendingMetalGRINRun

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV, "METAMATERIALS_BENDING_METAL_GRIN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "bending_metal_grin"),
)
const MAX_WORKERS = parse(Int, get(ENV, "METAMATERIALS_BENDING_METAL_GRIN_WORKERS", "3"))
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "topology_bending_metal_grin_mesher.jl"))
using .TopologyBendingMetalGRINMesher
include(joinpath(@__DIR__, "ndt_benchmark.jl"))
using .NDTBenchmark

const CONFIG = BendingMetalGRINConfig()
const VARIANTS = BENDING_METAL_GRIN_VARIANTS
const COARSE_FREQUENCIES_HZ = collect(162.0e3:20.0e3:322.0e3)
const FINE_FREQUENCIES_HZ = collect(162.0e3:4.0e3:322.0e3)
const FINE_VARIANTS_FILE = joinpath(OUTPUT_ROOT, "fine_variants.txt")

is_solver_stage(stage) = startswith(stage, "coarse-worker-") ||
                         startswith(stage, "fine-worker-") ||
                         stage in ("coarse", "fine", "analyze-coarse", "analyze-fine", "chain-check")

if REQUESTED_STAGE in ("mesh", "chain-mesh", "mesh-convergence")
    using Gmsh: gmsh
elseif REQUESTED_STAGE in ("convert", "chain-convert", "convert-convergence")
    include(joinpath(@__DIR__, "step1b_convert_models.jl"))
elseif !isnothing(REQUESTED_STAGE) && is_solver_stage(REQUESTED_STAGE)
    include(joinpath(@__DIR__, "modal_harmonic_solver.jl"))
    include(joinpath(@__DIR__, "modal_projection.jl"))
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    using .ModalHarmonicElasticity
    using .ElasticModalProjection
    using .ImpulseRiskAnalysis
    using FFTW
    using JLD2
end

mesh_path(variant) = joinpath(OUTPUT_ROOT, "meshes", "mesh_$(variant).msh")
model_path(variant) = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
point_path(stage, variant, frequency_hz) = joinpath(
    OUTPUT_ROOT, String(stage), "$(variant)_$(round(Int, frequency_hz))hz.jld2",
)

function write_csv(path, rows)
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, c)) for c in columns), ','))
        end
    end
    println("[+] $path")
    path
end

required_center_delay_s() =
    (hypot(60.0e-3, 71.6e-3 / 2) - 60.0e-3) / 6122.102437409232

function approximate_geometry_metrics(config, variant)
    variant == :solid && return (
        ligament_mm=config.port_height_mm,
        mass_width_mm=config.period_mm,
        mass_height_mm=config.port_height_mm,
        corner_mass_mm=config.port_height_mm,
        elbow_y_mm=0.0,
        relative_density=1.0,
    )
    p = variant_parameters(config, variant)
    # Exact OCC area is evaluated once per state here.  It is only a G0
    # descriptor; no effective wave property is assigned from porosity.
    (; p..., relative_density=relative_areal_density(config, variant))
end

function run_g0_stage()
    rows = [begin
        g = approximate_geometry_metrics(CONFIG, variant)
        (
            variant,
            total_length_mm=total_length_mm(CONFIG),
            active_length_mm=CONFIG.period_count * CONFIG.period_mm,
            period_mm=CONFIG.period_mm,
            period_count=CONFIG.period_count,
            port_height_mm=CONFIG.port_height_mm,
            ligament_mm=g.ligament_mm,
            mass_width_mm=g.mass_width_mm,
            mass_height_mm=g.mass_height_mm,
            corner_mass_mm=g.corner_mass_mm,
            elbow_y_mm=g.elbow_y_mm,
            relative_areal_density=g.relative_density,
            minimum_feature_mm=minimum_feature_mm(CONFIG, variant),
            required_center_delay_us=required_center_delay_s() * 1e6,
            cell_to_p_wavelength_at_322khz=CONFIG.period_mm /
                (6122.102437409232 / 322e3 * 1e3),
            named_process="wire_EDM_Al6061_plate_stack",
            plate_stack_thickness_mm=CONFIG.plate_stack_thickness_mm,
        )
    end for variant in VARIANTS]
    write_csv(joinpath(OUTPUT_ROOT, "g0_manufacturing.csv"), rows)
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            build_bending_metal_grin_mesh(mesh_path(variant); config=CONFIG, variant)
            println("[+] $variant mesh -> $(mesh_path(variant))")
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
    right_modes = filter(m -> m.kind == :propagating && m.direction == :right, state.modes)
    [begin
        left_candidates = filter(
            m -> m.kind == :propagating && m.direction == :left && m.parity == rm.parity,
            state.modes,
        )
        lm = first(sort(left_candidates; by=m -> abs(m.wavenumber_per_m + rm.wavenumber_per_m)))
        projection = project_boundary_mode_pair(state, tag, rm, lm; quadrature_degree=6)
        (; right_mode=rm, left_mode=lm, projection...)
    end for rm in right_modes]
end

function solve_scattering_at(path, local_config, variant, frequency_hz)
    hconfig = ModalHarmonicElasticity.HarmonicElasticity.HarmonicConfig(
        density=2700.0,
        pressure_wave_speed=6122.102437409232,
        shear_wave_speed=3083.810277185563,
        rayleigh_alpha=0.0,
        rayleigh_beta=0.0,
        element_order=2,
        quadrature_degree=4,
    )
    state = solve_symmetry_reduced_modal_state(
        path, frequency_hz;
        config=hconfig,
        port_height_m=local_config.port_height_mm * 1e-3,
        port_element_count=110,
    )
    left = project_all(state, "Source")
    right = project_all(state, "Microphone")
    input_index = findfirst(pair -> pair.right_mode.parity == :symmetric, left)
    input_index === nothing && error("symmetric incident mode not found")
    li = left[input_index]
    ri = right[input_index]
    incident = abs2(li.right_amplitude)
    reflection = abs2(li.left_amplitude) / incident
    transmission = abs2(ri.right_amplitude) / incident
    conversion = sum(
        abs2(lp.left_amplitude) + abs2(rp.right_amplitude)
        for (index, (lp, rp)) in enumerate(zip(left, right)) if index != input_index
    ) / incident
    s21 = ri.right_amplitude / li.right_amplitude
    s11 = li.left_amplitude / li.right_amplitude
    (
        variant,
        frequency_hz=Float64(frequency_hz),
        T00=transmission,
        R00=reflection,
        C0=conversion,
        S11_complex=ComplexF64(s11),
        S21_complex=ComplexF64(s21),
        transmission_amplitude=ComplexF64(s21),
        propagating_mode_count=length(left),
        power_sum=transmission + reflection + conversion,
        balance_residual=state.metrics.balance_residual_w_per_m / incident,
    )
end

solve_scattering(variant, frequency_hz) =
    solve_scattering_at(model_path(variant), CONFIG, variant, frequency_hz)

function fine_variants()
    isfile(FINE_VARIANTS_FILE) || return (:solid, :bd_a, :bd_b, :bd_c, :bd_d)
    names = Symbol.(filter(!isempty, strip.(readlines(FINE_VARIANTS_FILE))))
    :solid in names || pushfirst!(names, :solid)
    Tuple(unique(names))
end

function stage_data(stage::Symbol)
    stage == :coarse && return VARIANTS, COARSE_FREQUENCIES_HZ
    stage == :fine && return fine_variants(), FINE_FREQUENCIES_HZ
    throw(ArgumentError("unknown stage $stage"))
end

specs(stage) = [(; variant, frequency_hz=f) for variant in stage_data(stage)[1]
                for f in stage_data(stage)[2]]

function solve_and_save(stage, spec)
    row = solve_scattering(spec.variant, spec.frequency_hz)
    path = point_path(stage, spec.variant, spec.frequency_hz)
    mkpath(dirname(path))
    JLD2.jldsave(path; row)
    println(
        "[+] $stage $(spec.variant) $(round(spec.frequency_hz / 1e3; digits=1)) kHz",
        ": T=", round(row.T00; digits=4),
        ", R=", round(row.R00; digits=4),
        ", C=", round(row.C0; digits=5),
        ", sum=", round(row.power_sum; digits=4),
    )
end

function run_worker(stage, worker_index)
    jobs = specs(stage)
    for index in worker_index:MAX_WORKERS:length(jobs)
        spec = jobs[index]
        path = point_path(stage, spec.variant, spec.frequency_hz)
        if isfile(path)
            stored = JLD2.load(path)["row"]
            :S11_complex in propertynames(stored) && continue
        end
        solve_and_save(stage, spec)
    end
end

function child_command(stage)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --stage=$stage`
end

function run_parallel_stage(stage::Symbol)
    @sync for worker in 1:MAX_WORKERS
        @async run(child_command("$(stage)-worker-$worker"))
    end
end

function unwrap_phase(values)
    out = Float64[first(values)]
    for raw in Iterators.drop(values, 1)
        push!(out, last(out) + mod(raw - last(out) + pi, 2pi) - pi)
    end
    out
end

function interp(x, y, q)
    q <= first(x) && return first(y)
    q >= last(x) && return last(y)
    ir = searchsortedfirst(x, q)
    il = ir - 1
    a = (q - x[il]) / (x[ir] - x[il])
    (1 - a) * y[il] + a * y[ir]
end

load_raw(stage) = [JLD2.load(point_path(stage, s.variant, s.frequency_hz))["row"]
                   for s in specs(stage)]

function calibrated_rows(raw, variant)
    rows = sort(filter(r -> r.variant == variant, raw); by=r -> r.frequency_hz)
    refs = Dict(r.frequency_hz => r for r in raw if r.variant == :solid)
    transfer = ComplexF64[
        r.transmission_amplitude / refs[r.frequency_hz].transmission_amplitude for r in rows
    ]
    phase = unwrap_phase(angle.(transfer))
    omega = 2pi .* getproperty.(rows, :frequency_hz)
    delay = similar(phase)
    delay[1] = -(phase[2] - phase[1]) / (omega[2] - omega[1])
    for i in 2:length(phase)-1
        delay[i] = -(phase[i + 1] - phase[i - 1]) / (omega[i + 1] - omega[i - 1])
    end
    delay[end] = -(phase[end] - phase[end - 1]) / (omega[end] - omega[end - 1])
    [merge(r, (
        relative_amplitude=abs(transfer[i]),
        relative_phase_rad=phase[i],
        excess_group_delay_s=delay[i],
        trusted_excess_group_delay_s=all(
            rows[j].T00 >= 0.5 for j in max(1, i - 1):min(length(rows), i + 1)
        ) ? delay[i] : NaN,
    )) for (i, r) in enumerate(rows)]
end

function pulse_metrics_for_variant(rows)
    spectrum = pulse_spectrum(PulseConfig(
        center_frequency_hz=242e3, cycles=5.0,
        samples_per_period=80, fft_length=65536,
    ))
    band = spectral_band(spectrum, -20.0)
    f = getproperty.(rows, :frequency_hz)
    a = getproperty.(rows, :relative_amplitude)
    p = getproperty.(rows, :relative_phase_rad)
    transfer = zeros(ComplexF64, length(spectrum.spectrum))
    refspec = zeros(ComplexF64, length(spectrum.spectrum))
    for i in band.lower_index:band.upper_index
        fq = spectrum.frequency_hz[i]
        transfer[i] = interp(f, a, fq) * cis(interp(f, p, fq))
        refspec[i] = spectrum.spectrum[i]
    end
    reference = irfft(refspec, length(spectrum.signal))
    candidate = irfft(refspec .* transfer, length(spectrum.signal))
    pm = pulse_metrics(candidate, reference, reference, spectrum.time_s[2] - spectrum.time_s[1])
    w = abs2.(spectrum.spectrum[band.lower_index:band.upper_index])
    sf = spectrum.frequency_hz[band.lower_index:band.upper_index]
    weighted(column) = sum(w[i] * interp(f, getproperty.(rows, column), sf[i])
                           for i in eachindex(w)) / sum(w)
    (
        relative_peak=pm.gain_peak,
        pulse_broadening_ratio=pm.broadening_ratio,
        pulse_correlation=pm.pulse_correlation,
        postcursor_ratio=pm.postcursor_ratio,
        weighted_transmission=weighted(:T00),
        weighted_reflection=weighted(:R00),
        weighted_mode_conversion=weighted(:C0),
    )
end

function summary_rows(stage, all_rows)
    variants, _ = stage_data(stage)
    [begin
        rows = sort(filter(r -> r.variant == variant, all_rows); by=r -> r.frequency_hz)
        pulse = pulse_metrics_for_variant(rows)
        f = getproperty.(rows, :frequency_hz)
        delay = interp(f, getproperty.(rows, :trusted_excess_group_delay_s), 242e3)
        phase = interp(f, getproperty.(rows, :relative_phase_rad), 242e3)
        p = approximate_geometry_metrics(CONFIG, variant)
        pass_transfer = pulse.weighted_transmission >= 0.80
        pass_pulse = pulse.pulse_broadening_ratio <= 1.25 &&
                     pulse.pulse_correlation >= 0.90 && pulse.postcursor_ratio <= 0.10
        (
            variant,
            ligament_mm=p.ligament_mm,
            mass_width_mm=p.mass_width_mm,
            mass_height_mm=p.mass_height_mm,
            relative_areal_density=p.relative_density,
            minimum_feature_mm=minimum_feature_mm(CONFIG, variant),
            min_T00=minimum(getproperty.(rows, :T00)),
            max_R00=maximum(getproperty.(rows, :R00)),
            weighted_transmission=pulse.weighted_transmission,
            weighted_reflection=pulse.weighted_reflection,
            weighted_mode_conversion=pulse.weighted_mode_conversion,
            relative_peak=pulse.relative_peak,
            pulse_broadening_ratio=pulse.pulse_broadening_ratio,
            pulse_correlation=pulse.pulse_correlation,
            postcursor_ratio=pulse.postcursor_ratio,
            carrier_relative_phase_rad=phase,
            delay_at_242_us=delay * 1e6,
            required_center_delay_us=required_center_delay_s() * 1e6,
            max_power_sum_error=maximum(abs.(getproperty.(rows, :power_sum) .- 1)),
            pass_transfer,
            pass_pulse,
            pass_local=pass_transfer && pass_pulse && isfinite(delay),
        )
    end for variant in variants]
end

function choose_fine(summary)
    candidates = filter(r -> r.variant != :solid &&
                             r.weighted_transmission >= 0.55 &&
                             r.pulse_correlation >= 0.75, summary)
    ranked = sort(candidates; by=r -> -abs(r.carrier_relative_phase_rad))
    chosen = [:solid]
    append!(chosen, getproperty.(first(ranked, min(5, length(ranked))), :variant))
    length(chosen) == 1 && append!(chosen, (:bd_a, :bd_b, :bd_c))
    open(FINE_VARIANTS_FILE, "w") do io
        foreach(v -> println(io, v), chosen)
    end
    println("[+] fine variants: ", join(chosen, ", "))
end

function write_benchmark(summary)
    byv = Dict(r.variant => r for r in summary)
    passing = filter(r -> r.variant != :solid && r.pass_local, summary)
    selected = isempty(passing) ? first(filter(r -> r.variant != :solid, summary)) :
               last(sort(passing; by=r -> r.delay_at_242_us))
    straight = byv[:solid]
    nan = NaN
    row(topology, reference, state, verdict) = (
        topology, case_id="none", reference, frequency_hz=242e3,
        input_energy_j_per_m=nan,
        echo_mf_gain_bare=nan, echo_energy_gain_bare=nan,
        echo_to_ringing_db=nan, equal_noise_snr_proxy_db=nan,
        dead_zone_depth_mm=nan, focus_peak_gain_bare=nan,
        focus_energy_gain_bare=nan, fwhm_mm=nan, dof_minus3db_mm=nan,
        sidelobe_ratio=nan,
        pulse_broadening_ratio=state === nothing ? nan : state.pulse_broadening_ratio,
        pulse_correlation=state === nothing ? nan : state.pulse_correlation,
        postcursor_ratio=state === nothing ? nan : state.postcursor_ratio,
        weighted_transmission=state === nothing ? nan : state.weighted_transmission,
        weighted_reflection=state === nothing ? nan : state.weighted_reflection,
        weighted_mode_conversion=state === nothing ? nan : state.weighted_mode_conversion,
        minimum_feature_mm=state === nothing ? nan : state.minimum_feature_mm,
        verdict,
    )
    rows = [
        row("bending_metal_grin", "bare", nothing, "not_calculated_g1_gate"),
        row("bending_metal_grin", "straight_footprint", straight, "local_reference_only"),
        row("bending_metal_grin_$(selected.variant)", "straight_footprint", selected, "g1_pending_library_gate"),
    ]
    propertynames(first(rows)) == metric_columns() || error("shared metric columns changed")
    write_csv(joinpath(OUTPUT_ROOT, "benchmark_metrics.csv"), rows)
end

function run_analysis_stage(stage::Symbol)
    raw = load_raw(stage)
    all_rows = reduce(vcat, (calibrated_rows(raw, v) for v in stage_data(stage)[1]))
    write_csv(joinpath(OUTPUT_ROOT, "g1_$(stage)_frequency.csv"), all_rows)
    summary = summary_rows(stage, all_rows)
    write_csv(joinpath(OUTPUT_ROOT, "g1_$(stage)_summary.csv"), summary)
    stage == :coarse && choose_fine(summary)
    stage == :fine && write_benchmark(summary)
    for r in summary
        println(
            "[+] $(r.variant): Tbar=", round(r.weighted_transmission; digits=3),
            ", delay=", round(r.delay_at_242_us; digits=3), " us",
            ", Bt=", round(r.pulse_broadening_ratio; digits=3),
            ", rho=", round(r.pulse_correlation; digits=3),
            ", local=", r.pass_local,
        )
    end
end

const CHAIN_COUNTS = (1, 3, 6)
const CHAIN_FREQUENCIES_HZ = (182e3, 222e3, 242e3, 262e3, 282e3, 302e3, 322e3)

chain_config(count) = BendingMetalGRINConfig(period_count=count)
chain_mesh_path(count, variant) = joinpath(
    OUTPUT_ROOT, "chain", "meshes", "mesh_n$(count)_$(variant).msh",
)
chain_model_path(count, variant) = joinpath(
    OUTPUT_ROOT, "chain", "models", "model_n$(count)_$(variant).json",
)

function run_chain_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for count in CHAIN_COUNTS, variant in (:solid, :dense_c)
            build_bending_metal_grin_mesh(
                chain_mesh_path(count, variant);
                config=chain_config(count), variant,
            )
            println("[+] chain n=$count $variant mesh")
        end
    finally
        gmsh.finalize()
    end
end

function run_chain_convert_stage()
    for count in CHAIN_COUNTS, variant in (:solid, :dense_c)
        convert_mesh(
            chain_mesh_path(count, variant);
            output_dir=dirname(chain_model_path(count, variant)),
        )
    end
end

function unwrap_complex_phase(values)
    unwrap_phase(angle.(values))
end

function run_chain_check_stage()
    raw = NamedTuple[]
    for count in CHAIN_COUNTS, frequency_hz in CHAIN_FREQUENCIES_HZ
        config = chain_config(count)
        solid = solve_scattering_at(
            chain_model_path(count, :solid), config, :solid, frequency_hz,
        )
        dense = solve_scattering_at(
            chain_model_path(count, :dense_c), config, :dense_c, frequency_hz,
        )
        push!(raw, (
            cell_count=count,
            frequency_hz,
            solid_T=solid.T00,
            dense_T=dense.T00,
            dense_R=dense.R00,
            solid_s11=solid.S11_complex,
            solid_s21=solid.S21_complex,
            dense_s11=dense.S11_complex,
            dense_s21=dense.S21_complex,
            relative_s21=dense.S21_complex / solid.S21_complex,
            dense_power_sum=dense.power_sum,
        ))
        println(
            "[+] chain n=$count $(frequency_hz / 1e3) kHz: T=",
            round(dense.T00; digits=4), ", R=", round(dense.R00; digits=4),
        )
    end
    rows = NamedTuple[]
    for count in CHAIN_COUNTS
        selected = sort(filter(r -> r.cell_count == count, raw); by=r -> r.frequency_hz)
        phase = unwrap_complex_phase(getproperty.(selected, :relative_s21))
        for (index, r) in enumerate(selected)
            # The unit-cell transfer eigenphase is evaluated only for N=1.
            # `acos(trace(T)/2)` becomes complex in a stop band; its imaginary
            # part is retained rather than forcing a propagating branch.
            trace_half_dense = (1 - r.dense_s11^2 + r.dense_s21^2) / (2r.dense_s21)
            trace_half_solid = (1 - r.solid_s11^2 + r.solid_s21^2) / (2r.solid_s21)
            q_dense = acos(trace_half_dense)
            q_solid = acos(trace_half_solid)
            push!(rows, (
                cell_count=count,
                frequency_hz=r.frequency_hz,
                dense_T=r.dense_T,
                dense_R=r.dense_R,
                relative_phase_rad=phase[index],
                relative_phase_per_cell_rad=phase[index] / count,
                unit_transfer_q_real_rad=count == 1 ? real(q_dense) : NaN,
                unit_transfer_q_imag_rad=count == 1 ? imag(q_dense) : NaN,
                unit_transfer_excess_q_real_rad=count == 1 ? real(q_dense - q_solid) : NaN,
                unit_transfer_excess_q_imag_rad=count == 1 ? imag(q_dense - q_solid) : NaN,
                dense_power_sum=r.dense_power_sum,
            ))
        end
    end
    write_csv(joinpath(OUTPUT_ROOT, "g1_chain_bloch_consistency.csv"), rows)

    trusted = filter(r -> r.dense_T >= 0.5, rows)
    checks = NamedTuple[]
    for frequency_hz in CHAIN_FREQUENCIES_HZ
        points = filter(r -> r.frequency_hz == frequency_hz && r.dense_T >= 0.5, trusted)
        phases = getproperty.(points, :relative_phase_per_cell_rad)
        push!(checks, (
            frequency_hz,
            trusted_chain_count=length(points),
            per_cell_phase_min_rad=isempty(phases) ? NaN : minimum(phases),
            per_cell_phase_max_rad=isempty(phases) ? NaN : maximum(phases),
            per_cell_phase_spread_rad=length(phases) < 2 ? NaN : maximum(phases) - minimum(phases),
            finite_chain_consistent=length(phases) >= 2 && maximum(phases) - minimum(phases) <= 0.15,
        ))
    end
    write_csv(joinpath(OUTPUT_ROOT, "g1_chain_consistency_summary.csv"), checks)
end

function run_child(stage)
    stage == "g0" && return run_g0_stage()
    stage == "mesh" && return run_mesh_stage()
    stage == "convert" && return run_convert_stage()
    stage == "chain-mesh" && return run_chain_mesh_stage()
    stage == "chain-convert" && return run_chain_convert_stage()
    stage == "chain-check" && return run_chain_check_stage()
    stage in ("coarse", "fine") && return run_parallel_stage(Symbol(stage))
    if startswith(stage, "coarse-worker-") || startswith(stage, "fine-worker-")
        solve_stage = startswith(stage, "coarse") ? :coarse : :fine
        worker = parse(Int, split(stage, "worker-"; limit=2)[2])
        return run_worker(solve_stage, worker)
    end
    stage == "analyze-coarse" && return run_analysis_stage(:coarse)
    stage == "analyze-fine" && return run_analysis_stage(:fine)
    error("unknown bending-metal-GRIN stage: $stage")
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("g0", "mesh", "convert", "coarse", "analyze-coarse")
            println("\n=== Bending-metal-GRIN stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child(REQUESTED_STAGE)
    else
        error("usage: topology_bending_metal_grin_run.jl [--stage=STAGE]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
