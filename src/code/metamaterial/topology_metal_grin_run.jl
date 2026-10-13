module TopologyMetalGRINRun

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_METAL_GRIN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "metal_grin"),
)
const MAX_WORKERS = parse(Int, get(ENV, "METAMATERIALS_METAL_GRIN_WORKERS", "3"))
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing

include(joinpath(@__DIR__, "topology_metal_grin_mesher.jl"))
using .TopologyMetalGRINMesher
include(joinpath(@__DIR__, "ndt_benchmark.jl"))
using .NDTBenchmark

const CONFIG = MetalGRINCellConfig()
const VARIANTS = METAL_GRIN_VARIANTS
const COARSE_FREQUENCIES_HZ = collect(162.0e3:20.0e3:322.0e3)
const FINE_FREQUENCIES_HZ = collect(162.0e3:4.0e3:322.0e3)
# Fine follow-up is restricted to states that bracket the useful-delay region
# after the coarse screen.  The solid state is always retained for calibration.
const FINE_VARIANTS = (:solid, :neck_2p0, :neck_1p6, :neck_1p2)

is_solver_stage(stage) = startswith(stage, "coarse-worker-") ||
                         startswith(stage, "fine-worker-") ||
                         stage in ("coarse", "fine", "analyze-coarse", "analyze-fine")

if REQUESTED_STAGE == "mesh"
    using Gmsh: gmsh
elseif REQUESTED_STAGE == "convert"
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
    OUTPUT_ROOT,
    String(stage),
    "$(variant)_$(round(Int, frequency_hz))hz.jld2",
)

function write_csv(path, rows)
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, column)) for column in columns), ','))
        end
    end
    println("[+] $path")
    path
end

function quasi_static_area_surrogate(config, variant; samples=20001)
    x = range(0.0, config.period_mm; length=samples)
    # Sample inside an arbitrary active cell.  Endpoints have zero measure.
    height = [
        local_height_mm(config, variant, config.lead_length_mm + value)
        for value in x
    ]
    area_ratio = height ./ config.port_height_mm
    mean_area = sum(area_ratio) / length(area_ratio)
    mean_inverse_area = sum(inv, area_ratio) / length(area_ratio)
    speed_ratio = inv(sqrt(mean_area * mean_inverse_area))
    impedance_ratio = sqrt(mean_area / mean_inverse_area)
    (
        mean_area_ratio=mean_area,
        mean_inverse_area_ratio=mean_inverse_area,
        speed_ratio,
        impedance_ratio,
    )
end

function required_center_delay_s()
    radius_m = 71.6e-3 / 2
    focus_m = 60.0e-3
    speed_m_s = 6122.102437409232
    (hypot(focus_m, radius_m) - focus_m) / speed_m_s
end

function run_g0_stage()
    cp = 6122.102437409232
    rows = [begin
        surrogate = quasi_static_area_surrogate(CONFIG, variant)
        delay = CONFIG.active_length_mm * 1.0e-3 / cp *
                (inv(surrogate.speed_ratio) - 1)
        (
            variant,
            total_length_mm=total_length_mm(CONFIG),
            active_length_mm=CONFIG.active_length_mm,
            period_mm=CONFIG.period_mm,
            period_count=CONFIG.period_count,
            transition_fraction=CONFIG.transition_fraction,
            neck_height_mm=neck_height_mm(CONFIG, variant),
            minimum_feature_mm=minimum_feature_mm(CONFIG, variant),
            minimum_curvature_radius_mm=minimum_curvature_radius_mm(CONFIG, variant),
            quasi_static_speed_ratio=surrogate.speed_ratio,
            quasi_static_impedance_ratio=surrogate.impedance_ratio,
            quasi_static_excess_delay_us=delay * 1.0e6,
            required_center_delay_us=required_center_delay_s() * 1.0e6,
            cell_to_wavelength_at_322khz=CONFIG.period_mm /
                (cp / 322.0e3 * 1.0e3),
        )
    end for variant in VARIANTS]
    write_csv(joinpath(OUTPUT_ROOT, "g0_manufacturing_and_surrogate.csv"), rows)
end

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for variant in VARIANTS
            path = build_metal_grin_cell_mesh(mesh_path(variant); config=CONFIG, variant)
            println("[+] $variant mesh -> $path")
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
    [begin
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
    end for right_mode in right_modes]
end

function solve_scattering(variant, frequency_hz)
    harmonic_config = ModalHarmonicElasticity.HarmonicElasticity.HarmonicConfig(
        density=2700.0,
        pressure_wave_speed=6122.102437409232,
        shear_wave_speed=3083.810277185563,
        rayleigh_alpha=0.0,
        rayleigh_beta=0.0,
        element_order=2,
        quadrature_degree=4,
    )
    state = solve_symmetry_reduced_modal_state(
        model_path(variant),
        frequency_hz;
        config=harmonic_config,
        port_height_m=CONFIG.port_height_mm * 1.0e-3,
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
    conversion = sum(
        abs2(left_pair.left_amplitude) + abs2(right_pair.right_amplitude)
        for (index, (left_pair, right_pair)) in enumerate(zip(left, right))
        if index != input_index
    ) / incident
    transmission_amplitude = input_right.right_amplitude / input_left.right_amplitude
    reflection_amplitude = input_left.left_amplitude / input_left.right_amplitude
    (
        variant,
        frequency_hz=Float64(frequency_hz),
        T00=transmission,
        R00=reflection,
        C0=conversion,
        S11_complex=ComplexF64(reflection_amplitude),
        S21_complex=ComplexF64(transmission_amplitude),
        transmission_amplitude=ComplexF64(transmission_amplitude),
        propagating_mode_count=length(left),
        power_sum=transmission + reflection + conversion,
        balance_residual=state.metrics.balance_residual_w_per_m / incident,
    )
end

function stage_data(stage::Symbol)
    if stage == :coarse
        return VARIANTS, COARSE_FREQUENCIES_HZ
    elseif stage == :fine
        return FINE_VARIANTS, FINE_FREQUENCIES_HZ
    end
    throw(ArgumentError("unknown solve stage"))
end

function specs(stage::Symbol)
    variants, frequencies = stage_data(stage)
    [(; variant, frequency_hz) for variant in variants for frequency_hz in frequencies]
end

function solve_and_save(stage, spec)
    row = solve_scattering(spec.variant, spec.frequency_hz)
    path = point_path(stage, spec.variant, spec.frequency_hz)
    mkpath(dirname(path))
    JLD2.jldsave(path; row)
    println(
        "[+] $stage $(spec.variant) ", round(spec.frequency_hz / 1e3; digits=1),
        " kHz: T=", round(row.T00; digits=4),
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
    any(!isfile(point_path(stage, spec.variant, spec.frequency_hz)) for spec in specs(stage)) || return
    @sync for worker_index in 1:MAX_WORKERS
        @async run(child_command("$(stage)-worker-$worker_index"))
    end
end

function unwrap_phase(values)
    result = Float64[first(values)]
    for raw in Iterators.drop(values, 1)
        push!(result, last(result) + mod(raw - last(result) + pi, 2pi) - pi)
    end
    result
end

function linear_interpolate(x, values, query)
    query <= first(x) && return first(values)
    query >= last(x) && return last(values)
    right = searchsortedfirst(x, query)
    left = right - 1
    fraction = (query - x[left]) / (x[right] - x[left])
    (1 - fraction) * values[left] + fraction * values[right]
end

function load_raw_rows(stage::Symbol)
    [JLD2.load(point_path(stage, spec.variant, spec.frequency_hz))["row"] for spec in specs(stage)]
end

function calibrated_rows(raw, variant)
    selected = sort(filter(row -> row.variant == variant, raw); by=row -> row.frequency_hz)
    references = Dict(row.frequency_hz => row for row in raw if row.variant == :solid)
    transfer = ComplexF64[
        row.transmission_amplitude / references[row.frequency_hz].transmission_amplitude
        for row in selected
    ]
    phase = unwrap_phase(angle.(transfer))
    omega = 2pi .* getproperty.(selected, :frequency_hz)
    delay = similar(phase)
    delay[1] = -(phase[2] - phase[1]) / (omega[2] - omega[1])
    for index in 2:(length(phase) - 1)
        delay[index] = -(phase[index + 1] - phase[index - 1]) /
                       (omega[index + 1] - omega[index - 1])
    end
    delay[end] = -(phase[end] - phase[end - 1]) / (omega[end] - omega[end - 1])
    [merge(row, (
        relative_amplitude=abs(transfer[index]),
        relative_phase_rad=phase[index],
        excess_group_delay_s=delay[index],
        trusted_excess_group_delay_s=all(
            selected[neighbor].T00 >= 0.5
            for neighbor in max(1, index - 1):min(length(selected), index + 1)
        ) ? delay[index] : NaN,
    )) for (index, row) in enumerate(selected)]
end

function pulse_metrics_for_variant(rows)
    spectrum = pulse_spectrum(PulseConfig(
        center_frequency_hz=242.0e3,
        cycles=5.0,
        samples_per_period=80,
        fft_length=65536,
    ))
    band = spectral_band(spectrum, -20.0)
    frequencies = getproperty.(rows, :frequency_hz)
    amplitude = getproperty.(rows, :relative_amplitude)
    phase = getproperty.(rows, :relative_phase_rad)
    transfer = zeros(ComplexF64, length(spectrum.spectrum))
    reference_spectrum = zeros(ComplexF64, length(spectrum.spectrum))
    for index in band.lower_index:band.upper_index
        frequency = spectrum.frequency_hz[index]
        transfer[index] = linear_interpolate(frequencies, amplitude, frequency) *
                          cis(linear_interpolate(frequencies, phase, frequency))
        reference_spectrum[index] = spectrum.spectrum[index]
    end
    reference = irfft(reference_spectrum, length(spectrum.signal))
    candidate = irfft(reference_spectrum .* transfer, length(spectrum.signal))
    metrics = pulse_metrics(
        candidate,
        reference,
        reference,
        spectrum.time_s[2] - spectrum.time_s[1],
    )
    weights = abs2.(spectrum.spectrum[band.lower_index:band.upper_index])
    sample_frequencies = spectrum.frequency_hz[band.lower_index:band.upper_index]
    weighted_transmission = sum(
        weights[index] * linear_interpolate(
            frequencies,
            getproperty.(rows, :T00),
            sample_frequencies[index],
        ) for index in eachindex(weights)
    ) / sum(weights)
    weighted_reflection = sum(
        weights[index] * linear_interpolate(
            frequencies,
            getproperty.(rows, :R00),
            sample_frequencies[index],
        ) for index in eachindex(weights)
    ) / sum(weights)
    weighted_conversion = sum(
        weights[index] * linear_interpolate(
            frequencies,
            getproperty.(rows, :C0),
            sample_frequencies[index],
        ) for index in eachindex(weights)
    ) / sum(weights)
    (
        relative_peak=metrics.gain_peak,
        pulse_broadening_ratio=metrics.broadening_ratio,
        pulse_correlation=metrics.pulse_correlation,
        postcursor_ratio=metrics.postcursor_ratio,
        weighted_transmission,
        weighted_reflection,
        weighted_mode_conversion=weighted_conversion,
    )
end

function run_analysis_stage(stage::Symbol)
    raw = load_raw_rows(stage)
    variants, _ = stage_data(stage)
    all_rows = reduce(vcat, (calibrated_rows(raw, variant) for variant in variants))
    write_csv(joinpath(OUTPUT_ROOT, "g1_$(stage)_frequency.csv"), all_rows)
    summary = [begin
        current = sort(filter(row -> row.variant == variant, all_rows); by=row -> row.frequency_hz)
        pulse = pulse_metrics_for_variant(current)
        frequencies = getproperty.(current, :frequency_hz)
        delays = getproperty.(current, :trusted_excess_group_delay_s)
        delay_242 = linear_interpolate(frequencies, delays, 242.0e3)
        carrier_phase = linear_interpolate(
            frequencies,
            getproperty.(current, :relative_phase_rad),
            242.0e3,
        )
        required_phase = 2pi * 242.0e3 * required_center_delay_s()
        pass_pulse = pulse.pulse_broadening_ratio <= 1.25 &&
                     pulse.pulse_correlation >= 0.90 &&
                     pulse.postcursor_ratio <= 0.10
        pass_transfer = pulse.weighted_transmission >= 0.80 &&
                        pulse.weighted_reflection <= 0.20 &&
                        pulse.weighted_mode_conversion <= 0.05
        (
            variant,
            neck_height_mm=neck_height_mm(CONFIG, variant),
            minimum_feature_mm=minimum_feature_mm(CONFIG, variant),
            min_T00=minimum(getproperty.(current, :T00)),
            max_R00=maximum(getproperty.(current, :R00)),
            weighted_transmission=pulse.weighted_transmission,
            weighted_reflection=pulse.weighted_reflection,
            weighted_mode_conversion=pulse.weighted_mode_conversion,
            relative_peak=pulse.relative_peak,
            pulse_broadening_ratio=pulse.pulse_broadening_ratio,
            pulse_correlation=pulse.pulse_correlation,
            postcursor_ratio=pulse.postcursor_ratio,
            carrier_relative_phase_rad=carrier_phase,
            required_center_phase_rad=required_phase,
            carrier_phase_coverage=abs(carrier_phase) / required_phase,
            delay_at_242_us=delay_242 * 1.0e6,
            required_center_delay_us=required_center_delay_s() * 1.0e6,
            max_power_balance_error=maximum(abs.(getproperty.(current, :power_sum) .- 1)),
            pass_transfer,
            pass_pulse,
            pass_g1=pass_transfer && pass_pulse && isfinite(delay_242) && delay_242 >= 0,
        )
    end for variant in variants]
    write_csv(joinpath(OUTPUT_ROOT, "g1_$(stage)_summary.csv"), summary)
    if stage == :fine
        by_variant = Dict(row.variant => row for row in summary)
        selected = by_variant[:neck_2p0]
        straight = by_variant[:solid]
        nan = NaN
        benchmark_rows = [
            (
                topology="metal_grin_waisted_rings",
                case_id="none",
                reference="bare",
                frequency_hz=242.0e3,
                input_energy_j_per_m=nan,
                echo_mf_gain_bare=nan,
                echo_energy_gain_bare=nan,
                echo_to_ringing_db=nan,
                equal_noise_snr_proxy_db=nan,
                dead_zone_depth_mm=nan,
                focus_peak_gain_bare=nan,
                focus_energy_gain_bare=nan,
                fwhm_mm=nan,
                dof_minus3db_mm=nan,
                sidelobe_ratio=nan,
                pulse_broadening_ratio=nan,
                pulse_correlation=nan,
                postcursor_ratio=nan,
                weighted_transmission=nan,
                weighted_reflection=nan,
                weighted_mode_conversion=nan,
                minimum_feature_mm=nan,
                verdict="not_calculated_g1_stop",
            ),
            (
                topology="metal_grin_waisted_rings",
                case_id="none",
                reference="straight_footprint",
                frequency_hz=242.0e3,
                input_energy_j_per_m=nan,
                echo_mf_gain_bare=nan,
                echo_energy_gain_bare=nan,
                echo_to_ringing_db=nan,
                equal_noise_snr_proxy_db=nan,
                dead_zone_depth_mm=nan,
                focus_peak_gain_bare=nan,
                focus_energy_gain_bare=nan,
                fwhm_mm=nan,
                dof_minus3db_mm=nan,
                sidelobe_ratio=nan,
                pulse_broadening_ratio=straight.pulse_broadening_ratio,
                pulse_correlation=straight.pulse_correlation,
                postcursor_ratio=straight.postcursor_ratio,
                weighted_transmission=straight.weighted_transmission,
                weighted_reflection=straight.weighted_reflection,
                weighted_mode_conversion=straight.weighted_mode_conversion,
                minimum_feature_mm=straight.minimum_feature_mm,
                verdict="reference_only",
            ),
            (
                topology="metal_grin_waisted_rings_neck_2p0",
                case_id="none",
                reference="straight_footprint",
                frequency_hz=242.0e3,
                input_energy_j_per_m=nan,
                echo_mf_gain_bare=nan,
                echo_energy_gain_bare=nan,
                echo_to_ringing_db=nan,
                equal_noise_snr_proxy_db=nan,
                dead_zone_depth_mm=nan,
                focus_peak_gain_bare=nan,
                focus_energy_gain_bare=nan,
                fwhm_mm=nan,
                dof_minus3db_mm=nan,
                sidelobe_ratio=nan,
                pulse_broadening_ratio=selected.pulse_broadening_ratio,
                pulse_correlation=selected.pulse_correlation,
                postcursor_ratio=selected.postcursor_ratio,
                weighted_transmission=selected.weighted_transmission,
                weighted_reflection=selected.weighted_reflection,
                weighted_mode_conversion=selected.weighted_mode_conversion,
                minimum_feature_mm=selected.minimum_feature_mm,
                verdict="stop",
            ),
        ]
        propertynames(first(benchmark_rows)) == metric_columns() ||
            error("benchmark result columns do not match the shared contract")
        write_csv(joinpath(OUTPUT_ROOT, "benchmark_metrics.csv"), benchmark_rows)
    end
    for row in summary
        println(
            "[+] $(row.variant): Tbar=", round(row.weighted_transmission; digits=3),
            ", Rbar=", round(row.weighted_reflection; digits=3),
            ", delay=", round(row.delay_at_242_us; digits=3), " us",
            ", Bt=", round(row.pulse_broadening_ratio; digits=3),
            ", rho=", round(row.pulse_correlation; digits=3),
            ", G1=", row.pass_g1,
        )
    end
end

function run_child_stage(stage)
    if stage == "g0"
        run_g0_stage()
    elseif stage == "mesh"
        run_mesh_stage()
    elseif stage == "convert"
        run_convert_stage()
    elseif stage in ("coarse", "fine")
        run_parallel_stage(Symbol(stage))
    elseif startswith(stage, "coarse-worker-") || startswith(stage, "fine-worker-")
        solve_stage = startswith(stage, "coarse") ? :coarse : :fine
        worker_index = parse(Int, split(stage, "worker-"; limit=2)[2])
        1 <= worker_index <= MAX_WORKERS || error("invalid worker index")
        run_worker(solve_stage, worker_index)
    elseif stage == "analyze-coarse"
        run_analysis_stage(:coarse)
    elseif stage == "analyze-fine"
        run_analysis_stage(:fine)
    else
        error("unknown metal-GRIN stage: $stage")
    end
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("g0", "mesh", "convert", "coarse", "analyze-coarse")
            println("\n=== Metal-GRIN stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_child_stage(REQUESTED_STAGE)
    else
        error("usage: julia topology_metal_grin_run.jl [--stage=STAGE]")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
