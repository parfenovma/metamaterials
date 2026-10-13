module RunOperatorO6HLayeredContact

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o6g_motion_boundary.jl"))
using .RunOperatorO6GMotionBoundary

export layer_power_transmission,
       calibrated_two_way_contact,
       run_o6h_layered_contact

const G = RunOperatorO6GMotionBoundary
const M = G.M
const O6H_CONTACT_OUTPUT_ROOT = joinpath(
    M.PROJECT_ROOT, "results", "operator_space", "o6h_layered_contact_v0",
)
const PROBE_IMPEDANCE_RAYL = 17.0e6
const AL_DENSITY_KG_M3 = 2700.0
const AL_PRESSURE_SPEED_M_S = 6122.1024
const AL_IMPEDANCE_RAYL = AL_DENSITY_KG_M3 * AL_PRESSURE_SPEED_M_S
const COUPLANT_IMPEDANCE_RAYL = 1.5e6
const COUPLANT_SPEED_M_S = 1480.0
const CARRIER_FREQUENCY_HZ = 242.0e3
const COUPLANT_WAVELENGTH_M = COUPLANT_SPEED_M_S / CARRIER_FREQUENCY_HZ
const NOMINAL_THICKNESS_M = 0.01COUPLANT_WAVELENGTH_M
const THICKNESS_RMS_M = 0.005COUPLANT_WAVELENGTH_M
const IMPEDANCE_LOG_RMS = 0.10
const SPEED_LOG_RMS = 0.05
const CONTACT_COARSE_SCALES = (0.0, 0.5, 1.0, 2.0, 3.0, 4.0)
const BISECTION_STEPS = 12
const MAXIMUM_SEARCH_SCALE = 64.0

"""Power-wave transmission through one exact acoustic layer."""
function layer_power_transmission(
    frequency_hz,
    thickness_m,
    layer_impedance,
    layer_speed;
    input_impedance=PROBE_IMPEDANCE_RAYL,
    output_impedance=AL_IMPEDANCE_RAYL,
)
    phase = 2pi * frequency_hz * thickness_m / layer_speed
    cosine = cos(phase)
    sine = sin(phase)
    m11 = cosine
    m12 = im * layer_impedance * sine
    m21 = im * sine / layer_impedance
    m22 = cosine
    denominator = m11 + m12 / output_impedance +
        input_impedance * m21 + input_impedance * m22 / output_impedance
    pressure_transmission = 2 / denominator
    pressure_transmission * sqrt(input_impedance / output_impedance)
end

function calibrated_two_way_contact(
    frequency_hz,
    thickness_m,
    layer_impedance,
    layer_speed,
)
    actual = layer_power_transmission(
        frequency_hz, thickness_m, layer_impedance, layer_speed,
    )
    nominal = layer_power_transmission(
        frequency_hz,
        NOMINAL_THICKNESS_M,
        COUPLANT_IMPEDANCE_RAYL,
        COUPLANT_SPEED_M_S,
    )
    (actual / nominal)^2
end

function contact_matrix(config, thickness, impedance, speed)
    stop_count = length(thickness)
    transfer = Matrix{ComplexF64}(undef,
        stop_count, length(config.frequency_ratios),
    )
    for stop in 1:stop_count, (frequency_index, ratio) in
        enumerate(config.frequency_ratios)
        transfer[stop, frequency_index] = calibrated_two_way_contact(
            CARRIER_FREQUENCY_HZ * ratio,
            thickness[stop],
            impedance[stop],
            speed[stop],
        )
    end
    transfer
end

function contact_features(library, config, contact)
    stop_count, _, frequency_count, atom_count = size(library)
    size(contact) == (stop_count, frequency_count) ||
        throw(DimensionMismatch("one contact factor per stop and frequency"))
    raw = Matrix{ComplexF64}(undef, stop_count * frequency_count, atom_count)
    for stop in 1:stop_count, frequency in 1:frequency_count
        row = (stop - 1) * frequency_count + frequency
        for atom in 1:atom_count
            raw[row, atom] = contact[stop, frequency] *
                library[stop, stop, frequency, atom]
        end
    end
    transform = M.pulse_whitening_transform(config)
    features = similar(raw)
    for stop in 1:stop_count
        rows = ((stop - 1) * frequency_count + 1):(stop * frequency_count)
        features[rows, :] .= transform * view(raw, rows, :)
    end
    features
end

function unwrap_contact_phase(phases)
    result = Float64.(phases)
    for index in 2:length(result)
        difference = result[index] - result[index - 1]
        while difference > pi
            result[index] -= 2pi
            difference -= 2pi
        end
        while difference < -pi
            result[index] += 2pi
            difference += 2pi
        end
    end
    result
end

function affine_contact_residual(contact, frequency_ratios)
    design = hcat(ones(length(frequency_ratios)), frequency_ratios)
    residuals = Float64[]
    delays = Float64[]
    log_amplitudes = Float64[]
    for stop in axes(contact, 1)
        values = view(contact, stop, :)
        log_amplitude = mean(log.(max.(abs.(values), eps(Float64))))
        phase = unwrap_contact_phase(angle.(values))
        coefficients = design \ phase
        fitted = exp(log_amplitude) .* exp.(im .* (design * coefficients))
        push!(residuals, norm(values - fitted) / max(norm(values), eps(Float64)))
        push!(delays, -coefficients[2] / (2pi))
        push!(log_amplitudes, log_amplitude)
    end
    (; residuals, delays, log_amplitudes)
end

function reflect_nonnegative(value)
    abs(value)
end

function evaluate_contact_scale(setup, scale)
    config = setup.config
    hypotheses = setup.hypotheses
    frequency_count = length(config.frequency_ratios)
    physical = setup.saved_o6f
    thickness_standard = Float64.(physical["position_standard"])
    impedance_standard = Float64.(physical["amplitude_standard"])
    speed_standard = Float64.(physical["delay_standard"])
    signal_noise_rng = Xoshiro(M.O6F_SIGNAL_NOISE_SEED)
    h0_noise_rng = Xoshiro(M.O6F_H0_NOISE_SEED)
    counts = Dict{Tuple{Symbol,Float64,Float64},NamedTuple}()
    false_alarms = Dict(architecture.id => 0 for architecture in setup.architectures)
    trials = Dict(architecture.id => 0 for architecture in setup.architectures)
    fit_residuals = Float64[]
    fitted_delay_cycles = Float64[]
    fitted_log_amplitudes = Float64[]
    for (context_index, context) in enumerate(setup.contexts)
        thickness = [reflect_nonnegative(
            NOMINAL_THICKNESS_M + scale * THICKNESS_RMS_M *
                thickness_standard[stop, context_index],
        ) for stop in 1:config.element_count]
        impedance = COUPLANT_IMPEDANCE_RAYL .* exp.(
            scale * IMPEDANCE_LOG_RMS .* view(
                impedance_standard, :, context_index,
            ),
        )
        speed = COUPLANT_SPEED_M_S .* exp.(
            scale * SPEED_LOG_RMS .* view(speed_standard, :, context_index),
        )
        contact = contact_matrix(config, thickness, impedance, speed)
        fit = affine_contact_residual(contact, config.frequency_ratios)
        append!(fit_residuals, fit.residuals)
        append!(fitted_delay_cycles, fit.delays)
        append!(fitted_log_amplitudes, fit.log_amplitudes)
        signal_library = M.response_library(context.hypotheses, context.config)
        signal_full = contact_features(signal_library, context.config, contact)
        test_full = Dict(category => contact_features(
            setup.split_libraries[:test][category], config, contact,
        ) for category in (:backwall, :ringing, :material))
        signal_latent = M.complex_noise(
            signal_noise_rng, length(hypotheses), length(hypotheses),
        )
        h0_latent = M.complex_noise(
            h0_noise_rng, length(hypotheses), M.H0_TRIALS_PER_CONTEXT,
        )
        for architecture in setup.architectures
            id = architecture.id
            model = setup.models[id]
            indices = architecture.mechanical_indices
            signal_dictionary = M.select_stops(
                signal_full, indices, frequency_count,
            )
            clutter = hcat((
                setup.scales[category] .* M.select_stops(
                    test_full[category], indices, frequency_count,
                ) for category in (:backwall, :ringing, :material)
            )...)
            factor = M.covariance_factor(M.matched_output_covariance(
                model.precision_applied, clutter,
            ))
            batch = M.evaluate_with_common_noise(
                factor * signal_latent,
                M.matched_mean(model.precision_applied, signal_dictionary),
                model.normalization,
                model.threshold,
                setup.amplitudes,
                hypotheses,
                context.hypotheses,
            )
            M.add_batch!(counts, id, scale, batch)
            scores, _ = M.normalized_scores(
                factor * h0_latent, model.normalization,
            )
            false_alarms[id] += count(>(model.threshold), scores)
            trials[id] += length(scores)
        end
    end
    curves = M.aggregate_curve_rows(
        counts, setup.architectures, (scale,), setup.amplitudes,
    )
    rows = NamedTuple[]
    for architecture in setup.architectures
        id = architecture.id
        alpha90 = M.curve_alpha90(curves, id, scale)
        reference = M.curve_reference(
            curves, id, scale, setup.reference_alpha,
        )
        pfa_upper = M.wilson_upper_bound(false_alarms[id], trials[id])
        ratio = alpha90 / setup.baseline_alpha[id]
        gates = (
            pfa=pfa_upper <= M.PFA_GATE,
            alpha90=isfinite(ratio) && ratio <= M.ALPHA90_RATIO_GATE,
            pdl=reference.pdl >= M.PDL_REFERENCE_GATE,
        )
        push!(rows, (
            architecture_id=id,
            scale,
            thickness_rms_um=1e6scale * THICKNESS_RMS_M,
            impedance_log_rms=scale * IMPEDANCE_LOG_RMS,
            speed_log_rms=scale * SPEED_LOG_RMS,
            h0_trials=trials[id],
            false_alarms=false_alarms[id],
            observed_pfa=false_alarms[id] / trials[id],
            pfa_upper_95=pfa_upper,
            alpha90,
            nominal_replay_alpha90=setup.baseline_alpha[id],
            alpha90_ratio_to_nominal=ratio,
            pd_reference=reference.pd,
            pdl_reference=reference.pdl,
            median_affine_fit_residual=median(fit_residuals),
            q95_affine_fit_residual=quantile(fit_residuals, 0.95),
            fitted_log_amplitude_rms=std(fitted_log_amplitudes),
            fitted_delay_rms_cycles=std(fitted_delay_cycles),
            pfa_gate=gates.pfa,
            alpha90_gate=gates.alpha90,
            pdl_gate=gates.pdl,
            joint_pass=all(values(gates)),
        ))
    end
    rows, curves
end

function cached_evaluation!(cache, curve_cache, setup, scale)
    key = Float64(scale)
    if !haskey(cache, key)
        println("[>] O6h-B layered contact scale=$(round(scale; digits=6))")
        rows, curves = evaluate_contact_scale(setup, key)
        cache[key] = rows
        curve_cache[key] = curves
    end
    cache[key]
end

function architecture_row(rows, architecture_id)
    only(filter(row -> row.architecture_id == architecture_id, rows))
end

function bracket!(cache, curve_cache, setup, architecture_id)
    scales = Float64.(collect(CONTACT_COARSE_SCALES))
    rows = [architecture_row(cached_evaluation!(
        cache, curve_cache, setup, scale,
    ), architecture_id) for scale in scales]
    first(rows).joint_pass || error("nominal calibrated contact must pass")
    failure = findfirst(row -> !row.joint_pass, rows)
    while isnothing(failure)
        next_scale = 2last(scales)
        next_scale <= MAXIMUM_SEARCH_SCALE || error("contact boundary not found")
        push!(scales, next_scale)
        push!(rows, architecture_row(cached_evaluation!(
            cache, curve_cache, setup, next_scale,
        ), architecture_id))
        failure = last(rows).joint_pass ? nothing : length(rows)
    end
    upper = something(failure)
    (scales[upper - 1], scales[upper])
end

function bisect!(cache, curve_cache, setup, architecture_id, lower, upper)
    for _ in 1:BISECTION_STEPS
        middle = (lower + upper) / 2
        row = architecture_row(cached_evaluation!(
            cache, curve_cache, setup, middle,
        ), architecture_id)
        if row.joint_pass
            lower = middle
        else
            upper = middle
        end
    end
    passing = architecture_row(cache[lower], architecture_id)
    failing = architecture_row(cache[upper], architecture_id)
    failed_gates = join((String(name) for (name, value) in pairs((
        pfa=failing.pfa_gate,
        alpha90=failing.alpha90_gate,
        pdl=failing.pdl_gate,
    )) if !value), ";")
    (
        architecture_id,
        maximum_passing_scale=lower,
        first_failing_scale=upper,
        threshold_width=upper - lower,
        thickness_rms_um=passing.thickness_rms_um,
        impedance_log_rms=passing.impedance_log_rms,
        speed_log_rms=passing.speed_log_rms,
        passing_pfa_upper_95=passing.pfa_upper_95,
        passing_alpha90_ratio=passing.alpha90_ratio_to_nominal,
        passing_pdl_reference=passing.pdl_reference,
        q95_affine_fit_residual=passing.q95_affine_fit_residual,
        fitted_log_amplitude_rms=passing.fitted_log_amplitude_rms,
        fitted_delay_rms_cycles=passing.fitted_delay_rms_cycles,
        first_failed_gates=failed_gates,
    )
end

function run_o6h_layered_contact()
    mkpath(O6H_CONTACT_OUTPUT_ROOT)
    setup = G.setup_frozen_replay()
    cache = Dict{Float64,Vector{NamedTuple}}()
    curve_cache = Dict{Float64,Vector{NamedTuple}}()
    threshold_rows = NamedTuple[]
    for architecture_id in M.ARCHITECTURE_IDS
        lower, upper = bracket!(
            cache, curve_cache, setup, architecture_id,
        )
        push!(threshold_rows, bisect!(
            cache, curve_cache, setup, architecture_id, lower, upper,
        ))
    end
    all_rows = reduce(vcat, values(cache))
    all_curves = reduce(vcat, values(curve_cache))
    sort!(all_rows; by=row -> (String(row.architecture_id), row.scale))
    sort!(all_curves; by=row ->
        (String(row.architecture_id), row.scale, row.amplitude))
    sort!(threshold_rows; by=row -> String(row.architecture_id))
    M.write_rows(joinpath(O6H_CONTACT_OUTPUT_ROOT, "contact_boundary_evaluations.csv"), all_rows)
    M.write_rows(joinpath(O6H_CONTACT_OUTPUT_ROOT, "contact_boundary_curves.csv"), all_curves)
    M.write_rows(joinpath(O6H_CONTACT_OUTPUT_ROOT, "contact_tolerance_thresholds.csv"), threshold_rows)
    assumptions = [
        (parameter="probe face impedance", value=PROBE_IMPEDANCE_RAYL / 1e6, unit="MRayl"),
        (parameter="Al impedance", value=AL_IMPEDANCE_RAYL / 1e6, unit="MRayl"),
        (parameter="couplant impedance", value=COUPLANT_IMPEDANCE_RAYL / 1e6, unit="MRayl"),
        (parameter="couplant speed", value=COUPLANT_SPEED_M_S, unit="m/s"),
        (parameter="nominal thickness", value=1e6NOMINAL_THICKNESS_M, unit="um"),
    ]
    M.write_rows(joinpath(O6H_CONTACT_OUTPUT_ROOT, "contact_assumptions.csv"), assumptions)
    jldsave(
        joinpath(O6H_CONTACT_OUTPUT_ROOT, "o6h_layered_contact.jld2");
        format_version=1,
        source_o6f_path=abspath(G.O6G_MOTION_SOURCE),
        assumptions,
        coarse_scales=CONTACT_COARSE_SCALES,
        bisection_steps=BISECTION_STEPS,
        all_rows,
        all_curves,
        threshold_rows,
    )
    open(joinpath(O6H_CONTACT_OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6h-B exact layered-contact motion replay")
        println(io)
        println(io, "The frozen P12/P17 detector is driven by an exact normal-incidence thin-layer transfer matrix with all layer multiple reflections. The nominal 61.16 um couplant response is calibrated once; exact per-stop, per-frequency residual transfer multiplies signal and nuisance without an affine approximation.")
        println(io)
        println(io, "| Layout | Max passing scale | Thickness RMS | Z log-RMS | c log-RMS | q95 non-affine residual | First failed gate |")
        println(io, "|---|---:|---:|---:|---:|---:|---|")
        for row in threshold_rows
            @printf(io,
                "| `%s` | `%.5f` | `%.2f um` | `%.2f%%` | `%.2f%%` | `%.3f%%` | `%s` |\n",
                String(row.architecture_id), row.maximum_passing_scale,
                row.thickness_rms_um, 100row.impedance_log_rms,
                100row.speed_log_rms, 100row.q95_affine_fit_residual,
                row.first_failed_gates,
            )
        end
        println(io)
        println(io, "The fitted equivalent log-amplitude and delay RMS are reported in the threshold CSV to compare this physical layer with O6g. The calculation is exact for normal-incidence longitudinal plane waves, but omits oblique P/S conversion, lateral squeeze-out and finite-footprint lift-off; those define the justified 2D spot check.")
    end
    frontier_plot = plot(
        xlabel="physical layer variation scale",
        ylabel="alpha90 / nominal",
        title="O6h-B layered contact boundary",
        legend=:topleft,
    )
    for (id, marker, color) in zip(
        M.ARCHITECTURE_IDS, (:circle, :diamond), (:steelblue, :darkorange),
    )
        selected = sort(filter(row -> row.architecture_id == id, all_rows);
            by=row -> row.scale)
        plot!(frontier_plot, getproperty.(selected, :scale),
            getproperty.(selected, :alpha90_ratio_to_nominal);
            marker, color, linewidth=2, label=String(id))
    end
    hline!(frontier_plot, [M.ALPHA90_RATIO_GATE];
        linestyle=:dash, color=:black, label="1.10 gate")
    savefig(frontier_plot, joinpath(
        O6H_CONTACT_OUTPUT_ROOT, "layered_contact_frontier.png",
    ))
    for row in threshold_rows
        println("[+] $(row.architecture_id) layered-contact threshold=$(round(row.maximum_passing_scale; digits=5))")
    end
    println("[+] wrote O6h-B to $O6H_CONTACT_OUTPUT_ROOT")
    (; threshold_rows, all_rows)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6h_layered_contact()

end
