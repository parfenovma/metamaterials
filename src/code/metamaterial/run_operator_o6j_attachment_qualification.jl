module RunOperatorO6JAttachmentQualification

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "run_operator_o6i_mitigation_map.jl"))
using .RunOperatorO6IMitigationMap

export run_o6j_attachment_qualification

const I = RunOperatorO6IMitigationMap
const A = I.A
const C = I.C
const G = C.G
const M = I.M
const OUTPUT_ROOT = joinpath(
    M.PROJECT_ROOT, "results", "operator_space", "o6j_attachment_qualification_v0",
)
const CP_STANDARD_SEED = UInt64(0x06f0e117)
const ATTACHMENT_SCALES = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
const BISECTION_STEPS = 12
const MAXIMUM_SEARCH_SCALE = 64.0
const ENCODER_POSITION_RMS_LAMBDA = 0.05
const THICKNESS_LOG_RMS = 0.10
const IMPEDANCE_LOG_RMS = 0.05
const PRESSURE_SPEED_LOG_RMS = 0.02
const SHEAR_SPEED_LOG_RMS = 0.05
const HANDOFF_IMPEDANCE = 0.75sqrt(
    C.PROBE_IMPEDANCE_RAYL * C.AL_IMPEDANCE_RAYL,
)
const HANDOFF_PRESSURE_SPEED = C.AL_PRESSURE_SPEED_M_S
const HANDOFF_SHEAR_SPEED = 0.25HANDOFF_PRESSURE_SPEED
const HANDOFF_DENSITY = HANDOFF_IMPEDANCE / HANDOFF_PRESSURE_SPEED
const HANDOFF_THICKNESS_M = 0.02HANDOFF_PRESSURE_SPEED / C.CARRIER_FREQUENCY_HZ
const P12_ID = :motion_compact_center_p12_native

function handoff_material(
    impedance=HANDOFF_IMPEDANCE,
    pressure_speed=HANDOFF_PRESSURE_SPEED,
    shear_speed=HANDOFF_SHEAR_SPEED,
)
    A.solid_material(impedance / pressure_speed, pressure_speed, shear_speed)
end

function angular_corner_rows(config)
    rows = NamedTuple[]
    corner_index = 0
    for thickness_sign in (-1.0, 1.0), impedance_sign in (-1.0, 1.0),
        pressure_sign in (-1.0, 1.0), shear_sign in (-1.0, 1.0)
        corner_index += 1
        thickness = HANDOFF_THICKNESS_M * exp(
            thickness_sign * THICKNESS_LOG_RMS,
        )
        impedance = HANDOFF_IMPEDANCE * exp(
            impedance_sign * IMPEDANCE_LOG_RMS,
        )
        pressure_speed = HANDOFF_PRESSURE_SPEED * exp(
            pressure_sign * PRESSURE_SPEED_LOG_RMS,
        )
        shear_speed = HANDOFF_SHEAR_SPEED * exp(
            shear_sign * SHEAR_SPEED_LOG_RMS,
        )
        material = handoff_material(impedance, pressure_speed, shear_speed)
        local_rows = NamedTuple[]
        for ratio in config.frequency_ratios, angle in A.ANGLES_DEG
            frequency = C.CARRIER_FREQUENCY_HZ * ratio
            layer = I.solid_layer_scattering(
                frequency, angle, thickness, material,
            )
            direct = A.direct_interface_scattering(frequency, angle)
            push!(local_rows, (
                frequency_ratio=ratio,
                angle,
                power=layer.transmitted_power /
                    max(direct.transmitted_power, eps(Float64)),
                s_change=abs(
                    layer.transmitted_s_fraction - direct.transmitted_s_fraction,
                ),
                energy_error=abs(layer.energy_balance - 1),
            ))
        end
        qualification = filter(row -> row.angle <= 45.0, local_rows)
        minimum_power = minimum(getproperty.(qualification, :power))
        maximum_s_change = maximum(getproperty.(qualification, :s_change))
        maximum_energy_error = maximum(getproperty.(local_rows, :energy_error))
        gates = (
            power=minimum_power >= A.POWER_RATIO_GATE,
            mode_conversion=maximum_s_change <= A.S_FRACTION_CHANGE_GATE,
            energy=maximum_energy_error <= A.ENERGY_BALANCE_TOLERANCE,
        )
        push!(rows, (
            corner_index,
            thickness_sign,
            impedance_sign,
            pressure_speed_sign=pressure_sign,
            shear_speed_sign=shear_sign,
            thickness_um=1e6thickness,
            pressure_impedance_mrayl=impedance / 1e6,
            pressure_speed_m_s=pressure_speed,
            shear_speed_m_s=shear_speed,
            density_kg_m3=impedance / pressure_speed,
            minimum_layer_to_direct_power_0_45=minimum_power,
            maximum_s_fraction_change_0_45=maximum_s_change,
            minimum_layer_to_direct_power_0_55=
                minimum(getproperty.(local_rows, :power)),
            maximum_energy_balance_error=maximum_energy_error,
            power_gate=gates.power,
            mode_conversion_gate=gates.mode_conversion,
            energy_gate=gates.energy,
            joint_pass=all(values(gates)),
        ))
    end
    rows
end

function solid_contact_matrix(config, thickness, impedance, pressure_speed)
    transfer = Matrix{ComplexF64}(undef,
        length(thickness), length(config.frequency_ratios),
    )
    nominal_material = handoff_material()
    for stop in eachindex(thickness)
        material = handoff_material(
            impedance[stop], pressure_speed[stop],
            HANDOFF_SHEAR_SPEED,
        )
        for (frequency_index, ratio) in enumerate(config.frequency_ratios)
            frequency = C.CARRIER_FREQUENCY_HZ * ratio
            actual = I.solid_layer_scattering(
                frequency, 0.0, thickness[stop], material,
            ).transmitted_p
            nominal = I.solid_layer_scattering(
                frequency, 0.0, HANDOFF_THICKNESS_M, nominal_material,
            ).transmitted_p
            transfer[stop, frequency_index] = (actual / nominal)^2
        end
    end
    transfer
end

function scan_features_with_contact(hypotheses, config, positions, contact)
    stop_count = length(positions)
    frequency_count = length(config.frequency_ratios)
    size(contact) == (stop_count, frequency_count) ||
        throw(DimensionMismatch("contact matrix does not match scan"))
    raw = Matrix{ComplexF64}(undef,
        stop_count * frequency_count, length(hypotheses),
    )
    for (stop, position) in enumerate(positions),
        (frequency_index, frequency) in enumerate(config.frequency_ratios)
        row = (stop - 1) * frequency_count + frequency_index
        for (hypothesis_index, hypothesis) in enumerate(hypotheses)
            raw[row, hypothesis_index] = contact[stop, frequency_index] *
                M.monostatic_value(position, hypothesis, frequency, config)
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

function ringing_features_with_contact(library, config, positions, contact)
    native_positions = M.element_positions(config)
    frequency_count = length(config.frequency_ratios)
    atom_count = size(library, 4)
    raw = Matrix{ComplexF64}(undef,
        length(positions) * frequency_count, atom_count,
    )
    diagonal_values = Matrix{ComplexF64}(undef,
        length(native_positions), atom_count,
    )
    for frequency_index in 1:frequency_count
        for atom in 1:atom_count, index in eachindex(native_positions)
            diagonal_values[index, atom] = library[
                index, index, frequency_index, atom,
            ]
        end
        for (stop, position) in enumerate(positions), atom in 1:atom_count
            row = (stop - 1) * frequency_count + frequency_index
            raw[row, atom] = contact[stop, frequency_index] * M.linear_sample(
                native_positions, view(diagonal_values, :, atom), position,
            )
        end
    end
    transform = M.pulse_whitening_transform(config)
    features = similar(raw)
    for stop in eachindex(positions)
        rows = ((stop - 1) * frequency_count + 1):(stop * frequency_count)
        features[rows, :] .= transform * view(raw, rows, :)
    end
    features
end

function evaluate_attachment_scale(setup, scale, pressure_speed_standard)
    config = setup.config
    hypotheses = setup.hypotheses
    architecture = only(filter(row -> row.id == P12_ID, setup.architectures))
    model = setup.models[P12_ID]
    physical = setup.saved_o6f
    position_standard = Float64.(physical["position_standard"])
    thickness_standard = Float64.(physical["amplitude_standard"])
    impedance_standard = Float64.(physical["delay_standard"])
    signal_noise_rng = Xoshiro(M.O6F_SIGNAL_NOISE_SEED)
    h0_noise_rng = Xoshiro(M.O6F_H0_NOISE_SEED)
    counts = Dict{Tuple{Symbol,Float64,Float64},NamedTuple}()
    false_alarms = 0
    h0_trials = 0
    frequency_count = length(config.frequency_ratios)
    for (context_index, context) in enumerate(setup.contexts)
        positions = M.element_positions(config) .+
            ENCODER_POSITION_RMS_LAMBDA .* view(
                position_standard, :, context_index,
            )
        thickness = HANDOFF_THICKNESS_M .* exp.(
            scale * THICKNESS_LOG_RMS .* view(
                thickness_standard, :, context_index,
            ),
        )
        impedance = HANDOFF_IMPEDANCE .* exp.(
            scale * IMPEDANCE_LOG_RMS .* view(
                impedance_standard, :, context_index,
            ),
        )
        pressure_speed = HANDOFF_PRESSURE_SPEED .* exp.(
            scale * PRESSURE_SPEED_LOG_RMS .* view(
                pressure_speed_standard, :, context_index,
            ),
        )
        contact = solid_contact_matrix(
            config, thickness, impedance, pressure_speed,
        )
        signal_full = scan_features_with_contact(
            context.hypotheses, context.config, positions, contact,
        )
        material_full = scan_features_with_contact(
            M.material_clutter_hypotheses(M.MOTION_SPLIT_SEEDS[:nuisance_test]),
            config,
            positions,
            contact,
        )
        backwall_full = scan_features_with_contact(
            M.backwall_hypotheses(M.MOTION_SPLIT_SEEDS[:nuisance_test]),
            config,
            positions,
            contact,
        )
        ringing_full = ringing_features_with_contact(
            setup.split_libraries[:test][:ringing], config, positions, contact,
        )
        indices = architecture.mechanical_indices
        signal_dictionary = M.select_stops(
            signal_full, indices, frequency_count,
        )
        clutter = hcat(
            setup.scales[:backwall] .* M.select_stops(
                backwall_full, indices, frequency_count,
            ),
            setup.scales[:ringing] .* M.select_stops(
                ringing_full, indices, frequency_count,
            ),
            setup.scales[:material] .* M.select_stops(
                material_full, indices, frequency_count,
            ),
        )
        factor = M.covariance_factor(M.matched_output_covariance(
            model.precision_applied, clutter,
        ))
        signal_latent = M.complex_noise(
            signal_noise_rng, length(hypotheses), length(hypotheses),
        )
        h0_latent = M.complex_noise(
            h0_noise_rng, length(hypotheses), M.H0_TRIALS_PER_CONTEXT,
        )
        batch = M.evaluate_with_common_noise(
            factor * signal_latent,
            M.matched_mean(model.precision_applied, signal_dictionary),
            model.normalization,
            model.threshold,
            setup.amplitudes,
            hypotheses,
            context.hypotheses,
        )
        M.add_batch!(counts, P12_ID, scale, batch)
        scores, _ = M.normalized_scores(
            factor * h0_latent, model.normalization,
        )
        false_alarms += count(>(model.threshold), scores)
        h0_trials += length(scores)
    end
    curves = M.aggregate_curve_rows(
        counts, [architecture], (scale,), setup.amplitudes,
    )
    alpha90 = M.curve_alpha90(curves, P12_ID, scale)
    reference = M.curve_reference(
        curves, P12_ID, scale, setup.reference_alpha,
    )
    pfa_upper = M.wilson_upper_bound(false_alarms, h0_trials)
    ratio = alpha90 / setup.baseline_alpha[P12_ID]
    gates = (
        pfa=pfa_upper <= M.PFA_GATE,
        alpha90=isfinite(ratio) && ratio <= M.ALPHA90_RATIO_GATE,
        pdl=reference.pdl >= M.PDL_REFERENCE_GATE,
    )
    row = (
        architecture_id=P12_ID,
        attachment_scale=scale,
        encoder_position_rms_lambda=ENCODER_POSITION_RMS_LAMBDA,
        thickness_log_rms=scale * THICKNESS_LOG_RMS,
        impedance_log_rms=scale * IMPEDANCE_LOG_RMS,
        pressure_speed_log_rms=scale * PRESSURE_SPEED_LOG_RMS,
        h0_trials,
        false_alarms,
        observed_pfa=false_alarms / h0_trials,
        pfa_upper_95=pfa_upper,
        alpha90,
        zero_error_o6f_alpha90=setup.baseline_alpha[P12_ID],
        alpha90_ratio_to_o6f=ratio,
        pd_reference=reference.pd,
        pdl_reference=reference.pdl,
        pfa_gate=gates.pfa,
        alpha90_gate=gates.alpha90,
        pdl_gate=gates.pdl,
        joint_pass=all(values(gates)),
    )
    row, curves
end

function cached_evaluation!(cache, curve_cache, setup, scale, cp_standard)
    key = Float64(scale)
    if !haskey(cache, key)
        println("[>] O6j attachment scale=$(round(scale; digits=6))")
        row, curves = evaluate_attachment_scale(setup, key, cp_standard)
        cache[key] = row
        curve_cache[key] = curves
    end
    cache[key]
end

function attachment_threshold!(cache, curve_cache, setup, cp_standard)
    scales = Float64.(collect(ATTACHMENT_SCALES))
    rows = [cached_evaluation!(
        cache, curve_cache, setup, scale, cp_standard,
    ) for scale in scales]
    first(rows).joint_pass || error("fixed encoder baseline does not pass")
    failure = findfirst(row -> !row.joint_pass, rows)
    while isnothing(failure)
        next_scale = 2last(scales)
        next_scale <= MAXIMUM_SEARCH_SCALE || error("attachment boundary not found")
        push!(scales, next_scale)
        push!(rows, cached_evaluation!(
            cache, curve_cache, setup, next_scale, cp_standard,
        ))
        failure = last(rows).joint_pass ? nothing : length(rows)
    end
    upper_index = something(failure)
    lower = scales[upper_index - 1]
    upper = scales[upper_index]
    for _ in 1:BISECTION_STEPS
        middle = (lower + upper) / 2
        row = cached_evaluation!(
            cache, curve_cache, setup, middle, cp_standard,
        )
        if row.joint_pass
            lower = middle
        else
            upper = middle
        end
    end
    passing = cache[lower]
    failing = cache[upper]
    failed_gates = join((String(name) for (name, value) in pairs((
        pfa=failing.pfa_gate,
        alpha90=failing.alpha90_gate,
        pdl=failing.pdl_gate,
    )) if !value), ";")
    (
        architecture_id=P12_ID,
        maximum_passing_attachment_scale=lower,
        first_failing_attachment_scale=upper,
        threshold_width=upper - lower,
        encoder_position_rms_lambda=ENCODER_POSITION_RMS_LAMBDA,
        thickness_log_rms=passing.thickness_log_rms,
        impedance_log_rms=passing.impedance_log_rms,
        pressure_speed_log_rms=passing.pressure_speed_log_rms,
        passing_pfa_upper_95=passing.pfa_upper_95,
        passing_alpha90_ratio=passing.alpha90_ratio_to_o6f,
        passing_pdl_reference=passing.pdl_reference,
        first_failed_gates=failed_gates,
    )
end

function run_o6j_attachment_qualification()
    mkpath(OUTPUT_ROOT)
    setup = G.setup_frozen_replay()
    corners = angular_corner_rows(setup.config)
    angular_passed = all(getproperty.(corners, :joint_pass))
    rng = Xoshiro(CP_STANDARD_SEED)
    cp_standard = randn(
        rng, setup.config.element_count, M.MOTION_TEST_CONTEXTS,
    )
    cache = Dict{Float64,NamedTuple}()
    curve_cache = Dict{Float64,Vector{NamedTuple}}()
    threshold = attachment_threshold!(
        cache, curve_cache, setup, cp_standard,
    )
    nominal = cache[1.0]
    detector_passed = nominal.joint_pass
    recommendation = angular_passed && detector_passed
    evaluations = sort(collect(values(cache)); by=row -> row.attachment_scale)
    curves = reduce(vcat, values(curve_cache))
    sort!(curves; by=row -> (row.scale, row.amplitude))
    M.write_rows(joinpath(OUTPUT_ROOT, "angular_attachment_corners.csv"), corners)
    M.write_rows(joinpath(OUTPUT_ROOT, "attachment_detector_evaluations.csv"), evaluations)
    M.write_rows(joinpath(OUTPUT_ROOT, "attachment_detector_curves.csv"), curves)
    M.write_rows(joinpath(OUTPUT_ROOT, "attachment_threshold.csv"), [threshold])
    M.write_rows(joinpath(OUTPUT_ROOT, "frozen_solid_contact.csv"), [(
        pressure_impedance_mrayl=HANDOFF_IMPEDANCE / 1e6,
        pressure_speed_m_s=HANDOFF_PRESSURE_SPEED,
        shear_speed_m_s=HANDOFF_SHEAR_SPEED,
        density_kg_m3=HANDOFF_DENSITY,
        thickness_um=1e6HANDOFF_THICKNESS_M,
    )])
    jldsave(
        joinpath(OUTPUT_ROOT, "o6j_attachment_qualification.jld2");
        format_version=1,
        source_o6f_path=abspath(G.O6G_MOTION_SOURCE),
        corners,
        evaluations,
        curves,
        threshold,
        cp_standard_seed=CP_STANDARD_SEED,
        cp_standard,
        angular_passed,
        detector_passed,
        recommendation,
    )
    open(joinpath(OUTPUT_ROOT, "README.md"), "w") do io
        println(io, "# O6j P12 solid-contact attachment qualification")
        println(io)
        println(io, recommendation ?
            "**CONDITIONAL P12 COUPON/PROTOTYPE RECOMMENDATION.**" :
            "**NO PROTOTYPE RECOMMENDATION: at least one frozen attachment gate failed.**")
        println(io)
        @printf(io,
            "All angular attachment corners pass: `%s` (`%d/%d`). Worst corner power `%.4f`, worst Delta S `%.4f`, maximum energy error `%.3e`.\n\n",
            angular_passed, count(getproperty.(corners, :joint_pass)), length(corners),
            minimum(getproperty.(corners, :minimum_layer_to_direct_power_0_45)),
            maximum(getproperty.(corners, :maximum_s_fraction_change_0_45)),
            maximum(getproperty.(corners, :maximum_energy_balance_error)),
        )
        @printf(io,
            "At the nominal attachment variation scale `1.0` plus fixed `0.05 lambda` encoder RMS, P12 has Pfa upper `%.3f%%`, alpha90 ratio `%.4f`, and reference Pdl `%.4f`: joint pass `%s`.\n\n",
            100nominal.pfa_upper_95, nominal.alpha90_ratio_to_o6f,
            nominal.pdl_reference, nominal.joint_pass,
        )
        @printf(io,
            "The finite maximum passing attachment scale is `%.5f`: thickness log-RMS `%.2f%%`, impedance log-RMS `%.2f%%`, cP log-RMS `%.2f%%`; first failed gate `%s`.\n\n",
            threshold.maximum_passing_attachment_scale,
            100threshold.thickness_log_rms,
            100threshold.impedance_log_rms,
            100threshold.pressure_speed_log_rms,
            threshold.first_failed_gates,
        )
        println(io, "The recommendation is for a material-identification and finite-footprint coupon, not a product or POD claim. The frozen effective layer is approximately 12.57 MRayl, cP 6122 m/s, cS 1531 m/s, rho 2054 kg/m3 and 506 um thick. Wear, debonding, squeeze-out edges and scan-time economics remain benchtop gates. q8 remains blocked on a concrete broadband/multiresonant matching topology.")
    end
    plot_object = plot(
        getproperty.(evaluations, :attachment_scale),
        getproperty.(evaluations, :alpha90_ratio_to_o6f);
        marker=:circle,
        linewidth=2,
        xlabel="attachment variation scale",
        ylabel="alpha90 / zero-error O6f",
        title="O6j P12 attachment qualification",
        label="P12 + fixed encoder error",
    )
    hline!(plot_object, [M.ALPHA90_RATIO_GATE];
        linestyle=:dash, color=:black, label="1.10 gate")
    vline!(plot_object, [threshold.maximum_passing_attachment_scale];
        linestyle=:dot, color=:red, label="attachment ceiling")
    savefig(plot_object, joinpath(OUTPUT_ROOT, "attachment_qualification.png"))
    println(recommendation ?
        "[+] O6j conditional P12 prototype recommendation" :
        "[-] O6j prototype recommendation withheld")
    println("[+] wrote O6j to $OUTPUT_ROOT")
    (; recommendation, corners, evaluations, threshold)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && run_o6j_attachment_qualification()

end
