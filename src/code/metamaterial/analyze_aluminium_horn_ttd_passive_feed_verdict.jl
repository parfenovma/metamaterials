module AnalyzeAluminiumHornTTDPassiveFeedVerdict

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_feed_verdict")
const SPLITTER_ROOT = joinpath(PROJECT_ROOT, "tmp")
const ARRAY_ROOT = joinpath(
    PROJECT_ROOT, "results", "aluminium_horn_ttd_p6_242khz",
    "fifteen_channel_diffuser",
)
const VERIFIED_PULSE_GAIN = 2.1260609847797616
const VERIFIED_CARRIER_UPPER_GAIN = 2.05240

using DelimitedFiles
using JLD2
using LinearAlgebra
using Plots

include(joinpath(@__DIR__, "design_aluminium_horn_ttd_passive_feed.jl"))
include(joinpath(@__DIR__, "analyze_aluminium_horn_ttd_power_weights.jl"))
using .DesignAluminiumHornTTDPassiveFeed

modal_path(root_name) = joinpath(
    SPLITTER_ROOT, root_name, "passive_feed_splitter_modal_f242000hz.jld2",
)

function optimistic_tree(eta_equal, eta_other)
    target = DesignAluminiumHornTTDPassiveFeed.selected_full_weights()
    root, _ = DesignAluminiumHornTTDPassiveFeed.balanced_tree(abs2.(target))
    leaf_power = zeros(Float64, length(target))
    function distribute(node, input_power)
        if DesignAluminiumHornTTDPassiveFeed.isleaf(node)
            leaf_power[only(node.leaves)] = input_power
            return
        end
        left_fraction = node.left.power / node.power
        small_fraction = min(left_fraction, 1 - left_fraction)
        efficiency = isapprox(small_fraction, 0.5; atol=1e-10) ?
                     eta_equal : eta_other
        distribute(node.left, input_power * efficiency * left_fraction)
        distribute(node.right, input_power * efficiency * (1 - left_fraction))
    end
    distribute(root, 1.0)
    target, leaf_power
end

function routed_carrier_gain(leaf_power)
    response = JLD2.load(joinpath(ARRAY_ROOT, "lens_half_response_matrix.jld2"))
    reference = JLD2.load(joinpath(ARRAY_ROOT, "uniform_half_harmonic.jld2"))
    focus_response = vec(response["focus_response_m"][1, :])
    power_matrix = AnalyseAluminiumHornTTDPowerWeights.active_power_matrix(
        response["source_normal_displacement_integral_m3"],
    )
    useful_efficiency = sum(leaf_power)
    raw_half_weights = Float64[
        sqrt(leaf_power[8]); sqrt.(leaf_power[7:-1:1])...
    ]
    uniform_power = abs(Float64(reference["active_input_power_w"]))
    raw_power = AnalyseAluminiumHornTTDPowerWeights.active_power(
        power_matrix, raw_half_weights,
    )
    weights = raw_half_weights .* sqrt(uniform_power * useful_efficiency / raw_power)
    gain = abs(sum(focus_response .* weights)) /
           abs(reference["focus_displacement_m"][1])
    gain, weights
end

function best_mmi_row()
    raw = readdlm(
        joinpath(SPLITTER_ROOT, "aluminium_horn_ttd_passive_mmi", "mmi_width_sweep.csv"),
        ',', Float64; skipstart=1,
    )
    index = argmax(raw[:, 4])
    (
        width_mm=raw[index, 1],
        mode_count=round(Int, raw[index, 2]),
        length_mm=raw[index, 3],
        coherence=raw[index, 4],
    )
end

function run()
    straight = JLD2.load(modal_path("aluminium_horn_ttd_passive_straight_control"))
    equal = JLD2.load(modal_path("aluminium_horn_ttd_passive_equal_long_splitter"))
    unequal = JLD2.load(modal_path("aluminium_horn_ttd_passive_calibrated_beat_splitter"))
    manifold = JLD2.load(joinpath(
        SPLITTER_ROOT, "aluminium_horn_ttd_passive_manifold",
        "passive_feed_manifold_f242000hz.jld2",
    ))
    mmi = best_mmi_row()
    eta_equal = Float64(equal["useful_fundamental_power_transmission"])
    eta_other = Float64(unequal["useful_fundamental_power_transmission"])
    target, leaf_power = optimistic_tree(eta_equal, eta_other)
    tree_efficiency = sum(leaf_power)
    required_efficiency = (2 / VERIFIED_PULSE_GAIN)^2
    pulse_power_only_upper = VERIFIED_PULSE_GAIN * sqrt(tree_efficiency)
    carrier_power_only_upper = VERIFIED_CARRIER_UPPER_GAIN * sqrt(tree_efficiency)
    routed_gain, routed_weights = routed_carrier_gain(leaf_power)
    passed = pulse_power_only_upper >= 2.0
    mkpath(OUTPUT_ROOT)

    open(joinpath(OUTPUT_ROOT, "passive_feed_candidate_summary.csv"), "w") do io
        println(io, "candidate,reflection_power_fraction,all_mode_power_transmission,useful_fundamental_power_transmission,amplitude_ratio,relative_phase_deg,target_coherence,carrier_or_reduced_gain,passed")
        println(io, join((
            "straight_modal_control", straight["reflection_power_fraction"],
            straight["all_mode_power_transmission"],
            straight["useful_fundamental_power_transmission"], 1.0, 0.0,
            1.0, NaN, true,
        ), ','))
        println(io, join((
            "equal_y_131mm", equal["reflection_power_fraction"],
            equal["all_mode_power_transmission"],
            equal["useful_fundamental_power_transmission"],
            equal["small_to_large_fundamental_amplitude_ratio"],
            equal["small_to_large_fundamental_phase_deg"], 1.0, NaN, false,
        ), ','))
        println(io, join((
            "calibrated_unequal_y_144.478mm", unequal["reflection_power_fraction"],
            unequal["all_mode_power_transmission"],
            unequal["useful_fundamental_power_transmission"],
            unequal["small_to_large_fundamental_amplitude_ratio"],
            unequal["small_to_large_fundamental_phase_deg"], NaN, NaN, false,
        ), ','))
        println(io, join((
            "global_15way_manifold", manifold["reflection_power_fraction"],
            manifold["all_mode_power_transmission"],
            manifold["useful_power_transmission"], NaN,
            manifold["maximum_fundamental_phase_error_deg"],
            manifold["target_complex_coherence"], NaN, false,
        ), ','))
        println(io, join((
            "mmi_reduced_ceiling", NaN, NaN, NaN, NaN, NaN,
            mmi.coherence, mmi.coherence, false,
        ), ','))
        println(io, join((
            "corporate_tree_optimistic", NaN, NaN, tree_efficiency, NaN, 0.0,
            NaN, routed_gain, passed,
        ), ','))
    end

    open(joinpath(OUTPUT_ROOT, "passive_feed_tree_upper_bound.csv"), "w") do io
        println(io, "equal_splitter_useful_efficiency,best_unequal_splitter_useful_efficiency,tree_useful_efficiency_upper_bound,minimum_tree_efficiency_for_pulse_gain_2,pulse_power_only_gain_upper_bound,carrier_power_only_gain_upper_bound,routed_profile_carrier_gain,passed")
        println(io, join((
            eta_equal, eta_other, tree_efficiency, required_efficiency,
            pulse_power_only_upper, carrier_power_only_upper, routed_gain, passed,
        ), ','))
    end

    centers_mm = collect(-57.4:8.2:57.4)
    target_relative = target ./ target[8]
    delivered_relative = sqrt.(leaf_power ./ leaf_power[8])
    open(joinpath(OUTPUT_ROOT, "passive_feed_leaf_profile.csv"), "w") do io
        println(io, "channel,y_mm,target_pressure_relative_center,optimistic_tree_pressure_relative_center,optimistic_tree_power_fraction")
        for index in eachindex(centers_mm)
            println(io, join((
                index, centers_mm[index], target_relative[index],
                delivered_relative[index], leaf_power[index],
            ), ','))
        end
    end

    profile_plot = plot(
        centers_mm, target_relative;
        marker=:circle, linewidth=2.5, label="target pulse30",
        xlabel="y, mm", ylabel="longitudinal amplitude / centre",
        title="Passive feed aperture profile",
        legend=:topright,
    )
    plot!(
        profile_plot, centers_mm, delivered_relative;
        marker=:diamond, linewidth=2.5, label="optimistic physical tree",
    )
    gain_values = [1.870, pulse_power_only_upper, VERIFIED_PULSE_GAIN]
    gain_plot = bar(
        ["unweighted\nlens", "passive tree\nupper", "segmented\ndrive"],
        gain_values;
        ylabel="G_peak", label="", title="Equal-input-power focal gain",
        ylim=(0.0, 2.25), color=[:steelblue, :darkorange, :seagreen],
    )
    hline!(gain_plot, [2.0]; color=:red, linestyle=:dash, linewidth=2, label="target 2x")
    for (index, value) in enumerate(gain_values)
        annotate!(gain_plot, index, value + 0.045, text(string(round(value; digits=3)), 10))
    end
    combined = plot(profile_plot, gain_plot; layout=(2, 1), size=(920, 820))
    savefig(combined, joinpath(OUTPUT_ROOT, "passive_feed_final_verdict.png"))

    verdict_path = joinpath(OUTPUT_ROOT, "passive_feed_verdict.txt")
    open(verdict_path, "w") do io
        println(io, passed ?
            "PASS: passive one-input feed retains enough useful power for G_peak >= 2." :
            "STOP: every tested passive one-input feed is below the G_peak = 2 gate.")
        println(io, "straight_modal_Tfund=$(straight["useful_fundamental_power_transmission"])")
        println(io, "equal_splitter_Tfund=$eta_equal")
        println(io, "unequal_splitter_Tfund=$eta_other")
        println(io, "unequal_ratio=$(unequal["small_to_large_fundamental_amplitude_ratio"])")
        println(io, "unequal_phase_deg=$(unequal["small_to_large_fundamental_phase_deg"])")
        println(io, "global_manifold_Tfund=$(manifold["useful_power_transmission"])")
        println(io, "global_manifold_coherence=$(manifold["target_complex_coherence"])")
        println(io, "mmi_reduced_best_coherence=$(mmi.coherence)")
        println(io, "tree_useful_efficiency_upper_bound=$tree_efficiency")
        println(io, "tree_required_efficiency=$required_efficiency")
        println(io, "pulse_power_only_gain_upper_bound=$pulse_power_only_upper")
        println(io, "carrier_power_only_gain_upper_bound=$carrier_power_only_upper")
        println(io, "routed_profile_carrier_gain=$routed_gain")
        println(io, "broadband_full_tree_solve=SKIPPED_BY_CARRIER_GATE")
    end
    JLD2.jldsave(
        joinpath(OUTPUT_ROOT, "passive_feed_verdict.jld2");
        format_version=1,
        eta_equal,
        eta_other,
        target_pressure_weights=target,
        optimistic_leaf_power=leaf_power,
        optimistic_routed_half_weights=routed_weights,
        tree_efficiency,
        required_efficiency,
        pulse_power_only_upper,
        carrier_power_only_upper,
        routed_gain,
        passed,
    )
    println("[+] passive tree efficiency upper=$tree_efficiency, required=$required_efficiency")
    println("[+] pulse power-only upper=$pulse_power_only_upper")
    println("[+] routed carrier gain=$routed_gain, passed=$passed")
    println("[+] $verdict_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
