module RunAluminiumHornTTDPassiveSplitter

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DESIGN_PATH = let provisional = joinpath(
        PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_feed", "passive_feed_design.jld2",
    )
    isfile(provisional) ? provisional : joinpath(
        PROJECT_ROOT, "results", "aluminium_horn_ttd_p6_242khz", "passive_feed",
        "design", "passive_feed_design.jld2",
    )
end

function argument(name, default=nothing)
    prefix = "--$name="
    match = findfirst(value -> startswith(value, prefix), ARGS)
    isnothing(match) ? default : split(ARGS[match], '='; limit=2)[2]
end

const STAGE = argument("stage", "analyze")
const VARIANT = argument("variant", "worst")
const FREQUENCY_HZ = parse(Float64, argument("frequency-hz", "242000"))
const FREQUENCIES_HZ = Float64[226.04e3, 242.0e3, 257.94e3]
const OUTPUT_ROOT = if VARIANT == "worst"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_splitter")
elseif VARIANT == "equal"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_equal_splitter")
elseif VARIANT == "equal_long"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_equal_long_splitter")
elseif VARIANT == "equal_compact"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_equal_compact_splitter")
elseif VARIANT == "equal_rounded"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_equal_rounded_splitter")
elseif VARIANT == "calibrated"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_calibrated_splitter")
elseif VARIANT == "calibrated_final"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_calibrated_final_splitter")
elseif VARIANT == "calibrated_beat"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_calibrated_beat_splitter")
elseif VARIANT == "calibrated_beat_long"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_calibrated_beat_long_splitter")
elseif VARIANT == "straight"
    joinpath(PROJECT_ROOT, "tmp", "aluminium_horn_ttd_passive_straight_control")
else
    error("unknown splitter variant: $VARIANT")
end

using JLD2
include(joinpath(@__DIR__, "passive_feed_splitter_mesher.jl"))
using .PassiveFeedSplitterMesher

const FEED_DESIGN = JLD2.load(DESIGN_PATH)
const CONFIG = if VARIANT == "worst"
    PassiveFeedSplitterConfig(
        small_power_fraction=Float64(FEED_DESIGN["worst_small_power_fraction"]),
        transformer_length_mm=Float64(FEED_DESIGN["common_transformer_length_mm"]),
    )
elseif VARIANT == "equal"
    PassiveFeedSplitterConfig(
        small_power_fraction=0.499999,
        transformer_length_mm=47.0,
    )
elseif VARIANT == "equal_long"
    PassiveFeedSplitterConfig(
        small_power_fraction=0.499999,
        transformer_length_mm=131.0,
    )
elseif VARIANT == "equal_compact"
    PassiveFeedSplitterConfig(
        small_power_fraction=0.499999,
        transformer_length_mm=131.0,
        output_center_offset_mm=4.0,
    )
elseif VARIANT == "equal_rounded"
    PassiveFeedSplitterConfig(
        small_power_fraction=0.499999,
        transformer_length_mm=131.0,
        output_center_offset_mm=4.5,
        splitter_tip_radius_mm=0.75,
        profile_samples=600,
    )
elseif VARIANT == "calibrated"
    PassiveFeedSplitterConfig(
        small_power_fraction=0.2183563961252619,
        transformer_length_mm=131.0,
    )
elseif VARIANT == "calibrated_final"
    PassiveFeedSplitterConfig(
        small_power_fraction=0.1955777738982686,
        transformer_length_mm=131.0,
    )
elseif VARIANT == "calibrated_beat"
    PassiveFeedSplitterConfig(
        small_power_fraction=0.1955777738982686,
        transformer_length_mm=144.4780677722518,
    )
elseif VARIANT == "calibrated_beat_long"
    PassiveFeedSplitterConfig(
        small_power_fraction=0.1955777738982686,
        transformer_length_mm=209.9212009074102,
    )
else
    PassiveFeedSplitterConfig(
        small_power_fraction=0.499999,
        transformer_length_mm=47.0,
    )
end
const CASCADE_DEPTH = VARIANT == "worst" ?
    Int(FEED_DESIGN["maximum_tree_depth"]) : VARIANT == "straight" ? 1 : 9

if STAGE == "mesh"
    using Gmsh: gmsh
elseif STAGE == "solve"
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
    include(joinpath(@__DIR__, "passive_feed_splitter_harmonic_solver.jl"))
    using .SinusoidalMaterialLens
    using .PassiveFeedSplitterHarmonicSolver
end

mesh_path() = joinpath(OUTPUT_ROOT, "passive_feed_splitter.msh")
result_path(frequency_hz=FREQUENCY_HZ) = joinpath(
    OUTPUT_ROOT, "passive_feed_splitter_f$(round(Int, frequency_hz))hz.jld2",
)

function run_mesh_stage()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        result = VARIANT == "straight" ?
            build_passive_feed_straight_control_mesh(mesh_path(); config=CONFIG) :
            build_passive_feed_splitter_mesh(mesh_path(); config=CONFIG)
        mkpath(OUTPUT_ROOT)
        open(joinpath(OUTPUT_ROOT, "passive_feed_splitter_mesh_summary.csv"), "w") do io
            println(io, "nodes,elements,input_width_mm,small_junction_width_mm,large_junction_width_mm,output_width_mm,transformer_length_mm,output_center_offset_mm,target_output_amplitude_ratio,small_adiabatic_parameter,large_adiabatic_parameter")
            println(io, join((
                result.node_count, result.element_count, CONFIG.input_width_mm,
                small_junction_width_mm(CONFIG), large_junction_width_mm(CONFIG),
                CONFIG.output_width_mm, CONFIG.transformer_length_mm,
                CONFIG.output_center_offset_mm, target_output_amplitude_ratio(CONFIG),
                adiabatic_parameter(CONFIG, small_junction_width_mm(CONFIG)),
                adiabatic_parameter(CONFIG, large_junction_width_mm(CONFIG)),
            ), ','))
        end
    finally
        gmsh.finalize()
    end
end

function run_solve_stage()
    result = solve_passive_feed_splitter_harmonic(
        mesh_path(), aluminium_6061();
        config=PassiveFeedSplitterHarmonicConfig(frequency_hz=FREQUENCY_HZ),
    )
    mkpath(OUTPUT_ROOT)
    JLD2.jldsave(
        result_path();
        format_version=1,
        frequency_hz=FREQUENCY_HZ,
        splitter_variant=VARIANT,
        target_output_amplitude_ratio=target_output_amplitude_ratio(CONFIG),
        small_displacement_m=result.small_displacement_m,
        large_displacement_m=result.large_displacement_m,
        active_input_power_w_per_m=result.active_input_power_w_per_m,
        small_absorbed_power_w_per_m=result.small_absorbed_power_w_per_m,
        large_absorbed_power_w_per_m=result.large_absorbed_power_w_per_m,
        small_modal_power_proxy_w_per_m=result.small_modal_power_proxy_w_per_m,
        large_modal_power_proxy_w_per_m=result.large_modal_power_proxy_w_per_m,
    )
    println("[+] splitter $(FREQUENCY_HZ / 1e3) kHz: small/large |ux|=" *
            "$(abs(result.small_displacement_m[1] / result.large_displacement_m[1]))")
    println("[+] input/output absorbed power=$(result.active_input_power_w_per_m) / " *
            "$(result.small_absorbed_power_w_per_m + result.large_absorbed_power_w_per_m) W/m")
    println("[+] output fundamental-mode proxy=" *
            "$(result.small_modal_power_proxy_w_per_m + result.large_modal_power_proxy_w_per_m) W/m")
    println("[+] $(result_path())")
end

wrap_phase_deg(value) = rad2deg(atan(sin(angle(value)), cos(angle(value))))

function run_analyze_stage()
    data = JLD2.load.(result_path.(FREQUENCIES_HZ))
    small = ComplexF64[item["small_displacement_m"][1] for item in data]
    large = ComplexF64[item["large_displacement_m"][1] for item in data]
    small_y = ComplexF64[item["small_displacement_m"][2] for item in data]
    large_y = ComplexF64[item["large_displacement_m"][2] for item in data]
    input_power = Float64[item["active_input_power_w_per_m"] for item in data]
    output_absorbed_power = Float64[
        item["small_absorbed_power_w_per_m"] +
        item["large_absorbed_power_w_per_m"] for item in data
    ]
    output_modal_power = Float64[
        item["small_modal_power_proxy_w_per_m"] +
        item["large_modal_power_proxy_w_per_m"] for item in data
    ]
    ratio = abs.(small ./ large)
    phase_deg = wrap_phase_deg.(small ./ large)
    transmission = output_absorbed_power ./ input_power
    fundamental_mode_fraction = output_modal_power ./ output_absorbed_power
    transverse = (abs2.(small_y) .+ abs2.(large_y)) ./
                 (abs2.(small) .+ abs2.(large))
    target = target_output_amplitude_ratio(CONFIG)
    ratio_error = ratio ./ target .- 1
    total_efficiency_requirement = Float64(
        FEED_DESIGN["minimum_total_feed_power_efficiency_for_gain_2"],
    )
    per_stage_efficiency_requirement = total_efficiency_requirement^(
        1.0 / CASCADE_DEPTH
    )
    passed = maximum(abs, ratio_error) <= 0.05 &&
             maximum(abs, phase_deg) <= 3.0 &&
             minimum(transmission) >= per_stage_efficiency_requirement &&
             maximum(transverse) <= 0.10

    rows_path = joinpath(OUTPUT_ROOT, "passive_feed_splitter_response.csv")
    open(rows_path, "w") do io
        println(io, "frequency_hz,target_amplitude_ratio,measured_amplitude_ratio,amplitude_ratio_relative_error,small_large_phase_deg,input_power_w_per_m,output_absorbed_power_w_per_m,output_modal_power_proxy_w_per_m,power_transmission,fundamental_mode_power_fraction,transverse_displacement_ratio")
        for index in eachindex(FREQUENCIES_HZ)
            println(io, join((
                FREQUENCIES_HZ[index], target, ratio[index], ratio_error[index],
                phase_deg[index], input_power[index], output_absorbed_power[index],
                output_modal_power[index], transmission[index],
                fundamental_mode_fraction[index], transverse[index],
            ), ','))
        end
    end
    summary_path = joinpath(OUTPUT_ROOT, "passive_feed_splitter_summary.csv")
    open(summary_path, "w") do io
        println(io, "target_amplitude_ratio,maximum_amplitude_ratio_error,maximum_phase_error_deg,minimum_power_transmission,minimum_fundamental_mode_power_fraction,maximum_transverse_displacement_ratio,total_feed_efficiency_requirement,per_stage_efficiency_requirement,maximum_tree_depth,passed")
        println(io, join((
            target, maximum(abs, ratio_error), maximum(abs, phase_deg),
            minimum(transmission), minimum(fundamental_mode_fraction),
            maximum(transverse), total_efficiency_requirement,
            per_stage_efficiency_requirement,
            CASCADE_DEPTH, passed,
        ), ','))
    end
    verdict_path = joinpath(OUTPUT_ROOT, "passive_feed_splitter_verdict.txt")
    open(verdict_path, "w") do io
        println(io, passed ?
            "PASS: worst corporate-feed splitter is phase-neutral and efficient enough." :
            "STOP: worst corporate-feed splitter fails amplitude/phase/efficiency gate.")
        println(io, "target_ratio=$target")
        println(io, "measured_ratio=$(join(ratio, ';'))")
        println(io, "phase_deg=$(join(phase_deg, ';'))")
        println(io, "transmission=$(join(transmission, ';'))")
    end
    println("[+] splitter ratios=$ratio, target=$target")
    println("[+] splitter phases=$phase_deg deg")
    println("[+] splitter power transmission=$transmission")
    println("[+] per-stage requirement=$per_stage_efficiency_requirement, passed=$passed")
    println("[+] $summary_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    STAGE == "mesh" ? run_mesh_stage() :
    STAGE == "solve" ? run_solve_stage() :
    STAGE == "analyze" ? run_analyze_stage() :
    error("unknown stage: $STAGE")
end

end
