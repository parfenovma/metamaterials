module HornBinaryLens

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PILOT_ROOT = joinpath(PROJECT_ROOT, "tmp", "horn_point_radiator_pilot")
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_HORN_LENS_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "horn_binary_lens_242khz"),
)
const REQUESTED_STAGE = length(ARGS) == 1 && startswith(ARGS[1], "--stage=") ?
                        split(ARGS[1], '='; limit=2)[2] : nothing
const FEM_ORDER = parse(Int, get(ENV, "METAMATERIALS_HORN_LENS_ORDER", "1"))
const FINE_MESH = lowercase(get(ENV, "METAMATERIALS_HORN_LENS_FINE_MESH", "false")) == "true"

include(joinpath(@__DIR__, "lens_design.jl"))
using .LensDesign
include(joinpath(@__DIR__, "point_radiator_aperture_mesher.jl"))
using .PointRadiatorApertureMesher

if REQUESTED_STAGE in ("design", "solve", "plot")
    using JLD2
end
if REQUESTED_STAGE == "design"
    include(joinpath(@__DIR__, "spectral_analysis.jl"))
    using .SpectralAnalysis
end
if REQUESTED_STAGE == "solve"
    include(joinpath(@__DIR__, "sinusoidal_material_lens.jl"))
    include(joinpath(@__DIR__, "measured_aperture_harmonic_solver.jl"))
    using .SinusoidalMaterialLens
    using .MeasuredApertureHarmonicSolver
end
if REQUESTED_STAGE in ("design", "plot")
    include(joinpath(@__DIR__, "impulse_risk_analysis.jl"))
    using .ImpulseRiskAnalysis
end
if REQUESTED_STAGE == "plot"
    ENV["GKSwstype"] = "100"
    using Plots
end

const APERTURE_CONFIG = PointRadiatorApertureConfig()
const LENS_CONFIG = LensConfig(
    frequency_hz=242.0e3,
    longitudinal_speed_m_per_s=2340.0,
    element_count=15,
    element_width_mm=1.6,
    slot_width_mm=3.2,
    focal_distance_mm=35.0,
    aperture_quadrature_points=7,
)

mesh_path() = joinpath(OUTPUT_ROOT, "mesh_point_radiator_aperture.msh")
selection_path() = joinpath(OUTPUT_ROOT, "horn_binary_selection.csv")
design_path() = joinpath(OUTPUT_ROOT, "horn_binary_design.jld2")
result_path(role) = joinpath(OUTPUT_ROOT, "$(role)_fem_order$(FEM_ORDER).jld2")

function child_command(stage)
    julia = joinpath(Sys.BINDIR, Base.julia_exename())
    `$julia --startup-file=no --project=$PROJECT_ROOT $(@__FILE__) --stage=$stage`
end

function probe_transfer(case_name, probe_name)
    data = JLD2.load(joinpath(PILOT_ROOT, "signals", "$(case_name).jld2"))
    index = findfirst(==(probe_name), data["probe_names"])
    isnothing(index) && error("probe $probe_name not found in $case_name")
    signal = vec(data["probe_velocity_x_m_per_s"][index, :])
    analyze_transfer(
        data["time_s"],
        data["source_drive_mpa"],
        signal;
        config=SpectrumConfig(
            window=:rectangular,
            zero_padding_factor=16,
            input_floor_relative=1e-3,
        ),
    )
end

function measured_delay_coefficient()
    device = probe_transfer("collector_device", "near_radiator")
    gentle = probe_transfer("collector_gentle", "near_radiator")
    value_at_frequency(relative_transfer(gentle, device), LENS_CONFIG.frequency_hz)
end

function uniform_selection(entry)
    centers = lens_centers_mm(LENS_CONFIG)
    LensSelection(centers, fill(entry, length(centers)), 0.0 + 0.0im)
end

function transverse_metrics(selection)
    y_mm = collect(-20.0:0.1:20.0)
    amplitude = abs.([
        field_at(LENS_CONFIG.focal_distance_mm, y, selection; config=LENS_CONFIG)
        for y in y_mm
    ])
    width = contiguous_width(y_mm, amplitude)
    y_mm, amplitude, width
end

function run_design_stage()
    mkpath(OUTPUT_ROOT)
    measured = measured_delay_coefficient()
    delay = measured.amplitude * cis(measured.phase_rad)
    entries = [
        LibraryEntry("device", 0, 26.0, 0.0, 0.0, 1.0 + 0.0im),
        LibraryEntry("gentle", 0, 30.42, 0.0, 0.0, delay),
    ]
    selection = select_lens(entries; config=LENS_CONFIG)
    uniform = uniform_selection(first(entries))
    focus = abs(field_at(35.0, 0.0, selection; config=LENS_CONFIG))
    uniform_focus = abs(field_at(35.0, 0.0, uniform; config=LENS_CONFIG))
    y_mm, transverse, width = transverse_metrics(selection)
    _, uniform_transverse, uniform_width = transverse_metrics(uniform)
    weights = ComplexF64[getproperty(entry, :transfer) for entry in selection.entries]
    states = String[getproperty(entry, :case_id) for entry in selection.entries]

    open(selection_path(), "w") do io
        println(io, "element_index,center_y_mm,state,weight_real,weight_imag,weight_amplitude,weight_phase_deg")
        for index in eachindex(weights)
            println(io, join((
                index,
                selection.centers_mm[index],
                states[index],
                real(weights[index]),
                imag(weights[index]),
                abs(weights[index]),
                rad2deg(angle(weights[index])),
            ), ','))
        end
    end
    jldsave(
        design_path();
        format_version=1,
        frequency_hz=LENS_CONFIG.frequency_hz,
        measured_frequency_hz=measured.frequency_hz,
        delay_amplitude=measured.amplitude,
        delay_phase_rad=measured.phase_rad,
        delay_group_delay_s=measured.group_delay_s,
        centers_mm=selection.centers_mm,
        states,
        weights,
        reduced_focus_amplitude=focus,
        reduced_uniform_focus_amplitude=uniform_focus,
        reduced_focus_gain=focus / uniform_focus,
        reduced_transverse_fwhm_mm=width.width,
        reduced_uniform_transverse_fwhm_mm=uniform_width.width,
        transverse_y_mm=y_mm,
        reduced_transverse_amplitude=transverse,
        reduced_uniform_transverse_amplitude=uniform_transverse,
    )
    println("[+] measured gentle/device=$(measured.amplitude) ∠ $(rad2deg(measured.phase_rad)) deg")
    println("[+] reduced gain=$(focus / uniform_focus), FWHM=$(width.width) mm")
    println("[+] $(selection_path())")
end

function run_mesh_stage()
    use_fine = FEM_ORDER == 1 || FINE_MESH
    build_point_radiator_aperture_mesh(
        mesh_path();
        config=APERTURE_CONFIG,
        size_radiator_mm=use_fine ? 0.25 : 0.40,
        size_focus_mm=use_fine ? 0.55 : 0.85,
        size_max_mm=use_fine ? 0.90 : 1.35,
    )
end

function selection_weights()
    rows = readlines(selection_path())
    ComplexF64[
        let values=split(line, ',')
            complex(parse(Float64, values[4]), parse(Float64, values[5]))
        end
        for line in Iterators.drop(rows, 1) if !isempty(strip(line))
    ]
end

function solve_role(role)
    weights = role == :lens ? selection_weights() :
              role == :uniform ? ones(ComplexF64, APERTURE_CONFIG.element_count) :
              error("unknown role: $role")
    result = solve_measured_aperture_harmonic(
        mesh_path(),
        weights,
        photopolymer();
        config=MeasuredApertureHarmonicConfig(
            frequency_hz=242.0e3,
            focal_distance_mm=35.0,
            element_order=FEM_ORDER,
            quadrature_degree=2FEM_ORDER,
            scan_step_mm=0.5,
        ),
    )
    jldsave(
        result_path(role);
        format_version=1,
        role=String(role),
        element_order=FEM_ORDER,
        weights,
        focus_displacement=result.focus_displacement,
        focus_longitudinal_amplitude_m=result.focus_longitudinal_amplitude_m,
        focus_total_amplitude_m=result.focus_total_amplitude_m,
        scan_x_mm=result.scan_x_mm,
        scan_y_mm=result.scan_y_mm,
        scan_ux_m=result.scan_ux_m,
        scan_uy_m=result.scan_uy_m,
        scan_total_amplitude_m=result.scan_total_amplitude_m,
    )
    println("[+] $role |ux(focus)|=$(result.focus_longitudinal_amplitude_m) m")
end

function run_solve_stage()
    solve_role(:lens)
    solve_role(:uniform)
end

function fem_metrics(data)
    x_mm = data["scan_x_mm"]
    y_mm = data["scan_y_mm"]
    ux = abs.(data["scan_ux_m"])
    focus_x = argmin(abs.(x_mm .- 35.0))
    axis_y = argmin(abs.(y_mm))
    transverse = ux[:, focus_x]
    axial = ux[axis_y, :]
    transverse_width = contiguous_width(y_mm, transverse)
    target_index = argmin(abs.(x_mm .- 35.0))
    target_threshold = axial[target_index] / sqrt(2.0)
    left, right = target_index, target_index
    while left > firstindex(axial) && axial[left - 1] >= target_threshold
        left -= 1
    end
    while right < lastindex(axial) && axial[right + 1] >= target_threshold
        right += 1
    end
    local_peak_offset = argmax(axial[left:right]) - 1
    local_peak_index = left + local_peak_offset
    center_peak = argmax(transverse)
    left_minimum = center_peak
    while left_minimum > firstindex(transverse) + 1
        left_minimum -= 1
        transverse[left_minimum] <= transverse[left_minimum - 1] &&
            transverse[left_minimum] <= transverse[left_minimum + 1] && break
    end
    right_minimum = center_peak
    while right_minimum < lastindex(transverse) - 1
        right_minimum += 1
        transverse[right_minimum] <= transverse[right_minimum - 1] &&
            transverse[right_minimum] <= transverse[right_minimum + 1] && break
    end
    outside_main_lobe = vcat(
        collect(firstindex(transverse):(left_minimum - 1)),
        collect((right_minimum + 1):lastindex(transverse)),
    )
    sidelobe_ratio = isempty(outside_main_lobe) ? 0.0 :
                     maximum(transverse[outside_main_lobe]) / transverse_width.peak_amplitude
    (
        transverse,
        axial,
        transverse_fwhm_mm=transverse_width.width,
        axial_dof_mm=x_mm[right] - x_mm[left],
        peak_x_mm=x_mm[local_peak_index],
        peak_y_mm=transverse_width.peak_coordinate,
        sidelobe_amplitude_ratio=sidelobe_ratio,
    )
end

function run_plot_stage()
    lens = JLD2.load(result_path(:lens))
    uniform = JLD2.load(result_path(:uniform))
    lens_metrics = fem_metrics(lens)
    uniform_metrics = fem_metrics(uniform)
    gain = lens["focus_longitudinal_amplitude_m"] / uniform["focus_longitudinal_amplitude_m"]

    maps = Any[]
    for (data, title) in ((lens, "measured binary lens"), (uniform, "uniform point aperture"))
        panel = heatmap(
            data["scan_x_mm"], data["scan_y_mm"], data["scan_total_amplitude_m"] .* 1e9;
            xlabel="x, mm", ylabel="y, mm", title=title,
            color=:viridis, colorbar_title="|u|, nm", aspect_ratio=:equal, grid=false,
        )
        scatter!(panel, [35.0], [0.0]; marker=:xcross, color=:white,
                 markerstrokewidth=2, label=false)
        push!(maps, panel)
    end
    map_path = joinpath(OUTPUT_ROOT, "horn_binary_lens_fem_maps_order$(FEM_ORDER).png")
    savefig(plot(maps...; layout=(2, 1), size=(1050, 920), left_margin=7Plots.mm), map_path)

    profile = plot(
        lens["scan_y_mm"], lens_metrics.transverse .* 1e9;
        linewidth=2.5, label="binary lens", xlabel="y at x=35 mm, mm",
        ylabel="|ux|, nm", title="FEM transverse focus profile", gridalpha=0.25,
    )
    plot!(profile, uniform["scan_y_mm"], uniform_metrics.transverse .* 1e9;
          linewidth=2, label="uniform aperture")
    profile_path = joinpath(OUTPUT_ROOT, "horn_binary_lens_profiles_order$(FEM_ORDER).png")
    savefig(profile, profile_path)

    summary_path = joinpath(OUTPUT_ROOT, "horn_binary_lens_fem_summary_order$(FEM_ORDER).csv")
    open(summary_path, "w") do io
        println(io, "focus_ux_m,uniform_focus_ux_m,focus_gain,transverse_fwhm_mm,uniform_transverse_fwhm_mm,axial_dof_mm,peak_x_mm,peak_y_mm,sidelobe_amplitude_ratio")
        println(io, join((
            lens["focus_longitudinal_amplitude_m"],
            uniform["focus_longitudinal_amplitude_m"],
            gain,
            lens_metrics.transverse_fwhm_mm,
            uniform_metrics.transverse_fwhm_mm,
            lens_metrics.axial_dof_mm,
            lens_metrics.peak_x_mm,
            lens_metrics.peak_y_mm,
            lens_metrics.sidelobe_amplitude_ratio,
        ), ','))
    end
    println("[+] FEM gain=$gain, transverse FWHM=$(lens_metrics.transverse_fwhm_mm) mm")
    println("[+] $map_path")
    println("[+] $summary_path")
end

function run_stage(stage)
    stage == "design" && return run_design_stage()
    stage == "mesh" && return run_mesh_stage()
    stage == "solve" && return run_solve_stage()
    stage == "plot" && return run_plot_stage()
    error("unknown horn-binary-lens stage: $stage")
end

function main(args=ARGS)
    if isempty(args)
        for stage in ("design", "mesh", "solve", "plot")
            println("\n=== Horn binary lens stage: $stage ===")
            run(child_command(stage))
        end
    elseif !isnothing(REQUESTED_STAGE)
        run_stage(REQUESTED_STAGE)
    else
        error("usage: julia run_horn_binary_lens.jl [--stage=...] ")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
