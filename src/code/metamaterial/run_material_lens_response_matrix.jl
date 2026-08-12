module RunMaterialLensResponseMatrix

ENV["GKSwstype"] = "100"

using Gridap: Point
using JLD2
using LinearAlgebra
using Plots
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_RESPONSE_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "aluminium_material_lens_response_matrix_242khz"),
)
const DESIGN_ROOT = get(
    ENV,
    "METAMATERIALS_MATERIAL_LENS_RESPONSE_DESIGN",
    joinpath(PROJECT_ROOT, "tmp", "sinusoidal_material_lens_5cycle_242khz"),
)
const FREQUENCY_HZ = 242.0e3
const JOBS = parse(Int, get(ENV, "METAMATERIALS_MATERIAL_LENS_RESPONSE_JOBS", "2"))
const ROLES = (:selected, :uniform)
const GROUP_COUNT = 8
const CENTERS_MM = collect(0.0:8.2:57.4)
const PREOUT_X_MM = 33.0

ENV["METAMATERIALS_MATERIAL_LENS_DESIGN"] = DESIGN_ROOT
ENV["METAMATERIALS_MATERIAL_LENS_FEM_OUTPUT"] = OUTPUT_ROOT
ENV["METAMATERIALS_MATERIAL_LENS_FEM_ORDER"] = "1"
ENV["METAMATERIALS_MATERIAL_LENS_H_LENS_MM"] = "0.42"
ENV["METAMATERIALS_MATERIAL_LENS_H_OUTPUT_MM"] = "1.55"
ENV["METAMATERIALS_MATERIAL_LENS_SCAN"] = "false"

include(joinpath(@__DIR__, "run_material_lens_fem.jl"))
using .RunMaterialLensFEM

const Solver = RunMaterialLensFEM.MaterialLensHarmonicSolver

group_strip_indices(group) = group == 1 ? [8] : [8 - (group - 1), 8 + (group - 1)]
group_source_tags(group) = ["SourceStrip$(index)" for index in group_strip_indices(group)]

basis_path(role, group) = joinpath(
    OUTPUT_ROOT,
    "basis",
    "$(role)_group$(group)_$(round(Int, FREQUENCY_HZ))hz.jld2",
)

function build_meshes()
    for role in ROLES
        RunMaterialLensFEM.build_case_mesh(:aluminium, role)
    end
end

function solve_one(role, group)
    role in ROLES || throw(ArgumentError("unknown role"))
    1 <= group <= GROUP_COUNT || throw(ArgumentError("group outside 1:8"))
    path = basis_path(role, group)
    isfile(path) && return println("[=] Existing basis: $path")
    mesh_info = RunMaterialLensFEM.build_case_mesh(:aluminium, role)
    mesh_path = mesh_info isa AbstractString ? mesh_info : mesh_info.mesh_path
    material = RunMaterialLensFEM.aluminium_6061()
    solved = Solver.solve_material_lens_harmonic(
        mesh_path,
        material,
        material;
        config=Solver.MaterialLensHarmonicConfig(
            frequency_hz=FREQUENCY_HZ,
            focus_x_mm=88.0,
            element_order=1,
            quadrature_degree=2,
            scan_before_focus_mm=0.0,
            scan_after_focus_mm=0.0,
            scan_transverse_half_width_mm=0.0,
            scan_step_mm=2.0,
            source_excitation_tags=group_source_tags(group),
        ),
    )
    preout_ux_m = ComplexF64[
        solved.displacement(Point(PREOUT_X_MM * 1.0e-3, y_mm * 1.0e-3))[1]
        for y_mm in CENTERS_MM
    ]
    mkpath(dirname(path))
    jldsave(
        path;
        format_version=1,
        role=String(role),
        group,
        strip_indices=group_strip_indices(group),
        source_tags=group_source_tags(group),
        frequency_hz=FREQUENCY_HZ,
        focus_ux_m=solved.focus_displacement[1],
        preout_x_mm=PREOUT_X_MM,
        preout_y_mm=CENTERS_MM,
        preout_ux_m,
        mesh_path=abspath(mesh_path),
    )
    println("[+] $role group $group: |focus ux|=$(abs(solved.focus_displacement[1])) m")
    println("[+] $path")
end

function child_command(role, group)
    julia = Base.julia_cmd()
    `$julia --startup-file=no --project=$PROJECT_ROOT --threads=1 $(@__FILE__) --solve-one=$(role):$(group)`
end

function solve_all()
    build_meshes()
    pending = [
        (; role, group)
        for role in ROLES for group in 1:GROUP_COUNT
        if !isfile(basis_path(role, group))
    ]
    isempty(pending) && return println("[=] All strip-group bases exist")
    jobs = min(max(JOBS, 1), length(pending), Sys.CPU_THREADS)
    semaphore = Base.Semaphore(jobs)
    println("=== Aluminium strip-group response matrix: $(length(pending)) solves, $jobs processes ===")
    @sync for spec in pending
        @async begin
            Base.acquire(semaphore)
            try
                run(child_command(spec.role, spec.group))
            finally
                Base.release(semaphore)
            end
        end
    end
end

function wrap_phase(value)
    mod(value + pi, 2pi) - pi
end

function load_role(role)
    data = [load(basis_path(role, group)) for group in 1:GROUP_COUNT]
    focus = ComplexF64[item["focus_ux_m"] for item in data]
    preout = hcat((ComplexF64.(item["preout_ux_m"]) for item in data)...)
    (; role, focus, preout)
end

function offdiagonal_ratio(matrix)
    diagonal_norm = norm(diag(matrix))
    offdiagonal = copy(matrix)
    for index in axes(offdiagonal, 1)
        offdiagonal[index, index] = 0
    end
    norm(offdiagonal) / diagonal_norm
end

function write_channel_csv(selected, uniform)
    common_phase = angle(sum(selected.focus))
    path = joinpath(OUTPUT_ROOT, "aluminium_strip_group_focus_contributions.csv")
    open(path, "w") do io
        println(io, "group,abs_y_mm,strip_indices,selected_focus_ux_real_m,selected_focus_ux_imag_m,selected_focus_ux_abs_m,selected_phase_error_deg,uniform_focus_ux_real_m,uniform_focus_ux_imag_m,uniform_focus_ux_abs_m,device_amplitude_ratio,device_phase_deg")
        for group in 1:GROUP_COUNT
            device = selected.focus[group] / uniform.focus[group]
            println(io, join((
                group, CENTERS_MM[group], join(group_strip_indices(group), ';'),
                real(selected.focus[group]), imag(selected.focus[group]), abs(selected.focus[group]),
                rad2deg(wrap_phase(angle(selected.focus[group]) - common_phase)),
                real(uniform.focus[group]), imag(uniform.focus[group]), abs(uniform.focus[group]),
                abs(device), rad2deg(angle(device)),
            ), ','))
        end
    end
    path
end

function make_plot(selected, uniform, metrics)
    common_phase = angle(sum(selected.focus))
    selected_error = rad2deg.([wrap_phase(angle(value) - common_phase) for value in selected.focus])
    uniform_error = rad2deg.([wrap_phase(angle(value) - angle(sum(uniform.focus))) for value in uniform.focus])
    amplitude = plot(
        CENTERS_MM, abs.(selected.focus) .* 1.0e9;
        xlabel="|y|, mm", ylabel="group |uₓ(focus)|, nm",
        title="Symmetric group contributions", marker=:circle, lw=2.5, label="lens",
    )
    plot!(amplitude, CENTERS_MM, abs.(uniform.focus) .* 1.0e9;
          marker=:square, lw=2.5, ls=:dash, label="straight")

    phase = plot(
        CENTERS_MM, selected_error;
        xlabel="|y|, mm", ylabel="phase error to total, deg",
        title="Focal phase coherence", marker=:circle, lw=2.5, label="lens",
    )
    plot!(phase, CENTERS_MM, uniform_error;
          marker=:square, lw=2.5, ls=:dash, label="straight")
    hline!(phase, [0.0]; color=:black, ls=:dot, label=false)

    device = selected.focus ./ uniform.focus
    transfer = plot(
        CENTERS_MM, abs.(device);
        xlabel="|y|, mm", ylabel="|C_lens/C_straight|",
        title="Per-group device transfer", marker=:circle, lw=2.5,
        label="amplitude",
    )
    transfer_phase = twinx(transfer)
    plot!(transfer_phase, CENTERS_MM, rad2deg.(angle.(device));
          ylabel="phase, deg", marker=:square, color=:darkorange, lw=2,
          label="phase")

    gains = bar(
        ["current", "phase-aligned", "target"],
        [metrics.current_gain, metrics.phase_aligned_gain, 2.0];
        ylabel="focus gain over straight", title="Carrier gain ceiling",
        color=[:steelblue, :goldenrod, :gray60], label=false,
    )

    path = joinpath(OUTPUT_ROOT, "aluminium_strip_group_response.png")
    savefig(plot(
        amplitude, phase, transfer, gains;
        layout=(2, 2), size=(1400, 950), margin=5Plots.mm,
    ), path)
    path
end

function analyze()
    all(isfile(basis_path(role, group)) for role in ROLES for group in 1:GROUP_COUNT) ||
        error("response basis is incomplete")
    selected = load_role(:selected)
    uniform = load_role(:uniform)
    selected_total = sum(selected.focus)
    uniform_total = sum(uniform.focus)
    current_gain = abs(selected_total) / abs(uniform_total)
    phase_aligned_gain = sum(abs, selected.focus) / abs(uniform_total)
    coherence_efficiency = abs(selected_total) / sum(abs, selected.focus)
    common_phase = angle(selected_total)
    phase_errors_deg = rad2deg.([
        wrap_phase(angle(value) - common_phase) for value in selected.focus
    ])
    metrics = (
        frequency_hz=FREQUENCY_HZ,
        selected_focus_ux_m=abs(selected_total),
        uniform_focus_ux_m=abs(uniform_total),
        current_gain,
        phase_aligned_gain,
        coherence_efficiency,
        maximum_abs_phase_error_deg=maximum(abs, phase_errors_deg),
        selected_preout_offdiagonal_l2_over_diagonal=offdiagonal_ratio(selected.preout),
        uniform_preout_offdiagonal_l2_over_diagonal=offdiagonal_ratio(uniform.preout),
    )
    summary_path = joinpath(OUTPUT_ROOT, "aluminium_strip_group_response_summary.csv")
    open(summary_path, "w") do io
        columns = propertynames(metrics)
        println(io, join(string.(columns), ','))
        println(io, join((getproperty(metrics, column) for column in columns), ','))
    end
    channel_path = write_channel_csv(selected, uniform)
    plot_path = make_plot(selected, uniform, metrics)
    matrix_path = joinpath(OUTPUT_ROOT, "aluminium_strip_group_response_matrix.jld2")
    jldsave(matrix_path; format_version=1, selected, uniform, metrics, centers_mm=CENTERS_MM)
    println("[+] Current carrier gain: $(metrics.current_gain)")
    println("[+] Phase-aligned same-amplitude ceiling: $(metrics.phase_aligned_gain)")
    println("[+] Coherence efficiency: $(metrics.coherence_efficiency)")
    println("[+] Max focal phase error: $(metrics.maximum_abs_phase_error_deg) deg")
    println("[+] Summary: $summary_path")
    println("[+] Channels: $channel_path")
    println("[+] Plot: $plot_path")
    println("[+] Matrix: $matrix_path")
end

function parse_solve_one(argument)
    value = split(argument, '='; limit=2)[2]
    role_text, group_text = split(value, ':'; limit=2)
    Symbol(role_text), parse(Int, group_text)
end

function main(args=ARGS)
    if isempty(args) || args == ["--all"]
        solve_all()
        analyze()
    elseif args == ["--mesh"]
        build_meshes()
    elseif args == ["--solve"]
        solve_all()
    elseif args == ["--analyze"]
        analyze()
    elseif length(args) == 1 && startswith(only(args), "--solve-one=")
        role, group = parse_solve_one(only(args))
        solve_one(role, group)
    else
        error("usage: run_material_lens_response_matrix.jl --mesh|--solve|--analyze|--all")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
