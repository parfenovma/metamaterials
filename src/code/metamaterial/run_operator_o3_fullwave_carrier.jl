module RunOperatorO3FullWaveCarrier

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots
using Printf

include(joinpath(@__DIR__, "ndt_benchmark.jl"))
include(joinpath(@__DIR__, "ndt_element_array_mesher.jl"))
include(joinpath(@__DIR__, "ndt_fullwave_harmonic_solver.jl"))

using .NDTBenchmark
using .NDTElementArrayMesher
using .NDTFullWaveHarmonicSolver

export project_operator,
       architecture_features,
       relative_residual,
       build_mesh_pair,
       run_campaign

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_O3_FULLWAVE_OUTPUT",
    joinpath(
        PROJECT_ROOT,
        "results",
        "operator_space",
        "fullwave_carrier_commonmesh_h1p0_v0",
    ),
)
const STATIC_FRONTIER_PATH = joinpath(
    PROJECT_ROOT, "results", "operator_space", "static_frontier_v0", "static_frontier.jld2",
)
const ELEMENT_COUNT = 17
const CARRIER_HZ = 242.0e3
const STATIC_K4Q2_ID = :trained_k4_q2_unconstrained
const STATIC_K3Q4_ID = :trained_k3_q4_reflection_balanced
const MOTION_P12_INDICES = collect(3:14)
const MOTION_P17_INDICES = collect(1:17)

function stand_config()
    NDTStandConfig(
        frequency_hz=CARRIER_HZ,
        mesh_size_mm=parse(Float64, get(ENV, "METAMATERIALS_O3_F0_MESH_MM", "1.0")),
        defect_mesh_size_mm=parse(
            Float64, get(ENV, "METAMATERIALS_O3_F0_DEFECT_MESH_MM", "0.25"),
        ),
    )
end

mesh_path(case_id) = joinpath(OUTPUT_ROOT, "meshes", "$(case_id).msh")
response_path(case_id) = joinpath(OUTPUT_ROOT, "responses", "$(case_id)_242000hz.jld2")

function build_mesh_pair(; force=false)
    config = stand_config()
    defects = Dict(defect.id => defect for defect in required_defects(config))
    free_path = mesh_path(:none)
    void_path = mesh_path(:void_r2_focus)
    if force || !isfile(free_path) || !isfile(void_path)
        return build_ndt_element_array_mesh_pair(
            free_path,
            void_path;
            config,
            defect=defects[:void_r2_focus],
            element_count=ELEMENT_COUNT,
        )
    end
    (
        free_mesh_path=abspath(free_path),
        void_mesh_path=abspath(void_path),
        defect_id=:void_r2_focus,
        element_count=ELEMENT_COUNT,
        common_triangulation=true,
    )
end

function solve_case(case_id; force=false)
    path = response_path(case_id)
    if isfile(path) && !force
        println("[=] Existing F0 response: $path")
        return path
    end
    mkpath(dirname(path))
    stand = stand_config()
    solver_config = NDTFullWaveHarmonicConfig(
        frequency_hz=stand.frequency_hz,
        density_kg_m3=stand.density_kg_m3,
        pressure_wave_speed_m_s=stand.pressure_wave_speed_m_s,
        shear_wave_speed_m_s=stand.shear_wave_speed_m_s,
        element_count=ELEMENT_COUNT,
    )
    println("[+] F0 solve $case_id: shared matrix + $ELEMENT_COUNT right-hand sides")
    result = solve_ndt_element_response_matrix(mesh_path(case_id); config=solver_config)
    jldsave(
        path;
        format_version=1,
        case_id,
        mesh_path=abspath(mesh_path(case_id)),
        frequency_hz=result.frequency_hz,
        normal_response_m2=result.normal_response_m2,
        tangent_response_m2=result.tangent_response_m2,
        element_lengths_m=result.element_lengths_m,
        matrix_size=result.matrix_size,
        matrix_nonzeros=result.matrix_nonzeros,
        stand_config=(
            specimen_length_mm=stand.specimen_length_mm,
            specimen_height_mm=stand.specimen_height_mm,
            nominal_focus_depth_mm=stand.nominal_focus_depth_mm,
            mesh_size_mm=stand.mesh_size_mm,
            defect_mesh_size_mm=stand.defect_mesh_size_mm,
        ),
        solver_config=(
            density_kg_m3=solver_config.density_kg_m3,
            pressure_wave_speed_m_s=solver_config.pressure_wave_speed_m_s,
            shear_wave_speed_m_s=solver_config.shear_wave_speed_m_s,
            rayleigh_alpha_s_inv=solver_config.rayleigh_alpha_s_inv,
            rayleigh_beta_s=solver_config.rayleigh_beta_s,
            element_count=solver_config.element_count,
            element_order=solver_config.element_order,
            quadrature_degree=solver_config.quadrature_degree,
        ),
    )
    println("[+] Saved F0 response: $path")
    path
end

positive_floor(value::Real) = max(Float64(value), floatmin(Float64))
relative_residual(first, second) = norm(first - second) / positive_floor(norm(first))

function mirror_matrix(matrix)
    matrix[end:-1:1, end:-1:1]
end

function project_operator(matrix, definition)
    if definition.kind == :a_scan
        port = ones(Float64, size(matrix, 1)) ./ sqrt(size(matrix, 1))
        return ComplexF64[only(transpose(port) * matrix * port)]
    elseif definition.kind == :static
        scale = inv(sqrt(definition.electronic_channels * definition.physical_states))
        return reduce(vcat, (
            vec(transpose(state) * matrix * state) .* scale
            for state in definition.states
        ))
    elseif definition.kind == :motion
        scale = inv(sqrt(length(definition.indices)))
        return ComplexF64[
            matrix[index, index] * scale for index in definition.indices
        ]
    elseif definition.kind == :fmc
        return vec(matrix) ./ sqrt(size(matrix, 1))
    end
    throw(ArgumentError("unknown architecture kind: $(definition.kind)"))
end

function architecture_features(free_matrix, defect_matrix, definition)
    free = project_operator(free_matrix, definition)
    defect = project_operator(defect_matrix, definition)
    (; free, defect, paired=defect - free)
end

function frozen_definitions()
    isfile(STATIC_FRONTIER_PATH) || error("static frontier not found: $STATIC_FRONTIER_PATH")
    states = JLD2.load(STATIC_FRONTIER_PATH, "architecture_states")
    [
        (
            id=:uniform_a_scan, kind=:a_scan, electronic_channels=1,
            physical_states=1, transmit_events=1, indices=Int[], states=(),
        ),
        (
            id=STATIC_K4Q2_ID, kind=:static, electronic_channels=4,
            physical_states=2, transmit_events=8, indices=Int[],
            states=states[STATIC_K4Q2_ID],
        ),
        (
            id=STATIC_K3Q4_ID, kind=:static, electronic_channels=3,
            physical_states=4, transmit_events=12, indices=Int[],
            states=states[STATIC_K3Q4_ID],
        ),
        (
            id=:motion_compact_center_p12_native, kind=:motion,
            electronic_channels=1, physical_states=1, transmit_events=12,
            indices=MOTION_P12_INDICES, states=(),
        ),
        (
            id=:motion_uniform_full_p17_a8, kind=:motion,
            electronic_channels=1, physical_states=1, transmit_events=17,
            indices=MOTION_P17_INDICES, states=(),
        ),
        (
            id=:paut_fmc_oracle, kind=:fmc, electronic_channels=17,
            physical_states=1, transmit_events=17, indices=Int[], states=(),
        ),
    ]
end

function finite_complex(array)
    all(isfinite, real.(array)) && all(isfinite, imag.(array))
end

function mesh_physical_names(path)
    names = Set{String}()
    for line in eachline(path)
        line == "\$EndPhysicalNames" && break
        matched = match(r"^\d+\s+\d+\s+\"(.+)\"$", line)
        isnothing(matched) || push!(names, matched.captures[1])
    end
    names
end

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

function analyze()
    free_data = JLD2.load(response_path(:none))
    defect_data = JLD2.load(response_path(:void_r2_focus))
    free = ComplexF64.(free_data["normal_response_m2"])
    defect = ComplexF64.(defect_data["normal_response_m2"])
    free_tangent = ComplexF64.(free_data["tangent_response_m2"])
    defect_tangent = ComplexF64.(defect_data["tangent_response_m2"])
    delta = defect - free
    delta_tangent = defect_tangent - free_tangent

    diagnostics = [
        (
            case_id=:none,
            reciprocity_residual=relative_residual(free, transpose(free)),
            mirror_residual=relative_residual(free, mirror_matrix(free)),
            tangent_over_normal=norm(free_tangent) / positive_floor(norm(free)),
        ),
        (
            case_id=:void_r2_focus,
            reciprocity_residual=relative_residual(defect, transpose(defect)),
            mirror_residual=relative_residual(defect, mirror_matrix(defect)),
            tangent_over_normal=norm(defect_tangent) / positive_floor(norm(defect)),
        ),
        (
            case_id=:paired_scattering,
            reciprocity_residual=relative_residual(delta, transpose(delta)),
            mirror_residual=relative_residual(delta, mirror_matrix(delta)),
            tangent_over_normal=norm(delta_tangent) / positive_floor(norm(delta)),
        ),
    ]

    definitions = frozen_definitions()
    feature_map = Dict(
        definition.id => architecture_features(free, defect, definition)
        for definition in definitions
    )
    a_scan = feature_map[:uniform_a_scan]
    fmc = feature_map[:paut_fmc_oracle]
    a_pair_norm = norm(a_scan.paired)
    a_contrast = a_pair_norm / positive_floor(norm(a_scan.free))
    fmc_pair_norm = norm(fmc.paired)
    architecture_rows = NamedTuple[]
    for definition in definitions
        features = feature_map[definition.id]
        free_norm = norm(features.free)
        defect_norm = norm(features.defect)
        paired_norm = norm(features.paired)
        contrast = paired_norm / positive_floor(free_norm)
        push!(architecture_rows, (
            architecture_id=definition.id,
            kind=definition.kind,
            electronic_channels=definition.electronic_channels,
            physical_states=definition.physical_states,
            transmit_events=definition.transmit_events,
            complex_observables=length(features.paired),
            free_response_norm=free_norm,
            defect_response_norm=defect_norm,
            paired_scattering_norm=paired_norm,
            paired_to_free_db=20log10(positive_floor(contrast)),
            paired_gain_vs_a_scan_db=20log10(positive_floor(paired_norm / positive_floor(a_pair_norm))),
            contrast_gain_vs_a_scan_db=20log10(positive_floor(contrast / positive_floor(a_contrast))),
            paired_norm_vs_fmc=paired_norm / positive_floor(fmc_pair_norm),
        ))
    end

    max_reciprocity = maximum(getproperty.(diagnostics, :reciprocity_residual))
    max_mirror = maximum(getproperty.(diagnostics, :mirror_residual))
    free_names = mesh_physical_names(mesh_path(:none))
    defect_names = mesh_physical_names(mesh_path(:void_r2_focus))
    required_element_names = Set("SourceElement$(index)" for index in 1:ELEMENT_COUNT)
    all_element_tags_present =
        issubset(required_element_names, free_names) &&
        issubset(required_element_names, defect_names)
    free_backwall_distinct = all(names ->
        "Backwall" in names && "RadiationBoundary" in names,
        (free_names, defect_names),
    )
    common_physics =
        free_data["frequency_hz"] == defect_data["frequency_hz"] &&
        free_data["solver_config"] == defect_data["solver_config"] &&
        isapprox(
            free_data["element_lengths_m"], defect_data["element_lengths_m"];
            rtol=1.0e-12, atol=1.0e-15,
        )
    gate_rows = [
        (gate=:all_element_tags_present, limit="true", value=string(all_element_tags_present), passed=all_element_tags_present),
        (gate=:free_backwall_distinct, limit="true", value=string(free_backwall_distinct), passed=free_backwall_distinct),
        (gate=:maximum_reciprocity_residual, limit="1e-6", value=string(max_reciprocity), passed=max_reciprocity <= 1.0e-6),
        (gate=:maximum_mirror_residual, limit="0.03", value=string(max_mirror), passed=max_mirror <= 0.03),
        (gate=:paired_normal_scattering_nonzero, limit=">0", value=string(norm(delta)), passed=finite_complex(delta) && norm(delta) > 0),
        (gate=:tangential_response_finite, limit="true", value=string(finite_complex(delta_tangent)), passed=finite_complex(delta_tangent)),
        (gate=:common_physics_for_architectures, limit="true", value=string(common_physics), passed=common_physics),
    ]
    gate_passed = all(getproperty.(gate_rows, :passed))

    mkpath(OUTPUT_ROOT)
    diagnostic_path = write_namedtuple_csv(joinpath(OUTPUT_ROOT, "matrix_diagnostics.csv"), diagnostics)
    architecture_path = write_namedtuple_csv(joinpath(OUTPUT_ROOT, "architecture_projection.csv"), architecture_rows)
    gate_path = write_namedtuple_csv(joinpath(OUTPUT_ROOT, "f0_gates.csv"), gate_rows)
    result_path = joinpath(OUTPUT_ROOT, "fullwave_carrier_summary.jld2")
    jldsave(
        result_path;
        format_version=1,
        contract="docs/operator_o3_fullwave_contract.md",
        gate_passed,
        diagnostics,
        architecture_rows,
        gate_rows,
        free_normal_response_m2=free,
        defect_normal_response_m2=defect,
        paired_normal_response_m2=delta,
        paired_tangent_response_m2=delta_tangent,
    )

    labels = replace.(string.(getproperty.(architecture_rows, :architecture_id)), '_' => ' ')
    contrast_values = getproperty.(architecture_rows, :paired_to_free_db)
    gain_values = getproperty.(architecture_rows, :paired_gain_vs_a_scan_db)
    contrast_plot = bar(
        labels, contrast_values;
        ylabel="20 log10 ||Delta y||/||y0||, dB",
        title="F0 paired carrier contrast",
        legend=false,
        xrotation=25,
    )
    gain_plot = bar(
        labels, gain_values;
        ylabel="paired norm gain vs A-scan, dB",
        title="Frozen operator capture",
        legend=false,
        xrotation=25,
    )
    plot_path = joinpath(OUTPUT_ROOT, "fullwave_carrier_architectures.png")
    savefig(plot(contrast_plot, gain_plot; layout=(2, 1), size=(1200, 1050), margin=6Plots.mm), plot_path)

    readme_path = joinpath(OUTPUT_ROOT, "README.md")
    open(readme_path, "w") do io
        println(io, "# O3 F0: paired full-wave carrier operator")
        println(io)
        println(io, "> F0 physics gate: **$(gate_passed ? "PASS" : "FAIL")**. This is not a transient O3/POD result.")
        println(io)
        println(io, "One 2D full-elastic FEM operator at 242 kHz was factored per paired mesh. All six architectures below are projections of the same 17 x 17 element response.")
        println(io)
        println(io, "## Matrix diagnostics")
        println(io)
        println(io, "| Case | Reciprocity | Mirror | Tangential / normal |")
        println(io, "|---|---:|---:|---:|")
        for row in diagnostics
            println(io, @sprintf("| `%s` | %.3e | %.3e | %.3e |", row.case_id, row.reciprocity_residual, row.mirror_residual, row.tangent_over_normal))
        end
        println(io)
        println(io, "## Frozen architecture projections")
        println(io)
        println(io, "| Architecture | E | Observables | Paired/free, dB | Paired gain vs A, dB | Paired norm / FMC |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for row in architecture_rows
            println(io, @sprintf("| `%s` | %d | %d | %.2f | %.2f | %.3f |", row.architecture_id, row.transmit_events, row.complex_observables, row.paired_to_free_db, row.paired_gain_vs_a_scan_db, row.paired_norm_vs_fmc))
        end
        println(io)
        println(io, "`Delta = void - free` is used here only to diagnose physical operator transfer. F1 must return to raw paired pulse-echo trials and may not give the detector a free reference.")
        println(io)
        println(io, "A single centred circular target at one carrier is close to an even rank-one task and is unusually favourable to the uniform A-scan. This table is therefore a transfer diagnostic, not an architecture ranking.")
        println(io)
        println(io, "Artifacts: `matrix_diagnostics.csv`, `architecture_projection.csv`, `f0_gates.csv`, `fullwave_carrier_summary.jld2`, and `fullwave_carrier_architectures.png`.")
    end
    println("[+] F0 gate: $(gate_passed ? "PASS" : "FAIL")")
    println("[+] $readme_path")
    (; gate_passed, diagnostics, architecture_rows, gate_rows, result_path, plot_path, readme_path, diagnostic_path, architecture_path, gate_path)
end

function run_campaign(; force_mesh=false, force_solve=false)
    mkpath(OUTPUT_ROOT)
    build_mesh_pair(; force=force_mesh)
    solve_case(:none; force=force_solve)
    solve_case(:void_r2_focus; force=force_solve)
    analyze()
end

function main(args=ARGS)
    if "--build-only" in args
        build_mesh_pair(; force="--force" in args)
    elseif "--analyze-only" in args
        analyze()
    else
        run_campaign(
            ; force_mesh="--force-mesh" in args,
            force_solve="--force-solve" in args,
        )
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

end
