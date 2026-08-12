module BaselineR0PortCheck

ENV["GKSwstype"] = "100"

using Gmsh: gmsh
using Plots
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_R0_PORT_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "baseline_r0_port_check"),
)
const FREQUENCY_HZ = parse(Float64, get(ENV, "METAMATERIALS_PORT_FREQUENCY_HZ", "242000"))
const LEAD_LENGTHS_MM = [6.0, 12.0, 18.0, 24.0]
const NARROW_LENGTH_MM = 18.0

include(joinpath(@__DIR__, "baseline_r0_mesher.jl"))
using .BaselineR0Mesher
include(joinpath(@__DIR__, "step1b_convert_models.jl"))
include(joinpath(@__DIR__, "modal_harmonic_solver.jl"))
using .ModalHarmonicElasticity
include(joinpath(@__DIR__, "modal_projection.jl"))
using .ElasticModalProjection

slug(value) = replace(@sprintf("%.1f", value), "." => "p")
mesh_path(lead_mm) = joinpath(OUTPUT_ROOT, "mesh_baseline_R0_lead$(slug(lead_mm)).msh")
model_path(lead_mm) = joinpath(OUTPUT_ROOT, "model_baseline_R0_lead$(slug(lead_mm)).json")

function ensure_models()
    mkpath(OUTPUT_ROOT)
    missing = filter(lead_mm -> !isfile(model_path(lead_mm)), LEAD_LENGTHS_MM)
    isempty(missing) && return
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for lead_mm in missing
            build_baseline_r0_mesh(
                mesh_path(lead_mm);
                config=BaselineR0Config(
                    narrow_length_mm=NARROW_LENGTH_MM,
                    lead_length_mm=lead_mm,
                ),
            )
        end
    finally
        gmsh.finalize()
    end
    for lead_mm in missing
        convert_mesh(mesh_path(lead_mm); output_dir=OUTPUT_ROOT)
    end
end

function project_all(state, tag)
    right_modes = filter(
        mode -> mode.kind == :propagating && mode.direction == :right,
        state.modes,
    )
    [
        begin
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
        end
        for right_mode in right_modes
    ]
end

function solve_lead(lead_mm)
    config = ModalHarmonicElasticity.HarmonicElasticity.HarmonicConfig(
        rayleigh_alpha=0.0,
        rayleigh_beta=0.0,
        element_order=2,
        quadrature_degree=4,
    )
    state = solve_symmetry_reduced_modal_state(
        model_path(lead_mm),
        FREQUENCY_HZ;
        config,
    )
    left = project_all(state, "Source")
    right = project_all(state, "Microphone")
    input_index = findfirst(pair -> pair.right_mode.parity == :symmetric, left)
    input_index === nothing && error("symmetric m0 was not found")
    input_left = left[input_index]
    input_right = right[input_index]
    incident = abs2(input_left.right_amplitude)
    R00 = abs2(input_left.left_amplitude) / incident
    T00 = abs2(input_right.right_amplitude) / incident
    right_incoming = abs2(input_right.left_amplitude) / incident
    C0 = sum(
        abs2(left_pair.left_amplitude) + abs2(right_pair.right_amplitude)
        for (index, (left_pair, right_pair)) in enumerate(zip(left, right))
        if index != input_index
    ) / incident
    total_length_m = (NARROW_LENGTH_MM + 2lead_mm) * 1e-3
    t = input_right.right_amplitude / input_left.right_amplitude
    deembedded_t = t * exp(im * real(state.right_mode.wavenumber_per_m) * total_length_m)
    (
        lead_mm,
        total_length_mm=NARROW_LENGTH_MM + 2lead_mm,
        incident,
        R00,
        T00,
        C0,
        right_incoming,
        phase_rad=angle(deembedded_t),
        amplitude=abs(t),
    )
end

function write_summary(path, rows)
    open(path, "w") do io
        println(io, "lead_mm,total_length_mm,incident_power,R00,T00,C0,right_incoming,deembedded_phase_rad,amplitude")
        for row in rows
            @printf(
                io,
                "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                row.lead_mm,
                row.total_length_mm,
                row.incident,
                row.R00,
                row.T00,
                row.C0,
                row.right_incoming,
                row.phase_rad,
                row.amplitude,
            )
        end
    end
end

function plot_summary(path, rows)
    lead = getproperty.(rows, :lead_mm)
    top = plot(
        lead,
        getproperty.(rows, :T00);
        marker=:circle,
        linewidth=2,
        label="T00",
        xlabel="full-height lead, mm",
        ylabel="power",
        title="Absolute baseline_R0 modal scattering",
    )
    plot!(top, lead, getproperty.(rows, :R00); marker=:circle, linewidth=2, label="R00")
    bottom = plot(
        lead,
        getproperty.(rows, :phase_rad);
        marker=:circle,
        linewidth=2,
        label="de-embedded arg(t)",
        xlabel="full-height lead, mm",
        ylabel="phase, rad",
    )
    savefig(plot(top, bottom; layout=(2, 1), size=(850, 700)), path)
end

function run()
    ensure_models()
    rows = NamedTuple[]
    for lead_mm in LEAD_LENGTHS_MM
        println("[>] solving baseline_R0 with $(lead_mm) mm leads")
        push!(rows, solve_lead(lead_mm))
    end
    csv_path = joinpath(OUTPUT_ROOT, "baseline_R0_lead_check.csv")
    figure_path = joinpath(OUTPUT_ROOT, "baseline_R0_lead_check.png")
    write_summary(csv_path, rows)
    plot_summary(figure_path, rows)
    for row in rows
        @printf(
            "[+] lead=%4.1f mm: T00=%.5f, R00=%.5f, C0=%.3g, phase=%+.4f rad\n",
            row.lead_mm,
            row.T00,
            row.R00,
            row.C0,
            row.phase_rad,
        )
    end
    println("[+] $csv_path")
    println("[+] $figure_path")
    rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
