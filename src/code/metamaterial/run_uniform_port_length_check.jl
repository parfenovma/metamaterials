module UniformPortLengthCheck

ENV["GKSwstype"] = "100"

using Gmsh: gmsh
using Plots
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_PORT_LENGTH_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "uniform_port_length_check"),
)
const FREQUENCY_HZ = parse(Float64, get(ENV, "METAMATERIALS_PORT_FREQUENCY_HZ", "242000"))
const LENGTHS_MM = [12.0, 24.0, 36.0]

include(joinpath(@__DIR__, "uniform_reference_mesher.jl"))
using .UniformReferenceMesher
include(joinpath(@__DIR__, "step1b_convert_models.jl"))
include(joinpath(@__DIR__, "modal_harmonic_solver.jl"))
using .ModalHarmonicElasticity
include(joinpath(@__DIR__, "modal_projection.jl"))
using .ElasticModalProjection

slug(length_mm) = replace(@sprintf("%.1f", length_mm), "." => "p")
mesh_path(length_mm) = joinpath(OUTPUT_ROOT, "mesh_reference_uniform_L$(slug(length_mm)).msh")
model_path(length_mm) = joinpath(OUTPUT_ROOT, "model_reference_uniform_L$(slug(length_mm)).json")

function ensure_models()
    mkpath(OUTPUT_ROOT)
    missing_lengths = filter(length_mm -> !isfile(model_path(length_mm)), LENGTHS_MM)
    isempty(missing_lengths) && return
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        for length_mm in missing_lengths
            build_uniform_reference_mesh(
                mesh_path(length_mm);
                config=UniformReferenceConfig(length_mm=length_mm),
                size_min_mm=0.12,
                size_max_mm=0.45,
            )
        end
    finally
        gmsh.finalize()
    end
    for length_mm in missing_lengths
        convert_mesh(mesh_path(length_mm); output_dir=OUTPUT_ROOT)
    end
end

function solve_length(length_mm)
    config = ModalHarmonicElasticity.HarmonicElasticity.HarmonicConfig(
        rayleigh_alpha=0.0,
        rayleigh_beta=0.0,
        element_order=2,
        quadrature_degree=4,
    )
    state = solve_symmetry_reduced_modal_state(
        model_path(length_mm),
        FREQUENCY_HZ;
        config,
    )
    left = project_boundary_mode_pair(
        state,
        "Source",
        state.right_mode,
        state.left_mode,
        quadrature_degree=6,
    )
    right = project_boundary_mode_pair(
        state,
        "Microphone",
        state.right_mode,
        state.left_mode,
        quadrature_degree=6,
    )
    incident = abs2(left.right_amplitude)
    reflection = abs2(left.left_amplitude) / incident
    transmission = abs2(right.right_amplitude) / incident
    right_incoming = abs2(right.left_amplitude) / incident
    r = left.left_amplitude / left.right_amplitude
    t = right.right_amplitude / left.right_amplitude
    k = real(state.right_mode.wavenumber_per_m)
    deembedded_t = t * exp(im * k * length_mm * 1e-3)
    (
        length_mm,
        k_per_m=k,
        incident_power=incident,
        R00=reflection,
        T00=transmission,
        right_incoming,
        r,
        t,
        deembedded_phase_rad=angle(deembedded_t),
        deembedded_amplitude=abs(deembedded_t),
    )
end

function write_summary(path, rows)
    open(path, "w") do io
        println(io, "length_mm,k_per_m,incident_power,R00,T00,right_incoming,deembedded_phase_rad,deembedded_amplitude")
        for row in rows
            @printf(
                io,
                "%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",
                row.length_mm,
                row.k_per_m,
                row.incident_power,
                row.R00,
                row.T00,
                row.right_incoming,
                row.deembedded_phase_rad,
                row.deembedded_amplitude,
            )
        end
    end
end

function plot_summary(path, rows)
    lengths = getproperty.(rows, :length_mm)
    power_panel = plot(
        lengths,
        getproperty.(rows, :R00);
        marker=:circle,
        linewidth=2,
        label="R00",
        xlabel="uniform length, mm",
        ylabel="power fraction",
        yscale=:log10,
        title="Symmetry-reduced modal boundary",
    )
    plot!(power_panel, lengths, getproperty.(rows, :right_incoming); marker=:circle, linewidth=2, label="right incoming")
    phase_panel = plot(
        lengths,
        getproperty.(rows, :deembedded_phase_rad);
        marker=:circle,
        linewidth=2,
        label="de-embedded arg(t)",
        xlabel="uniform length, mm",
        ylabel="phase, rad",
    )
    figure = plot(power_panel, phase_panel; layout=(2, 1), size=(850, 700))
    savefig(figure, path)
end

function run()
    ensure_models()
    rows = NamedTuple[]
    for length_mm in LENGTHS_MM
        println("[>] solving uniform length $(length_mm) mm")
        push!(rows, solve_length(length_mm))
    end
    csv_path = joinpath(OUTPUT_ROOT, "uniform_length_check.csv")
    figure_path = joinpath(OUTPUT_ROOT, "uniform_length_check.png")
    write_summary(csv_path, rows)
    plot_summary(figure_path, rows)
    for row in rows
        @printf(
            "[+] L=%4.1f mm: R00=%.3g, T00=%.7f, right-in=%.3g, deembedded phase=%+.3g rad\n",
            row.length_mm,
            row.R00,
            row.T00,
            row.right_incoming,
            row.deembedded_phase_rad,
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
