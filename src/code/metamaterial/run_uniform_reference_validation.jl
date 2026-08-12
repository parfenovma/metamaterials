module UniformReferenceValidation

using Gmsh: gmsh
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV,
    "METAMATERIALS_PORT_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "port_mode_calibration"),
)
const FREQUENCY_HZ = parse(Float64, get(ENV, "METAMATERIALS_PORT_FREQUENCY_HZ", "242000"))
const PORT_ELEMENT_COUNT = parse(Int, get(ENV, "METAMATERIALS_PORT_ELEMENTS", "80"))

include(joinpath(@__DIR__, "uniform_reference_mesher.jl"))
using .UniformReferenceMesher
include(joinpath(@__DIR__, "step1b_convert_models.jl"))
include(joinpath(@__DIR__, "harmonic_solver.jl"))
using .HarmonicElasticity
include(joinpath(@__DIR__, "port_mode_solver.jl"))
using .ElasticPortModes
include(joinpath(@__DIR__, "modal_projection.jl"))
using .ElasticModalProjection
include(joinpath(@__DIR__, "modal_harmonic_solver.jl"))
using .ModalHarmonicElasticity

mesh_path() = joinpath(OUTPUT_ROOT, "mesh_reference_uniform.msh")
model_path() = joinpath(OUTPUT_ROOT, "model_reference_uniform.json")
summary_path() = joinpath(OUTPUT_ROOT, "uniform_modal_validation.csv")

function ensure_model()
    isfile(model_path()) && return model_path()
    mkpath(OUTPUT_ROOT)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try
        build_uniform_reference_mesh(mesh_path())
    finally
        gmsh.finalize()
    end
    convert_mesh(mesh_path(); output_dir=OUTPUT_ROOT)
end

function projected_pairs(state, modes, tag)
    right_modes = propagating_modes(modes; direction=:right)
    [
        begin
            left_mode = counterpropagating_mode(right_mode, modes)
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

function write_rows(path, left, right, source_power)
    open(path, "w") do io
        println(io, "k_per_m,parity,left_right_power,left_left_power,right_right_power,right_left_power")
        for (left_pair, right_pair) in zip(left, right)
            @printf(
                io,
                "%.12g,%s,%.12g,%.12g,%.12g,%.12g\n",
                real(left_pair.right_mode.wavenumber_per_m),
                left_pair.right_mode.parity,
                abs2(left_pair.right_amplitude),
                abs2(left_pair.left_amplitude),
                abs2(right_pair.right_amplitude),
                abs2(right_pair.left_amplitude),
            )
        end
        left_net = sum(abs2(pair.right_amplitude) - abs2(pair.left_amplitude) for pair in left)
        right_net = sum(abs2(pair.right_amplitude) - abs2(pair.left_amplitude) for pair in right)
        println(io, "# source_power_w_per_m,$source_power")
        println(io, "# left_modal_net_w_per_m,$left_net")
        println(io, "# right_modal_net_w_per_m,$right_net")
        println(io, "# modal_net_ratio,$(right_net / left_net)")
    end
end

function run()
    ensure_model()
    modes = solve_port_modes(PortModeConfig(
        frequency_hz=FREQUENCY_HZ,
        element_count=PORT_ELEMENT_COUNT,
    ))
    modal_config = ModalHarmonicElasticity.HarmonicElasticity.HarmonicConfig(
        rayleigh_alpha=0.0,
        rayleigh_beta=0.0,
        element_order=2,
        quadrature_degree=4,
    )
    state = solve_symmetry_reduced_modal_state(
        model_path(),
        FREQUENCY_HZ;
        config=modal_config,
        port_element_count=PORT_ELEMENT_COUNT,
    )
    modes = state.modes
    left = projected_pairs(state, modes, "Source")
    right = projected_pairs(state, modes, "Microphone")
    write_rows(summary_path(), left, right, state.metrics.source_work_w_per_m)

    input_index = findfirst(pair -> pair.right_mode.parity == :symmetric, left)
    input_index === nothing && error("symmetric input mode was not found")
    input_left = left[input_index]
    input_right = right[input_index]
    incident = abs2(input_left.right_amplitude)
    reflected = abs2(input_left.left_amplitude)
    transmitted = abs2(input_right.right_amplitude)
    right_incoming = abs2(input_right.left_amplitude)
    converted = sum(
        abs2(left_pair.left_amplitude) + abs2(right_pair.right_amplitude)
        for (index, (left_pair, right_pair)) in enumerate(zip(left, right))
        if index != input_index
    )
    left_net = sum(abs2(pair.right_amplitude) - abs2(pair.left_amplitude) for pair in left)
    right_net = sum(abs2(pair.right_amplitude) - abs2(pair.left_amplitude) for pair in right)
    @printf("[+] m0 incident / reflected: %.6f / %.6f W/m\n", incident, reflected)
    @printf("[+] m0 transmitted / right-incoming: %.6f / %.6f W/m\n", transmitted, right_incoming)
    @printf("[+] R00 / T00 / C0: %.8g / %.8g / %.8g\n", reflected / incident, transmitted / incident, converted / incident)
    @printf("[+] all-mode net flux ratio: %.8f\n", right_net / left_net)
    @printf("[+] right-boundary contamination / incident: %.6f\n", right_incoming / incident)
    println("[+] $(summary_path())")
    (; modes, left, right, state, left_net, right_net)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end

end
