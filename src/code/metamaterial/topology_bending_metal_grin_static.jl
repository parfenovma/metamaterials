module TopologyBendingMetalGRINStatic

using Gridap

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUTPUT_ROOT = get(
    ENV, "METAMATERIALS_BENDING_METAL_GRIN_OUTPUT",
    joinpath(PROJECT_ROOT, "tmp", "topology_benchmark", "bending_metal_grin"),
)

include(joinpath(@__DIR__, "topology_bending_metal_grin_mesher.jl"))
using .TopologyBendingMetalGRINMesher

function solve_static_stiffness(model_path, length_m; imposed_displacement_m=1e-6)
    model = DiscreteModelFromFile(model_path)
    rho = 2700.0
    cp = 6122.102437409232
    cs = 3083.810277185563
    mu = rho * cs^2
    lambda = rho * cp^2 - 2mu
    identity_tensor = one(TensorValue{2, 2, ComplexF64})
    sigma(strain) = lambda * tr(strain) * identity_tensor + 2mu * strain
    reference = ReferenceFE(lagrangian, VectorValue{2, Float64}, 2)
    test = TestFESpace(
        model, reference; conformity=:H1, vector_type=Vector{ComplexF64},
        dirichlet_tags=["Source", "Microphone"],
    )
    zero_u(_) = VectorValue(0.0 + 0im, 0.0 + 0im)
    right_u(_) = VectorValue(ComplexF64(imposed_displacement_m), 0.0 + 0im)
    trial = TrialFESpace(test, [zero_u, right_u])
    domain = Triangulation(model)
    measure = Measure(domain, 4)
    a(u, v) = ∫(sigma(ε(u)) ⊙ ε(v))measure
    zero_vector = VectorValue(0.0 + 0im, 0.0 + 0im)
    b(v) = ∫(zero_vector ⋅ v)measure
    displacement = solve(AffineFEOperator(a, b, trial, test))

    coordinate = get_physical_coordinate(domain)
    lifting = (point -> VectorValue(ComplexF64(point[1] / length_m), 0.0 + 0im)) ∘ coordinate
    lifting_strain = TensorValue(
        ComplexF64(1 / length_m), 0.0 + 0im,
        0.0 + 0im, 0.0 + 0im,
    )
    reaction = sum(∫(sigma(ε(displacement)) ⊙ lifting_strain)measure)
    real(reaction / imposed_displacement_m)
end

function write_csv(path, rows)
    mkpath(dirname(path))
    columns = propertynames(first(rows))
    open(path, "w") do io
        println(io, join(string.(columns), ','))
        for row in rows
            println(io, join((string(getproperty(row, c)) for c in columns), ','))
        end
    end
    println("[+] $path")
end

function main()
    config = BendingMetalGRINConfig()
    stiffness = Dict{Symbol, Float64}()
    density = Dict{Symbol, Float64}()
    g0_lines = readlines(joinpath(OUTPUT_ROOT, "g0_manufacturing.csv"))
    header = split(first(g0_lines), ',')
    variant_col = findfirst(==("variant"), header)
    density_col = findfirst(==("relative_areal_density"), header)
    for line in Iterators.drop(g0_lines, 1)
        fields = split(line, ',')
        density[Symbol(fields[variant_col])] = parse(Float64, fields[density_col])
    end
    for variant in BENDING_METAL_GRIN_VARIANTS
        path = joinpath(OUTPUT_ROOT, "models", "model_$(variant).json")
        stiffness[variant] = solve_static_stiffness(path, total_length_mm(config) * 1e-3)
        println("[+] static $variant: K=", stiffness[variant])
    end
    solid_k = stiffness[:solid]
    rows = [(
        variant,
        relative_areal_density=density[variant],
        static_axial_stiffness_n_per_m2=stiffness[variant],
        static_stiffness_ratio=stiffness[variant] / solid_k,
        quasistatic_speed_ratio=sqrt((stiffness[variant] / solid_k) / density[variant]),
        quasistatic_impedance_ratio=sqrt((stiffness[variant] / solid_k) * density[variant]),
        ligament_mm=variant_parameters(config, variant).ligament_mm,
        mass_height_mm=variant_parameters(config, variant).mass_height_mm,
    ) for variant in BENDING_METAL_GRIN_VARIANTS]
    write_csv(joinpath(OUTPUT_ROOT, "g0_static_stiffness_density.csv"), rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
