module MaterialLensResponseMatrixTests

using Test

include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial", "sinusoidal_material_lens.jl",
))
using .SinusoidalMaterialLens
include(joinpath(
    @__DIR__, "..", "src", "code", "metamaterial",
    "material_lens_harmonic_solver.jl",
))
using .MaterialLensHarmonicSolver

@testset "material-lens response-matrix configuration" begin
    config = MaterialLensResponseConfig()
    @test MaterialLensHarmonicSolver.validate(config) === config
    @test config.source_count == 15
    @test length(config.output_y_mm) == 15
    @test config.output_y_mm ≈ reverse(-config.output_y_mm)
    @test_throws ArgumentError MaterialLensHarmonicSolver.validate(
        MaterialLensResponseConfig(source_count=0),
    )
    @test_throws ArgumentError MaterialLensHarmonicSolver.validate(
        MaterialLensResponseConfig(output_y_mm=Float64[]),
    )
end

end
