using Test

@testset "STG metamaterial" begin
    include("staggered_waveguide_geometry_tests.jl")
    include("step3_integration.jl")
    include("stg_packet_tests.jl")
    include("multimode_prepared_tests.jl")
end
