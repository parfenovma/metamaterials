using Test, LinearAlgebra
include(joinpath(@__DIR__,"..","src","code","metamaterial","analyze_stg_completion.jl"))
using .AnalyzeSTGCompletion

@testset "Power-normalized pulse reconstruction" begin
    f=collect(20e3:2.5e3:360e3)
    H=hcat(0.6exp.(-2pi*im.*f.*8e-6),0.3exp.(-2pi*im.*f.*13e-6))
    R=hcat(sqrt(.1)exp.(-2pi*im.*f.*5e-6),zeros(ComplexF64,length(f)))
    s=(;f,H,R,V=copy(H),D=fill(.45,length(f)),labels=["S0","A0"])
    p=AnalyzeSTGCompletion.packet(s,140e3)
    @test p.coverage>.9999
    @test p.energy ≈ [.36,.09] atol=1e-12
    @test sum(p.reflected_energy)≈.1 atol=1e-12
    @test p.dissipated≈.45 atol=1e-12
    @test sum(p.energy)+sum(p.reflected_energy)+p.dissipated≈1. atol=1e-12
    # A positive propagation delay must move the signal to later time under
    # exp(+i omega t), not create an apparent time advance.
    @test p.output[:,1]≈.6circshift(p.incoming,80) atol=1e-12
    @test p.output[:,2]≈.3circshift(p.incoming,130) atol=1e-12
end
