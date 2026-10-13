using Test
include(joinpath(@__DIR__,"..","src","code","metamaterial","staggered_waveguide_mesher.jl"))
using .StaggeredWaveguideMesher

@testset "STG continuation preserves the complete core" begin
    b0,t0=extended_boundaries()
    for lead in (7.,14.)
        b,t=extended_boundaries(;lead_mm=lead)
        @test all(isapprox(x-lead, p[1]; atol=1e-12) && y == p[2]
                  for ((x,y),p) in zip(b[2:end-1],b0))
        @test all(isapprox(x-lead, p[1]; atol=1e-12) && y == p[2]
                  for ((x,y),p) in zip(t[2:end-1],t0))
        @test first(b)==(0.,0.)
        @test last(b)==(17+2lead,0.)
        @test first(t)==(17+2lead,7.)
        @test last(t)==(0.,7.)
    end
    @test_throws ArgumentError extended_boundaries(;lead_mm=-1.)
end

@testset "Changing N preserves period, gap and end transitions" begin
    profiles=StaggeredWaveguideMesher.MetamaterialProfiles
    for n in (1,2,4)
        g=fixed_period_geometry(n);p=profiles.SinusoidalProfile(2.5;periods=n)
        @test profiles.active_span_mm(p,g)/n ≈ 6.48
        @test profiles.half_period_shift_mm(p,g) ≈ 3.24
        @test g.samples==100n
        @test profiles.minimum_gap_mm(p,g) ≈ 4.5
        for x in (0.4,1.2,2.5,3.64)
            @test profiles.depth_at_x(p,g,x) ≈ profiles.depth_at_x(profiles.SinusoidalProfile(2.5;periods=2),profiles.GeometryConfig(),x)
        end
        bottom,top=profiles.generate_boundary_points(p,g);poly=vcat(bottom,top)
        area=abs(sum(poly[i][1]*poly[mod1(i+1,length(poly))][2]-poly[mod1(i+1,length(poly))][1]*poly[i][2] for i in eachindex(poly)))/2
        @test area ≈ 7g.length_mm-2.5n*6.48
    end
end
