using Test
include(joinpath(@__DIR__,"..","src","code","metamaterial","multimode_port_harmonic_solver.jl"))
using .MultimodePortHarmonicSolver, .SinusoidalMaterialLens

@testset "Prepared sweep agrees with direct volume assembly" begin
    path=joinpath(@__DIR__,"..","results","stg_modal_completion","models","model_stg_n2_lead14.json")
    @test isfile(path)
    cfg=MultimodePortConfig(frequency_hz=140e3,element_order=2,quadrature_degree=6,
        port_element_count=120,evanescent_modes_per_port=12,damping_domain_tag="Core")
    mat=ElasticMaterial(name="control",density_kg_m3=1210.,pressure_wave_speed_m_s=2340.,
        shear_wave_speed_m_s=1170.,rayleigh_alpha_s_inv=79560.,rayleigh_beta_s=2.5e-9)
    direct=solve_multimode_port_harmonic(path,["Microphone"],mat;config=cfg)
    prepared=prepare_multimode_model(path,mat;config=cfg)
    reused=solve_multimode_port_harmonic(path,["Microphone"],mat;config=cfg,prepared)
    for field in (:incident_power_w_per_m,:reflected_power_w_per_m,:internal_dissipated_power_w_per_m)
        @test getproperty(reused,field)≈getproperty(direct,field) rtol=1e-8
    end
    @test reused.transmitted_power_w_per_m≈direct.transmitted_power_w_per_m rtol=1e-8
end
