module ElasticMaterials

export ElasticMaterial, photopolymer, aluminium_6061

"""Parameters of a homogeneous isotropic elastic material."""
Base.@kwdef struct ElasticMaterial
    name::String
    density_kg_m3::Float64
    pressure_wave_speed_m_s::Float64
    shear_wave_speed_m_s::Float64
    rayleigh_alpha_s_inv::Float64 = 0.0
    rayleigh_beta_s::Float64 = 0.0
end

"""
Nominal photopolymer used by the existing STG calculations.

These numbers describe the current numerical model; they are not a claim
about every printable resin.
"""
photopolymer() = ElasticMaterial(
    name="project_photopolymer",
    density_kg_m3=1210.0,
    pressure_wave_speed_m_s=2340.0,
    shear_wave_speed_m_s=1170.0,
    rayleigh_alpha_s_inv=79560.0,
    rayleigh_beta_s=2.5e-9,
)

"""Nominal aluminium 6061 used as the inspected specimen."""
aluminium_6061() = ElasticMaterial(
    name="aluminium_6061",
    density_kg_m3=2700.0,
    pressure_wave_speed_m_s=6122.102437409232,
    shear_wave_speed_m_s=3083.810277185563,
)

end
