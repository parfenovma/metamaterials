module StaggeredWaveguidePortAudit

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUT = joinpath(ROOT,"results","staggered_waveguide_port_audit")
const STAGE = isempty(ARGS) ? "--pilot" : first(ARGS)
slug(x) = replace(string(Float64(x)),"."=>"p")
stem(role,lead,fine=false) = "$(role)_lead$(slug(lead))$(fine ? "_fine" : "")"
meshpath(role,lead,fine=false) = joinpath(OUT,"meshes","mesh_$(stem(role,lead,fine)).msh")
modelpath(role,lead,fine=false) = joinpath(OUT,"models","model_$(stem(role,lead,fine)).json")

if STAGE=="--mesh"
    using Gmsh: gmsh
    include(joinpath(@__DIR__,"staggered_waveguide_mesher.jl"))
elseif STAGE=="--convert"
    include(joinpath(@__DIR__,"step1b_convert_models.jl"))
else
    using JLD2, Gridap, LinearAlgebra
    include(joinpath(@__DIR__,"multimode_port_harmonic_solver.jl"))
    using .MultimodePortHarmonicSolver
    using .SinusoidalMaterialLens: ElasticMaterial
    using .ElasticPortModes
    include(joinpath(@__DIR__,"harmonic_solver.jl"))
    using .HarmonicElasticity
end

function meshes()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal",0)
    try
        for fine in (false,true), role in ("straight","stg"), lead in (0.,7.,14.)
            fine && lead != 14 && continue
            path=meshpath(role,lead,fine)
            isfile(path) && continue
            gmsh.clear();gmsh.model.add(stem(role,lead,fine))
            StaggeredWaveguideMesher.build_staggered_waveguide_mesh(path;
                lead_mm=lead,straight=role=="straight",
                size_min_mm=fine ? 0.10 : 0.15,size_max_mm=fine ? 0.40 : 0.60)
            println("[+] ",path)
        end
    finally
        gmsh.finalize()
    end
end

function convert()
    for path in sort(readdir(joinpath(OUT,"meshes");join=true))
        endswith(path,".msh") || continue
        dest=joinpath(OUT,"models",model_filename(path))
        isfile(dest) || convert_mesh(path;output_dir=joinpath(OUT,"models"))
    end
end

function write_rows(path,rows)
    isempty(rows) && return
    mkpath(dirname(path))
    open(path,"w") do io
        println(io,join(string.(propertynames(first(rows))),','))
        for row in rows;println(io,join(string.(values(row)),','));end
    end
end

const BRANCH = Dict{Float64,Any}()
function incident(f)
    get!(BRANCH,Float64(f)) do
        # Follow S0 from low frequency rather than switching at an extra-mode cut-on.
        frequencies=sort(unique(vcat(collect(20e3:5e3:f),[Float64(f)])))
        continuation=track_fundamental_branch(
            PortModeConfig(height_m=7e-3,element_count=120),frequencies)
        minimum(continuation.overlaps)>0.8 || error("ambiguous incident-mode continuation")
        last(continuation.modes)
    end
end

function modal_case(role,lead,f;loss=false,evanescent=8,fine=false)
    id="$(stem(role,lead,fine))_f$(round(Int,f))_$(loss ? "lossy" : "lossless")_ev$(evanescent)"
    path=joinpath(OUT,"cases",id*".jld2")
    isfile(path) && return JLD2.load(path,"row")
    mkpath(dirname(path))
    material=ElasticMaterial(name="project_photopolymer",density_kg_m3=1210.,
        pressure_wave_speed_m_s=2340.,shear_wave_speed_m_s=1170.,
        rayleigh_alpha_s_inv=loss ? 79560. : 0.,rayleigh_beta_s=loss ? 2.5e-9 : 0.)
    println("[>] ",id);flush(stdout)
    result=solve_multimode_port_harmonic(modelpath(role,lead,fine),["Microphone"],material;
        config=MultimodePortConfig(frequency_hz=f,element_order=2,quadrature_degree=6,
            port_element_count=120,evanescent_modes_per_port=evanescent,damping_domain_tag="Core"),
        incident_mode=incident(f))
    left=result.source_amplitudes;right=only(result.output_amplitudes)
    k0=real(result.incident_mode.wavenumber_per_m)
    index=argmin([abs(real(m.wavenumber_per_m)-k0) for m in left])
    a0=left[index].right_amplitude
    input=result.incident_power_w_per_m
    r=result.reflected_power_w_per_m/input
    t=only(result.transmitted_power_w_per_m)/input
    d=result.internal_dissipated_power_w_per_m/input
    incoming=only(result.incoming_output_power_w_per_m)/input
    t00=abs2(right[index].right_amplitude)/input
    r00=abs2(left[index].left_amplitude)/input
    transmission=right[index].right_amplitude/a0*exp(2im*k0*lead*1e-3)
    row=(geometry=role,lead_mm=lead,frequency_hz=Float64(f),lossy=loss,fine=fine,
        evanescent=evanescent,mode_count=length(left),incident_power=input,
        unwanted_incident_fraction=1-abs2(a0)/input,R_all=r,T_all=t,D=d,
        R00=r00,T00=t00,T_converted=t-t00,R_converted=r-r00,
        right_incoming_fraction=incoming,balance_residual=1+incoming-r-t-d,
        t00_core_real=real(transmission),t00_core_imag=imag(transmission),
        t00_core_phase=angle(transmission),k0_per_m=k0,
        max_port_gram_condition=maximum(result.port_gram_condition_numbers))
    mode_rows=[(mode=i,parity=string(m.parity),k_per_m=real(m.wavenumber_per_m),
        incoming_fraction=abs2(left[i].right_amplitude)/input,
        reflected_fraction=abs2(left[i].left_amplitude)/input,
        transmitted_fraction=abs2(m.right_amplitude)/input,
        t_core_real=real(m.right_amplitude/a0*exp(im*(k0+real(m.wavenumber_per_m))*lead*1e-3)),
        t_core_imag=imag(m.right_amplitude/a0*exp(im*(k0+real(m.wavenumber_per_m))*lead*1e-3)))
        for (i,m) in enumerate(right)]
    write_rows(joinpath(OUT,"cases",id*"_modes.csv"),mode_rows)
    ports=Dict(tag=>ElasticModalProjection.sample_boundary_field(result,tag)
               for tag in ("Source","Microphone"))
    jldsave(path;row,mode_rows,ports,source_model=modelpath(role,lead,fine),
        core_length_mm=17.,height_mm=7.,amplitude_mm=role=="stg" ? 2.5 : 0.,periods=2,
        input_mode="S0 continuation from 20 kHz",damping_region="Core only; extensions lossless")
    if lead==14 && !fine && evanescent==8 && f in (140e3,220e3)
        grid=Triangulation(result.model)
        data=only(Gridap.Visualization.visualization_data(grid,"field";order=1,
            cellfields=Dict("u"=>result.displacement)))
        coords=collect(Gridap.Geometry.get_node_coordinates(data.grid));u=collect(data.nodaldata["u"])
        jldsave(joinpath(OUT,"cases",id*"_field.jld2");
            x_m=[p[1] for p in coords],y_m=[p[2] for p in coords],
            ux_m=[p[1] for p in u],uy_m=[p[2] for p in u],
            triangles=[Int.(collect(ids)) for ids in Gridap.Geometry.get_cell_node_ids(data.grid)])
    end
    println("[+] R/T/D=",(r,t,d)," converted=",t-t00," balance=",row.balance_residual)
    flush(stdout)
    row
end

function local_case(role,lead)
    path=joinpath(OUT,"cases","local_$(stem(role,lead)).jld2")
    isfile(path) && return JLD2.load(path,"row")
    println("[>] local absorber ",role," lead=",lead);flush(stdout)
    p=solve_harmonic_point(modelpath(role,lead),140e3;
        config=HarmonicConfig(element_order=2,quadrature_degree=4))
    row=(geometry=role,lead_mm=lead,frequency_hz=140e3,
        traction_gain=abs(p.right_traction_pa)/1e6,mean_output_velocity=abs(2pi*140e3*p.right_displacement),
        source_power=p.source_power_w_per_m,right_absorbed_power=p.right_absorbed_power_w_per_m,
        left_absorbed_power=p.left_absorbed_power_w_per_m,dissipated_power=p.internal_dissipated_power_w_per_m)
    jldsave(path;row,damping_region="Entire domain, reproducing the old local-absorber model")
    row
end

function main()
    if STAGE=="--mesh";meshes()
    elseif STAGE=="--convert";convert()
    elseif STAGE=="--pilot"
        rows=[modal_case(role,lead,140e3) for role in ("straight","stg") for lead in (0.,7.,14.)]
        write_rows(joinpath(OUT,"lead_pilot.csv"),rows)
    elseif STAGE=="--refine"
        rows=[modal_case(role,14.,f;evanescent=ev,fine=fine)
              for f in (140e3,220e3) for role in ("straight","stg")
              for (ev,fine) in ((8,false),(16,false),(16,true))]
        append!(rows,[modal_case("stg",14.,f;loss=true,evanescent=16,fine=true)
                      for f in (140e3,220e3)])
        write_rows(joinpath(OUT,"refinement.csv"),rows)
    elseif STAGE=="--band"
        rows=[modal_case(role,14.,f;loss) for loss in (false,true)
              for f in (120e3,130e3,140e3,150e3,160e3,220e3) for role in ("straight","stg")]
        write_rows(joinpath(OUT,"frequency_scan.csv"),rows)
    elseif STAGE=="--local"
        rows=[local_case(role,lead) for role in ("straight","stg") for lead in (0.,7.,14.)]
        write_rows(joinpath(OUT,"local_absorber_leads.csv"),rows)
    else
        error("expected --mesh, --convert, --pilot, --refine, --band or --local")
    end
end

if abspath(PROGRAM_FILE)==@__FILE__;main();end
end
