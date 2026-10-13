module STGCompletion

using JLD2, LinearAlgebra, SHA
const ROOT=normpath(joinpath(@__DIR__,"..","..",".."))
const OUT=joinpath(ROOT,"results","stg_modal_completion")
const STAGE=isempty(ARGS) ? "--modes" : ARGS[1]
const FREQ=collect(20e3:2.5e3:360e3)
const LEAD=14.0

if STAGE=="--mesh"
    using Gmsh:gmsh
    include(joinpath(@__DIR__,"staggered_waveguide_mesher.jl"))
    using .StaggeredWaveguideMesher
elseif STAGE=="--convert"
    include(joinpath(@__DIR__,"step1b_convert_models.jl"))
else
    using Gridap
    include(joinpath(@__DIR__,"multimode_port_harmonic_solver.jl"))
    using .MultimodePortHarmonicSolver, .ElasticPortModes
    using .SinusoidalMaterialLens:ElasticMaterial
end

core_length(n)=0.8+(n+0.5)*6.48
stem(n,role;fine=false,lead=LEAD)="$(role)_n$(n)_lead$(round(Int,lead))$(fine ? "_fine" : "")"
modelpath(n,role;kwargs...)=joinpath(OUT,"models","model_$(stem(n,role;kwargs...)).json")
function material(loss)
    ElasticMaterial(name="project_photopolymer",density_kg_m3=1210.,pressure_wave_speed_m_s=2340.,
        shear_wave_speed_m_s=1170.,rayleigh_alpha_s_inv=79560.0*loss,rayleigh_beta_s=2.5e-9*loss)
end
function write_rows(path,rows)
    mkpath(dirname(path)); isempty(rows) && return
    open(path,"w") do io
        println(io,join(string.(propertynames(first(rows))),','))
        for r in rows;println(io,join(string.(values(r)),','));end
    end
end

function meshes()
    gmsh.initialize();gmsh.option.setNumber("General.Terminal",0)
    try
        specs=vcat([(n,role,false,LEAD) for n in (1,2,4) for role in ("stg","straight")],
            [(2,role,true,LEAD) for role in ("stg","straight")],[(4,"stg",true,LEAD)],[(2,"stg",false,28.)])
        for (n,role,fine,lead) in specs
            path=joinpath(OUT,"meshes","mesh_$(stem(n,role;fine,lead)).msh")
            isfile(path) && continue
            gmsh.clear();gmsh.model.add(stem(n,role;fine,lead))
            build_staggered_waveguide_mesh(path;lead_mm=lead,straight=role=="straight",periods=n,
                geometry=fixed_period_geometry(n),size_min_mm=fine ? 0.075 : 0.15,
                size_max_mm=fine ? 0.30 : 0.60)
            println("[+] mesh ",path);flush(stdout)
        end
    finally;gmsh.finalize();end
end

function rotate_mode(m,z)
    PortMode(m.gamma_per_m,m.wavenumber_per_m,m.kind,m.direction,m.power_w_per_m,
        m.relative_residual,m.parity_score,m.parity,m.p_fraction,m.s_fraction,
        m.axial_displacement_fraction,m.y_m,m.displacement_x*z,m.displacement_y*z,
        m.traction_xx_pa*z,m.traction_xy_pa*z)
end
vector(m)=vcat(m.displacement_x,m.displacement_y)
function modes_file(elements=120)
    joinpath(OUT,"port_modes_$(elements).jld2")
end

"""Track identities and fix the otherwise arbitrary eigenvector phase before IFFT."""
function build_modes(frequencies=FREQ;elements=120,save=true)
    banks=Any[];previous=Dict{String,Any}();nextlabel=Dict(:symmetric=>0,:antisymmetric=>0)
    diagnostics=NamedTuple[]
    for f in frequencies
        allm=solve_port_modes(PortModeConfig(height_m=7e-3,element_count=elements,frequency_hz=f))
        right=filter(m->m.kind==:propagating && m.direction==:right,allm)
        labels=fill("",length(right)); used=Set{String}(); overlaps=ones(length(right))
        pairs=sort([(mode_overlap(old,m),label,i) for (label,old) in previous
                    for (i,m) in enumerate(right) if old.parity==m.parity];by=first,rev=true)
        for (overlap,label,i) in pairs
            (overlap<0.55 || label in used || !isempty(labels[i])) && continue
            labels[i]=label;push!(used,label);overlaps[i]=overlap
        end
        canonical=PortMode[]
        for (i,m) in enumerate(right)
            if isempty(labels[i])
                prefix=m.parity==:symmetric ? "S" : "A"
                labels[i]=prefix*string(nextlabel[m.parity]);nextlabel[m.parity]+=1
                # At a branch birth, choose a nonzero longitudinal component.
                v=m.displacement_x; pivot=v[argmax(abs.(v))]
                phase=conj(pivot)/abs(pivot)
            else
                overlap=dot(vector(previous[labels[i]]),vector(m))
                phase=conj(overlap)/abs(overlap)
            end
            push!(canonical,rotate_mode(m,phase))
        end
        # A right/left pair shares the real-time conjugation convention.
        left=PortMode[]
        for m in canonical
            raw=ElasticModalProjection.counterpropagating_mode(m,allm)
            overlap=dot(conj.(vector(m)),vector(raw))
            push!(left,rotate_mode(raw,conj(overlap)/abs(overlap)))
        end
        kept=vcat(canonical,left)
        for direction in (:right,:left)
            ev=sort(filter(m->m.kind==:evanescent && m.direction==direction,allm);by=m->abs(real(m.gamma_per_m)))
            append!(kept,first(ev,min(24,length(ev))))
        end
        push!(banks,(frequency=f,modes=kept,right=canonical,labels,overlaps))
        for (m,label,overlap) in zip(canonical,labels,overlaps)
            push!(diagnostics,(frequency_hz=f,branch=label,k_per_m=real(m.wavenumber_per_m),
                parity=string(m.parity),overlap,axial_fraction=m.axial_displacement_fraction))
        end
        previous=Dict(label=>m for (label,m) in zip(labels,canonical))
    end
    if save
        mkpath(OUT);jldsave(modes_file(elements);banks,frequencies,phase_convention="parallel transport; conjugate counterpropagating partner")
        write_rows(joinpath(OUT,"mode_branches_$(elements).csv"),diagnostics)
    end
    banks
end

const PREPARED=Dict{Any,Any}()
function solve_case(n,role,loss,bank;fine=false,lead=LEAD,ev=12,save_field=false)
    f=bank.frequency
    id="$(stem(n,role;fine,lead))_f$(round(Int,f))_loss$(loss)_ev$(ev)_p$(length(first(bank.right).y_m)-1)"
    path=joinpath(OUT,"cases",id*".jld2")
    isfile(path) && return JLD2.load(path,"row")
    mat=material(loss)
    cfg=MultimodePortConfig(frequency_hz=f,element_order=2,quadrature_degree=6,
        port_element_count=length(first(bank.right).y_m)-1,evanescent_modes_per_port=ev,damping_domain_tag="Core")
    mesh=modelpath(n,role;fine,lead)
    prep=get!(PREPARED,(n,role,fine,lead)) do
        prepare_multimode_model(mesh,mat;config=cfg)
    end
    incident_index=findfirst(==("S0"),bank.labels); isnothing(incident_index) && error("S0 lost")
    result=solve_multimode_port_harmonic(mesh,["Microphone"],mat;config=cfg,
        prepared=prep,supplied_modes=bank.modes,incident_mode=bank.right[incident_index])
    left=result.source_amplitudes;right=only(result.output_amplitudes)
    a=left[incident_index].right_amplitude;k0=real(bank.right[incident_index].wavenumber_per_m)
    k=real.(getproperty.(bank.right,:wavenumber_per_m))
    # Reference planes are at the original core faces, excluding added lead phases.
    transmission=[m.right_amplitude/a for m in right].*exp.(im.*(k0.+k).*lead*1e-3)
    reflection=[m.left_amplitude/a for m in left].*exp.(im.*(k0.+k).*lead*1e-3)
    inc=result.incident_power_w_per_m
    T=only(result.transmitted_power_w_per_m)/inc;R=result.reflected_power_w_per_m/inc
    D=result.internal_dissipated_power_w_per_m/inc
    row=(n=n,geometry=role,frequency_hz=f,loss=Float64(loss),fine=fine,lead_mm=lead,
        T=T,R=R,D=D,T00=abs2(transmission[incident_index]),Tother=T-abs2(transmission[incident_index]),
        balance=1+only(result.incoming_output_power_w_per_m)/inc-T-R-D,
        incident_power=inc,unwanted_incident=1-abs2(a)/inc,
        incoming_right=only(result.incoming_output_power_w_per_m)/inc)
    mkpath(dirname(path))
    # Physical, branchwise velocity response also avoids any eigenvector-gauge
    # dependence of a scalar displacement receiver reconstruction.
    mean_vx=[im*2pi*f*sum((m.displacement_x[1:end-1]+m.displacement_x[2:end]).*diff(m.y_m))/2/7e-3 for m in bank.right]
    jldsave(path;row,labels=bank.labels,transmission,reflection,k,mean_vx,
        model_sha256=bytes2hex(sha256(read(mesh))),core_length_mm=core_length(n))
    if save_field
        grid=Triangulation(result.model)
        data=only(Gridap.Visualization.visualization_data(grid,"field";order=1,cellfields=Dict("u"=>result.displacement)))
        coords=collect(Gridap.Geometry.get_node_coordinates(data.grid));u=collect(data.nodaldata["u"])
        jldsave(joinpath(OUT,"cases",id*"_field.jld2");x_m=[p[1] for p in coords],y_m=[p[2] for p in coords],
            ux_m=[p[1] for p in u],uy_m=[p[2] for p in u],
            triangles=[Int.(collect(ids)) for ids in Gridap.Geometry.get_cell_node_ids(data.grid)])
    end
    println("[+] ",id," T=",round(T;digits=5)," R=",round(R;digits=5)," D=",round(D;digits=5)," err=",round(row.balance;digits=6));flush(stdout)
    row
end

function sweep(group)
    banks=JLD2.load(modes_file(),"banks");rows=NamedTuple[]
    for b in banks
        if group==0
            for loss in (0.,1.)
                push!(rows,solve_case(2,"stg",loss,b;save_field=b.frequency in (120e3,130e3,140e3,150e3,160e3,220e3)))
            end
        else
            push!(rows,solve_case(2,"straight",1.,b))
            if 80e3<=b.frequency<=280e3 && iszero(mod(b.frequency,5e3))
                for n in (1,4);push!(rows,solve_case(n,"stg",0.,b));end
            end
            if b.frequency in (130e3,140e3,160e3,220e3)
                for n in (1,4);push!(rows,solve_case(n,"stg",1.,b));end
                for loss in (0.5,2.);push!(rows,solve_case(2,"stg",loss,b));end
            end
        end
    end
    write_rows(joinpath(OUT,"sweep_group$(group).csv"),rows)
end

function main()
    if STAGE=="--mesh";meshes()
    elseif STAGE=="--convert"
        for p in sort(readdir(joinpath(OUT,"meshes");join=true))
            endswith(p,".msh") || continue
            isfile(joinpath(OUT,"models",model_filename(p))) || convert_mesh(p;output_dir=joinpath(OUT,"models"))
        end
    elseif STAGE=="--modes";build_modes()
    elseif STAGE=="--sweep";sweep(parse(Int,ARGS[2]))
    elseif STAGE=="--pilot"
        b=only(filter(b->b.frequency==140e3,JLD2.load(modes_file(),"banks")))
        solve_case(2,"stg",0.,b);solve_case(2,"straight",0.,b)
    else;error("unknown stage $STAGE");end
end
if abspath(PROGRAM_FILE)==@__FILE__;main();end
end
