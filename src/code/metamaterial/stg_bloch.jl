module STGBloch
using LinearAlgebra, SparseArrays, ArnoldiMethod, JLD2
include(joinpath(@__DIR__,"port_mode_solver.jl"))
using .ElasticPortModes
const OUT=normpath(joinpath(@__DIR__,"..","..","..","results","stg_modal_completion"))

"""Independent Q1 plane-strain FE model of one continuously undulating cell.

Both walls are vertical translates: y_bottom=A(1-cos(2πx/a))/2,
y_top=y_bottom+g. No straight connector is part of the Bloch cell.
"""
function assemble_cell(;nx=64,ny=40,amplitude=2.5e-3,period=6.48e-3,gap=4.5e-3)
    rho=1210.;mu=rho*1170^2;lambda=rho*2340^2-2mu
    C=[lambda+2mu lambda 0.;lambda lambda+2mu 0.;0. 0. mu]
    node(i,j)=1+i+(nx+1)*j
    coords=zeros((nx+1)*(ny+1),2)
    for j in 0:ny,i in 0:nx
        x=period*i/nx;coords[node(i,j),:]=[x,amplitude*(1-cos(2pi*x/period))/2+gap*j/ny]
    end
    I=Int[];J=Int[];vk=Float64[];vm=Float64[]
    for j in 0:ny-1,i in 0:nx-1
        nodes=[node(i,j),node(i+1,j),node(i+1,j+1),node(i,j+1)]
        xy=coords[nodes,:];dofs=vec(permutedims(hcat(2nodes.-1,2nodes)))
        ke=zeros(8,8);me=zeros(8,8)
        for s in (-1/sqrt(3),1/sqrt(3)),t in (-1/sqrt(3),1/sqrt(3))
            N=[(1-s)*(1-t),(1+s)*(1-t),(1+s)*(1+t),(1-s)*(1+t)]/4
            dN=[-(1-t) -(1-s);(1-t) -(1+s);(1+t) (1+s);-(1+t) (1-s)]/4
            jac=xy'*dN;grad=dN/jac;weight=det(jac)
            weight>0 || error("inverted cell")
            B=zeros(3,8);Nv=zeros(2,8)
            for a in 1:4
                B[1,2a-1]=grad[a,1];B[2,2a]=grad[a,2]
                B[3,2a-1]=grad[a,2];B[3,2a]=grad[a,1]
                Nv[1,2a-1]=N[a];Nv[2,2a]=N[a]
            end
            ke.+=B'*C*B*weight;me.+=rho*(Nv'*Nv)*weight
        end
        for b in 1:8,a in 1:8
            push!(I,dofs[a]);push!(J,dofs[b]);push!(vk,ke[a,b]);push!(vm,me[a,b])
        end
    end
    ndof=2size(coords,1)
    K=sparse(I,J,vk,ndof,ndof);M=sparse(I,J,vm,ndof,ndof)
    (;K,M,nx,ny,period,gap,amplitude)
end

struct ShiftInverse{F,M}
    factor::F
    mass::M
    work::Vector{ComplexF64}
end
Base.size(a::ShiftInverse)=size(a.mass)
Base.size(a::ShiftInverse,i::Int)=size(a.mass,i)
Base.eltype(::ShiftInverse)=ComplexF64
function LinearAlgebra.mul!(y,a::ShiftInverse,x)
    mul!(a.work,a.mass,x);ldiv!(y,a.factor,a.work);y
end

function bands(cell,q;nev=28)
    nx=cell.nx;ny=cell.ny;rows=Int[];cols=Int[];values=ComplexF64[]
    for j in 0:ny,i in 0:nx,d in 1:2
        push!(rows,2*(i+(nx+1)*j)+d);push!(cols,2*(mod(i,nx)+nx*j)+d)
        push!(values,i==nx ? exp(im*q) : 1.)
    end
    P=sparse(rows,cols,values,2*(nx+1)*(ny+1),2*nx*(ny+1))
    K=P'*cell.K*P;M=P'*cell.M*P
    shift=(2pi*1e3)^2
    op=ShiftInverse(lu(K+shift*M),M,zeros(ComplexF64,size(K,1)))
    decomp,history=partialschur(op;nev,which=:LM,tol=1e-9,restarts=100)
    history.converged || error("Bloch eigensolver did not converge at q=$q")
    vals,vecs=partialeigen(decomp)
    omega2=real.(1 ./vals).-shift
    frequencies=sqrt.(max.(omega2,0.))/(2pi)
    order=sortperm(frequencies);frequencies=frequencies[order];vecs=vecs[:,order];omega2=omega2[order]
    residuals=[norm(K*v-w*M*v)/(norm(K*v)+abs(w)*norm(M*v)+eps()) for (v,w) in zip(eachcol(vecs),omega2)]
    maximum(residuals[frequencies.>1e3])<1e-6 || error("large eigenpair residual")
    (;frequencies,residuals)
end

"""Complex Bloch multipliers from exact dynamic condensation of this FE cell."""
function complex_bloch(cell,frequency)
    nx=cell.nx;ny=cell.ny;omega=2pi*frequency
    left=[2*(nx+1)*j+d for j in 0:ny for d in 1:2]
    right=[2*(nx+(nx+1)*j)+d for j in 0:ny for d in 1:2]
    boundary=vcat(left,right);interior=setdiff(collect(axes(cell.K,1)),boundary)
    D=cell.K-omega^2*cell.M
    S=Matrix(D[boundary,boundary])-D[boundary,interior]*(lu(D[interior,interior])\Matrix(D[interior,boundary]))
    n=length(left);LL=S[1:n,1:n];LR=S[1:n,n+1:end];RL=S[n+1:end,1:n];RR=S[n+1:end,n+1:end]
    A0=RL;A1=LL+RR;A2=LR
    scale=maximum(norm.((A0,A1,A2)));A0=A0/scale;A1=A1/scale;A2=A2/scale
    eye=Matrix{ComplexF64}(I,n,n);z=zeros(ComplexF64,n,n)
    eig=eigen([z eye;-A0 -A1],[eye z;z A2])
    rows=NamedTuple[]
    for (j,multiplier) in enumerate(eig.values)
        isfinite(multiplier) && 1e-7<abs(multiplier)<1e7 || continue
        u=eig.vectors[1:n,j]
        residual=norm(A0*u+multiplier*A1*u+multiplier^2*A2*u)/
            (norm(A0*u)+abs(multiplier)*norm(A1*u)+abs2(multiplier)*norm(A2*u))
        residual<1e-7 || continue
        force=LL*u+LR*(multiplier*u)
        power=.5real(dot(im*omega*u,force))
        propagating=abs(abs(multiplier)-1)<1e-5
        forward=propagating ? power>0 : abs(multiplier)<1
        forward || continue
        push!(rows,(frequency_hz=frequency,lambda_real=real(multiplier),lambda_imag=imag(multiplier),
            magnitude=abs(multiplier),propagating,attenuation_per_m=propagating ? 0. : -log(abs(multiplier))/cell.period,
            phase_per_cell=angle(multiplier),residual))
    end
    sort(rows;by=r->r.attenuation_per_m)
end

function run_complex()
    allrows=NamedTuple[];checks=NamedTuple[]
    for (nx,ny) in ((64,40),(96,60))
        cell=assemble_cell(;nx,ny)
        for f in (80e3,130e3,140e3,160e3,200e3,220e3)
            r=complex_bloch(cell,f);isempty(r) && error("no forward Bloch multipliers")
            append!(allrows,[merge((nx=nx,ny=ny,),x) for x in r])
            println("[+] complex Bloch nx=$nx f=$f propagating=",count(x->x.propagating,r)," min kappa=",first(r).attenuation_per_m);flush(stdout)
        end
        # Independent straight-strip modal solver predicts every propagating multiplier.
        straight=assemble_cell(;nx,ny,amplitude=0.)
        for f in (140e3,220e3)
            r=complex_bloch(straight,f)
            modes=solve_port_modes(PortModeConfig(height_m=4.5e-3,frequency_hz=f,element_count=160))
            for m in filter(m->m.kind==:propagating && m.direction==:right,modes)
                exact=exp(-im*real(m.wavenumber_per_m)*straight.period)
                phase_error=minimum(abs(complex(x.lambda_real,x.lambda_imag)-exact) for x in r if x.propagating)
                push!(checks,(nx,frequency_hz=f,k_per_m=real(m.wavenumber_per_m),multiplier_error=phase_error))
            end
        end
    end
    maximum(getproperty.(checks,:multiplier_error))<.02 || error("complex Bloch straight benchmark failed")
    jldsave(joinpath(OUT,"complex_bloch.jld2");allrows,checks)
    open(joinpath(OUT,"complex_bloch.csv"),"w") do io
        println(io,join(string.(propertynames(first(allrows))),','))
        for r in allrows;println(io,join(values(r),','));end
    end
end

function run()
    mkpath(OUT);cell=assemble_cell();q=collect(range(0,pi;length=65))
    basepath=joinpath(OUT,"bloch_base.jld2")
    if isfile(basepath)
        frequencies=JLD2.load(basepath,"frequencies")
    else
        results=map(q) do value
            b=bands(cell,value);println("[+] Bloch q/pi=",round(value/pi;digits=3));flush(stdout);b
        end
        frequencies=hcat(getproperty.(results,:frequencies)...)
        jldsave(basepath;q,frequencies)
    end
    minimum(frequencies[end,:])>360e3 || error("insufficient number of bands")
    fine=assemble_cell(nx=96,ny=60)
    checks=NamedTuple[]
    fine_frequencies=similar(frequencies)
    for i in eachindex(q)
        f=bands(fine,q[i]).frequencies
        fine_frequencies[:,i]=f
        println("[+] fine Bloch q/pi=",round(q[i]/pi;digits=3));flush(stdout)
        for band in eachindex(f)
            20e3<f[band]<360e3 || continue
            push!(checks,(q_over_pi=q[i]/pi,band,base_hz=frequencies[band,i],fine_hz=f[band],relative_error=abs(frequencies[band,i]/f[band]-1)))
        end
    end
    straight=assemble_cell(amplitude=0.);analytic_checks=NamedTuple[]
    for f in (140e3,220e3)
        modes=solve_port_modes(PortModeConfig(height_m=4.5e-3,frequency_hz=f,element_count=160))
        for mode in filter(m->m.kind==:propagating && m.direction==:right,modes)
            phase=acos(clamp(cos(real(mode.wavenumber_per_m)*cell.period),-1,1))
            b=bands(straight,phase).frequencies;nearest=b[argmin(abs.(b.-f))]
            push!(analytic_checks,(frequency_hz=f,k_per_m=real(mode.wavenumber_per_m),bloch_hz=nearest,relative_error=abs(nearest/f-1)))
        end
    end
    zero=bands(straight,0.).frequencies
    maximum(zero[1:2])<10 || error("rigid-body modes missing")
    maximum(getproperty.(checks,:relative_error))<.01 || error("Bloch mesh convergence exceeds 1%")
    maximum(getproperty.(analytic_checks,:relative_error))<.005 || error("straight-strip dispersion benchmark failed")
    jldsave(joinpath(OUT,"bloch.jld2");q,frequencies,fine_frequencies,checks,analytic_checks,period_m=cell.period,gap_m=cell.gap)
    open(joinpath(OUT,"bloch_checks.csv"),"w") do io
        println(io,"q_over_pi,band,base_hz,fine_hz,relative_error")
        for r in checks;println(io,join(values(r),','));end
    end
    println("[+] Bloch: mesh error=",maximum(getproperty.(checks,:relative_error))," straight benchmark=",maximum(getproperty.(analytic_checks,:relative_error)))
end
if abspath(PROGRAM_FILE)==@__FILE__
    "--complex" in ARGS ? run_complex() : run()
end
end
