module AnalyzeSTGCompletion
ENV["GKSwstype"]="100"
using JLD2, FFTW, LinearAlgebra, Plots
include(joinpath(@__DIR__,"port_mode_solver.jl"))
include(joinpath(@__DIR__,"wavelet_analysis.jl"))
using .WaveletAnalysis
const OUT=normpath(joinpath(@__DIR__,"..","..","..","results","stg_modal_completion"))

function write_rows(path,rows)
    isempty(rows) && return
    open(path,"w") do io
        println(io,join(string.(propertynames(first(rows))),','))
        for r in rows;println(io,join(string.(values(r)),','));end
    end
end
function load_cases()
    records=[]
    for p in readdir(joinpath(OUT,"cases");join=true)
        endswith(p,".jld2") && !endswith(p,"_field.jld2") || continue
        d=JLD2.load(p);push!(records,(path=p,data=d,row=d["row"]))
    end
    records
end
function select(records,n,role,loss)
    sort(filter(d->d.row.n==n && d.row.geometry==role && d.row.loss==loss &&
        !d.row.fine && d.row.lead_mm==14 && occursin("_ev12_p120",d.path),records);by=d->d.row.frequency_hz)
end
function interpolate(x,y,q)
    (q<first(x) || q>last(x)) && return zero(eltype(y))
    i=searchsortedlast(x,q);i==length(x) && return y[end]
    a=(q-x[i])/(x[i+1]-x[i]);(1-a)*y[i]+a*y[i+1]
end
function spectra(records;stride=1,df=2500.)
    selected=records[1:stride:end]
    f=[d.row.frequency_hz for d in selected]
    labels=sort(unique(vcat([d.data["labels"] for d in selected]...)))
    H=zeros(ComplexF64,length(f),length(labels));R=copy(H);V=copy(H)
    for (i,d) in enumerate(selected),(j,label) in enumerate(labels)
        k=findfirst(==(label),d.data["labels"]);isnothing(k) && continue
        H[i,j]=d.data["transmission"][k];R[i,j]=d.data["reflection"][k]
        V[i,j]=H[i,j]*d.data["mean_vx"][k]
    end
    (;f,labels,H,R,V,D=[d.row.D for d in selected])
end
function envelope(x)
    n=length(x);a=zeros(ComplexF64,n);f=rfft(x)
    a[1]=f[1];a[2:n÷2]=2f[2:n÷2];a[n÷2+1]=f[end]
    abs.(ifft(a))
end

function packet(s,carrier;df=2500.)
    dt=1e-7;n=round(Int,1/(df*dt));t=collect(0:n-1).*dt;f=collect(0:n÷2).*df
    start=50e-6;duration=4/carrier
    input=[0<=v-start<=duration ? sin(2pi*carrier*(v-start))*(1-cos(2pi*(v-start)/duration))/2 : 0. for v in t]
    original=rfft(input);a=copy(original);mask=(f.>=20e3).&(f.<=360e3);a[.!mask].=0
    coverage=sum(abs2,a)/sum(abs2,original);incoming=irfft(a,n)
    H=hcat([[interpolate(s.f,s.H[:,j],v) for v in f] for j in eachindex(s.labels)]...)
    R=hcat([[interpolate(s.f,s.R[:,j],v) for v in f] for j in eachindex(s.labels)]...)
    V=hcat([[interpolate(s.f,s.V[:,j],v) for v in f] for j in eachindex(s.labels)]...)
    output=hcat([irfft(a.*H[:,j],n) for j in eachindex(s.labels)]...)
    reflected=hcat([irfft(a.*R[:,j],n) for j in eachindex(s.labels)]...)
    scalar_vx=irfft(a.*vec(sum(V;dims=2)),n)
    normin=sum(abs2,incoming);energy=vec(sum(abs2,output;dims=1))./normin
    reflected_energy=vec(sum(abs2,reflected;dims=1))./normin
    dissipated=sum(abs2.(a).*[interpolate(s.f,s.D,v) for v in f])/sum(abs2,a)
    edge=(t.<20e-6).|(t.>maximum(t)-30e-6)
    edge_fraction=sum(abs2,output[edge,:])/max(sum(abs2,output),eps())
    centroids=[sum(t.*envelope(output[:,j]).^2)/max(sum(envelope(output[:,j]).^2),eps())-start for j in eachindex(s.labels)]
    (;t=t.-start,incoming,output,reflected,scalar_vx,labels=s.labels,energy,reflected_energy,
        dissipated,coverage,edge_fraction,centroids)
end

function cwt_plot(p,carrier,loss)
    order=sortperm(p.energy;rev=true)
    active=filter(i->p.energy[i]>0.001,order)
    selected=first(active,min(5,length(active)))
    signals=vcat([p.incoming],[p.output[:,i] for i in selected])
    titles=vcat(["Incident S0"],["$(p.labels[i]): $(round(100p.energy[i];digits=2))% of incident packet energy" for i in selected])
    cfg=MorletConfig(frequencies_hz=collect(40e3:5e3:360e3))
    transforms=[continuous_wavelet_transform(p.t,s;config=cfg) for s in signals]
    normalization=maximum(abs2.(transforms[1].coefficients))
    panels=[]
    for (c,title) in zip(transforms,titles)
        power=10log10.(abs2.(c.coefficients)./normalization.+1e-12)
        power[.!c.cone_of_influence_valid].=NaN
        push!(panels,heatmap(c.time_s.*1e6,c.frequency_hz./1e3,power;clims=(-45,0),
            color=:inferno,title,xlims=(-15,170),ylims=(40,360),ylabel="kHz",xlabel="Time from pulse start, us",
            colorbar_title="dB",colorbar=true,legend=false))
    end
    fig=plot(panels...;layout=(length(panels),1),size=(1000,210length(panels)),margin=4Plots.mm,
        plot_title="STG N=2 | $(round(Int,carrier/1000)) kHz Hann pulse | loss scale $loss | same color scale")
    savefig(fig,joinpath(OUT,"modal_scalograms_$(round(Int,carrier/1000))k_loss$(loss).png"))
end

function analyze_packets(records)
    summary=NamedTuple[];mode_rows=NamedTuple[]
    for loss in (0.,1.),carrier in (140e3,220e3)
        cases=select(records,2,"stg",loss)
        length(cases)>=137 || error("incomplete STG frequency sweep")
        s=spectra(cases);frequency_step=minimum(diff(s.f))
        all(isapprox.(diff(s.f),frequency_step)) || error("incomplete uniform frequency sweep")
        first(s.f)==20e3 && last(s.f)==360e3 || error("wrong spectral coverage")
        p=packet(s,carrier;df=frequency_step)
        coarse=packet(spectra(cases;stride=2),carrier;df=frequency_step)
        waveform_error=norm(p.output-coarse.output)/norm(p.output)
        energy_error=abs(sum(p.energy)-sum(coarse.energy))
        if loss==1
            ref=packet(spectra(select(records,2,"straight",1.)),carrier;df=frequency_step)
            reference_transmission=sum(ref.energy)
        else;reference_transmission=1.;end
        antisymmetric=sum(p.energy[startswith.(p.labels,"A")])
        push!(summary,(carrier_hz=carrier,loss,frequency_step_hz=frequency_step,transmitted=sum(p.energy),reflected=sum(p.reflected_energy),
            dissipated=p.dissipated,balance=1-sum(p.energy)-sum(p.reflected_energy)-p.dissipated,
            straight_transmitted=reference_transmission,relative_transmission=sum(p.energy)/reference_transmission,
            antisymmetric_transmitted=antisymmetric,spectral_coverage=p.coverage,
            coarse_fine_waveform_error=waveform_error,coarse_fine_energy_error=energy_error,edge_fraction=p.edge_fraction))
        for (i,label) in enumerate(p.labels)
            push!(mode_rows,(carrier_hz=carrier,loss,branch=label,T=p.energy[i],R=p.reflected_energy[i],centroid_us=p.centroids[i]*1e6))
        end
        jldsave(joinpath(OUT,"packet_$(round(Int,carrier/1000))k_loss$(loss).jld2");p,coarse,summary=last(summary))
        cwt_plot(p,carrier,loss)
        println(last(summary));flush(stdout)
    end
    write_rows(joinpath(OUT,"packet_summary.csv"),summary);write_rows(joinpath(OUT,"packet_modes.csv"),mode_rows)
    summary
end

function packet_comparison(records)
    panels=[]
    for carrier in (140e3,220e3)
        s=spectra(select(records,2,"stg",1.));df=minimum(diff(s.f))
        p=packet(s,carrier;df);r=packet(spectra(select(records,2,"straight",1.)),carrier;df)
        scale=maximum(envelope(p.incoming))
        plotpanel=plot(;xlabel="Time from pulse start, us",ylabel="Modal envelope / incident peak",
            title="$(round(Int,carrier/1000)) kHz four-cycle pulse; nominal damping",xlims=(-5,120),
            legend=:outertop,legend_columns=3)
        vspan!(plotpanel,[0,4/carrier*1e6];color=:gray,alpha=.08,label=false)
        for i in first(sortperm(p.energy;rev=true),4)
            plot!(plotpanel,p.t.*1e6,envelope(p.output[:,i])./scale;label="STG $(p.labels[i])")
        end
        i=findfirst(==("S0"),r.labels)
        plot!(plotpanel,r.t.*1e6,envelope(r.output[:,i])./scale;color=:black,linestyle=:dash,label="Straight S0")
        push!(panels,plotpanel)
    end
    savefig(plot(panels...;layout=(2,1),size=(1050,700),margin=5Plots.mm),joinpath(OUT,"pulse_envelopes.png"))
end

function overview(records)
    default(fontfamily="Helvetica",titlefontsize=11,guidefontsize=10,tickfontsize=9,legendfontsize=9,linewidth=2)
    bloch=JLD2.load(joinpath(OUT,"bloch.jld2"));f=get(bloch,"fine_frequencies",bloch["frequencies"])
    q=bloch["q"];gaps=[(maximum(f[i,:]),minimum(f[i+1,:])) for i in 1:size(f,1)-1 if minimum(f[i+1,:])-maximum(f[i,:])>1e3 && minimum(f[i+1,:])<360e3]
    p1=plot(;xlabel="Bloch phase / pi",ylabel="Frequency, kHz",title="Infinite periodic STG interior",legend=false,ylims=(20,300))
    for band in axes(f,1);plot!(p1,q./pi,f[band,:]./1e3;color=:navy,linewidth=1);end
    for (a,b) in gaps;hspan!(p1,[a,b]./1e3;color=:orange,alpha=.2);end
    p2=plot(;xlabel="Frequency, kHz",ylabel="Transmitted power",title="Same period and end transitions; no damping",ylims=(0,1),xlims=(80,280),legend=:outertop,legend_columns=3)
    for (a,b) in gaps;vspan!(p2,[a,b]./1e3;color=:orange,alpha=.15,label=false);end
    for n in (1,2,4)
        r=select(records,n,"stg",0.);plot!(p2,[d.row.frequency_hz/1000 for d in r],[d.row.T for d in r];label="N=$n")
    end
    p3=plot(;xlabel="Frequency, kHz",ylabel="Power fraction",title="Original N=2 STG; nominal damping",ylims=(0,1),xlims=(60,300),legend=:outertop,legend_columns=4)
    r=select(records,2,"stg",1.)
    for (key,label) in ((:T,"T"),(:R,"R"),(:D,"D"));plot!(p3,[d.row.frequency_hz/1000 for d in r],[getproperty(d.row,key) for d in r];label);end
    ref=select(records,2,"straight",1.);plot!(p3,[d.row.frequency_hz/1000 for d in ref],[d.row.T for d in ref];label="Straight T",color=:black,linestyle=:dash)
    p4=plot(;xlabel="Frequency, kHz",ylabel="Transmitted power",title="N=2: modal conversion without damping",ylims=(0,1),xlims=(60,300),legend=:outertop,legend_columns=3)
    r=select(records,2,"stg",0.)
    for (key,label) in ((:T,"All modes"),(:T00,"S0"),(:Tother,"Other modes"));plot!(p4,[d.row.frequency_hz/1000 for d in r],[getproperty(d.row,key) for d in r];label);end
    savefig(plot(p1,p2,p3,p4;layout=(2,2),size=(1200,850),margin=5Plots.mm),joinpath(OUT,"overview.png"))
    write_rows(joinpath(OUT,"bloch_gaps.csv"),[(lower_hz=a,upper_hz=b) for (a,b) in gaps])
end

function main()
    records=load_cases();overview(records);summary=analyze_packets(records);packet_comparison(records)
    checks=[(check="spectral_coverage",value=minimum(r.spectral_coverage for r in summary),limit=.999,passed=all(r.spectral_coverage>=.999 for r in summary)),
        (check="frequency_step_waveform",value=maximum(r.coarse_fine_waveform_error for r in summary),limit=.05,passed=all(r.coarse_fine_waveform_error<.05 for r in summary)),
        (check="frequency_step_packet_energy",value=maximum(r.coarse_fine_energy_error for r in summary),limit=.01,passed=all(r.coarse_fine_energy_error<.01 for r in summary)),
        (check="packet_power_balance",value=maximum(abs(r.balance) for r in summary),limit=.01,passed=all(abs(r.balance)<.01 for r in summary)),
        (check="time_window_edge_energy",value=maximum(r.edge_fraction for r in summary),limit=.005,passed=all(r.edge_fraction<.005 for r in summary))]
    write_rows(joinpath(OUT,"packet_checks.csv"),checks)
    all(r.passed for r in checks) || error("pulse reconstruction did not meet frozen checks")
end
if abspath(PROGRAM_FILE)==@__FILE__;main();end
end
