ENV["GKSwstype"]="100"
using JLD2, Plots
include(joinpath(@__DIR__,"profiles.jl"))
using .MetamaterialProfiles
const OUT=normpath(joinpath(@__DIR__,"..","..","..","results","stg_modal_completion"))

function raster(d)
    x=d["x_m"].*1e3;y=d["y_m"].*1e3;ux=d["ux_m"];uy=d["uy_m"]
    xx=range(0,maximum(x);step=.04);yy=range(0,7;step=.04)
    z=fill(NaN,length(yy),length(xx))
    for ids in d["triangles"]
        a,b,c=ids
        denominator=(y[b]-y[c])*(x[a]-x[c])+(x[c]-x[b])*(y[a]-y[c])
        abs(denominator)>eps() || continue
        for j in searchsortedfirst(yy,minimum(y[ids])):searchsortedlast(yy,maximum(y[ids])),
            i in searchsortedfirst(xx,minimum(x[ids])):searchsortedlast(xx,maximum(x[ids]))
            wa=((y[b]-y[c])*(xx[i]-x[c])+(x[c]-x[b])*(yy[j]-y[c]))/denominator
            wb=((y[c]-y[a])*(xx[i]-x[c])+(x[a]-x[c])*(yy[j]-y[c]))/denominator
            wc=1-wa-wb;min(wa,wb,wc)>=-1e-9 || continue
            vx=wa*ux[a]+wb*ux[b]+wc*ux[c];vy=wa*uy[a]+wb*uy[b]+wc*uy[c]
            z[j,i]=sqrt(abs2(vx)+abs2(vy))*1e9
        end
    end
    xx,yy,z
end

function main()
    fs=(140000,160000,200000,220000)
    data=[JLD2.load(joinpath(OUT,"cases","stg_n4_lead14_fine_f$(f)_loss0.0_ev24_p240_field.jld2")) for f in fs]
    limit=maximum(maximum(sqrt.(abs2.(d["ux_m"])+abs2.(d["uy_m"]))) for d in data)*1e9
    g=GeometryConfig(length_mm=.8+4.5*6.48,samples=400)
    b,t=generate_boundary_points(SinusoidalProfile(2.5;periods=4),g)
    outline=vcat([(0.,0.)],[(x+14,y) for (x,y) in b],[(g.length_mm+28,0.),(g.length_mm+28,7.)],
        [(x+14,y) for (x,y) in t],[(0.,7.),(0.,0.)])
    panels=[]
    for (d,f) in zip(data,fs)
        xx,yy,z=raster(d)
        p=heatmap(xx,yy,z;clims=(0,limit),color=:viridis,aspect_ratio=:equal,
            title="$(f÷1000) kHz",xlabel="x, mm",ylabel="y, mm",yticks=[0,3.5,7],
            legend=false,colorbar=true,colorbar_title="nm")
        plot!(p,first.(outline),last.(outline);color=:black,linewidth=.5,label=false)
        vline!(p,[14,14+g.length_mm];color=:gray,linestyle=:dash,linewidth=.8,label=false)
        push!(panels,p)
    end
    savefig(plot(panels...;layout=(4,1),size=(1200,1000),margin=4Plots.mm,
        plot_title="N=4 STG | no damping | displacement amplitude norm | incident S0 approximately 1 W/m"),joinpath(OUT,"fields_n4.png"))
end
if abspath(PROGRAM_FILE)==@__FILE__;main();end
