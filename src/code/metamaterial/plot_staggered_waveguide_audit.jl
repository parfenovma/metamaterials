ENV["GKSwstype"] = "100"
using Plots, JLD2
include(joinpath(@__DIR__, "profiles.jl"))
using .MetamaterialProfiles

const AUDIT = normpath(joinpath(@__DIR__, "..", "..", "..", "results", "staggered_waveguide_port_audit"))
function read_table(name)
    lines = readlines(joinpath(AUDIT, name))
    keys = Symbol.(split(first(lines), ','))
    value(s) = s == "true" ? true : s == "false" ? false : something(tryparse(Float64, s), s)
    [NamedTuple{Tuple(keys)}(Tuple(value.(split(line, ',')))) for line in lines[2:end]]
end

# Interpolate the saved nodal complex field on its own triangles; do not fill
# the notches with an unconstrained Delaunay triangulation of all sample points.
function field_raster(data, component)
    x=data["x_m"].*1e3; y=data["y_m"].*1e3; u=data[component]
    xx=range(0,45;step=0.04); yy=range(0,7;step=0.04)
    z=fill(NaN,length(yy),length(xx))
    for ids in data["triangles"]
        a,b,c=ids
        denominator=(y[b]-y[c])*(x[a]-x[c])+(x[c]-x[b])*(y[a]-y[c])
        abs(denominator)>eps() || continue
        for j in searchsortedfirst(yy,minimum(y[ids])):searchsortedlast(yy,maximum(y[ids])),
            i in searchsortedfirst(xx,minimum(x[ids])):searchsortedlast(xx,maximum(x[ids]))
            wa=((y[b]-y[c])*(xx[i]-x[c])+(x[c]-x[b])*(yy[j]-y[c]))/denominator
            wb=((y[c]-y[a])*(xx[i]-x[c])+(x[a]-x[c])*(yy[j]-y[c]))/denominator
            wc=1-wa-wb
            min(wa,wb,wc)>=-1e-9 || continue
            z[j,i]=abs(wa*u[a]+wb*u[b]+wc*u[c])*1e9
        end
    end
    xx,yy,z
end

function audit_plots()
    rows = read_table("frequency_scan.csv")
    selected(role, loss) = sort(filter(r -> r.geometry == role && r.lossy == loss, rows); by=r -> r.frequency_hz)
    default(fontfamily="Helvetica", framestyle=:box, linewidth=2, markersize=5,
            titlefontsize=11, guidefontsize=10, tickfontsize=9, legendfontsize=9)
    bottom, top = generate_boundary_points(SinusoidalProfile(2.5; periods=2), GeometryConfig())
    points = vcat([(0.,0.)], [(x+14,y) for (x,y) in bottom], [(45.,0.),(45.,7.)],
                  [(x+14,y) for (x,y) in top], [(0.,7.)])
    geometry = plot(Shape(first.(points), last.(points)); fillcolor=:lightsteelblue,
        linecolor=:navy, label=false, aspect_ratio=:equal, xlims=(-1,46), ylims=(-1,9),
        xlabel="x, mm", ylabel="y, mm", title="Same STG core + straight port extensions")
    vline!(geometry,[14,31]; color=:gray, linestyle=:dash, label=false)
    annotate!(geometry,22.5,8.2,text("17 mm core",9))
    annotate!(geometry,7,8.2,text("14 mm",9))
    annotate!(geometry,38,8.2,text("14 mm",9))
    function power_panel(loss)
        r = selected("stg",loss); f = getproperty.(r,:frequency_hz)./1000
        p = plot(; title=loss ? "STG: nominal core damping" : "STG: no material damping",
            xlabel="Frequency, kHz", ylabel="Fraction of incident power", ylims=(-0.02,1.03),legend=:outertop,
            legend_columns=loss ? 4 : 3)
        for (field,label,color) in ((:T_all,"T (all modes)",:royalblue),(:R_all,"R (all modes)",:darkorange),(:D,"D",:gray))
            !loss && field==:D && continue
            plot!(p,f,getproperty.(r,field);marker=:circle,color,label,linestyle=:dot)
        end
        ref = selected("straight",loss)
        plot!(p,f,getproperty.(ref,:T_all);marker=:diamond,color=:black,linestyle=:dash,label="Straight T")
        p
    end
    mode = plot(;title="Transmission: original S0 and converted modes",xlabel="Frequency, kHz",
        ylabel="Fraction of incident power",ylims=(-0.02,1.03),legend=:outertop,legend_columns=2)
    for loss in (false,true)
        r=selected("stg",loss);f=getproperty.(r,:frequency_hz)./1000
        for (field,label,color) in ((:T00,"S0",:royalblue),(:T_converted,"Other modes",:darkorange))
            plot!(mode,f,getproperty.(r,field);marker=loss ? :diamond : :circle,
                linestyle=loss ? :dash : :dot,color,label=label*(loss ? ", damped" : ", lossless"))
        end
    end
    fig=plot(geometry,power_panel(false),power_panel(true),mode;layout=(2,2),size=(1150,750),
        margin=5Plots.mm,plot_title="Asymmetric STG | one incident S0 mode | six sampled frequencies")
    savefig(fig,joinpath(AUDIT,"overview.png"))

    fields=[JLD2.load(joinpath(AUDIT,"cases","stg_lead14p0_f$(f)_lossless_ev8_field.jld2")) for f in (140000,220000)]
    panels=[]
    for (data,f) in zip(fields,(140,220)), component in ("ux_m","uy_m")
        maximum_nm=maximum(maximum(abs.(d[component])) for d in fields)*1e9
        xx,yy,z=field_raster(data,component)
        p=heatmap(xx,yy,z;color=:viridis,clims=(0,maximum_nm),label=false,
            aspect_ratio=:equal,xlabel="x, mm",ylabel="y, mm",xlims=(0,45),ylims=(0,7),
            yticks=[0,3.5,7],colorbar_ticks=[0,round(maximum_nm/2;digits=1),round(maximum_nm;digits=1)],
            title="$f kHz: |$(component[1:2])|, nm")
        plot!(p,first.(vcat(points,[first(points)])),last.(vcat(points,[first(points)]));color=:black,linewidth=0.6,label=false)
        push!(panels,p)
    end
    fig=plot(panels...;layout=(4,1),size=(1150,980),margin=5Plots.mm,
        plot_title="Distributed displacement | lossless STG | incident S0 approximately 1 W/m")
    savefig(fig,joinpath(AUDIT,"fields.png"))
end

if abspath(PROGRAM_FILE)==@__FILE__;audit_plots();end
