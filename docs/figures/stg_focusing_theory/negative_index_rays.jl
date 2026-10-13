# Geometrical-ray illustration only. Homogeneous isotropic n = -1 slab,
# matched to n = +1 exterior. It does not calculate evanescent waves or STG.
ENV["GKSwstype"]="100"
using Plots
const OUT=@__DIR__

function main()
    a=1.;d=2.;ximage=2d-a
    panels=[]
    default(fontfamily="DejaVu Sans",guidefontsize=10,tickfontsize=9,titlefontsize=12)
    for negative in (false,true)
        p=plot(;title=negative ? "Пластина n = −1: лучи сходятся" : "Пластина n = +1: лучи расходятся",
            xlabel="Продольная координата, усл. ед.",ylabel="Поперечная координата",
            xlims=(-1.4,3.8),ylims=(-2.3,2.3),legend=false,aspect_ratio=:equal,
            grid=false,framestyle=:box)
        vspan!(p,[0,d];color=negative ? "#d2ebfa" : "#e5e9ed",alpha=.7)
        for slope in (-.5,-.25,0.,.25,.5)
            # Energy rays travel toward increasing x in both cases.
            y0=a*slope
            yd=y0+(negative ? -1 : 1)*d*slope
            x=[-a,0.,d,3.8];y=[0.,y0,yd,yd+(3.8-d)*slope]
            plot!(p,x,y;color="#1c6b93",linewidth=1.7)
            if negative
                @assert isapprox(y0-a*slope,0.;atol=1e-12)
                @assert isapprox(yd+(ximage-d)*slope,0.;atol=1e-12)
            end
        end
        scatter!(p,[-a],[0.];color="#cf6037",markersize=5)
        annotate!(p,-a,-.45,text("Источник",9))
        if negative
            scatter!(p,[ximage],[0.];color="#cf6037",markersize=5)
            annotate!(p,ximage,-.45,text("Изображение",9))
        end
        annotate!(p,1.,2.,text("Плоские границы",9))
        push!(panels,p)
    end
    fig=plot(panels...;layout=(1,2),size=(1160,490),margin=8Plots.mm)
    savefig(fig,joinpath(OUT,"negative_index_rays.png"))
    savefig(fig,joinpath(OUT,"negative_index_rays.svg"))
    println("Ray intersections verified; no wave-resolution claim.")
end
main()
