# Analytic paraxial illustration, not an STG or defect FEM calculation.
ENV["GKSwstype"]="100"
using Plots
const OUT=@__DIR__
const EPSILON=0.05 # wavelength * focal_distance / aperture_width^2
# Julia sinc(q) = sin(pi*q)/(pi*q). The echo PSF squares the one-way
# Fresnel field and retains the residual complex quadratic phase.
echo_psf(q)=sinc(q)^2*exp(-2pi*im*EPSILON*q^2)

function main()
    q=collect(range(-2.6,2.6;length=3001))
    default(fontfamily="DejaVu Sans",linewidth=2.5,guidefontsize=10,tickfontsize=9)
    panels=[]
    for d in (nothing,0.45,1.4)
        positions=isnothing(d) ? [0.] : [-d/2,d/2]
        title=isnothing(d) ? "Один точечный дефект" : "Два дефекта, d = $d ℓ"
        p=plot(;title,xlabel="Положение фокуса / ℓ",xlims=(-2.6,2.6),ylims=(0,2.05),
            ylabel=isempty(panels) ? "Амплитуда эха, усл. ед." : "",legend=false,
            gridalpha=.15,framestyle=:box,titlefontsize=12)
        components=[echo_psf.(q.-x) for x in positions]
        if length(components)>1
            for c in components;plot!(p,q,abs.(c);color="#a3adb7",linestyle=:dash,linewidth=1.2);end
        end
        plot!(p,q,abs.(sum(components));color="#175b8f")
        for x in positions
            plot!(p,[x,x],[0,.15];color="#cf6037")
            scatter!(p,[x],[.15];color="#cf6037",marker=:dtriangle,markersize=4)
        end
        push!(panels,p)
    end
    fig=plot(panels...;layout=(1,3),size=(1320,410),margin=8Plots.mm)
    savefig(fig,joinpath(OUT,"echo_psf.png"))
    savefig(fig,joinpath(OUT,"echo_psf.svg"))
    lo=0.;hi=1.
    for _ in 1:60
        mid=(lo+hi)/2
        sinc(mid)^2>.5 ? (lo=mid) : (hi=mid)
    end
    width=lo+hi
    @assert abs(width-.88589294)<1e-7
    @assert abs(echo_psf(0)-1)<1e-12 && abs(echo_psf(1))<1e-12
    for (d,ismaximum) in ((.45,true),(1.4,false))
        f(x)=abs(echo_psf(x-d/2)+echo_psf(x+d/2))
        curvature=f(.001)+f(-.001)-2f(0)
        @assert ismaximum ? curvature<0 : curvature>0
    end
    println("Analytic checks passed; pulse-echo half-amplitude width = $width ell")
end
main()
