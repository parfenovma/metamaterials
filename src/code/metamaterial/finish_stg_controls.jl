include(joinpath(@__DIR__,"run_stg_completion.jl"))
using .STGCompletion
using JLD2, LinearAlgebra

function resolved_bank(base;elements=240)
    E=STGCompletion.ElasticPortModes
    allm=E.solve_port_modes(E.PortModeConfig(height_m=7e-3,element_count=elements,frequency_hz=base.frequency))
    right=filter(m->m.kind==:propagating && m.direction==:right,allm)
    labels=String[];canonical=E.PortMode[]
    for m in right
        candidates=findall(old->old.parity==m.parity,base.right)
        i=candidates[argmin([abs(base.right[j].wavenumber_per_m-m.wavenumber_per_m) for j in candidates])]
        push!(labels,base.labels[i]);old=base.right[i]
        interp=STGCompletion.ElasticModalProjection.linear_interpolate
        oldvec=vcat([interp(old.y_m,old.displacement_x,y) for y in m.y_m],
                    [interp(old.y_m,old.displacement_y,y) for y in m.y_m])
        phase=dot(oldvec,STGCompletion.vector(m))
        push!(canonical,STGCompletion.rotate_mode(m,conj(phase)/abs(phase)))
    end
    length(unique(labels))==length(labels) || error("ambiguous resolution mapping")
    left=E.PortMode[]
    for m in canonical
        raw=STGCompletion.ElasticModalProjection.counterpropagating_mode(m,allm)
        phase=dot(conj.(STGCompletion.vector(m)),STGCompletion.vector(raw))
        push!(left,STGCompletion.rotate_mode(raw,conj(phase)/abs(phase)))
    end
    ev=filter(m->m.kind==:evanescent,allm)
    (;frequency=base.frequency,modes=vcat(canonical,left,ev),right=canonical,labels)
end

function finish_controls()
    banks=JLD2.load(STGCompletion.modes_file(),"banks");rows=NamedTuple[]
    for f in (130e3,140e3,160e3,200e3,220e3)
        b=only(filter(b->b.frequency==f,banks))
        for n in (1,4)
            push!(rows,STGCompletion.solve_case(n,"straight",1.,b))
            f==200e3 && push!(rows,STGCompletion.solve_case(n,"stg",1.,b))
        end
        if f in (140e3,160e3,200e3,220e3)
            fine=resolved_bank(b)
            push!(rows,STGCompletion.solve_case(4,"stg",0.,fine;fine=true,ev=24,save_field=true))
            f in (160e3,200e3) && push!(rows,STGCompletion.solve_case(4,"stg",1.,fine;fine=true,ev=24))
        end
    end
    STGCompletion.write_rows(joinpath(STGCompletion.OUT,"length_controls.csv"),rows)
end
if abspath(PROGRAM_FILE)==@__FILE__;finish_controls();end
