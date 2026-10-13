include(joinpath(@__DIR__,"run_stg_completion.jl"))
using .STGCompletion
using JLD2, LinearAlgebra
const EXTRA=joinpath(STGCompletion.OUT,"midpoint_modes.jld2")

function midpoint_modes(base=JLD2.load(STGCompletion.modes_file(),"banks");indices=1:length(base)-1,destination=EXTRA)
    banks=[]
    E=STGCompletion.ElasticPortModes
    for i in indices
        f=(base[i].frequency+base[i+1].frequency)/2;reference=base[i+1]
        modes=E.solve_port_modes(E.PortModeConfig(height_m=7e-3,element_count=120,frequency_hz=f))
        right=filter(m->m.kind==:propagating && m.direction==:right,modes)
        labels=fill("",length(right));used=Set{String}();overlaps=zeros(length(right));canonical=copy(right)
        pairs=sort([(E.mode_overlap(old,m),label,j,old) for (old,label) in zip(reference.right,reference.labels)
                    for (j,m) in enumerate(right) if old.parity==m.parity];by=first,rev=true)
        for (overlap,label,j,old) in pairs
            (label in used || !isempty(labels[j])) && continue
            overlap>.98 || error("ambiguous refined mode at $f")
            labels[j]=label;push!(used,label);overlaps[j]=overlap
            phase=dot(STGCompletion.vector(old),STGCompletion.vector(right[j]))
            canonical[j]=STGCompletion.rotate_mode(right[j],conj(phase)/abs(phase))
        end
        any(isempty,labels) && error("unmatched refined mode")
        left=E.PortMode[]
        for m in canonical
            raw=STGCompletion.ElasticModalProjection.counterpropagating_mode(m,modes)
            phase=dot(conj.(STGCompletion.vector(m)),STGCompletion.vector(raw))
            push!(left,STGCompletion.rotate_mode(raw,conj(phase)/abs(phase)))
        end
        kept=vcat(canonical,left)
        for direction in (:right,:left)
            ev=sort(filter(m->m.kind==:evanescent && m.direction==direction,modes);by=m->abs(real(m.gamma_per_m)))
            append!(kept,first(ev,min(24,length(ev))))
        end
        push!(banks,(frequency=f,modes=kept,right=canonical,labels,overlaps))
    end
    jldsave(destination;banks)
    banks
end
function solve_midpoints(group)
    banks=JLD2.load(EXTRA,"banks");rows=NamedTuple[]
    for b in banks
        if group==0
            for loss in (0.,1.);push!(rows,STGCompletion.solve_case(2,"stg",loss,b));end
        else;push!(rows,STGCompletion.solve_case(2,"straight",1.,b));end
    end
    STGCompletion.write_rows(joinpath(STGCompletion.OUT,"refined_group$(group).csv"),rows)
end
function second_refinement(group)
    base=sort(vcat(JLD2.load(STGCompletion.modes_file(),"banks"),JLD2.load(EXTRA,"banks"));by=b->b.frequency)
    destination=joinpath(STGCompletion.OUT,"second_modes_group$(group).jld2")
    banks=isfile(destination) ? JLD2.load(destination,"banks") :
        midpoint_modes(base;indices=(group+1):2:(length(base)-1),destination)
    rows=[STGCompletion.solve_case(2,"stg",0.,b) for b in banks]
    STGCompletion.write_rows(joinpath(STGCompletion.OUT,"second_refined_group$(group).csv"),rows)
end
if abspath(PROGRAM_FILE)==@__FILE__
    if ARGS[1]=="--modes";midpoint_modes()
    elseif ARGS[1]=="--level2";second_refinement(parse(Int,ARGS[2]))
    else;solve_midpoints(parse(Int,ARGS[2]));end
end
