include(joinpath(@__DIR__,"run_stg_completion.jl"))
using .STGCompletion
using JLD2

function validate_completion()
    finebanks=STGCompletion.build_modes(collect(20e3:5e3:260e3);elements=240,save=false)
    basebanks=JLD2.load(STGCompletion.modes_file(),"banks")
    results=NamedTuple[]
    for f in (140e3,160e3,200e3,220e3,260e3),loss in (0.,1.)
        b=only(filter(b->b.frequency==f,finebanks))
        row=STGCompletion.solve_case(2,"stg",loss,b;fine=true,ev=24)
        push!(results,merge((check="mesh_port_basis",),row))
    end
    for f in (140e3,220e3),loss in (0.,1.)
        b=only(filter(b->b.frequency==f,finebanks))
        row=STGCompletion.solve_case(2,"straight",loss,b;fine=true,ev=24)
        push!(results,merge((check="straight_control",),row))
    end
    for f in (140e3,152.5e3,220e3),loss in (0.,1.)
        b=only(filter(b->b.frequency==f,basebanks))
        row=STGCompletion.solve_case(2,"stg",loss,b;lead=28.)
        push!(results,merge((check="lead_length",),row))
    end
    STGCompletion.write_rows(joinpath(STGCompletion.OUT,"validation_cases.csv"),results)
end
if abspath(PROGRAM_FILE)==@__FILE__;validate_completion();end
