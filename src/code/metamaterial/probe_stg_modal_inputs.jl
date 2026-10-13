include(joinpath(@__DIR__,"run_stg_completion.jl"))
using .STGCompletion
using JLD2, LinearAlgebra

function probe_inputs()
    banks=JLD2.load(STGCompletion.modes_file(),"banks");rows=NamedTuple[];checks=NamedTuple[]
    M=STGCompletion.MultimodePortHarmonicSolver
    for n in (2,4),f in (200e3,220e3)
        b=only(filter(b->b.frequency==f,banks));mat=STGCompletion.material(0.)
        cfg=M.MultimodePortConfig(frequency_hz=f,element_order=2,quadrature_degree=6,
            port_element_count=120,evanescent_modes_per_port=12,damping_domain_tag="Core")
        mesh=STGCompletion.modelpath(n,"stg");prep=M.prepare_multimode_model(mesh,mat;config=cfg)
        count=length(b.right);A=zeros(ComplexF64,count,count);B=copy(A);C=copy(A)
        for (j,mode) in enumerate(b.right)
            result=M.solve_multimode_port_harmonic(mesh,["Microphone"],mat;
                config=cfg,prepared=prep,supplied_modes=b.modes,incident_mode=mode)
            A[:,j]=getproperty.(result.source_amplitudes,:right_amplitude)
            B[:,j]=getproperty.(result.source_amplitudes,:left_amplitude)
            C[:,j]=getproperty.(only(result.output_amplitudes),:right_amplitude)
        end
        R=B/A;T=C/A
        power_error=opnorm(T'*T+R'*R-I)
        reciprocity_error=norm(R-transpose(R))/norm(R)
        push!(checks,(n,frequency_hz=f,power_error,reflection_reciprocity_error=reciprocity_error,
            maximum_transmission=maximum(svdvals(T))^2))
        for (j,label) in enumerate(b.labels)
            push!(rows,(n,frequency_hz=f,incident_branch=label,T=sum(abs2,T[:,j]),R=sum(abs2,R[:,j])))
        end
        jldsave(joinpath(STGCompletion.OUT,"modal_input_matrix_n$(n)_f$(round(Int,f)).jld2");R,T,labels=b.labels,checks=last(checks))
        println(last(checks));flush(stdout)
    end
    STGCompletion.write_rows(joinpath(STGCompletion.OUT,"modal_input_scan.csv"),rows)
    STGCompletion.write_rows(joinpath(STGCompletion.OUT,"modal_matrix_checks.csv"),checks)
    all(r.power_error<.01 && r.reflection_reciprocity_error<.01 for r in checks) || error("modal matrix checks failed")
end
if abspath(PROGRAM_FILE)==@__FILE__;probe_inputs();end
