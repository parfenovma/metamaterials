using Gridap
using Gridap.ODEs
using GridapGmsh
using JLD2


const FREQ = 100e3
const MESH_FILE = "Acoustic_Metalens_F65.msh"
const SAVE_VTK = true
const VTK_DIR = "2_vtks_metalens"
const SIG_DIR = "3_signals_metalens"

mkpath(VTK_DIR)
mkpath(SIG_DIR)

function run_metalens_simulation()
    if !isfile(MESH_FILE)
        error("No $MESH_FILE")
    end

    println("=== Metalens (F=$(FREQ/1000) kHz) ===")
    model = GmshDiscreteModel(MESH_FILE)

    rho_solid = 1150.0; cp_solid = 2340.0; cs_solid = 1170.0        
    mu_solid = rho_solid * cs_solid^2; lam_solid = rho_solid * cp_solid^2 - 2*mu_solid 
    alpha_damp = 0.0; beta_damp = 2.0e-8 

    function pzt_signal(t)
        duration = 4.0 / FREQ
        t < duration ? 0.5 * (1.0 - cos(2.0 * pi * t / duration)) * sin(2.0 * pi * FREQ * t) : 0.0
    end

    reffe_vec = ReferenceFE(lagrangian, VectorValue{2, Float64}, 1)
    V0 = TestFESpace(model, reffe_vec, conformity=:H1) 
    U = TransientTrialFESpace(V0)

    degree = 2
    Ω = Triangulation(model); dΩ = Measure(Ω, degree)
    
    Γ_in = BoundaryTriangulation(model, tags=["Source"]); dΓ_in = Measure(Γ_in, degree)
    Γ_rad = BoundaryTriangulation(model, tags=["RadiationBoundary"]); dΓ_rad = Measure(Γ_rad, degree)

    n_in = get_normal_vector(Γ_in)
    n_rad = get_normal_vector(Γ_rad)

    function absorbing_traction(∂tu, n)
        ∂tu_n = (∂tu ⋅ n) * n
        ∂tu_t = ∂tu - ∂tu_n
        return rho_solid * cp_solid * ∂tu_n + rho_solid * cs_solid * ∂tu_t
    end

    σ(ε) = lam_solid * tr(ε) * one(ε) + 2.0 * mu_solid * ε
    P_amplitude = 1e6 

    res(t, u, v) = ∫( rho_solid * ∂tt(u) ⋅ v + alpha_damp * rho_solid * ∂t(u) ⋅ v + 
                      beta_damp * (σ∘(ε(∂t(u))) ⊙ ε(v)) + σ∘(ε(u)) ⊙ ε(v) )dΩ +
                   ∫( absorbing_traction(∂t(u), n_rad) ⋅ v )dΓ_rad + # <--- Излучение наружу!
                   ∫( absorbing_traction(∂t(u), n_in) ⋅ v )dΓ_in -
                   ∫( (P_amplitude * pzt_signal(t)) * (n_in ⋅ v) )dΓ_in

    jac(t, u, du, v) = ∫( σ∘(ε(du)) ⊙ ε(v) )dΩ
    jac_t(t, u, dut, v) = ∫( alpha_damp * rho_solid * dut ⋅ v + beta_damp * (σ∘(ε(dut)) ⊙ ε(v)) )dΩ +
                          ∫( absorbing_traction(dut, n_rad) ⋅ v )dΓ_rad + 
                          ∫( absorbing_traction(dut, n_in) ⋅ v )dΓ_in
    jac_tt(t, u, dutt, v) = ∫( rho_solid * dutt ⋅ v )dΩ

    t1 = 140.0e-6
    dt = (1.0 / FREQ) / 30.0  
    
    op = TransientFEOperator(res, (jac, jac_t, jac_tt), U, V0)
    
    U_at_t0 = U(0.0)
    uh0 = interpolate_everywhere(VectorValue(0.0, 0.0), U_at_t0)
    vh0 = interpolate_everywhere(VectorValue(0.0, 0.0), U_at_t0)

    nonlinear_solver = NLSolver(show_trace=false, method=:newton) 
    ode_solver = Newmark(nonlinear_solver, dt, 0.5, 0.25)
    
    println("  -> Сборка матриц и запуск Time Loop (dt=$(round(dt*1e6, digits=3)) мкс, tf=$(t1*1e6) мкс)...")
    sol_t = solve(ode_solver, op, 0.0, t1, (uh0, vh0))

    time_history = Float64[]
    focal_signal_MPa = Float64[]

    focal_point = Point(47.0 * 1e-3, 0.0) 

    pvd = SAVE_VTK ? createpvd(joinpath(VTK_DIR, "Metalens_Focus_Animation")) : nothing
    step = 0
    
    for (tn, uh) in sol_t
        push!(time_history, tn)
        u_focal = uh(focal_point) 
        push!(focal_signal_MPa, norm(u_focal))
        
        if step % 20 == 0
            println("    [Evaluated] t = $(round(tn*1e6, digits=2)) microseconds $(round(t1*1e6, digits=0))")
        end
        
        if SAVE_VTK && step % 3 == 0
            pvd[tn] = createvtk(Ω, joinpath(VTK_DIR, "Metalens_$(step).vtu"), cellfields=["u"=>uh, "Stress"=>σ∘(ε(uh))])
        end
        step += 1
    end
    
    if SAVE_VTK; savepvd(pvd); end

    in_signal_history = (P_amplitude / 1e6) .* pzt_signal.(time_history)

    save_path = joinpath(SIG_DIR, "Metalens_Data.jld2")
    jldsave(save_path; 
        time = time_history, 
        signal_in = in_signal_history, 
        signal_focal = focal_signal_MPa
    )
    
    println("=== Done! ===")
    println("Animations: $VTK_DIR.")
    println("focus data: $save_path.")
end

run_metalens_simulation()
