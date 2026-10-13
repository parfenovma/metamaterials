module ExportO8MeshAnimations

ENV["GKSwstype"] = "100"

using JLD2
using LinearAlgebra
using Plots

include(joinpath(@__DIR__, "operator_o8_plane_strain_solver.jl"))

using .OperatorO8PlaneStrainSolver

export export_o8_mesh_animations

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const PACKAGE_ROOT = joinpath(PROJECT_ROOT, "manufacturing", "o8_pzt_coupon_v1")
const OUTPUT_ROOT = joinpath(PACKAGE_ROOT, "visuals")
const MESH_PATH = joinpath(PACKAGE_ROOT, "coupon_validation.msh")
const CARRIER_HZ = 242.0e3
const PRESSURE_SPEED_M_S = 6122.102437409232
const BLOCK_HEIGHT_M = 8PRESSURE_SPEED_M_S / CARRIER_HZ
const BLOCK_LENGTH_M = 2PRESSURE_SPEED_M_S / CARRIER_HZ
const EPSILON_0_F_M = 8.8541878128e-12
const PZT_CAPACITANCE_F = EPSILON_0_F_M * 1900.0 * 20.0e-3 * 8.0e-3 / 1.0e-3
const FREQUENCY_RATIOS = collect(range(0.8, 1.2; length=7))
const FREQUENCIES_HZ = CARRIER_HZ .* FREQUENCY_RATIOS
const STATE_NAMES = ("q00", "q10", "q01", "q11")
const PZT_INTERVALS_MM = (
    (-101.191775824946, -81.191775824946, :darkorange),
    (-30.072127123967, -10.072127123967, :dodgerblue),
    (10.072127123967, 30.072127123967, :dodgerblue),
    (81.191775824946, 101.191775824946, :darkorange),
)

function tuned_load(target_hz, target_q)
    omega = 2pi * target_hz
    inductance_h = inv(omega^2 * 2PZT_CAPACITANCE_F)
    resistance_ohm = omega * inductance_h / target_q
    O8Load(; inductance_h, resistance_ohm)
end

function continuous_field(field)
    coordinate_to_node = Dict{Tuple{Int,Int},Int}()
    old_to_new = Vector{Int}(undef, length(field.x_m))
    x_m = Float64[]
    y_m = Float64[]
    counts = Int[]
    ux = Matrix{ComplexF64}(undef, 0, size(field.ux, 2))
    uy = similar(ux)
    for old in eachindex(field.x_m)
        key = (
            round(Int, field.x_m[old] * 1.0e12),
            round(Int, field.y_m[old] * 1.0e12),
        )
        node = get(coordinate_to_node, key, 0)
        if node == 0
            node = length(x_m) + 1
            coordinate_to_node[key] = node
            push!(x_m, field.x_m[old])
            push!(y_m, field.y_m[old])
            push!(counts, 0)
            ux = vcat(ux, zeros(ComplexF64, 1, size(field.ux, 2)))
            uy = vcat(uy, zeros(ComplexF64, 1, size(field.uy, 2)))
        end
        old_to_new[old] = node
        counts[node] += 1
        ux[node, :] .+= field.ux[old, :]
        uy[node, :] .+= field.uy[old, :]
    end
    for node in eachindex(counts)
        ux[node, :] ./= counts[node]
        uy[node, :] ./= counts[node]
    end
    cells = [[old_to_new[node] for node in cell] for cell in field.cells]
    (; port=field.port, x_m, y_m, cells, ux, uy)
end

function solve_fields()
    pzt = O8PZTBoundary(
        thickness_m=1.0e-3,
        density_kg_m3=7800.0,
        c33e_pa=117.0e9,
        mechanical_quality_factor=80.0,
        coupling_squared=0.49^2,
        capacitance_per_patch_f=PZT_CAPACITANCE_F,
    )
    load_a = tuned_load(258.989e3, 6.382)
    load_b = tuned_load(253.411e3, 4.232)
    results = NamedTuple[]
    for (index, frequency_hz) in enumerate(FREQUENCIES_HZ)
        println("[>] animation field $index/$(length(FREQUENCIES_HZ)): $(frequency_hz / 1e3) kHz")
        config = O8PlaneStrainConfig(
            ; frequency_hz, aperture_height_m=BLOCK_HEIGHT_M,
            block_length_m=BLOCK_LENGTH_M,
        )
        raw = solve_o8_frequency(
            MESH_PATH; mechanism=:pzt, config, pzt, load_a, load_b,
            visualization_port=1,
        )
        push!(results, merge(raw, (field_data=continuous_field(raw.field_data),)))
    end
    results
end

function mesh_lines(field)
    horizontal = Float64[]
    vertical = Float64[]
    for cell in field.cells
        for node in (cell[1], cell[2], cell[3], cell[1])
            push!(horizontal, 1e3field.y_m[node])
            push!(vertical, 1e3field.x_m[node])
        end
        push!(horizontal, NaN)
        push!(vertical, NaN)
    end
    horizontal, vertical
end

function add_pzt_marks!(figure)
    for (lower, upper, color) in PZT_INTERVALS_MM
        plot!(figure, [lower, upper], [0.0, 0.0]; color, linewidth=8, label=false)
    end
    figure
end

function plot_mesh(field)
    horizontal, vertical = mesh_lines(field)
    figure = plot(
        horizontal, vertical;
        color=:black, alpha=0.22, linewidth=0.32, label=false,
        xlabel="aperture y / mm", ylabel="propagation x / mm",
        title="O8 coupon: 3,073 nodes / 5,894 triangles (+250 boundary elements)",
        xlims=(-104, 104), ylims=(-2, 52.5),
        aspect_ratio=:equal, size=(1200, 390),
        background_color=:white, foreground_color=:black,
        grid=false, framestyle=:box,
    )
    add_pzt_marks!(figure)
    annotate!(figure, -91, -1.2, text("B", 9, :darkorange))
    annotate!(figure, -20, -1.2, text("A", 9, :dodgerblue))
    annotate!(figure, 20, -1.2, text("A", 9, :dodgerblue))
    annotate!(figure, 91, -1.2, text("B", 9, :darkorange))
    savefig(figure, joinpath(OUTPUT_ROOT, "coupon_mesh.png"))
end

function field_panel(field, state, phase; title=STATE_NAMES[state])
    ux = field.ux[:, state]
    uy = field.uy[:, state]
    common = maximum(sqrt.(abs2.(field.ux) .+ abs2.(field.uy)))
    instantaneous_x = real.(ux .* cis(phase)) ./ common
    instantaneous_y = real.(uy .* cis(phase)) ./ common
    deformation_mm = 1.1
    horizontal = 1e3 .* field.y_m .+ deformation_mm .* instantaneous_y
    vertical = 1e3 .* field.x_m .+ deformation_mm .* instantaneous_x
    figure = scatter(
        horizontal, vertical;
        marker_z=instantaneous_x, color=:balance, clims=(-1, 1),
        markersize=1.75, markerstrokewidth=0, colorbar=false, label=false,
        xlims=(-104, 104), ylims=(-2.5, 53), aspect_ratio=:equal,
        xlabel="y / mm", ylabel="x / mm", title,
        background_color=:white, foreground_color=:black,
        grid=false, framestyle=:box,
    )
    add_pzt_marks!(figure)
    figure
end

function carrier_animation(field)
    animation = @animate for phase in range(0, 2pi; length=41)[1:end-1]
        panels = [field_panel(field, state, phase) for state in 1:4]
        plot(
            panels...; layout=(2, 2), size=(1200, 700),
            plot_title="242 kHz harmonic field — port 1; displacement exaggerated 1.1 mm",
        )
    end
    gif(animation, joinpath(OUTPUT_ROOT, "four_state_carrier.gif"); fps=16)
    mp4(animation, joinpath(OUTPUT_ROOT, "four_state_carrier.mp4"); fps=24)
end

function pulse_field(results, state, time_s)
    weights = Float64[0.0, 0.25, 0.75, 1.0, 0.75, 0.25, 0.0]
    ux = zeros(ComplexF64, length(results[1].field_data.x_m))
    uy = similar(ux)
    for index in eachindex(results)
        phase = cis(-2pi * FREQUENCIES_HZ[index] * time_s)
        ux .+= weights[index] .* results[index].field_data.ux[:, state] .* phase
        uy .+= weights[index] .* results[index].field_data.uy[:, state] .* phase
    end
    ux, uy
end

function pulse_panel(results, state, time_s, common; title)
    field = results[4].field_data
    ux, uy = pulse_field(results, state, time_s)
    instantaneous_x = real.(ux) ./ common
    instantaneous_y = real.(uy) ./ common
    deformation_mm = 0.8
    horizontal = 1e3 .* field.y_m .+ deformation_mm .* instantaneous_y
    vertical = 1e3 .* field.x_m .+ deformation_mm .* instantaneous_x
    figure = scatter(
        horizontal, vertical;
        marker_z=instantaneous_x, color=:balance, clims=(-1, 1),
        markersize=2.0, markerstrokewidth=0, colorbar=false, label=false,
        xlims=(-104, 104), ylims=(-2.5, 53), aspect_ratio=:equal,
        xlabel="y / mm", ylabel="x / mm", title,
        background_color=:white, foreground_color=:black,
        grid=false, framestyle=:box,
    )
    add_pzt_marks!(figure)
    figure
end

function pulse_animation(results)
    probe_times = range(0.0, 62.0e-6; length=65)[1:end-1]
    common = maximum(
        maximum(abs, first(pulse_field(results, state, time_s)))
        for state in (1, 4) for time_s in probe_times
    )
    animation = @animate for time_s in probe_times
        q00 = pulse_panel(results, 1, time_s, common; title="q00 — both shunts open")
        q11 = pulse_panel(results, 4, time_s, common; title="q11 — both shunts closed")
        plot(
            q00, q11; layout=(2, 1), size=(1200, 650),
            plot_title="Seven-anchor wave-packet reconstruction; t=$(round(1e6time_s; digits=1)) us",
        )
    end
    gif(animation, joinpath(OUTPUT_ROOT, "q00_q11_wave_packet.gif"); fps=20)
    mp4(animation, joinpath(OUTPUT_ROOT, "q00_q11_wave_packet.mp4"); fps=30)
end

function js_int_array(values)
    "[" * join(values, ',') * "]"
end

function quantized_field_data(field)
    common = maximum(sqrt.(abs2.(field.ux) .+ abs2.(field.uy)))
    quantize(value) = round(Int, clamp(value / common, -1, 1) * 32767)
    x = round.(Int, 100_000 .* field.x_m)
    y = round.(Int, 100_000 .* field.y_m)
    triangles = reduce(vcat, field.cells)
    uxr = Int[]
    uxi = Int[]
    uyr = Int[]
    uyi = Int[]
    for state in 1:4
        append!(uxr, quantize.(real.(field.ux[:, state])))
        append!(uxi, quantize.(imag.(field.ux[:, state])))
        append!(uyr, quantize.(real.(field.uy[:, state])))
        append!(uyi, quantize.(imag.(field.uy[:, state])))
    end
    (; x, y, triangles, uxr, uxi, uyr, uyi)
end

const HTML_TEMPLATE = raw"""
<div id="o8-field-lab" aria-label="O8 finite-element mesh and harmonic displacement animation">
  <style>
    #o8-field-lab{color:var(--foreground);font:13px/1.35 ui-sans-serif,system-ui,sans-serif;max-width:920px}
    #o8-field-lab .o8-head{display:flex;align-items:center;gap:12px;flex-wrap:wrap;margin:0 0 8px}
    #o8-field-lab .o8-title{font-size:16px;font-weight:650;flex:1 0 100%}
    #o8-field-lab button{color:var(--foreground);background:transparent;border:1px solid var(--border);padding:5px 9px;border-radius:6px;cursor:pointer}
    #o8-field-lab button[aria-pressed="true"]{background:var(--accent);color:var(--accent-foreground);outline:2px solid currentColor;outline-offset:-2px}
    #o8-field-lab label{display:flex;align-items:center;gap:6px;white-space:nowrap}
    #o8-field-lab input[type="range"]{width:128px}
    #o8-field-lab .o8-stage{position:relative;border:1px solid var(--border);background:var(--background)}
    #o8-field-lab canvas{display:block;width:100%;height:auto;aspect-ratio:4/1.35}
    #o8-field-lab .o8-legend{display:flex;align-items:center;gap:8px;margin-top:7px;color:var(--muted-foreground)}
    #o8-field-lab .o8-bar{width:150px;flex:0 0 150px;height:9px;background:linear-gradient(90deg,#275dad,#f6f7f8,#bd3632);border:1px solid var(--border)}
    #o8-field-lab .o8-legend>span:not(.o8-meta){white-space:nowrap}
    #o8-field-lab .o8-meta{margin-left:auto}
    @media(max-width:560px){#o8-field-lab .o8-bar{flex-basis:86px;width:86px}#o8-field-lab .o8-meta{margin-left:0;width:100%}#o8-field-lab input[type="range"]{width:100px}}
  </style>
  <div class="o8-head">
    <span class="o8-title">242 kHz FEM field - source port 1</span>
    <button type="button" data-play aria-pressed="true">Pause</button>
    <button type="button" data-state="0" aria-pressed="true">q00</button>
    <button type="button" data-state="1" aria-pressed="false">q10</button>
    <button type="button" data-state="2" aria-pressed="false">q01</button>
    <button type="button" data-state="3" aria-pressed="false">q11</button>
    <label>phase <input data-phase type="range" min="0" max="628" value="0" step="1"></label>
    <label><input data-mesh type="checkbox" checked> mesh</label>
  </div>
  <div class="o8-stage"><canvas></canvas></div>
  <div class="o8-legend"><span>-ux</span><span class="o8-bar"></span><span>+ux</span><span class="o8-meta">aperture y horizontal / propagation x downward / deformation exaggerated</span></div>
  <script type="application/json" data-o8>__DATA__</script>
  <script>
  (()=>{
    const root=document.getElementById('o8-field-lab');
    const D=JSON.parse(root.querySelector('[data-o8]').textContent);
    const canvas=root.querySelector('canvas'),ctx=canvas.getContext('2d');
    const n=D.x.length, states=['q00','q10','q01','q11'];
    let state=0,phase=0,playing=true,last=0,showMesh=true;
    const xmin=Math.min(...D.x),xmax=Math.max(...D.x),ymin=Math.min(...D.y),ymax=Math.max(...D.y);
    function mix(a,b,t){return Math.round(a+(b-a)*t)}
    function color(v){v=Math.max(-1,Math.min(1,v));if(v<0){const t=v+1;return `rgb(${mix(39,246,t)},${mix(93,247,t)},${mix(173,248,t)})`}const t=v;return `rgb(${mix(246,189,t)},${mix(247,54,t)},${mix(248,50,t)})`}
    function resize(){const r=canvas.getBoundingClientRect(),d=Math.min(devicePixelRatio||1,2);canvas.width=Math.max(1,Math.round(r.width*d));canvas.height=Math.max(1,Math.round(r.width/3.25*d));draw()}
    function draw(){
      const W=canvas.width,H=canvas.height,pad=30*(devicePixelRatio||1),c=Math.cos(phase),s=Math.sin(phase),scale=32767;
      ctx.clearRect(0,0,W,H);ctx.fillStyle=getComputedStyle(root).getPropertyValue('--background')||'#fff';ctx.fillRect(0,0,W,H);
      const px=new Float32Array(n),py=new Float32Array(n),val=new Float32Array(n),off=state*n,def=0.9;
      for(let i=0;i<n;i++){
        const ux=(D.uxr[off+i]*c-D.uxi[off+i]*s)/scale,uy=(D.uyr[off+i]*c-D.uyi[off+i]*s)/scale;
        const xx=D.x[i]/100+def*ux,yy=D.y[i]/100+def*uy;
        px[i]=pad+(yy-ymin/100)/(ymax/100-ymin/100)*(W-2*pad);
        py[i]=pad+(xx-xmin/100)/(xmax/100-xmin/100)*(H-2*pad);val[i]=ux;
      }
      for(let k=0;k<D.tri.length;k+=3){const a=D.tri[k]-1,b=D.tri[k+1]-1,d=D.tri[k+2]-1;ctx.beginPath();ctx.moveTo(px[a],py[a]);ctx.lineTo(px[b],py[b]);ctx.lineTo(px[d],py[d]);ctx.closePath();ctx.fillStyle=color((val[a]+val[b]+val[d])/3);ctx.fill();if(showMesh){ctx.strokeStyle='rgba(20,24,30,.10)';ctx.lineWidth=.45;ctx.stroke()}}
      const patches=[[-101.192,-81.192,'#d87500'],[-30.072,-10.072,'#1777c8'],[10.072,30.072,'#1777c8'],[81.192,101.192,'#d87500']];
      ctx.lineWidth=7*(devicePixelRatio||1);for(const p of patches){ctx.strokeStyle=p[2];ctx.beginPath();ctx.moveTo(pad+(p[0]-ymin/100)/(ymax/100-ymin/100)*(W-2*pad),pad);ctx.lineTo(pad+(p[1]-ymin/100)/(ymax/100-ymin/100)*(W-2*pad),pad);ctx.stroke()}
      ctx.fillStyle=getComputedStyle(root).color;ctx.font=`${12*(devicePixelRatio||1)}px system-ui`;ctx.fillText('source / PZT boundary',pad,18*(devicePixelRatio||1));ctx.fillText('radiation face',pad,H-7*(devicePixelRatio||1));
    }
    function tick(t){if(playing){if(last)phase=(phase+(t-last)*0.0032)%(Math.PI*2);last=t;root.querySelector('[data-phase]').value=Math.round(phase*100);draw()}requestAnimationFrame(tick)}
    root.querySelector('[data-play]').addEventListener('click',e=>{playing=!playing;e.currentTarget.textContent=playing?'Pause':'Play';e.currentTarget.setAttribute('aria-pressed',String(playing));last=0});
    root.querySelectorAll('[data-state]').forEach(b=>b.addEventListener('click',()=>{state=+b.dataset.state;root.querySelectorAll('[data-state]').forEach(x=>x.setAttribute('aria-pressed',String(x===b)));draw()}));
    root.querySelector('[data-phase]').addEventListener('input',e=>{phase=+e.target.value/100;draw()});
    root.querySelector('[data-mesh]').addEventListener('change',e=>{showMesh=e.target.checked;draw()});
    new ResizeObserver(resize).observe(root.querySelector('.o8-stage'));requestAnimationFrame(tick);
  })();
  </script>
</div>
"""

function write_interactive(path, field)
    data = quantized_field_data(field)
    payload = "{" * join([
        "\"x\":" * js_int_array(data.x),
        "\"y\":" * js_int_array(data.y),
        "\"tri\":" * js_int_array(data.triangles),
        "\"uxr\":" * js_int_array(data.uxr),
        "\"uxi\":" * js_int_array(data.uxi),
        "\"uyr\":" * js_int_array(data.uyr),
        "\"uyi\":" * js_int_array(data.uyi),
    ], ',') * "}"
    html = replace(HTML_TEMPLATE, "__DATA__" => payload)
    open(path, "w") do io
        write(io, html)
    end
    path
end

function export_o8_mesh_animations(; visualization_path=nothing)
    mkpath(OUTPUT_ROOT)
    results = solve_fields()
    field = results[4].field_data
    plot_mesh(field)
    carrier_animation(field)
    pulse_animation(results)
    jldsave(joinpath(OUTPUT_ROOT, "field_bundle.jld2"); frequencies_hz=FREQUENCIES_HZ, results)
    if !isnothing(visualization_path)
        mkpath(dirname(visualization_path))
        write_interactive(visualization_path, field)
    end
    (; output_root=OUTPUT_ROOT, visualization_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    path = isempty(ARGS) ? nothing : abspath(first(ARGS))
    println(export_o8_mesh_animations(; visualization_path=path))
end

end
