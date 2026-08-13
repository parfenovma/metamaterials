# Hybrid TTD--GRIN robustness sweep

## Result

The refined coupled monostatic G3 advantage of the nominal hybrid over the
pure smooth TTD control is only `+0.9696 dB`.  The sweep therefore defines a
strict coupon qualification envelope rather than demonstrating a large
engineering margin.

Independent one-parameter break-even limits are:

| added imperfection | last passing value |
|---|---:|
| one-way excess insertion loss | `0.480 dB` |
| edge-to-centre profile delay error, either sign | `610 ns` |
| differential phase ripple, worst of four ripple phases | `26.5 deg` |
| absolute unwanted-mode power fraction | `13.70%` |
| extra ringing penalty | `17.55 dB` |

These limits are **not simultaneous allowances**.  Each independently spends
the same `0.9696 dB` hybrid-to-pure margin.  A proposed combined coupon gate is
therefore much tighter:

```text
one-way excess insertion loss       <= 0.10 dB
edge-to-centre delay error          <= +/-100 ns
differential phase ripple           <= 5 deg
absolute mode conversion            <= 5.0%
extra echo-window ringing penalty   <= 6 dB
```

The worst sign/ripple phase at this combined point retains `+0.534 dB` over
pure TTD, `+11.12 dB` echo-to-ringing, `B_t=1.017`, and pulse correlation
`0.9963`.  A relaxed screen (`0.15 dB`, `100 ns`, `10 deg`, `6%`, `6 dB`)
retains only `+0.243 dB`, so it is not a comfortable manufacturing target.

The most decisive number is the insertion-loss limit.  Because the insert is
traversed on transmit and receive, `0.48 dB` of additional one-way loss erases
approximately `0.96 dB` of echo advantage.  If a real coupon exceeds this
value, no refinement of the present nominal topology can recover its advantage
over pure TTD without changing the delay/power design.

## Model level

The starting waveform is not a scalar focus proxy.  It is the stored refined
2D plane-strain paired result

```text
Delta s(t) = s_void(t) - s_none(t)
```

recorded on the same transmitting aperture by the explicit coupled
`hybrid + Al specimen + r=2 mm void` G3 model.  The additional imperfections
are then applied as a reduced complex two-way transfer:

- insertion loss is specified per one-way traversal;
- delay error is a linear radial error from aperture edge to centre;
- phase ripple is differential across the aperture and is screened at four
  phase offsets;
- mode conversion removes useful-mode power relative to the already modelled
  nominal `3.602%` conversion;
- the original refined device ringing is retained rather than attenuated with
  the echo, and a separate positive penalty scales it conservatively.

The 15 channels use their actual radial delay fractions but equal aperture
weights.  The reduced transfer captures loss of coherence, two-way loss and
pulse distortion without pretending to predict a particular imperfect bond
geometry.  It does not replace measurement of complex coupon `S21/S11`, and it
does not alter the G3 dead-zone trace because the defect-free baseline is held
fixed.  A real joint that creates a new delayed reflection must be represented
by the measured ringing penalty or a new full FEM geometry.

## Actual refined mesh

The mesh renderer uses the exact
`mesh_refined_hybrid_void_r2_focus.msh` employed by G3:

| quantity | value |
|---|---:|
| nodes | `49,197` |
| triangles | `87,800` |
| specimen target mesh | `1.00 mm` |
| guide target mesh | `0.32 mm` |
| void target mesh | `0.22 mm` |
| displacement DOF | `98,394` |

The generated four-panel image shows the full stand, all 15 routed channels,
the graded Al--Cu blocks and the locally refined free void boundary.  It is
written to:

```text
tmp/topology_benchmark/hybrid_robustness_sweep/hybrid_refined_g3_mesh.png
```

## Reproduction

From this worktree root:

```bash
julia --startup-file=no --project=. \
  test/topology_hybrid_robustness_sweep.jl

julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_robustness_sweep.jl

julia --startup-file=no --project=. \
  src/code/metamaterial/plot_topology_hybrid_monostatic_mesh.jl
```

The sweep writes dense axis and two-parameter CSVs plus the robustness map to
`tmp/topology_benchmark/hybrid_robustness_sweep/`.  The compact qualification
tables are retained with this report.  The exact G3 source signals and meshes
remain under `tmp/topology_benchmark/hybrid_monostatic_g3/`.

## Decision

The branch remains `conditional`.  The numerical sweep makes the next
experimental decision unambiguous:

- continue to a multi-channel sector only if a representative Al--Cu--Al
  coupon passes the combined qualification point across `162.2--321.7 kHz`;
- treat `0.48 dB` one-way excess loss as a hard stop relative to pure TTD;
- do not spend compute on a composition sweep before complex coupon transfer
  and echo-window ringing are measured.
