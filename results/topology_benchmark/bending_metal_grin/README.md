# Bending-dominated all-metal GRIN cell

## Verdict

`stop` at G1 for the tested single-phase Al6061 corner-mass/re-entrant-ligament
connectivity at `242 kHz`, a `55 mm` device length and the shared `5-cycle
Hann` pulse.

The topology successfully separates two geometric controls: ligament
thickness/re-entrant offset primarily changes axial bending stiffness, while
the retained centre/corner blocks change inertia.  Static full-elastic solves
span stiffness ratios `0.188--0.707` and relative areal densities
`0.705--0.912`; this is not the periodically waisted-ring mechanism rejected
in the earlier metal-GRIN branch.

It nevertheless does not produce the required library of five broadband
states.  The three dense states recover useful transmission but fail the
pulse gate.  The best one, `dense_c`, has scalar-port fine-grid metrics

```text
weighted T                 0.93542
weighted R                 0.03556
carrier excess phase      -2.55657 rad
group delay at 242 kHz     3.41789 us
B_t                        1.44503       (required <= 1.25)
pulse correlation          0.94762
postcursor/main             0.08790
```

Its carrier phase is large enough to cover the requested `2.451 rad`, but it
is not a true-time-delay state: trusted group delay rises from `1.47 us` at
`162 kHz` to `12.02 us` at `282 kHz` and `25.44 us` at `294 kHz`.  A
multimode DtN check closes power to `1e-5` and gives `T/R=0.9725/0.0275` at
`242 kHz`, `0.3105/0.6895` at `302 kHz`, and approximately
`0.000043/0.999958` at `322 kHz`.  Thus the upper part of the required
`162.2--321.7 kHz` band crosses a structural stop-band edge.

The stronger bending states have coarse weighted transmission only
`0.003--0.283` and pulse broadening from `2.06` to hundreds.  No non-solid
state passes all local gates, versus five required states and at least
`1.5 us` combined admissible delay coverage.  G2 aperture, G3 paired defect
and G4 robustness calculations were therefore not run; their common metrics
remain `NaN` rather than being replaced by a focus or echo proxy.

This verdict closes this specific corner-mass/re-entrant plate connectivity
at the shared dimensions.  It does not close all metallic pentamode lattices,
dissimilar-metal inertia loading or a topology with a substantially smaller
cell pitch.

## Reproduction

The numerical source used for the reported calculations is commit
`af18508170f19683ff21f788c92c4ef6dff83564` on
`codex/topology-bending-metal-grin`.

```bash
julia --startup-file=no --project=. test/topology_bending_metal_grin.jl

julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=g0
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=mesh
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=convert
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=coarse
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=analyze-coarse
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=fine
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=analyze-fine

julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_static.jl

julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=chain-mesh
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=chain-convert
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_run.jl --stage=chain-check

julia --startup-file=no --project=. \
  src/code/metamaterial/topology_bending_metal_grin_multimode_check.jl

for stage in mesh convert solve; do
  julia --startup-file=no --project=. \
    src/code/metamaterial/topology_bending_metal_grin_convergence.jl --stage=$stage
done
```

Raw meshes, Gridap models, per-frequency JLD2 files and diagnostic tables are
under `tmp/topology_benchmark/bending_metal_grin/`.  The permanent common
metric table is `benchmark_metrics.csv` beside this report.  The fine sweep is
`162:4:322 kHz`; coarse screening is `162:20:322 kHz`.

The standalone topology suite passed `13/13`.  The full repository suite also
passed after copying its pre-existing ignored `gap_pilot_220khz` and
`period_length_pilot_220khz` fixtures from the main worktree into this
worktree's `tmp/`; no fixture is part of this branch.

## Research basis and topology choice

The design was bounded by primary sources rather than importing published
effective properties:

1. [Broadband Lamb wave trapping in cellular metamaterial plates with
   multiple local resonances](https://doi.org/10.1038/srep09376) distinguishes
   bending- and stretch-dominated cellular mechanisms and motivates using
   ligament flexure as a stiffness degree of freedom.
2. [Extreme stiffness hyperbolic elastic metamaterial for total transmission
   subwavelength imaging](https://doi.org/10.1038/srep24026) demonstrates that
   void-only aluminium cells normally reduce density and stiffness together,
   and that dynamic stiffness compensation is required for high transmission.
   Its resonant imaging cell is not treated as a broadband NDT material here.
3. [Broadband focusing of underwater sound using a transparent pentamode
   lens](https://doi.org/10.1121/1.4985195) demonstrates an aluminium plate-stack
   manufacturing route, but its water-matched constitutive target is not
   transferred to aluminium-to-aluminium bulk-wave coupling.

The primary topology is a periodic `5 mm` Al6061 cell with full-height
boundary nodes, one retained centre block, two staggered re-entrant ligament
paths per half-cell, and corner blocks at the elbows.  Eleven states were
screened: solid, three dense states and seven stronger bending states.  The
upper-band ratio is `a/lambda_P=0.263`, so this is not assumed to be deeply
homogenized; finite-chain and transfer-eigenphase checks are mandatory.

## Geometry and manufacturing route

| Quantity | Value |
|---|---:|
| Total device length | `55.0 mm` |
| Active lattice length | `45.0 mm` |
| Period / cells | `5.0 mm / 9` |
| Strip/radial-zone width | `7.0 mm` |
| Ligament library | `0.8--2.7 mm` |
| Boundary-node minimum | `0.9 mm` |
| Retained centre-block height | `2.8--5.4 mm` |
| Relative areal density | `0.705--0.912` |
| Nominal plate thickness | `1.5 mm` |
| Named process | through-cut micro wire EDM |

Prototype route: through-cut the repeated void pattern in Al6061 plates by
wire EDM, align plates on external tooling features, and clamp or diffusion
bond only outside the active strip before lapping the two end faces.  A full
aperture would require separate annular/polygonal radial cartridges or a 3-D
redesign; the present calculation does not claim that an axisymmetric stack
can be obtained by simply revolving the 2-D pattern.  Any bond that constrains
the active free surfaces invalidates this G1 model and must be included in a
new 3-D gate.

## G0: independent stiffness and inertia controls

Static quadratic plane-strain FEM gives:

| State | ligament, mm | relative density | static `K/K_Al` | static `c/c_Al` | static `Z/Z_Al` |
|---|---:|---:|---:|---:|---:|
| `dense_a` | 2.00 | 0.807 | 0.488 | 0.777 | 0.627 |
| `dense_b` | 2.35 | 0.856 | 0.562 | 0.810 | 0.693 |
| `dense_c` | 2.70 | 0.912 | 0.707 | 0.880 | 0.803 |
| `bd_a` | 1.70 | 0.783 | 0.229 | 0.540 | 0.423 |
| `bd_g` | 0.80 | 0.765 | 0.188 | 0.496 | 0.380 |

These are quasi-static descriptors, not assigned effective materials.  They
confirm that the family changes stiffness much more strongly than inertia,
but also expose the impedance collapse of the high-delay states.

## G1: full-elastic modal-port results

The finite cell uses quadratic displacement elements, lossless Al6061 and
physical free surfaces.  Fine meshes use `0.16 mm` near ligaments and
`0.34 mm` in bulk regions.  The incident field is the power-normalized
fundamental symmetric quasi-longitudinal strip mode.

The full fine summary for the three dense states is:

| State | weighted T | weighted R | delay at 242 kHz, us | `B_t` | correlation | postcursor | local gate |
|---|---:|---:|---:|---:|---:|---:|---|
| `dense_a` | 0.769 | 0.161 | 8.279 | 1.801 | 0.911 | 0.0967 | fail T/Bt |
| `dense_b` | 0.814 | 0.104 | 7.083 | 2.105 | 0.883 | 0.0992 | fail Bt/correlation |
| `dense_c` | 0.935 | 0.0356 | 3.418 | 1.445 | 0.948 | 0.0879 | fail Bt |

Mirror symmetry suppresses antisymmetric propagating conversion in this 2-D
model.  A multimode DtN with all propagating strip modes and eight slow
evanescent modes was used to audit absolute power at seven anchors.  The
carrier result changes by less than `1.2e-7` in T when the evanescent basis is
increased from 8 to 12 modes.  The scalar fine sweep is used for dense phase
and pulse sampling; absolute power conclusions use the multimode anchors,
because a rank-one termination absorbs unresolved evanescent content near the
stop band.

### Transfer eigenphase and finite-chain check

The one-cell reciprocal S matrix was converted to a transfer-matrix
eigenphase.  Independent `N=1,3,6` finite chains give per-cell excess-phase
spreads of only `0.025--0.045 rad` from `182` through `282 kHz`, confirming a
propagating branch there.  At `302 kHz`, only the one-cell chain retains
`T>=0.5`; at `322 kHz` none does, and the transfer eigenphase is complex.
This is consistent with the stop-band seen by the 9-cell device.

This is a transfer-matrix Bloch check, not a separate periodic-boundary FE
eigensolve.  It is sufficient to diagnose the finite-chain branch used in the
verdict, but a promoted topology would still require a full 3-D Floquet model.

### Mesh convergence

At `242 kHz`, multimode `dense_c` results are:

| near-ligament / bulk size, mm | T | R |
|---|---:|---:|
| `0.20 / 0.40` | 0.973130 | 0.026907 |
| `0.16 / 0.34` | 0.972510 | 0.027481 |
| `0.11 / 0.24` | 0.971735 | 0.028271 |

The fine-to-baseline change is `7.75e-4` in T, small compared with every gate
margin relevant to the STOP decision.

## G2--G4 and common metrics

The required five-state G1 library does not exist, so an aperture would have
to insert an ideal phase/delay or use failed low-transmission states.  G2 and
the paired `none / void_r2_focus` monostatic G3 were deliberately skipped.
Consequently no CNR, defect echo gain, ringing metric, dead-zone depth, focus
gain, FWHM or POD is claimed.

## Limits of interpretation

- This is a 2-D plane-strain strip/supercell, not a 3-D annulus or complete
  aperture.  Hoop stiffness, plate interfaces and circumferential modes are
  absent.
- The structure is lossless.  Real plate bonds and surface damage add loss;
  damping may soften resonant peaks but cannot restore energy in a stop band.
- The excitation is perfectly symmetric.  Manufacturing asymmetry can couple
  energy into odd/flexural and out-of-plane modes.
- The high-frequency edge is not deeply subwavelength.  Reported delay is a
  dispersive finite-chain property, not a frequency-independent effective
  speed.
- Static stiffness and areal density demonstrate two geometric controls, but
  do not prove independent arbitrary tuning of dynamic effective tensors.
- `dense_c` is worth retaining only as evidence about this trade-off, not as a
  candidate lens cell: it fails the frozen pulse gate before aperture or
  defect evaluation.
