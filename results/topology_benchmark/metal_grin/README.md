# Metal GRIN: periodically waisted aluminium rings

## Verdict

`stop` for this connectivity at `242 kHz`, a `55 mm` device depth and the
shared `5-cycle Hann` pulse.

The manufacturable state that passes the local transfer and pulse gates
(`neck_2p0`) has weighted transmission `0.9288`, weighted reflection `0.0712`,
`B_t=1.000`, and pulse correlation `0.9993`, but supplies only `0.399 us` of
excess group delay and `0.561 rad` of carrier-phase range.  The shared stand
requires `1.612 us` or `2.451 rad` between aperture edge and centre.  Deeper
waists obtain more carrier phase only by entering a strongly dispersive,
reflecting branch.  Consequently no admissible radial state library exists
from which to assemble the aperture, so G2 and G3 were not run and their
metrics remain `NaN` rather than being replaced by a focus or one-way proxy.

This result closes only the periodically waisted-ring connectivity.  It does
not close all-metal GRIN media, pentamode lattices, or dissimilar-metal cells.

## Reproduction

The numerical source used for the reported calculations is commit
`3768c71a1e26c2079c303732ab26b20ebbd43f98` on
`codex/topology-metal-grin`.

```bash
julia --startup-file=no --project=. test/topology_metal_grin.jl

julia --startup-file=no --project=. \
  src/code/metamaterial/topology_metal_grin_run.jl --stage=g0
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_metal_grin_run.jl --stage=mesh
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_metal_grin_run.jl --stage=convert
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_metal_grin_run.jl --stage=coarse
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_metal_grin_run.jl --stage=analyze-coarse
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_metal_grin_run.jl --stage=fine
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_metal_grin_run.jl --stage=analyze-fine
```

Raw meshes, Gridap models and per-frequency JLD2 files are intentionally kept
under `tmp/topology_benchmark/metal_grin/`.  The principal generated tables
are:

- `g0_manufacturing_and_surrogate.csv`;
- `g1_coarse_frequency.csv` and `g1_coarse_summary.csv`;
- `g1_fine_frequency.csv` and `g1_fine_summary.csv`;
- `benchmark_metrics.csv`, whose columns exactly match
  `NDTBenchmark.metric_columns()`.

The standalone topology tests passed `17/17`.  The full repository suite also
passed after restoring its pre-existing ignored `gap_pilot` and
`period_length_pilot` fixtures into this worktree's `tmp/` directory.

## Research basis and topology choice

Three primary sources bound the design rather than supplying assumed material
properties:

1. [Longitudinal wave propagation in a rod with variable cross-section](https://doi.org/10.1016/j.jsv.2013.09.010)
   derives transfer models for variable-area rods and shows why Poisson and
   shear effects can widen stop bands.  This motivates using the area model
   only as G0 and requiring full elastic FEM for G1.
2. [High-transmission acoustic self-focusing and directional cloaking in a graded perforated metal slab](https://doi.org/10.1038/s41598-017-04512-z)
   retrieves index and group velocity from band structure and adjusts density
   and group velocity together.  Its water-loaded honeycomb cell is not
   transferred directly to metal-to-metal NDT.
3. [Broadband focusing of underwater sound using a transparent pentamode lens](https://doi.org/10.1121/1.4985195)
   demonstrates a real aluminium plate-stack route and explicitly uses FEM
   homogenization plus Bloch analysis.  It remains a water-impedance-matched
   device with centimetre cells, not evidence that the same cell is matched to
   solid aluminium at `242 kHz`.

The selected main topology is therefore a simpler non-resonant screen: each
radial zone is a separate Al6061 annulus with circumferential grooves machined
from its inner and outer sides.  The aligned grooves leave a smooth axial
waist.  Waist height varies between rings, producing different effective
phase velocity while retaining a direct metal load path.  This is distinct
from the previously rejected equal-area aligned through-hole array.

The pentamode plate stack was retained only as a backup concept.  It was not
simulated here: scaling its published waterborne cells to the present
frequency changes both the impedance target and minimum ligament, so treating
its published effective properties as aluminium-to-aluminium properties would
be unjustified.

## G0: geometry and manufacturing route

| Parameter | Value |
|---|---:|
| Material | Al6061 surrogate |
| Total axial length | `55.0 mm` |
| Straight lead at each end | `5.5 mm` |
| Active length | `44.0 mm` |
| Axial period | `4.0 mm` |
| Period count | `11` |
| Local radial zone width | `3.2 mm` |
| Cosine transition per side | `1.0 mm` |
| Constant-waist plateau | `2.0 mm` |
| Neck library | `3.2, 2.8, 2.4, 2.0, 1.6, 1.2, 0.8 mm` |
| Nominal minimum feature | `0.8 mm` |
| Minimum geometric curvature radius | `0.169 mm` at the `0.8 mm` state |
| Cell pitch / P wavelength | `0.106` at `162 kHz`, `0.210` at `322 kHz` |

Proposed prototype route:

1. Turn each concentric annulus separately from Al6061 tube or billet.
2. Cut the opposed circumferential groove profiles by CNC turning or wire EDM;
   the `0.8 mm` remaining web is the limiting feature.
3. Finish both end faces in one lapping setup.
4. Nest the rings with non-load-bearing radial clearance and locate them only
   at removable input/output collars.
5. Couple the coplanar ring ends to the source and specimen.  Do not shrink-fit
   the full ring length: that would introduce radial coupling absent from G1.

The expected first prototype is BOM class `high`: approximately ten annular
zones plus a central core, two locating collars and couplant.  It requires no
dissimilar acoustic material, adhesive-filled internal interface, PBF powder
removal, or diffusion-bonded stack.  The `straight_footprint` control is the
same ring set with no waists.

The G0 area-homogenized model used

```math
\frac{c_{eff}}{c_{Al}} =
\frac{1}{\sqrt{\langle A/A_0\rangle\langle A_0/A\rangle}},
\qquad
\frac{Z_{eff}}{Z_{Al}} =
\sqrt{\frac{\langle A/A_0\rangle}{\langle A_0/A\rangle}}.
```

It is a low-frequency one-dimensional bound, not an assigned effective
material.  Its most extreme printable state predicted only `1.047 us` delay,
already below the required `1.612 us`; G1 was nevertheless run because the
upper band reaches `p/lambda=0.210` and dynamic effects could be material.

## G1: full-elastic port FEM

The G1 model is 2-D plane strain with quadratic displacement elements,
approximately `0.10 mm` mesh size near the waists, lossless Al6061, and a
power-normalized quasi-longitudinal input.  The `3.2 mm` straight ports have
one propagating symmetric mode throughout `162--322 kHz`.  Mirror symmetry
forbids odd conversion.  `S11` and `S21` are obtained by reciprocity projection
on both ports; the lossless power residual is below `1.3e-4` for the reported
fine states.

The fine sweep is `162:4:322 kHz`.  Weighted quantities use the shared
five-cycle Hann spectrum over its `-20 dB` band.

| State | min feature, mm | min `T00` | weighted `T` | weighted `R` | relative phase at 242 kHz, rad | phase coverage | group delay at 242 kHz, us | `B_t` | correlation | postcursor | G1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| solid | 3.2 | 0.99997 | 0.99998 | 0.00000 | 0.000 | 0.000 | -0.000 | 1.000 | 1.0000 | 0.0666 | reference |
| neck 2.0 | 2.0 | 0.85452 | 0.92877 | 0.07120 | -0.561 | 0.229 | 0.399 | 1.000 | 0.9993 | 0.0641 | pass locally |
| neck 1.6 | 1.6 | 0.70946 | 0.85193 | 0.14800 | -1.040 | 0.424 | -0.307 | 1.016 | 0.9969 | 0.0681 | fail delay sign/coverage |
| neck 1.2 | 1.2 | 0.51674 | 0.73478 | 0.26528 | -1.648 | 0.672 | -1.010 | 1.063 | 0.9886 | 0.0781 | fail transfer/delay |

The coarse screen also evaluated the limiting `0.8 mm` state.  It reaches
relative carrier phase `-3.209 rad`, enough to span the target phase in a
single-frequency sense, but weighted transmission falls to `0.615`, weighted
reflection rises to `0.385`, and the group delay is untrusted because its
finite-difference stencil crosses `T00<0.5`.  This is the rejected dynamic
branch, not a valid slow effective medium.

The complex carrier values for the selected `neck_2p0` state are

```text
S11 = -0.154969 - 0.185866im
S21 = -0.745505 + 0.620995im
|S11|^2 = 0.058562
|S21|^2 = 0.941413
```

`S21` here includes propagation through the complete `55 mm` cell.  Carrier
phase in the table is calibrated against the solid cell of identical length.

## G2--G4 and shared NDT metrics

The aperture gate was not entered.  The best locally admissible state covers
only `22.9%` of the required centre-to-edge phase, and the deepest state that
still passes weighted transfer covers `42.4%` while having anomalous negative
group delay.  Therefore an aperture assembled from the passing states cannot
realize the prescribed focus without silently inserting an ideal delay.

| Topology / reference | `echo_mf_gain_bare` | `echo_energy_gain_bare` | `echo_to_ringing_db` | focus gain | weighted `T` | weighted `R` | `B_t` | correlation | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `bare` | `NaN` | `NaN` | `NaN` | `NaN` | `NaN` | `NaN` | `NaN` | `NaN` | not calculated after G1 stop |
| `straight_footprint` | `NaN` | `NaN` | `NaN` | `NaN` | 0.99998 | 0.00000 | 1.000 | 1.0000 | local reference only |
| `metal_grin_neck_2p0` | `NaN` | `NaN` | `NaN` | `NaN` | 0.92877 | 0.07120 | 1.000 | 0.9993 | `stop` |

Accordingly, the required `none` / `void_r2_focus` monostatic pair was not
run.  No one-way harmonic field, equal-noise assumption, or area-model delay
is reported as defect echo, CNR, POD, dead-zone depth, or ringing metric.

## Limits of interpretation

- G1 is a 2-D plane-strain strip surrogate, not an axisymmetric or full 3-D
  ring calculation.  It omits hoop stiffness, circumferential modes,
  inter-ring coupling, common faceplates, radial clearance and couplant.
- Free strip sides approximate isolated inner/outer cylindrical surfaces only
  qualitatively.  Their curvature is absent.
- The calculation is lossless.  This is conservative for transmission but can
  underpredict how real joints suppress narrow resonances; loss cannot recover
  the missing phase range without reducing echo energy.
- The quasi-static area model is a bounded G0 surrogate.  All verdict numbers
  come from absolute modal-port FEM, not from fitted effective properties.
- A future revival needs a topology with an additional independent mass
  degree of freedom, so speed can decrease without the impedance collapse
  seen here.  A properly re-scaled pentamode or multi-metal mass/ligament cell
  would be a new topology and must restart at Bloch/port G1.
