# Symmetric freeform void/solid oracle near 315 kHz

## Verdict

This freeform branch is a **categorical no-go** for the tested curved-slot
grammar.  None of 24 coarse candidate/frequency cases passes the prototype or
strict gate.  Fine-mesh validation and geometric simplification confirm the
failure rather than rescuing it.

The least-bad single-frequency point is the segmented design `f008@305 kHz`:

| geometry | mesh | J1 | J2 equal | J2 dB | exit energy | focus error | FWHM at x=115 | prototype |
|---|---|---:|---:|---:|---:|---:|---:|:---:|
| `f008` | coarse | 0.1681 | 0.03987 | -27.99 | 0.2314 | 11 mm | 30.66 mm | no |
| `f008` | fine | 0.1815 | 0.04649 | -26.65 | 0.2315 | 11 mm | 31.21 mm | no |
| `f008_simplified` | fine | 0.1824 | 0.04695 | -26.57 | 0.2318 | 11 mm | 31.26 mm | no |

The fine mesh changes the result modestly but leaves the physical verdict
unchanged.  The simplified five-control-point arc segments reproduce the
smooth seven-control-point version almost exactly; simplification is therefore
validated only as a faithful representation of a failed topology.

## Search and robust result

The campaign contains eight geometries at 305, 315 and 321.7 kHz: 24 coarse
candidate FEM solves plus matched straight references.

- `f001`--`f006`: 3--5 continuous tapered B-spline slots spanning most of the
  aperture.  Variables are component count, end and centre width, curvature,
  transverse extent and width exponent.
- `f007`--`f008`: 15 or 20 short symmetric curved segments.  These were added
  after the continuous channels were identified as free-surface barriers, and
  use the remaining six cases in the 24-case cap.

Continuous channels are the worst family: the best transmission among them is
only 0.106 and most values are 0.01--0.08.  Segmenting the slots opens direct
longitudinal paths and improves the best transmission to 0.231, but this is
still far below the prototype threshold of 0.70.

The required robust metric is

`robust_joint = minimum_f(min(J1(f),J2_equal(f)))`.

The largest sampled value is only 0.02539 for `f001`; that design has a worst
exit-energy proxy of 0.00931 and a worst focus error of 22 mm.  `f008`, the
least-bad gate-oriented design, has `robust_joint = 0.00370`.  Robust ranking
cannot override explicit transmission and focus constraints.

There are:

- zero prototype-feasible coarse rows;
- zero strict-feasible coarse rows;
- zero feasible fine rows;
- no feasible Pareto set and no fabrication candidate.

## Physical interpretation

The result is not evidence that arbitrary topology optimization is useless.
It rejects this particular low-dimensional grammar.

A through-cut slot is a traction-free boundary.  Long continuous slots
intercept most direct paths and behave as barriers, causing reflection and mode
conversion rather than a weak effective-index perturbation.  Breaking them
into short curved segments helps transmission, but the 15--20 components still
form several serial scattering screens.  The field comparison shows that the
amplitude has already collapsed before the nominal target; the poor CNR is not
merely a misplaced but otherwise strong focus.

This also explains why the earlier isolated-ellipse/matched-GRIN families were
substantially better: they preserve many continuous aluminium load paths.  A
future true pixel/level-set oracle should strongly constrain longitudinal
connectivity or optimize material additions/impedance variations, rather than
starting from aperture-spanning void channels.

## Geometry and regularization

Every design is mirror-symmetric about `y=0`, binary aluminium/void and a
monolithic through-cut pattern.  The aluminium surface remains connected and
no void touches the device boundary.

Across the eight designs:

- void fraction: 6.46--11.77%;
- void components: 3--20;
- reported boundary perimeter: 309.7--610.9 mm;
- smallest feature: at least 1.0 mm for the primary DOE;
- smallest conservative solid web: 0.920 mm;
- prescribed feature/web requirement: 0.8 mm.

`f008` has a 6.456% void fraction, 20 components, 353.1 mm perimeter, 1.192 mm
minimum realised slot width and 0.920 mm conservative web.  Conservative
bounding-circle separation is used between diagonal segmented components.

The simplification keeps the same 20 arc components and quantized dimensions,
but reduces each segment boundary from seven to five spline control points.

## Benchmark contract

Every metric is normalized against the solid straight footprint at the same
frequency, mesh resolution and energy-normalized uniform traction.  `J1`, the
reciprocal Born `J2`, clutter mask, exit-energy proxy, nominal-plane FWHM and
fixed-window focus error are identical to the common campaign contract.

The prototype gate is `exit_energy_proxy >= 0.70` and
`focus_error_mm <= 5`; the strict gate is `exit_energy_proxy >= 0.75` and
`focus_error_mm <= 3`.  Pareto dominance is computed on the raw `(J1,J2)` pair,
never on a hidden weighted score.

The fine `f008` meshes contain about 61k nodes and 120k 2D elements.  These are
linear steady-state carrier solves, not transient echo or ringing predictions.

## Artifacts

- `data/campaign_metrics.csv`: 24 coarse rows plus the two fine validation
  rows, using the required campaign-column prefix.
- `data/doe_manifest.csv`: grammar parameters, void fraction, perimeter,
  component count and manufacturing checks.
- `data/robust_summary.csv`: three-frequency robust metric and worst-case
  constraints.
- `data/fine_simplification_validation.csv`: coarse/fine/simplified comparison.
- `data/geometry_f008.csv` and `data/geometry_f008_simplified.csv`: sampled CAD
  boundaries.
- `figures/frequency_case_heatmap.png`: full 8-by-3 metric maps.
- `figures/J1_J2_pareto.png`: raw no-go Pareto plane.
- `figures/geometry_f008_fine.png`: least-bad segmented geometry.
- `figures/field_f008_f0305p0_fine.png`: fine field against the straight
  reference.
- `figures/geometry_f008_simplified_fine.png` and
  `figures/field_f008_simplified_f0305p0_fine_simplified.png`: simplified
  validation.
- `freeform_oracle_pipeline.jl`: reproducible geometry, meshing, FEM, metrics,
  plots and CLI.

## Reproduction

From the repository root:

```sh
julia --startup-file=no --project=. \
  results/freeform_oracle_p6_315khz/freeform_oracle_pipeline.jl --manifest
julia --startup-file=no --project=. \
  results/freeform_oracle_p6_315khz/freeform_oracle_pipeline.jl --sweep
julia --startup-file=no --project=. \
  results/freeform_oracle_p6_315khz/freeform_oracle_pipeline.jl --case f008 305 fine
julia --startup-file=no --project=. \
  results/freeform_oracle_p6_315khz/freeform_oracle_pipeline.jl \
  --case f008_simplified 305 fine
julia --startup-file=no --project=. \
  results/freeform_oracle_p6_315khz/freeform_oracle_pipeline.jl \
  --validate f008 305
```

The useful oracle output is the negative design rule: do not spend further FEM
budget on aperture-spanning or serial curved void screens.  Preserve continuous
longitudinal aluminium paths and use isolated, weak perturbations or a
connectivity-constrained material distribution in the next freeform search.
