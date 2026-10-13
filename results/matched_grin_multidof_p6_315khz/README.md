# Matched-GRIN topology campaign near 315 kHz

## Bottom line

The matched-GRIN hypothesis improved transmission and reduced the focal offset
relative to the earlier three-row ellipse plate, but it did **not** pass the
campaign's prototype gate after mesh refinement.  It is therefore a useful
near-miss, not a fabrication-ready winner.

The best gate-oriented design is `g001`: five longitudinal rows, circular
low-area outer matching rows and a progressively anisotropic three-row core.
At 315--316 kHz its amplitude and reduced-CNR gains are mesh-convergent and its
exit-energy proxy is very high, but the fine-grid local maximum is 6 mm
upstream of the nominal `(115,0) mm` target.  The prototype limit is 5 mm and
the strict limit is 3 mm.

| variant | f | mesh | J1 | J2 equal | J2 dB | J2 clutter-dominated | exit energy | FWHM at x=115 | focus error | prototype | strict |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|:---:|:---:|
| `g001` | 315 kHz | coarse | 1.0555 | 1.3295 | +2.474 | 1.5400 | 0.8895 | 23.00 mm | 5 mm | yes | no |
| `g001` | 315 kHz | fine | 1.0710 | 1.3533 | +2.628 | 1.5474 | 0.9612 | 26.77 mm | 6 mm | **no** | no |
| `g001` | 316 kHz | coarse | 1.0590 | 1.3320 | +2.490 | 1.5343 | 0.8944 | 24.19 mm | 5 mm | yes | no |
| `g001` | 316 kHz | fine | 1.0619 | 1.3110 | +2.352 | 1.4761 | 0.9714 | 28.20 mm | 6 mm | **no** | no |

There are two prototype-feasible coarse rows, zero prototype-feasible fine
rows and zero strict-feasible rows.  No fine-validated joint improvement is
claimed.

## Comparison with the previous `c005` ellipse plate

The fair same-frequency, fine-mesh comparison at 315 kHz is:

| topology | J1 | J2 equal | min(J1,J2) | exit energy | focus error |
|---|---:|---:|---:|---:|---:|
| ellipse `c005` | 1.2057 | 1.3859 | 1.2057 | 0.7090 | 11 mm |
| matched-GRIN `g001` | 1.0710 | 1.3533 | 1.0710 | 0.9612 | 6 mm |

Matched-GRIN sacrifices about 11% of the joint scalar gain while increasing
the exit-energy proxy by 35.6% relative (`0.9612/0.7090`) and nearly halving the
focus error.  This supports the physical decoupling idea: the outer matching
rows mainly recover transmitted energy, but the phase-forming core still needs
one more controlled degree of freedom to move the field maximum from x=109 to
x=115 mm.

## Raw Pareto oracle and why it is rejected

`g006@305 kHz` is the raw coarse Pareto point:

- `J1 = 1.2379`;
- `J2_equal = 1.4478` (`+3.214 dB`);
- `exit_energy_proxy = 0.9091`;
- `focus_error = 22.20 mm`, with the maximum at `(93,3) mm`.

It dominates the sampled raw `(J1,J2)` plane but fails both focus gates by a
large margin.  It is retained as an oracle showing that nominal-point gains can
coexist with a stronger unwanted upstream lobe; it must not enter the common
feasible Pareto set.

## Geometry and manufacturing checks

The 18-case primary DOE varies five independent design variables plus carrier
frequency:

1. 5 versus 7 longitudinal rows;
2. total void-area scale;
3. area fraction assigned to each circular outer matching row;
4. central-core ellipse aspect ratio;
5. transverse delay-profile exponent;
6. carrier frequency: 305, 315 and 321.7 kHz.

Two lower-void geometries and three local carrier points at 316--318 kHz bring
the coarse candidate total to 23.  Two same-frequency fine solutions validate
`g001` at 315 and 316 kHz.

All candidates are monolithic 2D through-cut patterns.  Geometry validation
uses exact gaps for aligned ellipses and a conservative bounding-circle gap
for every diagonal pair.  Across the campaign:

- specified minimum void dimension and solid web: 0.8 mm;
- smallest realised void dimension: 0.913 mm;
- smallest conservative solid web: 2.526 mm;
- `g001`: 45 holes, minimum void dimension 1.183 mm, minimum conservative web
  4.396 mm.

No manufacturing tolerance was relaxed to obtain a feasible mesh.

## Metric and reference contract

Every row is normalized to a solid straight footprint at the **same frequency,
mesh resolution and energy-normalized uniform traction**.  The solver, scan
grid and definitions of `J1`, reciprocal Born `J2`, clutter mask, exit-energy
proxy, nominal-plane FWHM and local-window focus error are exactly those used by
the ellipse campaign and the common
`topology_campaign_p6_315khz/benchmark_contract.md`.

The prototype gate is `exit_energy_proxy >= 0.70` and `focus_error <= 5 mm`.
The strict gate is `exit_energy_proxy >= 0.75` and `focus_error <= 3 mm`.
Pareto dominance is computed on raw `(J1,J2_equal)`; the scalar
`min(J1,J2_equal)` is only a ranking diagnostic.

The coarse lens meshes contain about 30--32k nodes and 59--61k 2D elements.
The fine `g001` mesh contains 62,769 nodes and 123,357 2D elements; its matching
straight mesh contains about 19.8k nodes and 39.7k total mesh elements.

These are linear steady-state carrier responses.  They do not predict group
delay, ringing or transient pulse-echo CNR.

## Artifacts

- `data/campaign_metrics.csv`: common campaign schema plus branch diagnostics;
  23 coarse and 2 fine rows, regenerated without duplicates.
- `data/doe_manifest.csv`: all geometry variables and conservative feature/web
  checks.
- `data/coarse_fine_convergence.csv`: paired 315/316 kHz convergence data.
- `data/geometry_g001.csv` and `data/zones_g001.csv`: fabrication coordinates
  and zone targets for the near-miss design.
- `figures/frequency_case_heatmap.png`: primary 6-by-3 DOE heatmaps.
- `figures/campaign_J1_J2_pareto.png`: raw Pareto and prototype-gate markers.
- `figures/g001_local_frequency.png`: local carrier/gate response.
- `figures/geometry_g001_fine_validated.png`: plate layout.
- `figures/field_g001_f0316p0_fine.png`: fine lens/straight field comparison.
- `figures/geometry_g006_coarse_best.png`: rejected raw-Pareto oracle geometry.
- `matched_grin_pipeline.jl`: geometry, meshing, FEM, metrics, plotting and CLI.

## Reproduction

From the repository root:

```sh
julia --startup-file=no --project=. \
  results/matched_grin_multidof_p6_315khz/matched_grin_pipeline.jl --manifest
julia --startup-file=no --project=. \
  results/matched_grin_multidof_p6_315khz/matched_grin_pipeline.jl --sweep coarse
julia --startup-file=no --project=. \
  results/matched_grin_multidof_p6_315khz/matched_grin_pipeline.jl --case g001 316 fine
julia --startup-file=no --project=. \
  results/matched_grin_multidof_p6_315khz/matched_grin_pipeline.jl --fine-report g001 316
julia --startup-file=no --project=. \
  results/matched_grin_multidof_p6_315khz/matched_grin_pipeline.jl --finalize
```

The next sensible refinement is not a wider blind sweep.  Add an independent
longitudinal phase-control variable to the core (for example, core-row x
positions or per-row area weights) and optimize directly against focus error
while holding the demonstrated high transmission.  That targets the remaining
6 mm offset without discarding the matching benefit.
