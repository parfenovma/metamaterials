# Passive standoff / known-depth reranking

## Bottom line

The 6--13 mm passive focal offsets are not all topology failures.  Re-evaluating
the saved fine fields at the actual local working depth converts five
frequency/topology combinations into strict candidates and one more into a
prototype candidate, without changing a hole or running a new FEM solve.

The strongest passive result is ellipse `c005` at 305 kHz with a 6 mm assembly
shift toward a fixed defect:

- target coordinate `x=109 mm`, or 54 mm beyond the plate exit;
- `J1=1.3642`;
- `J2_equal=1.7245` and `J2_clutter_dominated=1.6557`;
- exit-energy proxy `0.7581`;
- FWHM `25.81 mm`;
- local focus error `2 mm`;
- **strict gate passes**.

At 310 kHz, `c005` gives nearly the same outcome with a 9 mm shift:
`J1=1.3533`, `J2=1.7278`, exit energy `0.7806` and focus error `1 mm`.
Matched-GRIN `g001` and round-hole Fresnel `f006` also become strict candidates
at simple 5--6 mm longitudinal shifts.

## Critical validity condition

This result is a **coordinate translation**, not a new wave-path simulation.
It is valid only when:

1. the defect depth or inspection surface position is known accurately enough
   to set the plate-to-target distance;
2. the complete transducer/lens assembly is repositioned without changing the
   material, thickness, impedance or boundary conditions of the acoustic
   coupling path;
3. the same source normalization and carrier are retained.

If a physical spacer adds a new material layer, changes couplant thickness,
introduces extra interfaces, or alters clamping/contact, these saved fields are
not valid for it.  That configuration needs a new coupled-path FEM or
experiment.  Likewise, a fixed-focus standoff is unsuitable when defect depth
is unknown or varies over the scan.

## Selected standoffs

Positive assembly shift means moving the existing plate/transducer assembly
toward a defect fixed at the original nominal position.  Equivalently, it
reduces the plate-exit-to-target distance from the nominal 60 mm.

| topology | f | target x | assembly shift | J1 | J2 | exit | FWHM | error | result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| ellipse `c005` | 295.0 | 105 mm | +10 mm | 1.1376 | 1.3602 | 0.6301 | 25.05 mm | 22 mm | infeasible; competing lobe and low energy |
| ellipse `c005` | 300.0 | 110 mm | +5 mm | 1.2393 | 1.4987 | 0.6869 | 34.66 mm | 2 mm | focus aligned, energy fail |
| **ellipse `c005`** | **305.0** | **109 mm** | **+6 mm** | **1.3642** | **1.7245** | **0.7581** | **25.81 mm** | **2 mm** | **strict pass** |
| **ellipse `c005`** | **310.0** | **106 mm** | **+9 mm** | **1.3533** | **1.7278** | **0.7806** | **22.88 mm** | **1 mm** | **strict pass** |
| ellipse `c005` | 315.0 | 105 mm | +10 mm | 1.2793 | 1.5969 | 0.7090 | 22.14 mm | 1 mm | prototype pass |
| ellipse `c005` | 321.7 | 105 mm | +10 mm | 1.0818 | 1.2365 | 0.6175 | 21.97 mm | 3 mm | focus aligned, energy fail |
| **matched-GRIN `g001`** | **315.0** | **110 mm** | **+5 mm** | **1.1289** | **1.5343** | **0.9612** | **20.65 mm** | **1 mm** | **strict pass** |
| **matched-GRIN `g001`** | **316.0** | **110 mm** | **+5 mm** | **1.1508** | **1.5759** | **0.9714** | **20.32 mm** | **1 mm** | **strict pass** |
| **Fresnel `f006`** | **321.7** | **109 mm** | **+6 mm** | **1.0793** | **1.2715** | **0.8540** | **22.12 mm** | **3 mm** | **strict pass** |

The 295 kHz ellipse is a useful warning: moving the nominal coordinate toward
the earlier reported peak exposes a stronger upstream lobe inside the shifted
22 mm search window.  A spacer therefore cannot be inferred safely from the
old scalar `focus_error` alone; the full field must be reranked, as done here.

## Improvement relative to the original 115 mm target

| candidate | nominal J1 / J2 / error | selected J1 / J2 / error |
|---|---:|---:|
| c005, 305 kHz | 1.1929 / 1.2859 / 8 mm | **1.3642 / 1.7245 / 2 mm** |
| c005, 310 kHz | 1.1963 / 1.3116 / 10 mm | **1.3533 / 1.7278 / 1 mm** |
| c005, 315 kHz | 1.2057 / 1.3859 / 11 mm | **1.2793 / 1.5969 / 1 mm** |
| g001, 315 kHz | 1.0710 / 1.3533 / 6 mm | **1.1289 / 1.5343 / 1 mm** |
| g001, 316 kHz | 1.0619 / 1.3110 / 6 mm | **1.1508 / 1.5759 / 1 mm** |
| f006, 321.7 kHz | 1.0307 / 1.1512 / 9 mm | **1.0793 / 1.2715 / 3 mm** |

The straight reference is recomputed from the saved same-frequency fine field
at every target x.  J2's defect coordinate, axial clutter window and exclusion
radius are translated together.  These gains therefore do not arise from
normalizing a new target against the old 115 mm reference.

## Method

- Reused six `c005` fine fields from 295--321.7 kHz, `g001` fine fields at 315
  and 316 kHz, and `f006` at 321.7 kHz.
- Swept integer target coordinates from `x=100` through `120 mm`.
- At each coordinate recalculated J1, equal-noise and clutter-dominated J2,
  target-plane FWHM, and the strongest-point error in a translated
  `+/-22 mm` by `+/-30 mm` local window.
- Retained the original exit-energy proxy because longitudinal coordinate
  translation does not alter the saved exit line.
- Selected the maximum `min(J1,J2)` strict-feasible row; if none existed, the
  maximum prototype-feasible or focus-aligned energy-failing row.

These are still harmonic reciprocal-Born proxies.  Standoff reranking does not
predict transient ringing, group delay, real piezo bandwidth or pulse-echo
time separation.

## Artifacts

- `data/target_sweep.csv`: all 189 target/frequency/topology rows.
- `data/best_standoff.csv`: one selected longitudinal setting per saved field.
- `data/campaign_metrics.csv`: the same sweep using the common campaign prefix
  and topology labels `ellipse_standoff`, `matched_grin_standoff` and
  `fresnel_standoff`.
- `figures/ellipse_target_frequency_heatmap.png`: c005 J1/J2 depth-frequency
  map.
- `figures/focus_error_vs_target_x.png`: translated focus gates.
- `figures/J1_J2_vs_target_x.png`: same-x-normalized gain curves.
- `figures/standoff_J1_J2.png`: full reranked metric plane.
- `figures/best_standoff_shift.png`: selected mechanical shifts.
- `passive_standoff_rerank.jl`: reproducible field-only analysis.

## Reproduction

From the repository root:

```sh
julia --startup-file=no --project=. \
  results/passive_standoff_rerank_p6_315khz/passive_standoff_rerank.jl --report
```

The command only loads existing JLD2 field scans and writes CSV/plots.  It does
not assemble or solve a FEM system.
