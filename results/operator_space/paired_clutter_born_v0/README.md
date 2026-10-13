# Reduced paired-clutter operator benchmark O3 v0

## Outcome

The frozen four-channel, two-static-state entrant passes its intrinsic
detection, localization, false-alarm and PAUT-gap gates, but reduced `O3` as a
whole **does not pass** because the equal-time mechanical scan is better.

| architecture | channels | Tx events | positions / states | Pfa upper 95% | Pd at reference alpha | P(detect + localize) | alpha90 |
|---|---:|---:|---:|---:|---:|---:|---:|
| uniform A-scan | 1 | 1 | 1 position | 0.022% | 0.539 | 0.180 | 644.1 |
| equal-time mechanical scan | 1 | 8 | 8 positions | 0.926% | **0.990** | **0.983** | **177.8** |
| four-channel one-state control | 4 | 4 | 1 state | 2.408% | 0.953 | 0.936 | 233.4 |
| four-channel two-state entrant | 4 | 8 | 2 states | **0.745%** | 0.942 | 0.937 | 235.3 |
| 17-element FMC oracle | 17 | 17 | full aperture | 6.246% | 0.979 | 0.972 | 138.4 |

Lower `alpha90` is better. `reference alpha = 287.96` is frozen from the first
validation grid point where FMC reaches `Pd >= 0.95`; its absolute value is a
dimensionless Born amplitude and has no direct material calibration.

The entrant closes `80.8%` of the aggregate A-scan-to-FMC `alpha90` gap. Its
Q10 closure across all `100/100` held-out C3 contexts is `0.775`, above the
preregistered `0.50` gate. It also passes:

```text
Pfa upper 95%       0.00745 <= 0.0125
Pd(reference)       0.942   >= 0.90
Pdl(reference)      0.937   >= 0.80
Q90 location error  0.0948  <= 1.0 lambda_P0
```

It fails the final product comparator:

```text
alpha90(two static states) = 235.3
alpha90(mechanical scan)   = 177.8
ratio                      = 1.323
```

Thus the two-state operator requires about `32%` more scattering amplitude
than an eight-position, one-channel scan under the current equal-time
contract. No geometry handoff is authorized.

## Interpretation

The comparison separates three conclusions that would be conflated by a
focus-only study:

1. Relative to a stationary A-scan, the passive coded class is strongly
   informative: its `alpha90` is lower by a factor of `2.74` and localization
   is much better.
2. The second real symmetry state mainly stabilizes background rejection. The
   one-state and two-state `alpha90` values are almost identical, but the
   held-out false-alarm upper bound improves from `2.41%` to `0.75%`.
3. If equal-time mechanical motion is allowed, spatial translation is still
   the cheapest leading architecture in this reduced model.

The full FMC oracle has the best `alpha90`, but its frozen validation threshold
does not control false alarms on the independently generated test nuisance
subspace. This is not interpreted as PAUT being intrinsically poor. It shows
that high-dimensional FMC also pays a calibration/generalization cost; in this
campaign it remains an information oracle, not a deployable low-calibration
candidate.

## Frozen statistical protocol

- five architectures are assembled from one full virtual-aperture library;
- every sequential state/position shares a fixed total incident-energy budget;
- frequency features use the five-cycle Hann spectrum and correlated
  frequency covariance from O2c;
- GLRT templates and nuisance precision are trained only on the training
  library;
- thresholds use `20,000` independent validation no-defect trials;
- reported false alarms use `20,000` independent test no-defect trials;
- C3 contains `100` material/contact/calibration contexts and every one of the
  `136` frozen ROI positions at every amplitude;
- true signals cycle through isotropic, aligned anisotropic and 45-degree
  anisotropic P/S scattering tensors while the estimator retains its frozen
  isotropic templates;
- no defect-free reference is subtracted at inference.

Backwall, ringing and material nuisance atoms are generated in the full
17-element aperture before compression. Training, validation and test use
different spatial atoms. Their full-FMC mean powers relative to electronics
are fixed at `2`, `2` and `3`, respectively. These values define a
dimensionless stress test, not a material-specific clutter certification.

## Stage verdict

- `C0--C3` decision pipeline: **complete**;
- frozen two-state entrant intrinsic gates: **pass**;
- equal-time mechanical non-inferiority: **fail**;
- reduced `O3`: **fail**;
- full-wave O3: **not authorized**;
- `O4` transfer target and geometry handoff: **blocked**.

The result creates a product-level fork rather than a geometry task:

```text
motion allowed
    -> mechanical scan is the current reduced-model leader;

motion forbidden or strongly penalized
    -> explicitly revise O0 resource constraints, then evaluate the frozen
       static entrant under that new product contract;

static architecture must beat equal-time motion
    -> return to operator/code-state design with fresh train/validation/test
       splits; do not tune the frozen codebook on this test set.
```

## Reproduction

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_operator_o3_paired_clutter.jl
```

Primary artifacts:

- `architecture_metadata.csv` — channels, events, positions and energy rules;
- `nuisance_scales.csv` — frozen full-aperture clutter powers;
- `thresholds.csv` — validation thresholds and held-out Pfa;
- `detection_curves.csv` — C0--C3 aggregate `Pd/Pdl` curves;
- `c3_context_curves.csv` — per-context curves for audit;
- `c3_context_summary.csv` — per-context `alpha90` and reference metrics;
- `o3_summary.csv` — final gates;
- `o3_paired_clutter.jld2` — reproducible compact campaign bundle;
- `o3_detection_curves.png`, `o3_alpha90.png` — visual summaries.

## Scope limitation

This is a reciprocal elastic Born/reduced benchmark. It is not industrial POD,
does not model finite-crack contact or a manufactured front-end, and does not
replace the paired transient full-wave calculation required by the O3
contract. Its negative mechanical-comparator result is sufficient to stop the
current geometry handoff, but not to establish a universal no-go theorem for
passive coding.

