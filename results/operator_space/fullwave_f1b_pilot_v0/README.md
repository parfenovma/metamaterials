# O3 F1b finite-defect stochastic pilot

> This is a six-class full-wave pilot, not industrial POD or continuous-ROI validation.

Common electronics sigma: `2.487889052190143e-12`; FMC validation reference amplitude: `2.2`.

| Architecture | Pfa upper 95% | Pd ref | Pdl ref | Min class Pd | Worst class | alpha90 | A->FMC closure | Gate |
|---|---:|---:|---:|---:|---|---:|---:|---|
| `uniform_a_scan` | 1.361% | 1.000 | 0.994 | 1.000 | `void_r1_focus` | 0.635 | n/a | fail |
| `trained_k4_q2_unconstrained` | 1.532% | 0.858 | 0.851 | 0.183 | `void_r1_focus` | 2.702 | n/a | fail |
| `trained_k3_q4_reflection_balanced` | 1.756% | 0.849 | 0.838 | 0.237 | `void_r1_focus` | 2.713 | n/a | fail |
| `motion_compact_center_p12_native` | 1.113% | 0.581 | 0.562 | 0.058 | `void_r1_focus` | 5.046 | n/a | fail |
| `motion_uniform_full_p17_a8` | 1.393% | 0.603 | 0.584 | 0.071 | `void_r1_focus` | 5.067 | n/a | fail |
| `paut_fmc_oracle` | 1.575% | 0.978 | 0.976 | 0.892 | `void_r1_focus` | 1.523 | n/a | fail |

The A-scan -> FMC closure is undefined because the uniform A-scan alpha90 is lower than the FMC alpha90 on this six-class task. This is a failed PAUT positive-control condition, not evidence that A-scan dominates general NDT.

Thresholds and nuisance PCs were frozen on train/validation splits. The detector never receives a per-contact defect-free reference; `free` signals enter only the fixed calibration/nuisance model and raw H0 trials.
