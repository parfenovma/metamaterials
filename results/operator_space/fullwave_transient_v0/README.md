# O3 F1a: direct paired full-wave transient

> F1a physics gate: **PASS**. F1b stochastic detection is still required.

Five-cycle 242 kHz pulse, 0.138 us step, 510 samples, direct multi-RHS Newmark.

## Physics diagnostics

- reciprocity residual: `3.727e-14`
- paired mirror residual: `0.549%`
- precursor energy fraction: `0.000%`
- paired tangential/normal norm: `0.662`
- geometric void window: `18.70--41.12 us`
- backwall window: `40.57--61.75 us`

## Frozen pulse projections

| Architecture | E | Paired/free void, dB | Paired energy / FMC | Peak, us |
|---|---:|---:|---:|---:|
| `uniform_a_scan` | 1 | 15.20 | 7.119 | 31.22 |
| `trained_k4_q2_unconstrained` | 8 | -0.03 | 0.704 | 31.22 |
| `trained_k3_q4_reflection_balanced` | 12 | -5.77 | 0.350 | 31.91 |
| `motion_compact_center_p12_native` | 12 | -9.50 | 0.084 | 30.67 |
| `motion_uniform_full_p17_a8` | 17 | -10.18 | 0.070 | 30.53 |
| `paut_fmc_oracle` | 17 | -3.15 | 1.000 | 30.94 |

These are deterministic raw-pulse transfer metrics, not Pd/Pfa. F1b must add independent nuisance/noise splits and a finite defect ROI ensemble.
