# O3 F0: paired full-wave carrier operator

> F0 physics gate: **FAIL**. This is not a transient O3/POD result.

One 2D full-elastic FEM operator at 242 kHz was factored per paired mesh. All six architectures below are projections of the same 17 x 17 element response.

## Matrix diagnostics

| Case | Reciprocity | Mirror | Tangential / normal |
|---|---:|---:|---:|
| `none` | 2.060e-15 | 1.924e-03 | 3.336e-01 |
| `void_r2_focus` | 1.910e-15 | 3.123e-03 | 3.367e-01 |
| `paired_scattering` | 6.416e-14 | 1.018e-01 | 1.532e+00 |

## Frozen architecture projections

| Architecture | E | Observables | Paired/free, dB | Paired gain vs A, dB | FMC norm fraction |
|---|---:|---:|---:|---:|---:|
| `uniform_a_scan` | 1 | 1 | -23.17 | 0.00 | 1.770 |
| `trained_k4_q2_unconstrained` | 8 | 32 | -28.96 | -6.91 | 0.799 |
| `trained_k3_q4_reflection_balanced` | 12 | 36 | -32.63 | -10.37 | 0.537 |
| `motion_compact_center_p12_native` | 12 | 12 | -35.03 | -12.80 | 0.405 |
| `motion_uniform_full_p17_a8` | 17 | 17 | -36.13 | -13.92 | 0.357 |
| `paut_fmc_oracle` | 17 | 289 | -27.63 | -4.96 | 1.000 |

`Delta = void - free` is used here only to diagnose physical operator transfer. F1 must return to raw paired pulse-echo trials and may not give the detector a free reference.

Artifacts: `matrix_diagnostics.csv`, `architecture_projection.csv`, `f0_gates.csv`, `fullwave_carrier_summary.jld2`, and `fullwave_carrier_architectures.png`.
