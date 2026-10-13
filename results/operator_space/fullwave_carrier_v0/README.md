# O3 F0: paired full-wave carrier operator

> F0 physics gate: **FAIL**. This is not a transient O3/POD result.

One 2D full-elastic FEM operator at 242 kHz was factored per paired mesh. All six architectures below are projections of the same 17 x 17 element response.

## Matrix diagnostics

| Case | Reciprocity | Mirror | Tangential / normal |
|---|---:|---:|---:|
| `none` | 2.235e-15 | 6.289e-03 | 3.428e-01 |
| `void_r2_focus` | 4.143e-15 | 8.870e-03 | 3.365e-01 |
| `paired_scattering` | 9.761e-14 | 2.578e-01 | 2.477e+00 |

## Frozen architecture projections

| Architecture | E | Observables | Paired/free, dB | Paired gain vs A, dB | FMC norm fraction |
|---|---:|---:|---:|---:|---:|
| `uniform_a_scan` | 1 | 1 | -24.95 | 0.00 | 1.269 |
| `trained_k4_q2_unconstrained` | 8 | 32 | -28.77 | -4.68 | 0.741 |
| `trained_k3_q4_reflection_balanced` | 12 | 36 | -32.19 | -7.89 | 0.511 |
| `motion_compact_center_p12_native` | 12 | 12 | -34.05 | -9.88 | 0.407 |
| `motion_uniform_full_p17_a8` | 17 | 17 | -34.75 | -10.58 | 0.375 |
| `paut_fmc_oracle` | 17 | 289 | -26.75 | -2.07 | 1.000 |

`Delta = void - free` is used here only to diagnose physical operator transfer. F1 must return to raw paired pulse-echo trials and may not give the detector a free reference.

Artifacts: `matrix_diagnostics.csv`, `architecture_projection.csv`, `f0_gates.csv`, `fullwave_carrier_summary.jld2`, and `fullwave_carrier_architectures.png`.
