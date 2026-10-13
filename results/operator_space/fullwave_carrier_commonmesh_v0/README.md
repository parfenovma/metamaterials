# O3 F0: paired full-wave carrier operator

> F0 physics gate: **FAIL**. This is not a transient O3/POD result.

One 2D full-elastic FEM operator at 242 kHz was factored per paired mesh. All six architectures below are projections of the same 17 x 17 element response.

## Matrix diagnostics

| Case | Reciprocity | Mirror | Tangential / normal |
|---|---:|---:|---:|
| `none` | 1.861e-15 | 8.766e-03 | 3.352e-01 |
| `void_r2_focus` | 4.801e-15 | 8.870e-03 | 3.365e-01 |
| `paired_scattering` | 1.170e-13 | 4.695e-02 | 1.163e+00 |

## Frozen architecture projections

| Architecture | E | Observables | Paired/free, dB | Paired gain vs A, dB | FMC norm fraction |
|---|---:|---:|---:|---:|---:|
| `uniform_a_scan` | 1 | 1 | -21.03 | 0.00 | 2.291 |
| `trained_k4_q2_unconstrained` | 8 | 32 | -28.25 | -8.33 | 0.878 |
| `trained_k3_q4_reflection_balanced` | 12 | 36 | -32.38 | -12.28 | 0.557 |
| `motion_compact_center_p12_native` | 12 | 12 | -35.11 | -15.11 | 0.402 |
| `motion_uniform_full_p17_a8` | 17 | 17 | -36.09 | -16.10 | 0.359 |
| `paut_fmc_oracle` | 17 | 289 | -27.70 | -7.20 | 1.000 |

`Delta = void - free` is used here only to diagnose physical operator transfer. F1 must return to raw paired pulse-echo trials and may not give the detector a free reference.

Artifacts: `matrix_diagnostics.csv`, `architecture_projection.csv`, `f0_gates.csv`, `fullwave_carrier_summary.jld2`, and `fullwave_carrier_architectures.png`.
