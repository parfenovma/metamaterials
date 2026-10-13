# O3 F0: paired full-wave carrier operator

> F0 physics gate: **PASS**. This is not a transient O3/POD result.

One 2D full-elastic FEM operator at 242 kHz was factored per paired mesh. All six architectures below are projections of the same 17 x 17 element response.

## Matrix diagnostics

| Case | Reciprocity | Mirror | Tangential / normal |
|---|---:|---:|---:|
| `none` | 2.450e-15 | 3.185e-03 | 3.325e-01 |
| `void_r2_focus` | 2.070e-15 | 3.111e-03 | 3.367e-01 |
| `paired_scattering` | 5.848e-14 | 1.805e-02 | 1.250e+00 |

## Frozen architecture projections

| Architecture | E | Observables | Paired/free, dB | Paired gain vs A, dB | Paired norm / FMC |
|---|---:|---:|---:|---:|---:|
| `uniform_a_scan` | 1 | 1 | -21.43 | 0.00 | 2.250 |
| `trained_k4_q2_unconstrained` | 8 | 32 | -28.35 | -8.16 | 0.880 |
| `trained_k3_q4_reflection_balanced` | 12 | 36 | -32.48 | -12.08 | 0.560 |
| `motion_compact_center_p12_native` | 12 | 12 | -35.12 | -14.76 | 0.411 |
| `motion_uniform_full_p17_a8` | 17 | 17 | -36.22 | -15.88 | 0.362 |
| `paut_fmc_oracle` | 17 | 289 | -27.85 | -7.04 | 1.000 |

`Delta = void - free` is used here only to diagnose physical operator transfer. F1 must return to raw paired pulse-echo trials and may not give the detector a free reference.

A single centred circular target at one carrier is close to an even rank-one task and is unusually favourable to the uniform A-scan. This table is therefore a transfer diagnostic, not an architecture ranking.

Artifacts: `matrix_diagnostics.csv`, `architecture_projection.csv`, `f0_gates.csv`, `fullwave_carrier_summary.jld2`, and `fullwave_carrier_architectures.png`.
