# O3 F1c structural-clutter stochastic pilot

> Four actual full-wave physical contexts; nuisance/noise replication is stochastic and is not an industrial POD sample size.

Positive control: **PASS**. Common electronics sigma: `8.667116441957401e-13`; FMC validation reference amplitude: `3.0`.

FMC stable rank `4.269`, maximum template coherence `0.131`, alpha90 ratio FMC/A `0.145`, reference exact Pdl `0.980`, test Pfa upper `0.027%`.

| Architecture | Pfa upper 95% | Pd ref | exact Pdl ref | min Pd | alpha90 | localized alpha90 | A->FMC closure | Gate |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| `uniform_a_scan` | 0.027% | 0.369 | 0.156 | 0.000 | 14.817 | inf | 0.000 | fail |
| `trained_k4_q2_unconstrained` | 0.060% | 0.935 | 0.934 | 0.792 | 2.782 | 2.787 | 0.950 | pass |
| `trained_k3_q4_reflection_balanced` | 0.075% | 0.834 | 0.833 | 0.643 | 3.554 | 3.569 | 0.889 | fail |
| `motion_compact_center_p12_native` | 2.117% | 0.235 | 0.232 | 0.021 | 5.803 | 5.804 | 0.711 | fail |
| `motion_uniform_full_p17_a8` | 1.915% | 0.343 | 0.338 | 0.094 | 5.452 | 5.452 | 0.739 | fail |
| `paut_fmc_oracle` | 0.027% | 0.980 | 0.980 | 0.971 | 2.144 | 2.144 | 1.000 | fail |

All projectors/templates/nuisance bases use only `train_1` and `train_2`; threshold and reference amplitude use only `validation_1`; final metrics use only `test_1`.
