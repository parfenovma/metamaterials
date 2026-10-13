# F1d expanded physical-context calibration

> Five development/train physical contexts, two fresh validation contexts and two fresh final-test contexts. No per-contact H0 reference is supplied to the detector.

Positive control: **PASS**. Selected operator: **`o5a_dct_k3_q1`**.

| Architecture | Ch | Q | Events | Pfa upper | Pd@3 | Pdl@3 | Min Pd | alpha90 | Gap closure | Gate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `uniform_a_scan` | 1 | 1 | 1 | 0.027% | 0.537 | 0.221 | 0.000 | 9.697 | 0.000 | fail |
| `o5a_dct_k3_q1` | 3 | 1 | 3 | 0.089% | 0.982 | 0.976 | 0.957 | 2.427 | 0.910 | pass |
| `o5a_dct_k4_q1` | 4 | 1 | 4 | 0.312% | 0.993 | 0.983 | 0.988 | 1.911 | 0.974 | pass |
| `trained_k4_q2_unconstrained` | 4 | 2 | 8 | 1.361% | 0.980 | 0.979 | 0.947 | 2.120 | 0.948 | fail |
| `paut_fmc_oracle` | 17 | 1 | 17 | 1.026% | 0.989 | 0.983 | 0.980 | 1.706 | 1.000 | fail |

The old `validation_1`, `test_1` and failed `o5a_test_1` contexts are explicitly development data in this revised calibration. Final metrics use only `f1d_test_1/2`.
