# O5a fresh full-wave test

> Candidate set and priority were frozen on train/validation before `o5a_test_1` was solved.

Positive control: **FAIL**. Selected geometry-aware operator: **`none`**.

| Architecture | Channels | States | Events | Pfa upper | Pd@3 | Pdl@3 | Min Pd | alpha90 | Gap closure | Gate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `uniform_a_scan` | 1 | 1 | 1 | 0.851% | 0.192 | 0.147 | 0.000 | 23.879 | 0.000 | fail |
| `o5a_dct_k3_q1` | 3 | 1 | 3 | 43.395% | 0.922 | 0.893 | 0.620 | 2.863 | 0.969 | fail |
| `o5a_dct_k4_q1` | 4 | 1 | 4 | 13.146% | 0.989 | 0.977 | 0.969 | 2.477 | 0.987 | fail |
| `paut_fmc_oracle` | 17 | 1 | 17 | 1.841% | 0.974 | 0.972 | 0.948 | 2.192 | 1.000 | fail |
