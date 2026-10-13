# O5a smooth operator projection — validation freeze

> Train/validation only. No test response is loaded by this script.

| DCT modes M | State fidelity | Min principal cosine | Rank | Coherence | Pd@3 | Pdl@3 | Min Pd@3 | alpha90 val | Gate |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 4 | 0.028 | 0.236 | 4.124 | 0.178 | 0.983 | 0.958 | 0.938 | 2.270 | pass |
| 5 | 0.139 | 0.328 | 3.949 | 0.257 | 0.978 | 0.959 | 0.938 | 2.226 | pass |
| 6 | 0.174 | 0.346 | 3.710 | 0.298 | 0.984 | 0.966 | 0.958 | 2.192 | pass |
| 7 | 0.351 | 0.419 | 3.979 | 0.235 | 0.977 | 0.966 | 0.938 | 2.311 | pass |
| 8 | 0.368 | 0.464 | 3.881 | 0.262 | 0.967 | 0.944 | 0.917 | 2.642 | pass |
| 9 | 0.411 | 0.633 | 3.722 | 0.309 | 0.957 | 0.945 | 0.885 | 2.723 | pass |
| 10 | 0.685 | 0.749 | 3.678 | 0.334 | 0.950 | 0.941 | 0.865 | 2.799 | pass |
| 11 | 0.480 | 0.753 | 3.561 | 0.380 | 0.857 | 0.852 | 0.734 | 3.409 | fail |
| 12 | 0.715 | 0.764 | 3.556 | 0.381 | 0.878 | 0.871 | 0.755 | 3.239 | fail |

Frozen smallest passing projection: **M=4**. A fresh physical test context is required before geometry synthesis.
