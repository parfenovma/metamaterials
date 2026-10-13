# O6g-B frozen motion failure boundary

The exact O6f detector reconstruction, 60 signal contexts, physical standard draws and fresh-noise streams are reused. Every evaluated point contains 12,000 aggregate H0 trials; no layout or threshold was reoptimized.

| Layout | Scenario | Max passing scale | Position RMS | Contact log-amp RMS | Delay RMS | First failed gate |
|---|---|---:|---:|---:|---:|---|
| `motion_compact_center_p12_native` | `amplitude_only` | `5.52696` | `0.0000 lambda` | `27.63%` | `0.0000 cycles` | `pfa` |
| `motion_compact_center_p12_native` | `coupled` | `2.35522` | `0.1178 lambda` | `11.78%` | `0.0471 cycles` | `alpha90` |
| `motion_compact_center_p12_native` | `delay_only` | `3.14356` | `0.0000 lambda` | `0.00%` | `0.0629 cycles` | `pfa` |
| `motion_compact_center_p12_native` | `position_only` | `2.76808` | `0.1384 lambda` | `0.00%` | `0.0000 cycles` | `alpha90` |
| `motion_uniform_full_p17_a8` | `amplitude_only` | `4.53560` | `0.0000 lambda` | `22.68%` | `0.0000 cycles` | `pfa` |
| `motion_uniform_full_p17_a8` | `coupled` | `2.23883` | `0.1119 lambda` | `11.19%` | `0.0448 cycles` | `alpha90` |
| `motion_uniform_full_p17_a8` | `delay_only` | `2.69906` | `0.0000 lambda` | `0.00%` | `0.0540 cycles` | `pfa` |
| `motion_uniform_full_p17_a8` | `position_only` | `2.57225` | `0.1286 lambda` | `0.00%` | `0.0000 cycles` | `alpha90` |

These are finite boundaries of the reduced multiplicative/delay contact model, not real couplant specifications. O6g exhausts this surrogate; surviving motion risk is promoted to full-wave contact/lift-off in O6h.
