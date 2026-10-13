# O6d mechanically isolated three-bus bound

**FAIL: the low-complexity isolated three-bus face is stopped before geometry.**

The measured O6c differential residual must be reduced to `kappa <= 0.136698`: at least `86.33%` suppression (`-17.28 dB` residual amplitude relative to O6c). This is an operator requirement, not a kerf-width prediction.

At that boundary the active constraint is differential delay (`82.600 ns`);
amplitude RMS is already `1.236%`.

## Ideal electrical coding frontier

| Internal map | Active taps | Fidelity | Amp RMS | Useful | Gate |
|---|---:|---:|---:|---:|---|
| `sparse_1_bus` | 17 | 0.95519 | 25.55% | 0.91228 | FAIL |
| `sparse_2_bus` | 34 | 0.99919 | 4.04% | 0.99727 | FAIL |
| `full_3_bus_upper_bound` | 50 | 1.00000 | 0.00% | 1.00000 | PASS |

The best one-bus and two-bus maps use an ideal lossless public-port rotation, continuous signed tap weights, fixed column balancing, zero network delay and no network loss. Their misses therefore cannot be repaired by adding real passive parasitics. The SO(3) search is deterministic (`200000` Haar rotations plus local refinement) but not a global certificate.

The full map passes analytically, but needs `50` nonzero bus-to-column couplings and `18` distinct normalized magnitudes in its continuous form. The first passing quantization entry is `full_3_bus_8_levels_per_port` (`17` distinct normalized ratios across the map). This is the explicitly identified passive `3-to-17` beamformer fallback from O5, not a pass of the selected monolithic topology.

No kerf FEM or F1d replay is authorized. The next decision must explicitly compare accepting the full passive beamformer fabrication burden with returning to a different architecture-frontier point; it may not hide those taps inside the word `composite`.
