# O4 transfer target: smooth K3,Q1

> Final operator handoff after expanded F1d physical-context calibration.

The target has three fixed reciprocal modal ports and one passive state. The aperture maps are the first three orthonormal DCT modes: uniform even, first odd and second even. There is no element-wise phase control and no state switch.

```text
H(omega) = g(omega) exp(-i omega tau0) [S_DCT3 + E(omega)]
reverse transfer = transpose(H); sigma_max(H) <= 1
```

Qualified band: `193.6--290.4 kHz`. Common gain/delay are calibrated gauges; differential shape errors are acceptance variables.

## Hard acceptance numbers

- coefficient amplitude RMS <= 3%; differential delay RMS <= 82.6 ns (7.2 deg at carrier);
- transfer fidelity >= 0.98; port coherence <= 0.15; column imbalance <= 0.75 dB; reciprocity residual <= 3%;
- with no drive/noise compensation: one-way useful amplitude >= 0.8994, insertion loss <= 0.921 dB one way and <= 1.843 dB round trip;
- as-built calibrated H(omega) must replay the expanded F1d gates; coefficient tolerances alone are not a substitute.

The strict loss floor assumes equal device-side energy. Extra pulser headroom is allowed only when its cost and output-referred noise penalty are included explicitly in the resource vector.
