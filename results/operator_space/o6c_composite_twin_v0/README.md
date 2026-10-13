# O6c homogenized 1-3 composite finite twin

**FAIL: the fixed-thickness homogenized 1-3 mechanism is stopped.**

One frozen `65% ceramic` orthotropic effective medium was substituted into the unchanged `6.0 mm composite -> 0.10 mm epoxy -> 35 mm Al` O6b mesh. Its density is `5461.0 kg/m^3`, axial impedance `17.00 MRayl`, and tensor `(C11,C22,C12,C66)=(52.921, 17.465, 8.674, 4.391) GPa`. No composition, tensor, thickness or geometry sweep was performed.

## Worst dense-band metrics

- transfer fidelity: `0.993835`;
- maximum port coherence: `0.017475`;
- column imbalance: `0.1728 dB`;
- differential amplitude RMS: `9.024%`;
- differential-delay RMS: `640.985 ns`;
- minimum useful amplitude: `0.929702`;
- maximum sigma: `0.973904`;
- mechanical reciprocity residual: `5.940e-16`.

Compared with homogeneous O6b, the composite raises worst fidelity from
`0.97331` to `0.99383`, reduces differential amplitude RMS from `19.59%` to
`9.02%`, and reduces differential-delay RMS from `846 ns` to `641 ns`. The
direction is correct, but the remaining misses are still factors of `3.01`
and `7.76` beyond their gates. DCT2 is the worst individual port (`10.65%`,
`784 ns`).

The unchanged `6 mm` layer has a nominal axial half-wave frequency of `259.42 kHz`, inside the qualified band; that resonance was not retuned away. The common compliance and common delay remain calibrated gauges exactly as in O6b; FEM supplies the differential finite-aperture transfer and reciprocity.

Per the preregistered stopping rule, no fill-fraction/thickness sweep, explicit-cell geometry, electrical terminal model or F1d replay is authorized from this result.
