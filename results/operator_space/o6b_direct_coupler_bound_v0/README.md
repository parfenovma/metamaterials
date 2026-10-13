# O6b direct-coupler upper bound

**PASS: one finite mechanically coupled PZT/bond/Al twin is authorized.**

The reciprocal dense-band model uses exact plane-strain P/S loading of Al, a longitudinal PZT-5A mechanical port and a frozen `100 um` epoxy bond.

## Worst dense-band metrics

- normalized transfer fidelity: `0.99999309`;
- maximum port coherence: `2.400e-16`;
- maximum column imbalance: `0.07203 dB`;
- maximum differential amplitude RMS: `0.35276%`;
- maximum differential-delay RMS: `0.70457 ns`;
- minimum useful one-way amplitude: `0.90457124`;
- maximum singular value: `0.92416816`.

This pass removes the O6a lateral-guide mechanism but is only a necessary upper bound. It assumes an ideal distributed piezocomposite coupling profile and omits finite PZT lateral coupling, backing, electrode mass and dielectric loss. The next authorized calculation is exactly one finite 2D mechanical PZT/bond/Al twin with the same DCT0--2 profiles and bond thickness.
