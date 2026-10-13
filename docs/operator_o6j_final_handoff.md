# Final numerical architecture handoff after O6j

## Recommended next build

The only architecture that has passed the complete numerical sequence through
O6j is the compact P12 mechanical scan with a shear-transmitting solid contact:

```text
electronic channels       1
transmit events           12
mechanical positions      12 central stops
mechanical moves          11
aperture span             5.5 lambdaP
contact class             solid/dry, P and S transmitting
```

The frozen effective contact target is approximately:

```text
P-wave impedance          12.57 MRayl
cP                        6122 m/s
cS                        1531 m/s
density                   2054 kg/m3
thickness                 506 um
```

This is a property target, not a claim that a particular commercial film has
those values.

## Evidence

- The saved P12 detector, training nuisance precision and threshold are never
  retrained after architecture selection.
- Original signal/localization/H0 counts reproduce exactly.
- Original reduced coupled encoder/contact boundary: scale `2.35522`.
- Liquid thin-layer contact is rejected by exact 2D angular P/S transmission.
- A broad exact three-solid property map contains 122 passing points.
- All 16 one-sigma attachment corners of the frozen interior solid contact
  pass the angular gates.
- With fixed `0.05 lambda` encoder RMS and nominal solid-attachment variation,
  P12 has Pfa upper `0.441%`, alpha90 ratio `1.00874`, reference Pdl `0.9055`.

## Stopped alternatives

- P17 retains better nominal absolute sensitivity but fails the nominal
  physical liquid-layer Pfa gate and has smaller tolerance on every O6g ray.
- The q8 static network is feasible as an ideal reciprocal topology, but a
  single-resonance intrinsic-Q=35 BVD load needs effective Q `<=1.48148`.
  Resistive broadening would retain only `4.23%` useful power. q8 therefore
  remains blocked until a concrete lossless or multi-resonant matching
  topology is specified and loaded.
- Liquid couplant is stopped for the present ROI because its shear-free
  interfaces fail the exact angular power and mode-conversion gates.

## Meaning of “prototype”

The recommendation authorizes only material identification, a finite-footprint
model and a benchtop coupon. It is not a product design or POD certification.
The next evidence must measure contact repeatability, edge/squeeze-out fields,
wear/debonding, acquisition time and encoder performance on real hardware.
