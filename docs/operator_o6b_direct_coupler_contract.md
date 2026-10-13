# O6b direct load-local modal coupler contract

> Frozen before the dense-band evaluation. O6a rejected lateral propagation
> through a wide Al guide. O6b moves the three modal couplers to the loaded
> contact plane; it does not change the selected K3,Q1 measurement operator.

## 1. Topology class

```text
three electronic ports
    -> three fixed distributed coupling profiles DCT0 / DCT1 / DCT2
    -> thickness-mode piezocomposite source
    -> controlled thin bond
    -> Al specimen at the same contact plane.
```

The profiles are global even/odd/even functions. They may eventually be
implemented by patterned poling, effective `e33(y)`, 1-3 fill fraction or a
multilayer electrode stack, but O6b does not choose among those fabrication
mechanisms. A realization with 17 independently driven cells is forbidden.

There is no laterally guiding metal manifold. The coupling is local in the
depth direction, so the O6a strip cut-ons are removed by topology rather than
tuned away.

## 2. Frozen upper-bound model

The first gate is a reciprocal power-wave model for each lateral harmonic
`k_y=q*pi/W`, `q=0,1,2`, with:

```text
W                         71.6 mm
band                      193.6--290.4 kHz, 49 points
PZT mechanical surrogate rho=7750 kg/m3, cP=4350 m/s
bond surrogate            rho=1150 kg/m3, cP=2500 m/s
maximum controlled bond   100 um
Al 6061                   project nominal rho, cP, cS
```

The Al load is not reduced to `rho*cP`. Its normal surface impedance is
obtained from the exact outgoing plane-strain P/S potential pair under zero
tangential traction. The PZT and bond are an optimistic longitudinal
power-wave model: lateral piezo coupling, dielectric loss, electrode mass and
backing are not yet included.

For each `q`, the complex bond transfer is computed from its two-port matrix.
No per-frequency correction is allowed. The three DCT profiles are exactly
orthogonal in this infinite/local upper bound, so any later loss of isolation
belongs to the finite coupled twin rather than to the target definition.

## 3. Gates

The upper bound uses the unchanged O4 limits:

```text
normalized H fidelity           >= 0.98
maximum port coherence          <= 0.15
column imbalance                <= 0.75 dB
differential amplitude RMS      <= 3%
differential delay RMS          <= 82.6 ns
one-way useful amplitude        >= 0.8994
sigma_max(H)                    <= 1
reciprocity residual            <= 3%
```

Common gain and common delay remain calibrated gauges; the useful-amplitude
gate is evaluated before removing common insertion loss.

## 4. Decision rule

- **Fail:** stop the direct-coupler topology. Do not build piezo geometry.
- **Pass:** build one finite 2D mechanically coupled PZT/bond/Al twin with the
  same three source profiles and the same `100 um` bond. Its first purpose is
  to expose lateral PZT coupling, finite-aperture edge effects, loss and
  cross-port mixing.
- Only after the mechanical twin passes is an electrical terminal model and
  frozen F1d replay authorized.

The upper bound is deliberately necessary but not sufficient. Its pass is not
a prototype or fabrication claim.

## 5. Result

The 49-point upper bound passed, with minimum useful amplitude `0.90457`,
fidelity `0.999993`, differential amplitude RMS `0.353%` and differential
delay RMS `0.705 ns`.

The authorized finite homogeneous-PZT twin then failed: fidelity `0.97331`,
amplitude RMS `19.59%`, differential delay `846 ns` and useful amplitude
`0.88250`. A three-anchor refinement preserved worst fidelity
`0.97331 -> 0.97262` and amplitude RMS `19.59 -> 19.84%`.

The stopping rule was applied. No electrical terminal model and no F1d replay
were run. O6b rejects a homogeneous isotropic PZT slab, not the ideal
load-local bound. The next authorized mechanism is one homogenized laterally
decoupled `1-3` composite layer; thickness sweeps are not authorized.
