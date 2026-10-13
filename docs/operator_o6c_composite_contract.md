# O6c homogenized 1-3 composite contract

> Frozen before the O6c finite solve. O6b proved that ideal local coupling can
> satisfy O4 but a homogeneous isotropic PZT slab cannot. O6c tests one
> physically motivated laterally decoupled effective medium without changing
> thickness, bond or load geometry.

## 1. Fixed material point

The selected point is a `65%` ceramic random-fibre `1-3` composite. This is
not an optimization variable. Published experimental work reports about
`17 MRayl` acoustic impedance for such a composite and strongly attenuated
radial modes; classic 1-3 theory identifies reduced lateral clamping and
lower acoustic impedance as the relevant mechanism.

The project converts that one material point into a plane-strain orthotropic
mechanical surrogate:

```text
ceramic volume fraction       0.65
polymer volume fraction       0.35
effective density             rule of mixtures
axial stiffness C11           from Z_axial = 17 MRayl
lateral stiffness C22         harmonic mixture
cross stiffness C12           harmonic mixture of constituent lambda
shear stiffness C66           harmonic mixture
quality factor                35, frozen project surrogate
```

The constituent PZT and polymer constants are the same nominal mechanical
values already used in the repository. No fitted O4 data enter the tensor.

## 2. Unchanged geometry and excitation

```text
6.0 mm homogenized composite
    -> 0.10 mm epoxy bond
    -> 35 mm matched Al load
```

The mesh, 71.6-mm aperture, thirteen frequencies, three continuous DCT0--2
back-face tractions and 17-bin bonded-interface observer are unchanged from
O6b. The composite thickness is deliberately not retuned even though the
lower axial wave speed moves its nominal half-wave thickness resonance toward
the qualified band.

## 3. Gates

All O4 limits remain unchanged: fidelity `>=0.98`, coherence `<=0.15`,
imbalance `<=0.75 dB`, differential amplitude RMS `<=3%`, differential delay
RMS `<=82.6 ns`, useful amplitude `>=0.8994`, `sigma_max<=1`, reciprocity
residual `<=3%`.

The common power-wave scale is recomputed analytically from the frozen
composite impedance and the same bond. FEM supplies the finite-aperture
differential transfer. Common compliance and common delay remain calibrated
gauges.

## 4. Decision rule

- Failure stops the homogenized 1-3 mechanism at fixed thickness. No volume
  fraction, tensor or thickness sweep is allowed.
- Pass requires a three-anchor mesh check.
- Only after mesh-confirmed pass may one explicit pillar/kerf unit-cell be
  checked for homogenization validity and manufacturable pitch.
- Frozen F1d replay is still required before an electrical terminal model or
  prototype handoff.
