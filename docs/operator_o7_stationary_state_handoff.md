# O7 handoff: stationary two-bit global boundary

> Completed 15 August 2026.  Contract:
> [`operator_o7_global_state_contract.md`](operator_o7_global_state_contract.md).
> Numerical bundle:
> [`results/operator_space/o7_global_state_v0/`](../results/operator_space/o7_global_state_v0/README.md).

## Outcome

Both preregistered stationary two-bit mechanisms pass the reduced
coupled-modal gate on a fresh split.  This resolves the immediate question:
motion of the head or sample is not required *at the operator level*.  Four
global states can add information beyond the same three-port head observed at
seven frequencies in one fixed state.

| Architecture | K | Q | events | Pfa95 | alpha90 | Pd@FMC ref | Pdl@FMC ref | gain over matched K3Q1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| PZT two-bit | 3 | 4 | 12 | 1.134% | 285.15 | 0.987 | 0.851 | 17.84% |
| Mechanical two-bit | 3 | 4 | 12 | 0.796% | 270.67 | 0.996 | 0.950 | 34.16% |
| Frozen ideal K3Q4 | 3 | 4 | 12 | 0.229% | 259.44 | 0.979 | 0.976 | — |
| FMC | 17 | 1 | 17 | 0.116% | 220.19 | 0.984 | 0.982 | — |

The PZT candidate is a near-margin pass: `Pdl=0.85098` against a `0.85` gate
and `Pfa95=1.134%` against `1.2%`.  The mechanical candidate has visibly more
margin.  Neither is a product-level winner; both remain frozen entrants for
the next evidence level.

## What the frozen candidates look like

The aperture spans `-4...+4 lambdaP`.  Both heads are one mirror-symmetric
monolithic block with three reciprocal ports: two even and one odd.  Two
spatially separated symmetric control zones commute with reflection.  Their
binary settings produce `q00`, `q10`, `q01`, `q11` without moving the head on
the sample.

For the PZT entrant, zone centres and Gaussian reduced widths, normalized by
the aperture half-width, are:

```text
inner pair A    |x|=0.198, width=0.099
outer pair B    |x|=0.903, width=0.104
```

Thus the physical picture is four bonded patches: one mirrored inner pair and
one mirrored outer pair.  Each pair shares one passive shunt branch and one
switch bit.  With an assumed `4 nF` pair capacitance, the reduced resonances
translate to:

```text
pair A    258.99 kHz, Q=6.38, L~94.4 uH, series R~24.1 ohm
pair B    253.41 kHz, Q=4.23, L~98.6 uH, series R~37.1 ohm
```

These are order-of-magnitude circuit anchors, not a BVD-matched terminal
design.  Negative capacitance and cell-wise electronics are absent.

The mechanical entrant uses two symmetric stiffness/clamp zones at normalized
`|x|=0.366` and `0.902`.  It describes the boundary perturbation required from
an actuator; it does not yet specify a manufacturable latch.

The frozen centre-frequency aperture profiles are plotted in
[`frozen_state_aperture_profiles.png`](../results/operator_space/o7_global_state_v0/frozen_state_aperture_profiles.png).

## Positive controls

| Check | PZT | Mechanical | Gate |
|---|---:|---:|---:|
| maximum `sigma(W)` | 0.9800 | 0.9800 | <=1 |
| reflection commutator error | 2.48e-15 | 2.40e-15 | <=1e-10 |
| reciprocity error | 0 | 0 | <=1e-10 |
| minimum state-projector distance | 0.340 | 0.451 | >1e-4 |
| zero-coupling state collapse | 0 | 0 | <=1e-12 |

The four states are therefore not four labels on the same operator.  The
matched `K3Q1` controls also show that the gain is not explained by the seven
frequency samples alone.

## Evidence boundary

The model contains 11 reflection-resolved host modes, passive complex
feedback, a common contractive coupling scale, the same P/S sample propagator,
C3 clutter, pulse covariance and detector replay as the static frontier.  It
does **not** contain an explicit block geometry, piezoelectric constitutive
tensor, electrodes, adhesive, BVD terminals, switch parasitics, actuator
mechanics, finite contact footprint or thermal/manufacturing drift.

The physical states have a large subspace mismatch to the old random ideal
K3Q4 target (`~0.80` normalized projector cost).  Their pass therefore means
that the reduced physical family found a different useful four-state operator;
it is not evidence that the ideal random codebook has been synthesized.

## Next authorized stage: O8

Preserve both branches and do not reopen the O7 test split.

1. `O8-A — coupled PZT`: plane-strain monolithic block with the frozen inner
   and outer patch zones, explicit piezoelectric coupling, electrodes,
   adhesive and the two passive RLC/BVD loads above.  Extract the full
   frequency-dependent `17 x 3 x 4` transfer and replay the frozen detector.
2. `O8-B — mechanical boundary`: the same host and zones with explicit finite
   clamp/contact stiffness, actuator mass and loss.  Extract and replay the
   same operator.

No geometry optimization is allowed on the O8 final test.  If only one branch
can be afforded first, the mechanical branch has the larger O7 statistical
margin; the PZT branch remains the cleaner no-motion electronic architecture
but is near its preregistered Pdl/Pfa gates.

