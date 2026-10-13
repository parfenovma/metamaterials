# O6i mitigation-class contract

> Frozen after O6h. O6h blocked intrinsic-Q q8 without broadband matching and
> blocked liquid-coupled motion under 2D angular P/S physics. O6i may test one
> named mitigation class per branch; it may not reopen layout or detector
> optimization.

## O6i-A: q8 damping no-go and remaining matching requirement

The intrinsic normalized BVD load has `Q_i=35`; O6h requires effective
`Q_e<=1.481481` for 10 dB worst-band public return loss. For a series-resonant
load broadened only by added resistance,

```math
R_total/R_intrinsic = Q_i/Q_e,
eta_load = R_intrinsic/R_total = Q_e/Q_i.
```

The delivered useful-amplitude upper bound is `sqrt(eta_load)`. It is compared
with the frozen O4 useful-amplitude gate `0.8994`; matching-network and q8
losses are optimistically set to zero. Failure stops resistive damping as a
remedy. It does not stop lossless multi-resonant matching or a genuinely
multi-resonant piezo stack.

Classical [Fano theory](https://www.sciencedirect.com/science/article/pii/0016003250900068)
establishes integral bandwidth/match constraints for a
prescribed passive load and lossless matching network, but O6i does not apply
a generic single-pole formula to the current normalized series-BVD load. A
lossless-matching claim requires a concrete load rational function and network
topology at the next gate.

## O6i-B: shear-transmitting solid contact map

The exact O6h-C 2D boundary solver is generalized from a fluid layer to a
welded isotropic solid layer carrying P and SV waves in both directions. It
enforces continuity of both velocity components and both tractions at both
interfaces and includes all layer reflections.

The map is dimensionless and upper-level:

```text
ZP_layer / sqrt(ZP_probe ZP_sample) = 0.10, 0.25, 0.50, 0.75, 1.00
cP_layer / cP_Al                    = 0.25, 0.50, 0.75, 1.00
cS_layer / cP_layer                 = 0.10, 0.25, 0.50
d / lambdaP_layer                   = 0.0025, 0.005, 0.01, 0.02
```

Density is derived from `rho=ZP/cP`; only `800--8000 kg/m3` is retained. Each
candidate is checked at the seven operator frequencies and P-incidence angles
`0,15,30,45,55 deg` against the same direct welded interface.

Positive control: maximum lossless energy-balance error `<=1e-8`.
Qualification over `0--45 deg`:

```text
minimum layer/direct transmitted power >= 0.80
maximum change in transmitted S fraction <= 0.10
```

The output is a Pareto map over impedance ratio, wave-speed ratios, thickness,
density, minimum transmission and mode-conversion change. Passing proves that
a shear-supporting contact class exists; it does not yet prove repeatable
attachment or detector-level Pfa.

If the pass set is nonempty, a robust interior handoff is selected by a frozen
rule: require the stricter margins `power>=0.85` and `Delta S<=0.075`, exclude
the extreme impedance ratios `0.10` and `1.00`, then maximize physical layer
thickness. Ties minimize the absolute log-distance of density from
`sqrt(800*8000) kg/m3`, then use lexicographic parameter order. This favours a
resolvable layer with margin, not the thinnest numerical approximation to a
welded interface.

## Decision

- If resistive damping fails, q8 remains blocked on a specified lossless or
  multi-resonant matching topology.
- If no solid layer passes, the current moving-contact architecture stops.
- If solid layers pass, freeze one robust interior point (not the numerical
  boundary) for detector-level and attachment-tolerance replay in O6j.
