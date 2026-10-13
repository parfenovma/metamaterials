# Sinusoidal Al waveguides: carrier time-reversal audit

> Frozen before reading the audit outputs on 4 October 2026.

## Question

Does the frozen 15-strip sinusoidal aluminium lens create more useful spatial
channels than the matched straight-strip aperture, or does it improve its one
design focus by redistributing and phasing essentially the same channel set?

This audit is diagnostic only. It does not reopen the stopped gap sweep, tune
the geometry, or replace the existing broadband and defect results.

## Frozen evidence

The only input is the saved carrier response matrix in
`results/sinusoidal_material_transfer_p6_242khz/response_matrix/base/`.
It was produced by eight independent full-FEM source solves at `242 kHz` for
both the selected sinusoidal lens and its straight control.

Each column excites one mirror-symmetric source group. Each row samples
longitudinal displacement at one non-negative transverse coordinate,
`|y| = 0, 8.2, ..., 57.4 mm`, on the plane `x = 33 mm`, five millimetres after
the strips. The resulting `8 x 8` matrix represents only the reflection-even
sector. An off-axis time-reversal target therefore means the symmetric pair
`(+y,-y)`, not a signed point.

## Frozen calculations

For each geometry and normalization, the audit computes

```math
H=U\Sigma V^\dagger,
\qquad
K=HH^\dagger.
```

The singular spectrum, stable rank, entropy rank, condition number, observable
mode counts at `20/30/40/60 dB`, normalized time-reversal profiles, and
spatial-coherence matrix are reported.

Two source conventions are retained:

1. `as_solved`: the original `1 MPa` traction coefficient for each group;
2. `nominal_area`: columns for the seven two-strip groups are divided by
   `sqrt(2)`, a nominal equal-source-area energy proxy relative to the single
   centre strip.

The FEM files do not contain accepted input work for each basis solve, so
`nominal_area` must not be described as exact equal accepted energy. A physical
power-normalized channel comparison would require saving complex source work
or port power in a new solve. The main conclusion is considered robust only
if it has the same sign under both conventions.

## Interpretation boundary

- This is a single-frequency carrier test; it cannot measure temporal
  compression, modal lifetime, Q, or ringing.
- Eight point samples do not resolve sub-bin focal width.
- The current basis contains no odd source groups, so it cannot distinguish
  `+y` from `-y`.
- The time-reversal kernel describes what the stored multi-channel basis could
  do under phase-conjugate control. The existing bonded uniform transducer
  does not independently drive those eight groups.
- A broadband time-reversal claim requires an `8 x 8` response matrix over the
  already frozen 17-frequency grid, including source-work normalization.
