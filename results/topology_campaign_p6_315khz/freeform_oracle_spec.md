# Freeform oracle specification (wave 2)

This branch starts when a concurrency slot becomes available.  It is an upper
bound/search instrument, not automatically a fabrication candidate.

## Design domain

- Same aluminium device footprint, aperture, output domain and target as the
  common benchmark contract.
- Mirror symmetry about `y=0` is mandatory in the first search.
- Binary aluminium/void representation with through-cut voids.
- Initial void fraction bracket: 5--30%.

## Regularization

- Morphological/filter radius no smaller than 0.4 mm, yielding a nominal
  0.8 mm minimum feature/web after projection.
- Reject disconnected aluminium islands and enclosed features inaccessible to
  a through-cut manufacturing process.
- Penalize boundary perimeter and the number of disconnected void components;
  both quantities must be reported rather than hidden inside a score.

## Search hierarchy

1. Low-dimensional symmetric level-set/shape grammar: curved slots, circular
   holes and connected channels.  Use a space-filling DOE and a surrogate to
   avoid remeshing thousands of candidates.
2. Multiobjective screening at 305, 315 and 321.7 kHz using `(J1,J2_equal)`;
   transmission and focus error are explicit constraints.
3. Fine-mesh validation of the best feasible Pareto point.
4. Geometric simplification into arcs/holes, followed by a second fine solve.

## Robustness objective

In addition to single-frequency Pareto values, report

`robust_joint = minimum(min(J1(f),J2(f)) for f in [305,315,321.7])`.

This prevents a resonant needle from being mistaken for the campaign winner.

