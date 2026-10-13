# O6g joint-budget and failure-boundary contract

> Frozen after O6f and before O6g calculation. O6g is the final reduced-model
> stage. It may locate engineering ceilings, but it may not select a prototype
> or add another geometry optimization.

## O6g-A: end-to-end q8 passive network

The exact q8 factorization is expanded from the optimistic O6f output lattice
to the complete linear network:

```text
3 public ports
  -> signed 3 x 3 V' transform (3 Givens couplers)
  -> 3 singular attenuator/dump paths
  -> signed 3 x 3 R transform (3 Givens couplers)
  -> 17-output U lattice (45 Givens couplers)
  -> 17 loads
```

This produces `51` calibrated rotation couplers, three attenuator paths and
the required ideal sign inversions. The zero-parasitic network must reproduce
q8 to numerical precision. The singular gains are frozen; they are not
rebalanced after perturbation.

The calculation has three parts:

1. A two-dimensional deterministic frontier for identical insertion loss on
   all 51 couplers and identical excess loss on all three attenuator/dump
   paths. Both are swept from zero through failure; no scalar cost is formed.
2. A joint manufacturing ensemble with independent rotation-angle, coupler
   delay, attenuator log-amplitude and attenuator delay errors. Reference RMS
   values are `0.1 deg`, `2 ns`, `1%` and `2 ns`; `20,000` common draws and a
   one-sided Wilson lower-95 `>=99%` O4-pass gate are used.
3. The same ensemble is repeated at a preregistered loss allocation of half
   each zero-error scalar ceiling. If that allocation fails, the boundary is
   bisected along the ray from the origin. No post-error rebalance or
   contraction is allowed.

The loaded piezo impedance, electrical noise and component-specific FEM/RLC
remain outside this network abstraction. A pass promotes them to O6h; it is
not a fabrication pass.

## O6g-B: motion failure boundary and attribution

The O6f detector reconstruction, 60 signal contexts, physical standard-error
draws and fresh signal/H0 noise draws are reused exactly. Layouts, templates,
thresholds and nuisance training are frozen.

The coupled scale is extended through
`2.25, 2.5, 3.0, 3.5, 4.0`; if needed it doubles until the first joint gate
failure. The first failure is bracketed by twelve bisection steps using the
same standard draws. Every point retains `12,000` aggregate H0 trials.

At the coupled threshold, three one-at-a-time rays are calculated with only
one active error source:

```text
position RMS                 0.05 lambda * scale
contact log-amplitude RMS    5% * scale
contact delay RMS            0.02 cycles * scale
```

Each one-at-a-time ray is independently extended and bisected if it still
passes at the coupled threshold. This separates the limiting mechanism from
the arbitrary common-scale convention.

Frozen gates remain:

```text
Pfa Wilson upper 95%                  <= 1.25%
alpha90 / same-noise scale-zero replay <= 1.10
Pdl at frozen reference amplitude      >= 0.80
```

## O6g decision

- O6g-A reports a feasible joint budget or stops the q8 topology under the
  stated abstraction.
- O6g-B reports finite coupled and one-at-a-time tolerance ceilings. It does
  not claim real couplant robustness beyond the reduced multiplicative/delay
  model.
- No further Born or ideal-coupler parameter sweep is allowed after O6g.
  Surviving q8 risk moves to a loaded RLC/piezo network; surviving motion risk
  moves to a full-wave contact/lift-off gate.
