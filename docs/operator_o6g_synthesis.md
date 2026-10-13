# O6g synthesis and O6h promotion

O6g closes the final reduced-surrogate uncertainty on both surviving O6
branches. The q8 result now covers the full ideal linear network; the motion
result has finite coupled and one-at-a-time failure boundaries.

## q8 network

The complete factorization contains 51 signed Givens rotation couplers and
three singular attenuator/dump paths. It reconstructs q8 with relative error
`6.70e-16`. At half of each independent scalar loss ceiling, the network still
has a nonzero 99%-Wilson tolerance budget:

```text
loss per rotation coupler      0.022537 dB
excess loss per attenuator     0.456357 dB
coupler angle RMS              0.11472 deg
coupler delay RMS              2.29443 ns
attenuator log-amplitude RMS   1.14722%
attenuator delay RMS           2.29443 ns
```

This is a pass of a reciprocal matched-network abstraction. It does not yet
include reflection from capacitive/resonant piezo loads, matching-network
bandwidth or electrical-noise gain.

## Motion

The coupled encoder/contact failure boundary is scale `2.35522` for P12 and
`2.23883` for P17. The first failed gate is relative alpha90. One-at-a-time
boundaries show the ordering

```text
position error < contact delay < contact amplitude
```

for both layouts, while the coupled limit is lower than every scalar limit.
Attribution therefore rejects a scalar contact-quality specification.

## Decision

Both branches remain nondominated and both have positive reduced physical
margin. More reduced sweeps would only refine abstractions already shown to
pass. O6h must raise evidence:

- q8: reciprocal scattering completion plus loaded broadband piezo/BVD
  terminals and matching;
- motion: exact layered-contact/lift-off transfer including multiple
  reflections, followed by a 2D full-wave spot check only if it passes.

No prototype architecture is selected at O6g.
