# O6b finite PZT/bond/Al mechanical twin

**FAIL: the finite direct-coupler twin does not satisfy O4.**

One frozen `6.0 mm PZT-5A mechanical surrogate -> 0.10 mm epoxy -> 35 mm matched Al load` geometry was solved at 13 frequencies. The three right-hand sides are continuous DCT0--2 tractions on the common PZT back face; the 17-bin output is sampled at the bonded Al interface. No length, thickness or profile sweep was performed.

## Worst dense-band metrics

- transfer fidelity: `0.973307`;
- maximum port coherence: `0.044988`;
- column imbalance: `0.2305 dB`;
- differential amplitude RMS: `19.589%`;
- differential-delay RMS: `846.357 ns`;
- minimum useful amplitude: `0.882498`;
- maximum sigma: `0.935687`;
- mechanical reciprocity residual: `1.312e-15`.

The common power-wave scale is inherited from the preregistered O6b PZT/bond/Al bound and includes nominal PZT/bond propagation loss. The FEM contributes only the finite-aperture differential shape, phase and cross-port terms; common displacement compliance is a calibrated gauge.

This remains a mechanical distributed-port twin, not an electrical terminal model. A pass authorizes F1d replay first; electrode topology, capacitance and receive electronics remain blocked until that replay.

## Mesh check and decision

The three-anchor refinement from `PZT/Al h=0.60/1.00 mm` to
`0.45/0.75 mm` changed worst fidelity from `0.97331` to `0.97262` and worst
amplitude RMS from `19.59%` to `19.84%`. The failure is therefore not a mesh
artifact. Detailed values are stored in `convergence_anchors.csv` and
`convergence_verdict.md`.

The stopping rule was applied: no electrical model and no F1d replay were
run. The ideal load-local upper bound remains positive, so the next allowed
mechanism is one laterally decoupled homogenized `1-3` composite layer. A
thickness or geometry sweep of the homogeneous PZT slab is not authorized.
