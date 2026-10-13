# O6a modal gate: direct DCT3 transition plate

**FAIL: the uniform plate is rejected; no FEM geometry is authorized.**

Candidate: `71.6 mm x 3.000 mm` homogeneous Al 6061, plane strain, dense `49`-frequency grid from `193.6` to `290.4 kHz`.

The scalar `m=2` lower-band delay bound permits at most `4.413 mm`; the candidate length was frozen before the elastic solve.

## Worst dense-band metrics

- source-profile fidelity: `0.960987`;
- propagated fidelity: `0.981903`;
- maximum port coherence: `0.068964`;
- maximum column imbalance: `0.3662 dB`;
- maximum differential amplitude RMS: `16.2942%`;
- maximum differential-delay RMS: `4237.480 ns` (diagnostic; the local phase derivative is not mesh-converged at cut-on, but remains above `1 us` in the resolution check);
- minimum useful one-way amplitude: `0.937834`;
- maximum singular value after one global contraction: `0.987637941340`;
- global passive scale: `1.000000`.

Reciprocity residual is exactly zero in this homogeneous modal propagation model; the independent FEM solve must measure it numerically.

The model contains every propagating P/S-mixed free-strip branch. The number of open right-going modes changes from `14` to `21` across the band; the dense gate therefore includes the cut-on risk hidden by a three-mode scalar model.

## Robustness of the no-go

A targeted `120/160/200`-element cross-section check gives minimum source fidelity `0.96103/0.96109/0.96111` and differential amplitude RMS `15.72/15.26/15.27%`. The exact delay spike is sensitive to phase differentiation when a branch opens, so it is retained as a diagnostic and is not the sole rejection argument. The independently converged source-profile and amplitude gates already reject the candidate.

The stopping rule was applied: no 2D geometry and no F1d replay were run. This rejects propagation through the wide uniform free strip, not the frozen K3,Q1 measurement operator. The next allowed mechanism is a load-local direct modal transducer without a laterally guiding Al transition plate.
