# Seven-segment active-aperture upper bound at 315 kHz

## Bottom line

A seven-segment source adds a useful independent control direction.  On the
fine mesh at 315 kHz, the passive `c005` plate with uniform drive gives
`J1=1.2057`, `J2=1.3859`.  Equal-incident-energy complex weights, constrained
to amplitude dynamic range at most 4, raise this to `J1=1.3195`, `J2=1.6844`.
Moving along the complex-weight Pareto set can raise `J2` further to `1.7661`
while retaining `J1=1.2783`.  This is the requested partial decoupling: a
4.85% increase in J2 from the focus-matched active point costs 3.12% in J1.

An explicit target-lock re-optimization on the same response matrices also
removes the axial near-miss.  With mirror-symmetric weights (only four
independent channel settings), `c005` reaches `J1=1.2531`, `J2=1.6581`, exit
energy `0.7004` and focus error `2 mm`: the campaign prototype gate passes.
The corresponding symmetric straight aperture reaches `J1=1.1685`,
`J2=1.3587`, exit energy `1.0555` and focus error `3 mm`, so it also passes the
strict gate and is the simplest active prototype.

This branch is an **active upper-bound benchmark**, not a coupled
piezoelectric FEM.  The model applies seven independently weighted normal
tractions to equal source-boundary segments and superposes their linear elastic
responses.  It therefore answers how much controllable aperture weights could
buy before introducing uncertain PZT, bond-layer, electrode and electrical
impedance parameters.

## Fine-mesh results at 315 kHz

All weights satisfy

`sum(abs2(w)) = 7`,

so every row has the same prescribed incident-traction energy as uniform drive.
J1 and J2 are normalized to the fine solid-straight footprint with uniform
drive at the same frequency.

| footprint and drive | J1 | J2 equal | J2 clutter-dom. | exit energy | FWHM | focus error | fixed-weight robust min joint* | amplitude DR |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| straight, uniform | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 38.92 mm | 20.0 mm | 1.0000 | 1.00 |
| straight, amplitude-only | 1.1458 | 1.2626 | 1.2351 | 1.0313 | 38.03 mm | 14.0 mm | 1.0585 | 4.00 |
| straight, complex focus matched | 1.2006 | 1.3513 | 1.3047 | 0.9750 | 39.09 mm | 14.04 mm | 1.1016 | 3.19 |
| c005, uniform | 1.2057 | 1.3859 | 1.3496 | 0.7090 | 32.52 mm | 11.0 mm | 1.1527 | 1.00 |
| c005, amplitude-only | 1.2535 | 1.4711 | 1.4194 | 0.7547 | 33.00 mm | 11.0 mm | 1.2052 | 2.27 |
| c005, phase-only focus matched | 1.2684 | 1.5807 | 1.5644 | 0.7395 | 30.19 mm | 10.0 mm | 1.1165 | 1.00 |
| **c005, complex focus matched** | **1.3195** | **1.6844** | **1.6530** | **0.7695** | **30.76 mm** | **8.0 mm** | 1.1420 | 2.26 |
| c005, complex CNR end | 1.2783 | **1.7661** | 1.8625 | 0.7237 | 27.75 mm | 6.0 mm | 1.0155 | 3.50 |
| **c005, complex robust** | **1.2889** | **1.5783** | **1.5340** | **0.7684** | 32.21 mm | 9.0 mm | **1.2268** | 2.80 |
| **straight, symmetric target-lock (strict)** | 1.1685 | 1.3587 | 1.3548 | 1.0555 | 45.60 mm | **3.0 mm** | 1.0523 | 3.99 |
| **c005, symmetric target-lock (prototype)** | **1.2531** | **1.6581** | **1.7185** | **0.7004** | 30.07 mm | **2.0 mm** | 0.8821 | 3.96 |

\* Minimum of `min(J1,J2)` when the same 315-kHz complex weights are evaluated
on the coarse response matrices at 305, 315 and 321.7 kHz.  It is a harmonic
carrier robustness screen, not a transient or group-delay result.

Unconstrained Pareto weights improve both primary metrics but initially leave a
6--14 mm axial error.  Target-lock enforces that the target amplitude exceeds
every point outside a 5 mm radius, with exit energy at least 0.70 and DR at most
4.  It closes the focus gate, but the `c005` solution sits directly on the
energy constraint and its fixed-weight robust minimum is only 0.8821.  Retuning
the passive delay profile around 315 kHz therefore remains the preferable route
to a less narrowband prototype.

## What the independent controls do

- **Amplitude only** is the simplest hardware benchmark.  On `c005` it already
  adds 4.0% to J1 and 6.1% to J2 over uniform `c005`, while raising the exit
  proxy above 0.75.
- **Phase only** improves spatial rejection more strongly than amplitude only:
  J2 rises to 1.5807 and FWHM falls to 30.19 mm, with unit dynamic range.
- **Complex focus-matched weights** maximize target amplitude under the
  dynamic-range bound.  They give the best J1 and retain exit energy 0.7695.
- **Complex CNR-oriented weights** sacrifice a small amount of target amplitude
  to reduce the reduced-Born clutter term.  This is the clearest measured
  separation of the two objectives.
- **Robust complex weights** trade peak 315-kHz performance for the largest
  worst-frequency joint gain across 305--321.7 kHz.
- **Symmetric target-lock weights** reduce seven segment values to four
  independent settings and force the local peak inside the 5 mm focus gate.
  Straight is the simplest strict-feasible active prototype; `c005` is the
  prototype-feasible maximum-gain option.

The fine `c005` focus-matched weights from segment 1 (bottom) to 7 (top) are:

| segment | amplitude | phase (deg) |
|---:|---:|---:|
| 1 | 0.8748 | -12.71 |
| 2 | 0.7001 | -2.44 |
| 3 | 1.0057 | -43.72 |
| 4 | 1.5721 | 0.00 |
| 5 | 1.0095 | -43.67 |
| 6 | 0.6961 | -2.65 |
| 7 | 0.8707 | -12.37 |

The near mirror symmetry is an independent sanity check; symmetry was not
imposed during optimization.

For the prototype-gated `c005` solution, mirror symmetry is imposed explicitly:

| segment | amplitude | phase (deg) |
|---:|---:|---:|
| 1 / 7 | 0.8551 | -39.87 |
| 2 / 6 | 0.3693 | -46.96 |
| 3 / 5 | 1.2497 | -57.43 |
| 4 | 1.4634 | 0.00 |

## Frequency screen and convergence

Coarse response matrices were computed for both footprints at 305, 315 and
321.7 kHz.  The uniform `c005` results were:

| frequency | J1 | J2 | exit energy | FWHM | focus error |
|---:|---:|---:|---:|---:|---:|
| 305.0 kHz | 1.1729 | 1.2530 | 0.7565 | 37.43 mm | 7.07 mm |
| 315.0 kHz | 1.1671 | 1.2854 | 0.7344 | 34.85 mm | 12.0 mm |
| 321.7 kHz | 1.1527 | 1.3360 | 0.6330 | 29.43 mm | 12.0 mm |

Fine response matrices were then computed at 315 kHz: 39,656 algebraic DOF
for straight and 120,808 for `c005`.  Summing the seven fine basis fields with
uniform weights reproduces the earlier monolithic uniform-drive solutions with
relative L2 errors `9.55e-14` (straight) and `3.67e-13` (`c005`); the target
amplitude ratios are one to approximately 13 significant digits.

## Definitions and limitations

- Seven segments span the 71.6 mm aperture; each is 10.2286 mm high.
- Segment discontinuities are applied as piecewise-constant traction fields on
  the existing source boundary.  No remeshing or fictitious PZT constants are
  introduced.
- One FEM matrix is assembled and factored per footprint/frequency.  Seven
  source right-hand sides produce the response matrix; subsequent field and
  weight sweeps use exact linear superposition.
- J2 uses the campaign reciprocal-Born kernel `K=ux^2`, with clutter evaluated
  in the fixed spatial window.  It remains a reduced harmonic proxy, not a
  measured pulse-echo CNR.
- The side-lobe diagnostic is target amplitude divided by the largest amplitude
  at the nominal target plane with `|y| >= 12 mm`.  It is reported explicitly
  and is not hidden inside the optimization score.
- Hardware needs seven independently driven electrode/channel voltages and a
  phase-calibrated transfer function.  Real piezoelectric coupling, voltage and
  reactive-power limits, bond layers, ringing and channel crosstalk remain for
  a later coupled/transient model.

## Artifacts

- `data/campaign_metrics.csv`: campaign-contract-compatible combined coarse and
  fine table.  Required common columns come first; source constraints and
  robustness diagnostics follow.
- `data/optimized_weights_coarse.csv` and
  `data/optimized_weights_fine.csv`: complex weight for every segment and every
  candidate.
- `data/response_*.jld2`: reusable complex response matrices.
- `figures/pareto_315_fine.png`: active J1/J2 Pareto comparison.
- `figures/optimized_weights_fine.png`: amplitude/phase heatmaps.
- `figures/target_locked_weights_fine.png`: weights of the focus-gated
  prototype Pareto set.
- `figures/frequency_robustness_fine.png`: fixed-weight frequency screen.
- `figures/field_*_315_fine.png`: superposed fine field maps.

## Reproduction

```sh
julia --startup-file=no --project=. segmented_piezo_response.jl --basis straight 315 fine
julia --startup-file=no --project=. segmented_piezo_response.jl --basis c005 315 fine
julia --startup-file=no --project=. segmented_piezo_response.jl --report-coarse
julia --startup-file=no --project=. segmented_piezo_response.jl --report-fine
julia --startup-file=no --project=. segmented_piezo_response.jl --merge
```
