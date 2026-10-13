# O3 motion-allowed operator/resource frontier v0

## Verdict

Mechanical translation is a viable low-channel branch, but it is not a free
replacement for the frozen four-channel/two-state operator.  On a completely
fresh deterministic train/validation/test campaign:

- at the optimistic equal-event budget `T = Q = 8`, the best valid mechanical
  layout is the train-only D-optimal eight-position scan with
  `alpha90 = 193.58`; the frozen static operator has `alpha90 = 188.81` and is
  therefore about `2.5%` better in required scattering amplitude;
- the first tested mechanical budget that beats static in aggregate is
  `P = 12`: a compact `5.5 lambda_P0` scan reaches `alpha90 = 183.35`, about
  `2.9%` better than static, and wins in `40/60` paired test contexts;
- the full 17-position monostatic scan reaches `alpha90 = 159.73`, beats
  static in `60/60` contexts, and closes `93.1%` of the diagnostic
  A-scan-to-FMC gap, while using one electronic channel and only monostatic
  samples;
- FMC remains the information oracle at `alpha90 = 134.34`.

Thus the exact branch decision is conditional:

```text
motion allowed, but T <= 8 events
    -> frozen 4ch x 2-state remains the detection leader;

one channel is mandatory and 12--17 stops are affordable
    -> mechanical scan is a credible leading class;

17 stops are affordable
    -> motion clearly beats frozen static, but still does not reach FMC;

nonzero positioning overhead at a fixed T = 8 budget
    -> static wins by a wide margin in every tested time model.
```

This is an operator/resource verdict.  It does **not** authorize lens or probe
geometry synthesis.

## Selected held-out results

All entries below except the A-scan control their independently tested false
alarm upper bound to `<= 1.25%`.

| operator | K | Q / P | span, lambda | Pfa upper 95% | Pd at reference | Pdl at reference | alpha90 | median context ratio to static |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| stationary uniform A-scan | 1 | 1 | 0 | 6.609% (fail) | 0.559 | 0.107 | 504.87 | — |
| frozen 4ch x 2-state | 4 | 8 events | 0 | 1.080% | 0.901 | 0.890 | 188.81 | 1.000 |
| compact motion | 1 | 2 positions | 0.5 | 1.123% | 0.778 | 0.723 | 279.30 | — |
| aperture-sweep motion | 1 | 4 positions | 4.0 | 0.300% | 0.851 | 0.823 | 227.25 | — |
| train D-optimal motion | 1 | 8 positions | 3.5 | 0.484% | 0.895 | 0.891 | 193.58 | 1.019 |
| compact motion | 1 | 12 positions | 5.5 | 0.427% | 0.914 | 0.908 | 183.35 | 0.969 |
| uniform full motion | 1 | 17 positions | 8.0 | 0.630% | 0.951 | 0.948 | 159.73 | 0.830 |
| 17-element FMC oracle | 17 | 17 events | full | 0.884% | 0.974 | 0.971 | 134.34 | — |

`reference alpha = 189.34` is the first point where the new validation-only
FMC curve reaches `Pd >= 0.95`.  Lower `alpha90` is better.  The context ratio
is `alpha90_motion / alpha90_static`; below one favours motion.

The stationary A-scan threshold did not transfer to the held-out nuisance
subspace (`Pfa upper = 6.61%`).  Consequently the reported PAUT-gap closure is
a conservative diagnostic coordinate, not a formally equal-Pfa claim.  The
direct static-versus-motion conclusions above remain comparable because all
operators used in those conclusions pass the held-out Pfa gate.

## Layout map

The campaign freezes 35 operators before the signal test: three fixed
comparators and 32 unique mechanical layouts.  The mechanical set spans:

- `P = 1, 2, 3, 4, 6, 8, 12, 17` positions;
- full-aperture uniform placement;
- compact central placement on the native `0.5 lambda` grid;
- edge-biased Chebyshev placement;
- explicit `2/4/6/8 lambda` aperture sweeps at `P = 4` and `P = 8`;
- a nested D-optimal sequence selected only from the training dictionary and
  training nuisance subspace.

The train-only D-optimal criterion does not universally win held-out
`alpha90`.  It is best in the tested `P = 8` class, but compact placement is
better at `P = 12`.  This is useful negative evidence: a surrogate
log-determinant objective is not yet a sufficient scan-design target.

## Energy and resource contract

Every mechanical event is a collocated single-element Tx/Rx measurement.
Concatenated `P`-position features are divided by `sqrt(P)`, so total incident
energy is one for every scan.  The existing multiport and FMC dictionaries use
their corresponding unit-total-energy normalizations.

Every row carries the explicit vector

```text
(K electronic channels,
 Q transmit events,
 P mechanical positions,
 P-1 moves,
 aperture span,
 complex feature count,
 total incident energy)
```

No single mechanical-time constant is asserted.  Four dimensionless
sensitivity models are exported:

```text
electronics_only       T = Q
one_event_per_settle   T = Q + moves
settle_plus_travel     T = Q + moves + span/lambda
slow_motion            T = Q + 4*moves + 2*span/lambda
```

At the frozen static time budget `T <= 8`, the best Pfa-valid mechanical
choices are:

| time model | best allowed mechanical operator | T | alpha90 |
|---|---|---:|---:|
| electronics only | train D-optimal, P=8 | 8.0 | 193.58 |
| one event per settle | 4-position, 4-lambda aperture | 7.0 | 227.25 |
| settle plus travel | compact P=2 | 3.5 | 279.30 |
| slow motion | compact P=2 | 7.0 | 279.30 |

`resource_frontiers.csv` contains both a time-only Pareto flag and a
time-plus-channel Pareto flag.  It also gives gap closure per unit time and
per channel-time product.  These are sensitivity coordinates, not BOM prices;
real positioning latency and channel cost must eventually replace the
dimensionless weights.

## Fresh split contract

No signal, nuisance atom, threshold draw or context from
`paired_clutter_born_v0` is read.  The O2c codebook is frozen input, not test
data.  The common campaign base seed is `0x04f0a117`:

```text
nuisance train       base + 0x101
nuisance validation  base + 0x202
nuisance test        base + 0x303
signal validation    base + 0x404
signal test          base + 0x505
H0 validation        base + 0x606 + architecture offset
H0 test              base + 0x707 + architecture offset
```

The nuisance contract is unchanged for comparability: full-FMC relative mean
powers are backwall `2`, ringing `2`, and material clutter `3`, with the same
atom definitions and counts as reduced O3.  Thresholds use 10,000 validation
H0 trials; Pfa uses 10,000 independent test H0 trials.  The FMC reference uses
24 validation signal contexts.  Final curves use 60 new C3 contexts, all 136
ROI hypotheses, three cycling scattering tensors and common latent noise for
paired architecture comparisons.

## Reproduction and artifacts

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_operator_o3_motion_frontier.jl
```

- `split_contract.csv` — exact seeds and non-reuse declaration;
- `architecture_metadata.csv` — complete resource vectors and scan indices;
- `greedy_trace.csv` — training-only nested position selection;
- `nuisance_scales.csv` — preserved `2/2/3` full-FMC nuisance powers;
- `thresholds.csv` — validation threshold and independent Pfa audit;
- `detection_curves.csv` — aggregate fresh-test `Pd/Pdl` curves;
- `context_summary.csv` — paired context-level metrics;
- `motion_summary.csv` — architecture-level performance and resources;
- `pairwise_vs_static.csv` — held-out context ratios and win fractions;
- `resource_frontiers.csv` — time-only and time/channel Pareto analyses;
- `motion_frontier.jld2` — compact reproducible result bundle;
- `motion_alpha90_vs_positions.png`, `motion_resource_frontiers.png` — visual
  summaries.

## Limits

This remains a 2D reciprocal elastic Born surrogate with ideal point scan
positions, no couplant reset penalty, no encoder error, no continuous-motion
sampling, no probe footprint and no full-wave backwall interaction.  The
mechanical-time models are deliberately parametric.  The 12-position
crossover is a result inside this discrete virtual aperture, not an industrial
POD or a universal optimum.
