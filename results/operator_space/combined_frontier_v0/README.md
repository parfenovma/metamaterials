# Combined static/motion operator frontier v0

## Outcome

The motion-allowed and motion-forbidden branches were mapped independently on
the same fresh deterministic split contract. The union contains `51`
architecture points and does **not** have one universal winner. The leading
class changes with electronic-channel, event and motion budgets.

The main landmarks are:

| event budget | best Pfa-valid static | alpha90 | best Pfa-valid motion | alpha90 | full reference gate |
|---:|---|---:|---|---:|---|
| 2 | `K1,Q2` even state pair | 267.6 | compact `P2` | 279.3 | neither |
| 4 | `K1,Q4` static | 231.3 | four-position, `4 lambda` span | **227.3** | neither |
| 8 | `K4,Q2` unconstrained | **185.7** | train-D-optimal `P8` | 193.6 | static only |
| 9 | `K3,Q3` unconstrained | **183.9** | same `P8` | 193.6 | static only |
| 12 | `K3,Q4` balanced sectors | **168.6** | compact `P12` | 183.4 | both |
| 17 | `K3,Q4` balanced sectors | 168.6 | full monostatic `P17` | **159.7** | both |

Lower `alpha90` is better. “Full reference gate” here means held-out
`Pfa95 <= 1.25%`, `Pd(reference) >= 0.90` and
`P(detect and localize) >= 0.80`; each branch report retains its additional
context and localization audits.

The resulting design map is:

```text
few events, motion cheap, one channel mandatory
    -> motion and one-channel multistate coding are near peers;

8--12 events, 3--4 channels allowed
    -> static state diversity is the quality leader;

12 events, one channel mandatory
    -> compact motion is the first fully valid motion point;

17 positions and motion overhead acceptable
    -> full monostatic scan is the best low-channel reduced operator;

full electronic aperture acceptable
    -> FMC remains the information oracle.
```

## Strongest points in each branch

### Static

The strongest valid static endpoint is

```text
3 electronic channels x 4 balanced even/odd passive states
12 transmit events
alpha90 = 168.63
Pfa upper 95% = 0.862%
Pd(reference) = 0.937
Pdl(reference) = 0.928.
```

Two lower-overhead static co-leaders remain:

- `K3,Q3`, 9 events, `alpha90=183.95`;
- `K4,Q2`, 8 events, `alpha90=185.68`.

The apparent `K4,Q4` value `alpha90=147.97` is rejected because its held-out
false-alarm upper bound is `3.21%`. It is not part of the valid frontier.

### Motion

The first tested motion layout to pass the reference decision gate is compact
`P12`, with `alpha90=183.35`. The full `P17` monostatic scan reaches
`159.73`, beating the best tested static point by `5.3%`, but requires sixteen
moves and five more events.

At the important equal-event budget of eight, the best motion point has
`alpha90=193.58` and narrowly misses `Pd=0.90`, while static `K4,Q2` reaches
`185.68` and passes. Motion therefore does not dominate merely because it won
the earlier single-comparator O3 run.

## Why both branches remain open

The global resource vector is

```text
(electronic channels K,
 transmit events,
 passive states and switches,
 mechanical positions and moves,
 aperture travel,
 calibration coefficients,
 alpha90 and decision statistics).
```

Static `K3,Q4` and motion `P17` are mutually nondominated:

- static uses three electronic channels, four calibrated states and no motion;
- motion uses one electronic channel, seventeen events and sixteen moves;
- motion is `5.3%` better in `alpha90`;
- static is substantially better when moves/travel have nonzero time cost.

The exported sensitivity analysis therefore uses four time conventions rather
than inventing one BOM price:

```text
events_only
balanced_overhead
motion_expensive
state_expensive.
```

Under `events_only`, both branches retain Pareto points because motion trades
quality for one-channel electronics. Under balanced or motion-expensive time,
static state diversity occupies most of the high-quality low-time frontier.
Under state-expensive assumptions, motion gains territory. The branch cannot
be collapsed further without a real product constraint.

## Cross-campaign comparability

Both subcampaigns use:

- base seed `0x04f0a117` and the same derived split schedule;
- independent fresh train, validation and final-test data;
- `10,000` validation and `10,000` test H0 trials;
- `24/60` validation/test signal contexts;
- the same 15-point amplitude grid;
- full-FMC nuisance powers `backwall=2`, `ringing=2`, `material=3`;
- the same five-cycle pulse, covariance, ROI, perturbations and GLRT family;
- unit total incident energy across positions or passive states.

Independent implementations reproduce the common controls to about one
percent:

```text
FMC alpha90             134.34 motion / 133.82 static   (0.39%)
frozen K4,Q2 alpha90    188.81 motion / 190.20 static   (0.74%).
```

This residual Monte Carlo/interpolation variation is small relative to the
main branch differences. Exact paired claims should nevertheless use each
branch’s context-level tables, not subtract rounded headline values.

The fresh uniform A-scan fails held-out Pfa transfer in both campaigns.
Therefore A-scan-to-FMC closures remain diagnostic frozen-threshold proxies.
The direct static/motion landmarks above all state whether their own Pfa and
reference gates pass.

## Search coverage

The static branch evaluates `2,304` candidate codebooks, `192`
nuisance-aware finalists and one frozen winner in every one of sixteen
`(K=1..4,Q=1..4)` cells. It covers balanced reflection sectors, pure even,
pure odd and unconstrained real reciprocal port subspaces.

The motion branch freezes `32` unique layouts over
`P=1,2,3,4,6,8,12,17`, including full-uniform, compact, Chebyshev, aperture
sweeps and train-only D-optimal sequences.

This is a broad empirical map, not a proof that no unsampled codebook or scan
trajectory can improve it.

## Stage verdict

- parallel reduced static/motion mapping: **complete**;
- static intrinsic O3 entrants: **three valid points**;
- motion intrinsic O3 entrants: **valid from P12 upward in the tested map**;
- universal branch winner: **none**;
- full-wave paired confirmation: **pending**;
- O4 transfer target and physical geometry: **blocked**.

The correct next physics step, if both branches remain open, is a paired
full-wave comparison of no more than four frozen representatives:

```text
static K4,Q2  — minimum event count among valid static points
static K3,Q4  — best valid static alpha90
motion P12    — first valid motion crossover class
motion P17    — best tested low-channel alpha90.
```

That calculation would still compare operator classes. It would not yet
optimize a lens, switching network or scanner geometry.

## Reproduction and artifacts

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_operator_o3_motion_frontier.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/run_operator_o3_static_frontier.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/synthesize_operator_frontiers.jl
```

- `combined_architectures.csv` — normalized union of all points;
- `combined_resource_frontiers.csv` — four cost sensitivities and Pareto flags;
- `event_budget_envelope.csv` — best branch point at common event budgets;
- `combined_quality_frontier.png` — quality/event map;
- `combined_cost_sensitivity.png` — branch movement under overhead models.

Detailed search audits remain in the sibling `motion_frontier_v0/` and
`static_frontier_v0/` directories.

## Scope limitation

All results remain 2D reciprocal elastic Born surrogates. They omit real
couplant reset time, finite probe footprint, encoder error, switching
transients, transducer electromechanics, finite cracks and nonlinear contact.
They map which expensive full-wave comparisons are justified; they are not
industrial POD and do not authorize geometry synthesis.

