# Motion-forbidden static operator frontier O3 v0

## Outcome

The fresh reduced `C3` campaign leaves **three** static architectures through
all intrinsic gates.  The best valid point is not the four-channel entrant
from O2c, but a reflection-sector design with three electronic channels and
four sequential passive states:

| architecture | K | Q | Tx events | calibration coefficients | Pfa upper 95% | Pd(ref) | Pdl(ref) | alpha90 | frozen-threshold FMC-gap proxy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `trained_k3_q4_reflection_balanced` | 3 | 4 | 12 | 252 | **0.862%** | **0.937** | **0.928** | **168.6** | **0.903** |
| `trained_k3_q3_unconstrained` | 3 | 3 | 9 | 189 | 0.906% | 0.905 | 0.897 | 183.9 | 0.860 |
| `trained_k4_q2_unconstrained` | 4 | 2 | 8 | 224 | 1.167% | 0.901 | 0.889 | 185.7 | 0.856 |
| frozen O2c `K4,Q2` | 4 | 2 | 8 | 224 | 1.221% | 0.896 | 0.885 | 190.2 | 0.843 |
| fresh uniform A-scan | 1 | 1 | 1 | 7 | 6.814% | 0.560 | 0.109 | 493.0 | 0.000 |
| fresh FMC oracle | 17 | 1 | 17 | 2023 | 1.015% | 0.973 | 0.970 | 133.8 | 1.000 |

Lower `alpha90` is better.  The reference amplitude is `186.239`, selected
before final test as the first amplitude where FMC reaches `Pd >= 0.95` over
24 validation contexts.  The three passing designs simultaneously satisfy:

```text
test Pfa upper 95%                         <= 1.25%
Pd at the FMC-derived reference amplitude >= 0.90
P(detect + localize <= 0.5 lambda_P)       >= 0.80
Q90 localization error                    <= 1.0 lambda_P
Q10 context FMC-gap proxy                 >= 0.50
valid context-gap coverage                 >= 90%
```

The frozen O2c entrant narrowly controls false alarms but misses the detection
gate (`Pd=0.896`).  A newly selected design with the same `K=4,Q=2` resources
passes and reduces `alpha90` from `190.2` to `185.7`.  This is a modest but
real fresh-split gain; it is not a retuning of the old final test.

## What the frontier says

The resource vector is kept explicit instead of collapsed into an arbitrary
price score:

```text
K = simultaneous electronic channels
Q = sequential passive states
Tx events = K*Q
state switches = Q-1
calibration coefficients = 7*K^2*Q complex values
mechanical positions = 1 for every entrant
```

Every state uses a real orthonormal port matrix, hence is reciprocal and
contractive in the reduced model.  The incident-energy normalization is
`1/sqrt(K)` inside a transmit basis and `1/sqrt(Q)` across states.  More
states therefore add diversity, not fictitious energy.

The feasible nondominated set is broader than the three gate passers because
low-resource points can control false alarms without reaching the reference
`Pd/Pdl` gates.  Its main landmarks are:

```text
K1,Q2  alpha90=267.6   cheapest useful state diversity
K1,Q4  alpha90=231.3   one-channel multistate endpoint
K2,Q4  alpha90=199.8   two-channel endpoint
K3,Q3  alpha90=183.9   first full intrinsic pass
K3,Q4  alpha90=168.6   best full intrinsic pass
K4,Q2  alpha90=185.7   fewest events among full intrinsic passes
```

This is a genuine three-way engineering choice, not a single winner:

- `K3,Q4` minimizes scattering amplitude among statistically valid entrants;
- `K3,Q3` removes one state and 63 calibration coefficients for an `8.3%`
  `alpha90` penalty;
- `K4,Q2` uses the fewest Tx/state events and one state switch, at the price
  of one additional electronic channel and 35 more calibration coefficients
  than `K3,Q3`.

The apparent `K4,Q4` result (`alpha90=148.0`) is **not admissible**: its held-out
Pfa upper bound is `3.21%`.  It is evidence of background over-specialization,
not a better architecture.

## Symmetry finding

Four real static families were searched in every one of the 16 `(K,Q)` cells:

1. balanced even/odd reflection sectors;
2. even sector only;
3. odd sector only;
4. an unconstrained real orthonormal subspace.

There were 36 independent state sets per family and cell.  The best three by
nominal operator score advanced to the nuisance-aware validation proxy; only
one winner per `(K,Q)` was then frozen for final test.  Losing families never
saw the final-test split.

Balanced sectors won 7 cells, unconstrained states 6, even-only states 3, and
odd-only states none.  Pure even/odd designs retain unit nonlocal coherence in
many cells: respecting an irrep is not sufficient if the complementary irrep
is discarded.  At the useful high-information end, either mixed even/odd
sectors or unconstrained states are needed.  The best valid endpoint is the
analytically liftable balanced-sector class.

## Frozen split and statistical protocol

This campaign does not reuse the old `paired_clutter_born_v0` final test.
Its common base seed is `0x04f0a117`; all derived seeds are in
`seed_contract.csv`.  The exact counts and physics contract are in
`campaign_contract.csv`:

- independent train, threshold-validation and final-test nuisance libraries;
- the unchanged full-FMC relative powers `backwall=2`, `ringing=2`,
  `material=3`;
- unchanged atom definitions and counts: 17 backwall, 6 ringing and 48
  material atoms;
- `10,000` validation H0 and `10,000` independent test H0 trials;
- 24 signal-validation and 60 signal-test material/contact contexts;
- all 136 ROI hypotheses per context and the preregistered 15-point amplitude
  grid;
- isotropic, aligned anisotropic and 45-degree anisotropic P/S tensors;
- `cP/cS`, attenuation, registration, port-amplitude and port-delay mismatch;
- one frozen covariance-regularized matched-filter/GLRT pipeline for all
  architectures.

The A-scan threshold does not generalize to the fresh test nuisance subspace:
its Pfa rises to `6.4%` (upper bound `6.81%`).  Consequently the reported
A-scan-to-FMC closure is explicitly a **frozen-threshold proxy**, not a claim
that every row is measured at identical realized test Pfa.  This mismatch is
conservative for the static designs—the A-scan receives too many false alarms
and therefore an artificially favorable `alpha90`—but it prevents treating
the numerical closure as industrial POD.

## Stop rules and verdict

The following rules are now supported by the map:

1. Reject any attractive `alpha90` point with `Pfa95 > 1.25%`; this removes
   `K4,Q4`, `K3,Q2`, `K4,Q1` and `K3,Q1` from advancement.
2. Pareto-prune a state set if another frozen point has no larger `K`, `Q`,
   Tx-event count, calibration count or `alpha90`, with at least one strict
   improvement.
3. Do not enlarge a fixed-`K` state budget after two consecutive marginal
   `alpha90` improvements below `5%`, or after two consecutive Pfa failures.
4. Never replace a `(K,Q)` winner using these final-test results.  Any `Q>4`
   study requires a new preregistered split and starts from the validation
   rule, not from the observed test ranking.

The state-diversity frontier has **not** saturated inside `Q<=4`: the valid
`K3,Q4` point improves on `K3,Q3` by `8.3%`.  Therefore there is no empirical
no-go result for a fifth static state, but neither is such a state authorized
by this campaign.  The strongest present reduced-model target is:

```text
3 electronic channels x 4 passive reflection-balanced states
```

The lower-overhead alternatives `K3,Q3` and `K4,Q2` remain co-leaders for a
future product-cost decision.  This branch does **not** authorize geometry:
the result is still a reciprocal elastic Born surrogate and must first be
compared with the parallel motion-allowed branch and confirmed in paired
transient full-wave physics.

## Reproduction and artifacts

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_operator_o3_static_frontier.jl
```

Primary artifacts:

- `campaign_contract.csv`, `seed_contract.csv` — frozen common protocol;
- `search_stage1.csv` — all 2,304 nominal candidates;
- `validation_proxy.csv` — 192 nuisance-aware finalists;
- `family_winners.csv`, `cell_winners.csv` — pre-test selection audit;
- `architecture_metadata.csv` — energy, channels, states and calibration cost;
- `thresholds.csv` — validation thresholds and held-out Pfa;
- `detection_curves.csv`, `test_context_*.csv` — final signal results;
- `frontier_summary.csv` — gates and Pareto labels;
- `static_frontier.jld2` — compact reproducibility bundle;
- `static_resource_frontier.png`, `symmetry_family_validation.png` — plots.

## Scope limitation

This is a geometry-free 2D reciprocal elastic Born campaign.  It omits finite
crack contact, nonlinear scattering, transducer electromechanics, switching
transients and manufactured calibration drift.  It maps the static operator
class; it is not a POD certification and not a geometry design.
