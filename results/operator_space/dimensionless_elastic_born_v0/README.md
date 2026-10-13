# Dimensionless elastic Born operator pilot v0

## Outcome

`O1` passes. The common virtual-aperture response tensor is reciprocal to
`9.54e-17` relative error and exactly equivariant under reflection on the
chosen grid. The resulting `O2` pilot separates known-zone detection from
unknown-position localization before any physical geometry is synthesized.

The strongest architectural result is negative but useful: every symmetric
single-port front-end has exact `y <-> -y` ambiguity. A fixed focus strongly
raises response at its design point, but does not become a localization
operator.

| architecture | channels | Tx events | complex features | center detectability | minimum detectability | max coherence at separation `>= 1 lambda_P` | stable rank |
|---|---:|---:|---:|---:|---:|---:|---:|
| uniform A-scan | 1 | 1 | 7 | 0.601 | 0.489 | 1.000 | 3.070 |
| fixed reciprocal focus | 1 | 1 | 7 | **5.278** | 0.069 | 1.000 | 2.713 |
| even dispersive code | 1 | 1 | 7 | 0.250 | 0.080 | 1.000 | 3.134 |
| asymmetric dispersive code | 1 | 1 | 7 | 0.666 | 0.261 | 0.975 | 3.423 |
| two symmetry ports | 2 | 2 | 28 | 0.426 | 0.397 | 0.842 | 5.462 |
| four cosine ports | 4 | 4 | 112 | 0.903 | **0.820** | **0.732** | **7.140** |
| 17-element FMC oracle | 17 | 17 | 2023 | 1.534 | 1.123 | 0.159 | 29.804 |

All values use the `N0` unit-white-feature noise convention. They compare the
operators internally; they are not dimensional SNR or a hardware prediction.

## Known-zone detection versus localization

At the focus, the fixed reciprocal code has `8.79x` the N0 detectability of
the uniform A-scan. Across the full ROI its minimum detectability is only
`0.141x` the A-scan value, and mirror coherence is exactly one. The current
TTD branch therefore remains a credible `known-zone detection` point but is
not evidence for imaging an unknown defect position.

Breaking reflection symmetry with one dispersive port removes exact mirror
degeneracy, but only weakly: its worst nonlocal coherence is `0.975`. Two
even/odd symmetry ports reduce this to `0.842`, and four cosine ports to
`0.732`. The latter still remains far from the FMC nonlocal value `0.159`.

## Ideal linear compression ceiling

The SVD code is optimal for retained dictionary Frobenius energy and aggregate
Gaussian information under orthonormal linear measurement projections. It is
not yet optimized for worst-case detectability or mutual coherence.

| ideal measurements `K` | retained dictionary energy | stable rank | max pair coherence |
|---:|---:|---:|---:|
| 4 | 0.129 | 3.831 | 0.9999 |
| 8 | 0.244 | 7.271 | 0.9933 |
| 16 | 0.445 | 13.276 | 0.9850 |
| 32 | 0.729 | 21.731 | 0.9743 |
| 64 | 0.942 | 28.062 | 0.8301 |
| 128 | 0.9997 | 29.795 | 0.7327 |

The full dictionary has stable rank `29.8`: the selected dense ROI contains
many more than four aggregate measurement degrees of freedom. This does not
mean that 30 electronic channels are required. It means that a four-code
claim must exploit task priors, accept reduced field of view/resolution, use
additional frequency/time/position diversity, or explicitly demonstrate a
better task-driven codebook.

## P/S mode ablation

The shear contribution is not only a loss channel in this surrogate.

| architecture | P-only stable rank | P+S stable rank | ratio | information gain, nats |
|---|---:|---:|---:|---:|
| uniform A-scan | 2.884 | 3.070 | 1.065 | 1.326 |
| asymmetric one-port code | 2.951 | 3.423 | 1.160 | 1.721 |
| two ports | 4.989 | 5.462 | 1.095 | 2.343 |
| four ports | 6.674 | 7.140 | 1.070 | 6.039 |
| FMC oracle | 25.130 | 29.804 | **1.186** | **21.734** |

The fixed-focus stable rank decreases slightly when S content is added, even
though its aggregate information increases. Focusing and mode diversity are
therefore distinct objectives.

## O2b reciprocal codebook search

The first cosine ports were a reproducible symmetry baseline, not an optimum.
A deterministic random search evaluated `300` codebooks in each of four
families for both two and four ports:

```text
real / complex
x
reflection-sector preserving / unconstrained.
```

Every candidate has orthonormal spatial ports and uses the reciprocal
bilinear response `W^T H W`. The balanced selection score closes the ambiguity
gap to FMC while capping the reward for minimum detectability at the uniform
A-scan value.

| ports | family | minimum detectability | nonlocal coherence | stable rank | information, nats |
|---:|---|---:|---:|---:|---:|
| 2 | cosine baseline | 0.397 | 0.842 | 5.462 | 20.430 |
| 2 | best real sector | 0.321 | 0.661 | 8.359 | 23.568 |
| 2 | best complex sector | 0.280 | **0.615** | **8.566** | 25.444 |
| 4 | cosine baseline | **0.820** | 0.732 | 7.140 | 45.594 |
| 4 | best real sector | 0.587 | 0.362 | **13.337** | 49.943 |
| 4 | best complex sector | 0.563 | **0.343** | 12.259 | **52.686** |

The four-port sector codebooks roughly halve nonlocal ambiguity relative to
the cosine baseline and almost double its stable rank. The real-sector target
is especially attractive for the next constraint layer: it avoids arbitrary
carrier phases while retaining most of the gain.

No sampled two-port codebook simultaneously reached the uniform A-scan's
minimum detectability `0.489` and the four-port ambiguity level. The largest
two-port minimum was `0.475`, paired with nonlocal coherence `0.847`. This is
not a mathematical no-go because the search is finite and stochastic; it is a
reason to stop treating two ports as the default architecture.

The best families preserve the even/odd reflection sectors rather than using
fully unconstrained ports. Group structure is therefore useful in this pilot,
not merely a computational convenience.

## O2c: causal bandwidth and robustness

The `O2b` carrier codebooks were then lifted to finite-band transfer models,
weighted by a five-cycle Hann burst, whitened with correlated frequency
covariance and evaluated on `120` common perturbation realizations. The
ensemble varies `cP/cS` by `+/-1%`, attenuation by `+/-20%`, port amplitude by
`3% RMS`, port delay by `0.02` carrier cycles RMS and lateral registration by
`+/-0.10 lambda_P0`.

No single-state four-port target passed the complete preregistered gate. A
second search optimized another `1000` sector codebooks directly on the causal
pulse operator; its best nominal gate margin was still only `0.952`. The
failure is therefore treated as evidence against the tested one-state class,
not as a request for a different lens geometry.

The minimum successful extension keeps four electronic channels but uses two
sequential static real code states. Each state is reciprocal and contractive,
both preserve the even/odd reflection sectors, and their dictionaries are
stacked with a `1/sqrt(2)` factor so total incident energy is unchanged.

| O2c candidate | electronic channels | passive states | Q10 min / A-scan | Q90 nonlocal coherence | Q10 stable rank | Q10 calibration fidelity | pass |
|---|---:|---:|---:|---:|---:|---:|---|
| original real one-state | 4 | 1 | 0.843 | 0.568 | 7.364 | 0.909 | no |
| original causal complex one-state | 4 | 1 | 0.915 | 0.603 | 7.475 | 0.884 | no |
| optimized real two-state A | 4 | 2 | **1.120** | 0.493 | **10.074** | 0.908 | **yes** |
| optimized real two-state B | 4 | 2 | 0.995 | **0.489** | 9.846 | **0.910** | **yes** |

The two-state result closes `O2c` only for the reduced elastic operator. It
does not yet show probability of detection in structural clutter and does not
authorize geometry synthesis. Its architectural meaning is precise:

```text
4 electronic channels x 2 sequential passive real symmetry states
    -> current O3 entrant.
```

The causal complex-delay branch remains a useful upper-bound/control, but its
best tested two-state candidate misses the calibration-fidelity gate
(`0.884 < 0.90`) and is not advanced.

## Verdict

- `O1 — information ceiling`: **pass** for the dimensionless weak-scatterer
  task.
- `O2 — compression frontier`: **pass**. One-port symmetric localization is a
  no-go; two ports remain a lower-cost conditional branch.
- `O2c — causal robustness`: **pass** for the four-channel, two-state static
  real reflection-sector class. The tested one-state classes fail.
- `O3 — architecture frontier`: **current**. Structural clutter, ringing,
  backwall and defect/no-defect decision statistics are still absent.
- geometry handoff: **blocked by design**, not by a failed calculation.

The next permitted calculation is the common paired pulse-echo `O3` benchmark
of the two-state entrant against A-scan, mechanical scan and FMC/PAUT oracle.
No new FEM geometry is justified yet.

## Model and normalization

- `lambda_P0 = c_P = f0 = 1`;
- `c_S/c_P = 0.5`;
- aperture `D = 8 lambda_P0`, 17 virtual elements;
- ROI `x/lambda_P0 = 2.5--6.0`, `y/lambda_P0 = -2.0--2.0`;
- seven frequency anchors `f/f0 = 0.8--1.2`;
- reciprocal far-field P/S Green projectors with dimensionless attenuation;
- isotropic point-scatterer Born tensor `Q=I`;
- total transmit energy fixed across sequential port excitations.

The O2c layer includes pulse covariance and a calibration perturbation
ensemble but still omits finite-defect resonance, near-field corrections,
nonlinear crack contact, explicit electronics/ringing, bond layers and
material clutter. Plane strain/extruded interpretation uses `translation_kz`;
the pilot is not a sphere or finite penny-shaped-crack calculation.

## Reproduction

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_operator_space_pilot.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/optimize_operator_codebooks.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/run_operator_causal_robustness.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/optimize_operator_causal_codebooks.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/optimize_operator_multistate_codebooks.jl
```

Primary artifacts:

- `architecture_metrics.csv` — common operator and resource metrics;
- `ideal_k_frontier.csv` — SVD aggregate-information compression ceiling;
- `ambiguity_pairs.csv` — worst adjacent and nonlocal hypothesis pairs;
- `mode_ablation.csv` — P-only versus P+S comparison;
- `codebook_search.csv` and `codebook_selected.csv` — `O2b` candidate cloud
  and selected reciprocal targets;
- `optimized_codebooks.jld2` — exact selected port matrices;
- `o2c_summary.csv` and `o2c_ensemble.csv` — preregistered one-state causal
  gate and its `120` perturbation realizations;
- `o2c_codebook_finalists.csv` — direct causal one-state re-optimization;
- `o2c_multistate_summary.csv` and `o2c_multistate_ensemble.csv` — final
  two-state gate and ensemble;
- `o2c_multistate_codebooks.jld2` — exact O3 entrant state matrices;
- `singular_values.csv` — complete FMC singular spectrum;
- `operator_checks.csv` — reciprocity and reflection gates;
- `operator_dictionary.jld2` — compact reproducible operator bundle;
- `architecture_frontier.png`, `ideal_k_compression.png`,
  `mode_ablation.png`, `codebook_search.png`, `o2c_nominal.png`,
  `o2c_robustness.png`, `o2c_multistate_search.png` — visual summaries.
