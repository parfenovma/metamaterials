# O6f paired physical-risk contract

> Frozen before either O6f branch is calculated. O6e showed that q8 static
> coding and P12/P17 motion are nondominated only because each still lacks a
> different physical-risk coordinate. O6f computes both without choosing one.

## O6f-A: concrete reciprocal q8 lattice

The nominal `17 x 3` q8 contraction is factorized as

```text
B = U * Sigma * V'

3 public ports
  -> exact 3 x 3 reciprocal input rotation V'
  -> three singular-value attenuators with explicit dump ports
  -> exact 3 x 3 reciprocal rotation R
  -> 17-output adjacent-Givens lattice Q' whose first three inputs map to U
  -> 17 isolated piezo columns.
```

The rectangular QR elimination uses `16+15+14=45` adjacent two-mode lattice
couplers. The two small `3 x 3` rotations are counted as physical calibrated
hybrid blocks but treated as ideal in the first parasitic bound; all loss and
dispersion assigned to the 45-output lattice is therefore optimistic.

The zero-parasitic network must reproduce q8 to numerical precision. Then:

1. find the maximum identical insertion loss per Givens coupler with all
   other errors zero;
2. find the maximum identical group delay per coupler with all other errors
   zero;
3. use `20,000` deterministic draws of independent angle and group-delay
   errors on all 45 couplers. Reference component errors are `0.1 degree RMS`
   and `2 ns RMS`; a common scale is swept and bisected. The pass criterion is
   a one-sided Wilson lower-95 bound of at least `99%` joint O4 pass.

No post-error contraction or port rebalance is allowed: the synthesized
network itself must remain passive and meet the frozen transfer. Common input
hybrid loss, winding resistance, piezo capacitance and electrical noise remain
outside this optimistic first topology bound and must fit inside any residual
loss headroom.

## O6f-B: frozen-layout motion robustness

The exact compact P12 and uniform P17 index sets, nominal GLRT templates,
training nuisance subspace and validation thresholds from
`motion_frontier_v0` are frozen. Layouts and thresholds are not reoptimized.

For every physical realization, each stop receives fixed independent errors:

```text
encoder position RMS                 0.05 lambda_P0 * scale
two-way contact log-amplitude RMS    5% * scale
two-way contact delay RMS            0.02 carrier cycles * scale
```

The common scale is `0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0`. Errors are fixed
across frequency and all defect hypotheses within a contact realization.
They modify both signal and the material/backwall/ringing nuisance observables.

The campaign reconstructs the frozen nominal detector, verifies its baseline
against the saved P12/P17 result, then uses the same 60 frozen signal contexts
with fresh O6f physical-error and electronic-noise seeds. Each scale receives
at least `10,000` aggregate H0 trials per layout. Gates are:

```text
Pfa Wilson upper 95%                  <= 1.25%
alpha90 / new nominal replay          <= 1.10
Pdl at the frozen reference amplitude >= 0.80
```

This is still a reduced Born physical-risk replay. It does not claim probe
wear, lift-off dynamics or full-wave surface-wave coupling.

## Joint decision

- A network branch passes O6f-A only as a topology feasibility bound; actual
  RLC/piezo impedance remains required before fabrication.
- A motion branch passes O6f-B only up to the largest tested physical scale
  satisfying all three gates.
- No scalar winner is formed across evidence levels. O6f identifies which
  physical tolerances are credible enough to justify the next, more expensive
  full-wave or benchtop calculation.
