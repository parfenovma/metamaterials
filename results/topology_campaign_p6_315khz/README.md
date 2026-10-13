# Multi-topology J1/J2 campaign at 300--322 kHz

## Outcome

The campaign found two effective ways to decouple the nominal-focus amplitude
`J1` from the reduced defect-CNR proxy `J2`:

1. align the passive lens working distance with its actual focal plane;
2. use a segmented source and explicitly constrain the target to remain the
   strongest local response.

The original fixed target at `x=115 mm` made several passive designs look like
near-misses.  Reusing their converged fine fields at the correct assembly
standoff produces strict-feasible passive candidates without changing the lens.

## Recommended prototypes

| priority | candidate | frequency/control | J1 | J2 | exit energy | focus error | gate |
|---:|---|---|---:|---:|---:|---:|---|
| 1 | ellipse `c005` + standoff | 310 kHz; assembly 9 mm toward known defect | 1.3533 | 1.7278 | 0.7806 | 1 mm | strict |
| 2 | segmented straight | 315 kHz; 7 segments, 4 mirrored channels | 1.1685 | 1.3587 | 1.0555 | 3 mm | strict |
| 3 | matched-GRIN `g001` + standoff | 316 kHz; assembly 5 mm toward defect | 1.1508 | 1.5759 | 0.9714 | 1 mm | strict |
| 4 | segmented `c005` target-lock | 315 kHz; 4 mirrored channels | 1.2531 | 1.6581 | 0.7004 | 2 mm | prototype |
| 5 | drilled Fresnel `f006` + standoff | 321.7 kHz; assembly shift 6 mm | 1.0793 | 1.2715 | 0.8540 | 3 mm | strict |

The best discovered strict point for a known target depth is passive `c005` at
310 kHz.  The simplest mechanically is the segmented straight aperture, because
it needs no passive lens.  The matched-GRIN is the conservative high-transmission
passive fallback.

## What was screened

- Ellipse baseline and local frequency sweep.
- Matched-GRIN: 5/7 rows with circular matching rows and elliptical core;
  23 coarse plus two fine solutions.
- Quantized Fresnel: 3--5 fill levels using round through-holes; 24 coarse plus
  one fine solution.
- Segmented source: seven source segments, two footprints, three coarse
  frequencies and fine 315 kHz response matrices; field superposition was
  verified to relative L2 error below `4e-13`.
- Freeform oracle: continuous and segmented curved slots; 24 coarse plus two
  fine solutions.
- Passive target/standoff sweep using existing fine complex fields; no new FEM.
- Active-straight broadband screen from 295 to 335 kHz and a sparse pulse proxy.

All common rows use a same-frequency straight reference and equal incident
energy.  The shared definitions and gates are in `benchmark_contract.md`.

## Important branch conclusions

### Matched-GRIN

At the fixed 115 mm target, fine `g001@315/316` is a near-miss: high transmission
`0.96--0.97`, useful `J2=1.31--1.35`, but 6 mm focal error.  A 5 mm assembly
shift makes `g001@316` strict-feasible and improves its re-centered metrics to
`J1=1.1508`, `J2=1.5759`.

### Fresnel round holes

The fine result is weaker than its coarse resonance suggested.  Re-centering
still yields a simple strict backup with 36 round holes, five diameters and a
2.61 mm minimum feature.  It is attractive for manufacturability, not maximum
gain.

### Segmented source

This branch gives genuine independent control of J1/J2.  Symmetric target-lock
reduces seven segments to four independent channels.  With `c005` it reaches
`J1=1.2531`, `J2=1.6581`, but is narrowband: the fixed-weight coarse minimum
joint gain over 305/315/321.7 kHz is `0.882`.

The lens-free segmented straight passes the strict fine gate at 315 kHz.  Its
295--335 kHz coarse screen is qualified rather than fully positive: fixed phase
and true-time-delay controls have robust minimum joint gains `0.969` and
`0.960`, and the focal position is not coarse/fine robust.  The fine anchor is
credible, but hardware must be frequency-calibrated.

### Freeform slots

Categorical no-go.  Aperture-spanning and serial traction-free slots behave as
barriers to the longitudinal wave.  The best fine result has only
`J1=0.1815`, `J2=0.0465`, transmission `0.2315`.  Future freeform searches must
enforce strong longitudinal aluminium corridors and use weak isolated
perturbations, not transverse slot screens.

## Limits of the result

- `J2` is a reciprocal Born/reduced-CNR spatial proxy, not a simulated echo from
  a contact defect.
- Standoff re-ranking is valid when defect depth is known and moving the whole
  assembly does not change the coupling stack.  It is not a universal solution
  for unknown-depth inspection.
- The current model does not include coupled PZT impedance, bond layers, `S11`,
  electrical matching, amplifier noise or nonlinear crack contact.
- Sparse pulse reconstruction covers only 25/41/61/81% of the relevant spectrum
  for 3/5/8/12-cycle Hann bursts.  Only the 12-cycle result is semi-informative;
  it is not a replacement for a broadband transient solve.
- The second-harmonic paths tested earlier were both unfavorable.  The useful
  300--322 kHz response is a linear carrier window, not harmonic gain.

## Experimental sequence

1. Machine `c005` and `g001`, retaining the solid straight specimen.
2. Sweep 295--325 kHz and mechanically scan standoff around the predicted
   shifts before introducing a defect.
3. Build the four-channel mirrored segmented source controller; begin with the
   straight specimen and calibrate phases/amplitudes from measured transfer
   functions.
4. Record reference/defect pairs and compute the actual matched-filter echo,
   ringing-window energy and CNR.
5. Only after these measurements, fit coupling/bond/PZT parameters and run the
   full transient defect/no-defect optimization.

## Campaign artifacts

- `data/combined_topology_metrics.csv`: every compatible campaign row.
- `data/fine_prototype_shortlist.csv`: fine prototype-feasible points.
- `data/fine_strict_shortlist.csv`: fine strict-feasible points.
- `data/prototype_recommendations.csv`: concise fabrication ranking.
- `figures/fine_prototype_shortlist.png`: common feasible J1/J2 comparison.
- `figures/common_J1_J2_pareto.png`: raw campaign Pareto, including infeasible
  diagnostic points.
- `figures/common_constraints.png`: transmission/focal-error gates.
- `topology_campaign.jl`: reproducible aggregator.

Branch reports:

- `../matched_grin_multidof_p6_315khz/README.md`
- `../fresnel_zone_multidof_p6_315khz/README.md`
- `../segmented_piezo_response_p6_315khz/README.md`
- `../freeform_oracle_p6_315khz/README.md`
- `../passive_standoff_rerank_p6_315khz/README.md`
- `../active_straight_broadband_p6_315khz/README.md`

