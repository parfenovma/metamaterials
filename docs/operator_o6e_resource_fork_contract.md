# O6e passive-beamformer versus mechanical-motion resource fork

> Frozen before the O6e calculations. O6d leaves two materially different
> realizations alive: accept a full passive `3-to-17` coding network, or return
> to the already tested one-channel mechanical-motion frontier. O6e maps both;
> it does not collapse them into a guessed dollar price.

## 1. Common invariant

Both branches retain the same aperture, total incident-energy normalization,
qualified defect task and PAUT/FMC reference logic. Resources remain a vector:

```text
electronic channels
transmit events
mechanical positions and moves
normalized time under four existing motion models
active transducer columns
passive bus-to-column taps and weight levels
calibration sites / coefficients
unmodelled contact and positioning burden
detection evidence level
```

Metrics from the expanded full-wave F1d campaign and the fresh reduced O3
motion campaign are not numerically pooled. Their normalized gap closures are
shown side by side with an explicit evidence label, not treated as one common
ROC experiment.

## 2. Passive-beamformer branch

The frozen nominal map is the first O6d-passing quantization:

```text
3 public reciprocal terminals
17 mechanically isolated output columns
50 active signed bus-to-column taps
8 nonzero magnitude levels per DCT port
17 distinct normalized ratios across the complete map
3 transmit events, no moving position
```

The branch receives three optimistic checks before geometry:

1. **Common loss ceiling.** Compute the maximum one-way insertion loss allowed
   by the frozen O4 useful-amplitude gate.
2. **Tap tolerance ensemble.** Apply fixed log-amplitude and delay errors to
   all 50 taps, shared across the 13 frequencies. A single worst-band
   contractive scale enforces passivity. `50,000` deterministic draws are used
   for each tolerance point with a fresh O6e seed. The frontier uses scale
   factors `0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5` times the O4 reference
   (`3%` log-amplitude RMS and `82.6 ns` delay RMS). A bisection reports the
   largest common scale whose one-sided Wilson `95%` lower bound on joint O4
   pass probability is at least `99%`.
3. **Single open-tap audit.** Every nonzero nominal tap is removed once. This
   is a fault diagnostic, not part of the manufacturing pass.

Network parasitic loss, bus capacitance, transformer leakage, dielectric noise
and column mechanics are otherwise absent. Thus a tolerance pass is necessary,
not a device pass.

The nominal eight-level state is replayed once through the already frozen F1d
contexts with the same random draws as ideal DCT3. Reusing this test is allowed
because O4 explicitly requires as-built transfer replay; it is not a new
architecture selection campaign.

## 3. Motion branch

No motion detection data are regenerated. O6e imports only the held-out,
Pfa-controlled entries already frozen in `motion_frontier_v0`:

- compact `P12`, the first aggregate crossover over the frozen static
  comparator;
- uniform full `P17`, the strongest tested monostatic scan.

Their four existing time models are retained verbatim. Encoder error, couplant
reset, probe footprint and continuous-motion effects stay explicit unknowns;
O6e does not invent favourable values for them.

## 4. Decision rule

- Beamformer geometry remains blocked unless nominal q8 passes frozen F1d,
  common loss headroom is positive, and the `99%` tolerance scale is reported.
- Motion remains an architecture candidate, not a physical prototype pass,
  until positioning/contact robustness exists.
- Pareto dominance is asserted only on commensurate resource coordinates.
  Detection results from different evidence levels cannot by themselves break
  a fabrication-versus-motion tie.
- The output must state the exact exchange: extra channels/taps/weight levels
  versus saved events/moves/time. No synthetic BOM or weighted scalar score is
  permitted at O6e.
