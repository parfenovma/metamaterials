# Preregistered coupon acceptance: BUILD or KILL

## Question

Does the manufactured stationary coupon create four repeatable measurement
operators with the state separation predicted by the CTS-3195HD model?

This experiment validates the four-state boundary mechanism only.  It does
not estimate flaw POD and does not validate the final three transducer ports.

## Bench setup

- rigid kinematic support touching only the long free edges;
- three identical broadband source transducers centred at `y=-50, 0, +50 mm`
  on the source face, sequentially pulsed rather than phase-driven;
- scanning laser Doppler vibrometer on the opposite face at the 17 coordinates
  in `receiver_scan_points.csv`;
- AWG/pulser and amplifier operated in a verified linear regime;
- digitizer with a shared trigger and complex response extraction;
- LCR/impedance analyser capable of measurements near `250 kHz`.

If no LDV is available, a miniature contact receiver may be repositioned over
the same 17-point template for this coupon test.  That substitution is a bench
measurement method, not permission to scan the eventual inspection head.

## Frozen acquisition

1. Before PZT bonding, record the bare-block reference with the final source
   fixture at seven frequencies `193.6, 209.733, 225.867, 242.0, 258.133,
   274.267, 290.4 kHz`.
2. After complete assembly, acquire states `q00`, `q10`, `q01`, `q11` in a
   newly randomized order for each repeat.
3. Use ten independent repeats.  Re-seat nothing between states or repeats.
4. At every state/frequency/source, record complex displacement at all 17
   receiver points.  Retain raw waveforms and calibration metadata.
5. Form the complex `17 x 3` transfer matrix `H_q(f)` for every state.
6. Normalize all matrices by one common scalar from the bare-block reference;
   never normalize states independently.

The subspace distance is

```text
d(H1,H2) = ||Q1*Q1' - Q2*Q2'||_F / sqrt(2K),  K=3
```

where `Q` contains the three left singular vectors of `H`.

## Gates

All gates must pass.

### Electrical

| Gate | Requirement |
|---|---:|
| capacitance match inside each mirrored pair | <=5% |
| tuned L_total error at target frequency | <=2% |
| tuned R_total error including parasitics | <=5% |
| relay-open added capacitance | <20 pF |
| state truth table and B-pair reversal | no discrepancy |

### Acoustic

| Gate | Requirement |
|---|---:|
| minimum amplitude / bare-block reference over all states and frequencies | >=0.70 |
| maximum amplitude / bare-block reference over all states and frequencies | <=1.50 |
| minimum pairwise state distance at 242 kHz using repeat means | >=0.05 |
| 95th percentile within-state repeat distance at 242 kHz | <=0.015 |
| between-state minimum / within-state 95th percentile | >=3.0 |
| minimum state-distance result reproduced in both halves of the 10 repeats | >=0.05 in each half |

Mirror error is reported as a diagnostic because fixture and LDV alignment
break ideal symmetry; it is not allowed to replace the state-distance gate.

## One permitted repair

If an electrical gate fails, or if the measured capacitance moved the shunt
resonance, recalculate L/R from the measured terminal data and perform one
complete retune.  Then discard the first acoustic dataset and rerun the entire
randomized acquisition with fresh state order.

No patch movement, parity change, selective repeat deletion, frequency-grid
change, or threshold change is permitted.

## Terminal verdict

- **BUILD:** all gates pass on the first dataset or after the single permitted
  electrical retune.  Preserve this geometry/parity for the stationary head.
- **KILL:** any gate remains failed.  Stop this four-state PZT topology; do not
  advance it to a detector prototype.

