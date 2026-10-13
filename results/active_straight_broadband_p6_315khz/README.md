# Broadband validation of the symmetric active-straight target lock

## Bottom line

The seven-segment symmetric straight aperture remains a **strict fine-mesh
prototype at its 315 kHz design point**, but its target lock is not broadband or
mesh-converged.  Converting the optimized carrier phases into physical
true-time delays (TTD) does not materially improve the 295--335 kHz response
over simply holding the 315 kHz phases fixed.

Fine 315 kHz anchor, normalized to the equal-energy uniform straight source:

| controller | J1 | J2 | J2 clutter-dom. | exit | FWHM | focus error | strict |
|---|---:|---:|---:|---:|---:|---:|---|
| fixed 315 kHz phase | 1.1685 | 1.3587 | 1.3548 | 1.0555 | 45.60 mm | 3 mm | **pass** |
| true-time delay | 1.1685 | 1.3587 | 1.3548 | 1.0555 | 45.60 mm | 3 mm | **pass** |

The values are identical at 315 kHz by construction: the TTD weights reproduce
the optimized complex weights up to one irrelevant common phase with relative
error `9.1e-17`.

On the ten-point coarse basis, however, the same hardware gives:

| controller | min J1 | min J2 | min joint | min exit | max focus error | phase-fit RMS | group-delay p-p |
|---|---:|---:|---:|---:|---:|---:|---:|
| fixed phase | **0.9847** | **0.9691** | **0.9691** | 0.8632 | 22.80 mm | **1.170°** | **2.164 μs** |
| TTD | 0.9781 | 0.9601 | 0.9601 | **0.8686** | 22.80 mm | 1.195° | 2.244 μs |

Both controllers retain simultaneous `J1>1` and `J2>1` through 325 kHz.  At
330 kHz fixed phase still gives `J1=1.0038`, `J2=1.0157`, whereas TTD gives
`J1=0.9996`, `J2=1.0108`.  At 335 kHz both fall below the straight reference.
The energy proxy remains comfortably above 0.86 across the full interval.

The blocking result is focal position.  Only the 295 kHz coarse point reaches
the prototype boundary (`focus_error=5 mm`); no coarse active point passes the
strict gate.  At 315 kHz the fine solution has 3 mm error, but the coarse
solution has 10 mm error together with lower `J1=1.1124`, `J2=1.2068` and
`exit=0.9518`.  The strict target lock must therefore not be treated as a
converged broadband fabrication result.

## Fixed hardware

The fine symmetric target-lock weights were converted using the harmonic
convention

`w_i(f) = a_i exp(-i 2π f τ_i)`.

A common delay/phase does not affect any field metric, so all delays were
shifted to be nonnegative:

| segment pair | amplitude | phase at 315 kHz | delay |
|---|---:|---:|---:|
| 1 / 7 | 0.4505 | +70.303° | 0 μs |
| 2 / 6 | 1.1677 | +39.215° | 0.2741 μs |
| 3 / 5 | 1.3691 | -4.303° | 0.6579 μs |
| 4 | 0.3434 | 0° | 0.6200 μs |

The amplitude dynamic range is 3.987 and `sum(abs2(w))=7` at every frequency.
Only four settings are independent because mirror symmetry is retained.

The maximum differential delay is only 0.658 μs.  Across ±20 kHz it changes a
channel phase by at most about 4.7°, explaining why TTD and frozen carrier
phases are almost indistinguishable over this modest fractional bandwidth.
Fixed phase was included as the intended worst-case benchmark; the calculation
rejects that assumption and finds it marginally better in worst-case J1/J2 and
phase smoothness.

## Frequency basis and phase smoothness

The coarse response basis contains 295, 300, 305, 310, 315, 318, 321.7, 325,
330 and 335 kHz.  Existing 305, 315 and 321.7 kHz matrices are read from the
completed segmented-aperture branch; seven new matrices are stored only here.
Each matrix uses the same straight mesh and seven equal source segments.

Unwrapped target phase is close to linear for both active controllers.  The
reported group delay is `-d phase/d omega`, evaluated by finite differences on
the irregular ten-point grid.  It is a useful smoothness screen, not a
high-resolution dispersion measurement.  The active phase residuals are about
half those of uniform drive (`2.140°` RMS and `4.129 μs` group-delay ripple),
but TTD does not improve them relative to fixed phase.

## Sparse pulse reconstruction

Focus envelopes were synthesized for 3-, 5-, 8- and 12-cycle analytic Hann
bursts.  Complex target transfer was interpolated linearly in log-amplitude and
unwrapped phase between the ten harmonic samples, then integrated only over
295--335 kHz.

| cycles | spectral energy covered | fixed-phase peak gain | TTD peak gain | fixed late energy | TTD late energy |
|---:|---:|---:|---:|---:|---:|
| 3 | 25.0% | 1.0614 | 1.0622 | 0.1611 | 0.1619 |
| 5 | 40.6% | 1.0618 | 1.0626 | 0.0526 | 0.0530 |
| 8 | 60.9% | 1.0626 | 1.0634 | 0.0238 | 0.0239 |
| 12 | 80.8% | 1.0644 | 1.0651 | 0.00505 | 0.00481 |

`late energy` is the fraction after one burst duration beyond the reconstructed
envelope peak.  At 12 cycles the active controllers reduce this band-limited
proxy relative to uniform (`0.00745`) while increasing temporal FWHM slightly
from 31.01 to about 31.45 μs.

These traces are deliberately labelled **sparse, band-truncated spectral
reconstructions**.  The 3- and 5-cycle results omit 75% and 59% of input
spectral energy and are not valid transient ringing predictions.  Even the
12-cycle reconstruction omits 19%.  A defensible short-pulse result requires a
much wider response basis plus the real transducer, bond layer and receiver
transfer functions.

## Interpretation

- The amplitude apodization, not TTD, supplies most of the observed broadband
  benefit over a uniform aperture.
- Frozen phase is sufficient for this 40 kHz interval and is simpler hardware.
- Transmission is robust; target location is not.
- Before prototype fabrication, the 315 kHz target-lock solution needs a fine
  frequency bracket or a convergence study at several nearby frequencies.
- A later controller should optimize target position and worst-frequency J1/J2
  simultaneously rather than converting a single-frequency optimum to delays
  after the fact.

## Artifacts

- `data/campaign_metrics.csv`: campaign-compatible harmonic rows, including
  all 30 coarse controller/frequency combinations and three fine anchors.
- `data/broadband_metrics.csv`: coarse harmonic metrics plus unwrapped phase,
  finite-difference group delay and phase-fit residual.
- `data/broadband_summary.csv`: robust minima and phase-smoothness summary.
- `data/hardware_delays.csv`: amplitudes, carrier phases and physical delays.
- `data/fine_anchor_metrics.csv`: exact 315 kHz fine strict-gate check.
- `data/pulse_metrics.csv` and `data/pulse_traces.csv`: sparse reconstruction
  metrics and normalized focus envelopes.
- `figures/broadband_harmonic_metrics.png`: J1/J2/energy/focus/FWHM/group delay.
- `figures/target_phase_group_delay.png`: phase linearity diagnostics.
- `figures/hardware_drive_phase.png`: fixed-phase and TTD channel laws.
- `figures/sparse_pulse_focus_envelopes.png`: 3/5/8/12-cycle traces.
- `figures/pulse_width_ringing_proxy.png`: time-width and late-energy proxies.

## Reproduction

```sh
julia --startup-file=no active_straight_broadband.jl --delays
julia --startup-file=no active_straight_broadband.jl --basis 295
julia --startup-file=no active_straight_broadband.jl --run-missing
julia --startup-file=no active_straight_broadband.jl --report
```

This is a linear prescribed-traction controller/aperture upper bound.  It is
not a coupled-PZT, electrical-impedance or nonlinear defect model.
