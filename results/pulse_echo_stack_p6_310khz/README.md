# Mechanically coupled PZT--bond--Al pulse-echo pilot

## Outcome

The first paired transient calculation reverses the ranking from the harmonic
focus/CNR-proxy campaign.  Under one common 310 kHz, five-cycle, equal-input-
energy pulse-echo benchmark, neither passive plate beats the solid straight
reference.

| topology | matched-filter gain vs straight | scattered peak gain | echo-energy gain | echo / defect-free ringing |
|---|---:|---:|---:|---:|
| solid straight | 1.000 | 1.000 | 1.000 | +47.48 dB |
| ellipse `c005` | 0.631 | 0.602 | 0.416 | +3.15 dB |
| matched-GRIN `g001` | 0.619 | 0.547 | 0.427 | -1.38 dB |

The paired defect echo arrives in the expected `~40--58 us` window.  For the
straight reference its peak mechanical receive proxy is `2.704e4 Pa`.  The
passive structures retain delayed internal motion after the transmit burst;
this reduces the defect echo and raises the baseline in the same receive
window.  In `g001`, the integrated defect-free ringing is larger than the
paired echo energy.

This is the separation the earlier `J1/J2` hierarchy was designed to reveal:
good one-way harmonic focusing does not guarantee good double-pass pulse-echo
detectability.

## Seven-segment active design-point check

The symmetric seven-channel straight aperture was then evaluated with the
fine harmonic target-lock weights at their actual design point: `315 kHz`,
target/defect at `x=115 mm`.  The prescribed channel weights have
`sum(abs2(w))=7` and amplitude dynamic range `3.987`; reported pulse-echo
metrics are additionally normalized by the mechanical work delivered to the
stack.

| transmit / receive | matched-filter gain | scattered-peak gain | echo-energy gain | echo/ringing | gain vs straight |
|---|---:|---:|---:|---:|---:|
| uniform / uniform | 1.0000 | 1.0000 | 1.0000 | 48.69 dB | 0 dB |
| active / uniform | 0.9974 | 1.0010 | 0.9855 | 51.73 dB | **+3.04 dB** |
| active / reciprocal matched | 0.9888 | 0.9944 | 0.9713 | 52.31 dB | **+3.62 dB** |

This is a genuine, if modest, decoupling of the two objectives.  The active
controller does **not** amplify the paired defect echo: matched-filter
amplitude remains within `1.1%` of straight and echo energy falls by `2.9%`.
Instead it suppresses deterministic defect-free ringing in the echo window.
The resulting `+3.62 dB` echo/ringing improvement corresponds to about `1.52x`
in an amplitude-style CNR proxy.  Roughly `+3.04 dB` comes from transmit
shaping; reciprocal receive combining adds another `+0.58 dB`.

An off-design control at `310 kHz / x=106 mm` stays essentially at parity
(`matched-filter gain 1.005`, echo-energy gain 1.032, echo/ringing gain
`-1.11 dB`).  This confirms that the CNR improvement is target/frequency
specific and should not yet be called broadband.

## Reflective-backwall flaw-detector stand

The design-point stack was extended to `90 us`, and the absorbing boundary at
the end of the `180 mm` aluminium coupon was replaced by a traction-free
reflective backwall.  The raw monostatic trace now contains the three landmarks
expected on an A-scan: the initial pulse, the `115 mm` defect echo and the
backwall echo at approximately `177 mm` after zero calibration.

| transmit / receive | echo-energy gain | echo/ringing | gain vs straight |
|---|---:|---:|---:|
| uniform / uniform | 1.0000 | 35.66 dB | 0 dB |
| active / uniform | 0.9856 | 39.37 dB | +3.71 dB |
| active / reciprocal matched | 0.9714 | 42.99 dB | **+7.33 dB** |

The backwall raises deterministic background in the defect gate and makes the
benchmark materially harder.  The active controller still keeps defect-echo
amplitude within `1.2%` of straight while suppressing that background.  This is
the first result in the branch that is directly visible on an instrument-like
trace without defect/reference subtraction.  It is not yet an experimental
CNR because electrical loading and material clutter remain absent.

## Common physical setup

- 2D plane strain.
- 310 kHz, five-cycle Hann burst, 70 us record.
- `6.0 mm` mechanically represented PZT-5A layer.
- `0.15 mm` epoxy bond.
- Al 6061 surrogate.
- circular traction-free void, radius `2 mm`, at `x=106 mm`, `y=0`.
- monostatic receive observable: mean PZT `sigma_xx`.
- defect signal: `s_void(t) - s_free(t)`.
- all reported gains rescaled to the straight specimen's delivered mechanical
  input energy.

The meshes are conformal across PZT/bond/Al interfaces and around every lens
void and defect.  The current pilot meshes contain approximately `6.9k--9.7k`
nodes and use 347 average-acceleration Newmark steps.  Constant `M`, `C` and
`K` matrices are assembled once; one symmetric factorization is reused for the
whole record.

## What this result does and does not prove

This is a useful transient gate, not yet a fabrication-grade coupled-PZT
prediction.

- PZT is presently an isotropic mechanical surrogate driven by prescribed
  back-face traction.  Piezoelectric `e`, dielectric `epsilon`, electrodes,
  cable and electrical load are not included yet.
- Mean PZT stress is a mechanical voltage proxy, not a computed terminal
  voltage.
- The coarse pilot grid resolves timing and the large ranking difference, but
  the winning result still needs a medium/fine convergence pair.
- Absorbing lateral and far boundaries suppress specimen clutter.  Therefore
  the straight specimen's `+47.5 dB` echo/ringing value is an optimistic upper
  bound, not an experimental CNR.
- There is no stochastic noise or material-clutter ensemble; the reported
  value is echo relative to deterministic defect-free ringing.
- The comparison deliberately uses one common frequency and defect depth.
  `g001` still needs a qualified local check at its prior `316 kHz / 110 mm`
  harmonic working point before final rejection.

## Decision

The solid/segmented-active route is now the leading physical branch.  Its
first useful optimization direction is ringing cancellation rather than
maximum echo amplitude.  Passive
`c005` and `g001` must not be fabricated solely on the earlier harmonic scores.
They remain conditional until:

1. a coarse local carrier/depth transient check is completed;
2. the result is repeated on a medium/fine mesh;
3. the PZT constitutive and electrical receive model is enabled.

## Reproduction

```bash
julia --startup-file=no results/pulse_echo_stack_p6_310khz/pulse_echo_stack.jl --pilot-all
julia --startup-file=no results/pulse_echo_stack_p6_310khz/pulse_echo_stack.jl --pilot-active-design
julia --startup-file=no results/pulse_echo_stack_p6_310khz/pulse_echo_stack.jl --pilot-flaw-detector
```

Individual meshes and solves can be regenerated with:

```bash
julia --startup-file=no results/pulse_echo_stack_p6_310khz/pulse_echo_stack.jl --mesh c005 defect
julia --startup-file=no results/pulse_echo_stack_p6_310khz/pulse_echo_stack.jl --solve c005 defect
```

## Artifacts

- `data/pulse_echo_metrics.csv`: common paired metrics.
- `data/pulse_echo_*_none_310khz.jld2`: defect-free transients.
- `data/pulse_echo_*_void_r2_310khz.jld2`: defect transients.
- `figures/pulse_echo_comparison.png`: common equal-energy comparison.
- `figures/pulse_echo_straight.png`, `pulse_echo_c005.png`,
  `pulse_echo_g001.png`: paired waveforms and defect subtraction.
- `data/segmented_pulse_echo_metrics_315khz_x115.csv`: design-point active
  transmit/receive comparison.
- `data/active_straight_weights.csv`: exact seven complex drive weights and
  segment bounds.
- `figures/segmented_pulse_echo_comparison_315khz_x115.png`: magnified
  amplitude trade-off and deterministic CNR-proxy gain.
- `data/segmented_pulse_echo_metrics_315khz_x115_backwall.csv`: paired metrics
  with a reflective `180 mm` backwall.
- `figures/*_315khz_x115_backwall.png`: raw and paired reflective-backwall
  traces and their comparison.
- `pulse_echo_stack.jl`: mesher, constant-matrix Newmark solver and analysis.
