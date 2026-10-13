# O8: book-guided time-reversal and modal audit

> Frozen before reading the audit outputs on 4 October 2026.

## Question

Does the already frozen, physically passing O8 PZT/RL head gain useful
independent spatial modes and a sharper virtual time-reversal focus when all
four shunt states are used, or does Q4 mainly recode the same q00 subspace?

This is a post-hoc diagnostic, not a new optimization stage and not a new POD
claim. No geometry, load, transfer normalization, detector threshold, or
previous O8 result may be changed from this audit.

## Fixed input and energy rule

The only physical input is the saved complex `17 x 3 x 4 x 7` PZT transfer in
`results/operator_space/o8_plane_strain_v0/pzt_physical_bundle.jld2`. The
inherited common O8 dimensional scale is retained; state-wise and
frequency-wise transfer normalization remain forbidden.

The seven frequency anchors use the normalized amplitude spectrum of the
project-standard five-cycle Hann burst. For a set of `Q` physical states,

```math
A_Q=\frac{1}{\sqrt Q}[w_f H_q(f)]_{q,f},
\qquad
K_Q=A_QA_Q^\dagger
   =\frac{1}{Q}\sum_{q,f}w_f^2H_q(f)H_q^\dagger(f).
```

The `1/sqrt(Q)` field factor fixes total incident energy rather than granting
Q4 four times the q00 energy. By reciprocity, each column of `K_Q` is the
virtual time-reversal field produced by recording at one output bin and
playing the three port signals back with conjugated phase over the band.

## Frozen comparisons and diagnostics

Only two architectures are compared:

- `q00`: one open-shunt state, three ports, seven frequencies;
- `q4`: all four states, three ports per state, the same seven frequencies.

Reported diagnostics are fixed in advance:

1. singular-value spectrum of `A_Q`;
2. stable rank and entropy rank;
3. observable spatial-mode count at 20, 30, 40, and 60 dB amplitude dynamic
   range, where mode `n` is observable when
   `20 log10(sigma_n/sigma_1) >= -DR`;
4. energy in the three added states outside the observable q00 subspace at
   20, 30, and 40 dB;
5. for all 17 focal bins, target-peak validity, contiguous half-maximum width,
   mean off-target focus contrast, and peak intensity outside that
   half-maximum interval;
6. normalized spatial-coherence matrices and the centre-bin focus profile.

There is no pass/fail gate. The audit explains the mechanism behind the
already completed detector comparison; it cannot retroactively select O8.

## Interpretation boundary

The stored transfer already includes the O8 material, radiation, and PZT
losses. Varying an arbitrary scalar attenuation would not test modal loss—it
would merely rescale every singular value. Seven frequency anchors are also
too sparse to estimate linewidth, modal Q, decay time, or ringing. Therefore
the dynamic-range sweep is explicitly an observability sensitivity analysis,
not a loss sweep. A genuine loss/ringing study requires a dense-band solve or
a measured impulse response and must be preregistered separately.

The stochastic noise and structured clutter already used by the O8 detector
replay are not resampled here. Relating the modal result to Pd is an
interpretive comparison with that frozen replay, not another detector test.
