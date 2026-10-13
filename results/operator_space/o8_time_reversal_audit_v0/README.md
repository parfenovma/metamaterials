# O8 book-guided time-reversal and modal audit

This is a post-hoc diagnostic of the frozen, physically passing O8 PZT/RL transfer. It does not tune the head and does not open a new detector split.

The pulse-weighted operator is

```math
A_Q=[w_f H_q(f)]_{q,f}/\sqrt{Q}, \qquad K_Q=A_QA_Q^\dagger.
```

The `1/sqrt(Q)` factor keeps total incident energy fixed. The seven weights are the normalized amplitude spectrum of the project-standard five-cycle Hann burst. No state-wise or frequency-wise transfer normalization is used; the inherited common O8 scale is `5.622408245140293e12`.

| architecture | stable rank | entropy rank | modes @20 dB | @30 dB | @40 dB | median FWHM / lambdaP | median focus contrast / dB |
|---|---:|---:|---:|---:|---:|---:|---:|
| q00 | 2.622 | 3.313 | 5 | 9 | 12 | 2.353 | 4.36 |
| q4 | 2.62 | 3.385 | 5 | 10 | 13 | 1.882 | 4.42 |

Q4 exposes 1 additional spatial mode(s) at 30 dB.
 The Q4 median FWHM change relative to q00 is -20.0%.
 Additional-state energy outside the q00 observable subspace is 1.09% / 0.4% / 0.22% at 20/30/40 dB.

## Scope boundary

This kernel tests spatial focusing and mode observability in the already lossy O8 transfer. Seven frequency anchors are insufficient to estimate modal Q, decay time, or ringing. Those require a dense-band solve or a measured impulse response. Dynamic-range thresholds here are sensitivity diagnostics, not a replacement for the existing stochastic O8 detector replay.

Frequencies (kHz): `193.6, 209.733, 225.867, 242.0, 258.133, 274.267, 290.4`.
Pulse weights: `0.24004, 0.357156, 0.446502, 0.480079, 0.446771, 0.35743, 0.24004`.
