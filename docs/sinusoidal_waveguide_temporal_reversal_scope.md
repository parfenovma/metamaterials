# Sinusoidal Al waveguides: scalar broadband time reversal

## Exploratory question

Given the already frozen 17-frequency full-FEM response from the common
uniform lens input to the design focus, does a reciprocal record-and-replay
experiment produce a cleaner temporal focus through the selected sinusoidal
waveguides than through the matched straight strips?

This is an exploratory follow-up to the carrier spatial-channel audit. It uses
no new FEM solve and does not alter the stopped sinusoidal geometry campaign.

## Calculation

Let `H(f)` be the saved scalar transfer from the common input to the focus and
`S(f)` the project-standard five-cycle Hann pulse. The received record is

```math
R(f)=H(f)S(f).
```

The record is time reversed, normalized to unit replay energy, and sent back
through the reciprocal channel:

```math
Z(f)=H(f)\frac{R^*(f)}{\|r(-t)\|_2}
    =\frac{|H(f)|^2S^*(f)}{\|r(-t)\|_2}.
```

The same procedure and energy normalization are applied to the selected
sinusoidal lens and its matched straight control. Reported metrics are the
analytic-envelope peak, half-power width, energy inside that width, and energy
and maximum envelope outside a centred window equal to one input-burst
duration (`5/f0`).

## Scope boundary

- This is a scalar SISO temporal matched-filter test. It does not add spatial
  channels or contradict the separate `8 x 8` carrier result.
- The transfer is the interpolation of 17 full-FEM anchors over
  `162.2--321.7 kHz`; it was previously grid-refined for the five-cycle pulse,
  but it cannot resolve narrow modal linewidths or estimate Q.
- The calculation assumes reciprocity and a transducer at the original focus
  capable of recording and replaying the required waveform.
- It is not a defect-detection or probability-of-detection result.
