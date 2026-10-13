# Sinusoidal aluminium waveguides: scalar broadband time reversal

## Result

The saved 17-anchor full-FEM channel was time reversed with equal replay energy. The sinusoidal channel gives a `1.331x` larger refocused envelope peak than the matched straight strips. Its half-power width is `7.076 microseconds` versus `7.283 microseconds`.

| geometry | equal-energy peak / m | half-power width / us | energy in half-power lobe | ringing energy outside 5-cycle window | far peak / dB |
|---|---:|---:|---:|---:|---:|
| sinusoidal | 2.6524515174909296e-9 | 7.076 | 79.76% | 1.95% | -18.6 |
| straight | 1.9921152896231442e-9 | 7.283 | 75.99% | 8.13% | -13.64 |

The sinusoidal path leaves `1.95%` of refocused energy outside the centred five-cycle window, compared with `8.13%` for the straight control. Thus its advantage in this SISO experiment is temporal matched-filter concentration, even though the separate carrier `8 x 8` audit found worse spatial-channel conditioning.

## Boundary

This is a reciprocal scalar record-and-replay calculation from the common uniform input to the design focus. It neither demonstrates added spatial modes nor estimates modal Q. The 17-point response was already validated for the five-cycle waveform, but a broadband spatial time-reversal claim still requires the full response matrix at every frequency and odd as well as even source groups.
