# Ellipse-plate carrier and second-harmonic screening

## Bottom line

The original 242 kHz value was a common project benchmark, not an optimum of
the ellipse plate.  A same-frequency comparison against a solid straight
footprint finds a repeatable linear pass/focusing window near 300--322 kHz for
case `c005` (`matching_aspect_xy=1.6`, `central_area_scale=0.2`,
`matching_offset_fraction=0.28`).  The best sampled fine-mesh point is 315 kHz:

| frequency | J1 | reduced-CNR J2 | J2 (dB) | exit-energy proxy | FWHM at nominal x | local peak error |
|---:|---:|---:|---:|---:|---:|---:|
| 315 kHz | 1.2057 | 1.3859 | +2.835 dB | 0.7090 | 32.52 mm | 11.0 mm upstream |

Thus frequency is a useful additional design variable.  The present geometry
does not yet give the requested broadband pulse-echo optimum: when the original
five-cycle 242 kHz Hann spectrum is used to weight the coarse sweep, all three
tested plates remain below the straight reference in both aggregate proxies.
The 315 kHz point is a narrowband, steady-state result and must not be read as a
transient CNR prediction.

## What was calculated

1. A nine-point coarse sweep from 162.2 to 321.7 kHz, spanning the project's
   five-cycle Hann `B-20 dB` band and explicitly sampling its 193.8--290.1 kHz
   `B-6 dB` interval.
2. Three representative geometries (`c001`, `r025`, `c005`) normalized to a
   solid straight footprint at every frequency.
3. Fine-mesh validation at 305 and 321.7 kHz for `c001` and `c005`.
4. A local fine sweep of `c005` at 295, 300, 305, 310, 315 and 321.7 kHz.
5. Reciprocal mixed-frequency proxies for nonlinear second-harmonic paths
   157.5 -> 315 kHz and 242 -> 484 kHz.

The fine lens meshes contain approximately 60k nodes and 121k elements; the
matching straight mesh contains approximately 20k nodes and 40k elements.

## Local fine sweep

| frequency (kHz) | J1 | J2 | J2 (dB) | min(J1,J2) | exit energy | FWHM (mm) | focus error (mm) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 295.0 | 0.9878 | 1.0216 | +0.186 | 0.9878 | 0.6301 | 39.54 | 7.62 |
| 300.0 | 1.1095 | 1.1811 | +1.445 | 1.1095 | 0.6869 | 38.01 | 7.00 |
| 305.0 | 1.1929 | 1.2859 | +2.184 | 1.1929 | 0.7581 | 37.88 | 8.00 |
| 310.0 | 1.1963 | 1.3116 | +2.356 | 1.1963 | 0.7806 | 36.66 | 10.00 |
| **315.0** | **1.2057** | **1.3859** | **+2.835** | **1.2057** | 0.7090 | 32.52 | 11.00 |
| 321.7 | 1.0454 | 1.1318 | +1.075 | 1.0454 | 0.6175 | 27.69 | 13.00 |

The sampled joint-gain window is about 300--321.7 kHz.  The local field maximum
moves upstream as frequency rises; at 315 kHz it is at `x=104 mm`, versus the
nominal `x=115 mm`.  A plate redesigned at the new carrier should therefore
retune the delay profile rather than merely drive the old geometry harder.

## Harmonics

The FEM model is linear.  Its "harmonic" solve means one sinusoidal steady-state
frequency; it cannot generate 2f or 3f.  To screen a nonlinear crack without
introducing an uncalibrated contact law, reciprocity was used:

`s_2f ~ u_f^2 * G_2f`,

so the relative defect-signal proxy is `g_f^2 * g_2f`.  A reduced mixed-frequency
CNR proxy was built from the same spatial kernel.

| transmit -> receive | g(f) | g(2f) | mixed signal | mixed signal (dB) | mixed CNR | mixed CNR (dB) |
|---:|---:|---:|---:|---:|---:|---:|
| 157.5 -> 315 kHz | 0.6975 | 1.2057 | 0.5865 | -4.635 | 0.7416 | -2.597 |
| 242 -> 484 kHz | 0.9393 | 0.3211 | 0.2833 | -10.955 | 0.3869 | -8.248 |

The ellipse plate therefore does not benefit from either tested second-harmonic
path.  The strong response near 315 kHz is a different *linear carrier window*,
not a harmonic of 242 kHz.

## Interpretation for an experiment

- A 310--315 kHz transducer or tunable drive is worth testing as a new
  fundamental carrier.
- The useful simulated window is much narrower than the `B-6 dB` width of a
  five-cycle Hann pulse centered at 315 kHz.  A longer/narrower burst would use
  the window more efficiently, but increases ringing and worsens axial/time
  separation.  Pulse length must therefore enter the next CNR optimization.
- Before fabrication, run a transient or complex-frequency reconstruction with
  the real piezo, matching layer, bond layer and receiver bandwidth.  That is
  needed for group delay, pulse broadening, ringing and true echo CNR.
- If nonlinear harmonic imaging is a primary goal, the topology needs two
  deliberately aligned passbands at f and 2f; this plate has only the useful
  upper band.

## Artifacts

- `figures/frequency_response.png`: full coarse carrier sweep.
- `figures/c005_local_frequency_fine.png`: local fine optimum.
- `figures/high_frequency_fine_validation.png`: coarse/fine comparison.
- `figures/mixed_second_harmonic_fine.png`: mixed f/2f proxies.
- `data/frequency_metrics_coarse.csv`: all coarse metrics.
- `data/c005_local_frequency_fine.csv`: local fine metrics.
- `data/high_frequency_fine_validation.csv`: convergence check.
- `data/mixed_second_harmonic_fine.csv`: mixed harmonic proxies.

## Reproduction

The pipeline reuses the existing geometry and meshes, changing only the drive
frequency.  Examples:

```sh
julia --startup-file=no ellipse_frequency_sweep.jl --case c005 315 fine
julia --startup-file=no ellipse_frequency_sweep.jl --straight 315 fine
julia --startup-file=no ellipse_frequency_sweep.jl --local-fine-report
julia --startup-file=no ellipse_frequency_sweep.jl --harmonic-report
```

