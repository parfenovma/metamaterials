# Sinusoidal lens material experiment at 242 kHz

## Question and objective

This experiment returns to the original sinusoidal strips and changes the
material interface.  It compares:

1. photopolymer elements bonded directly to an aluminium output plate;
2. the same elements coupled through an impedance-transforming layer;
3. aluminium elements bonded to aluminium.

The objective is the maximum longitudinal displacement amplitude at the target
point under the same `1 MPa` incident traction.  Transverse resolution and
sidelobe level are not penalized.  Every selected lens is compared with a
straight-strip aperture of the same width, length, materials and source load.

## Assumptions

The carrier is `242 kHz`, the source is a four-cycle Hann pulse in the reduced
search, and the aperture contains 15 strips.  Strip height is `7 mm`, slot
width is `1.2 mm`, and the nominal target lies `60 mm` into the aluminium after
the last material interface.

The photopolymer uses the existing project calibration:

```text
rho = 1210 kg/m^3, cP = 2340 m/s, cS = 1170 m/s
Rayleigh alpha = 79560 1/s, beta = 2.5e-9 s
```

These wave speeds imply a stiffer resin than the midpoint of the current
VeroClear tensile-modulus range.  They are retained to remain comparable with
the existing cell library, not presented as universal resin data.  The
[Stratasys VeroClear data sheet](https://www.stratasys.com/siteassets/materials/materials-catalog/polyjet-materials/veroclear/mds_pj_veroclear_0320a.pdf)
reports density `1.18--1.19 g/cm^3` and tensile modulus `2--3 GPa`.

Aluminium is a nominal 6061 material with `rho = 2700 kg/m^3`, `E = 68.3 GPa`
and assumed `nu = 0.33`, giving `cP = 6122.1 m/s` and `cS = 3083.8 m/s`.
Density and modulus follow the
[Kaiser Aluminum 6061 technical data](https://online.kaiseraluminum.com/depot/PublicProductInformation/Document/1025/Kaiser_Aluminum_6061_Rod_and_Bar.pdf).

The nominal matching material uses geometric means of density, P speed and S
speed.  Its normal-incidence P impedance is therefore the geometric mean of
the two bounding impedances.  The plane-wave quarter-wave thickness is
`3.910 mm`; this value is a starting hypothesis, not imposed on the FEM result.

## Search model

`sinusoidal_material_lens.jl` uses a segmented longitudinal transfer matrix for
each sinusoidal strip and a broadband scalar propagation kernel in aluminium.
All strips in a physical lens share one length and one output plane.  The
search varies common length, the period count (`1` or `2`) and the minimum gap.
It maximizes the carrier field, retains a shortlist, and re-ranks it by the
peak of the reconstructed four-cycle pulse.

The reduced model predicted the following pulse peaks:

| configuration | peak | gain over straight aperture | selected length |
|---|---:|---:|---:|
| polymer direct | 1.658 | 1.565 | 12 mm |
| polymer + nominal match | 2.212 | 1.556 | 12 mm |
| aluminium | 1.184 | 1.528 | 28 mm |

These amplitudes are reduced-model units.  They screen geometry but do not rank
the real material configurations reliably, because the transfer matrix does
not contain guided-mode conversion, oblique incidence or inter-strip coupling.

## Full-lens FEM result

The heterogeneous 2D elastic FEM uses conforming bonded interfaces and
separate `Lens`, `MatchingLayer` and `Aluminium` domains.  A linear fine-mesh
pilot gave:

| configuration | selected `|ux|` at target, m | straight control, m | gain |
|---|---:|---:|---:|
| polymer direct | 1.258e-8 | 8.296e-9 | 1.516 |
| polymer + 3.5 mm match | 1.158e-8 | 1.181e-8 | 0.980 |
| aluminium | 2.104e-8 | 1.288e-8 | 1.633 |

Thus the nominal quarter-wave layer improves the straight polymer-to-aluminium
aperture but invalidates the phase profile chosen by the reduced model.  It is
an impedance transformer, not a phase-neutral add-on.

The leading aluminium design was repeated with quadratic elements on a coarser
geometric mesh.  It retained the lead and focusing gain:

```text
selected |ux(target)| = 2.2639e-8 m
straight |ux(target)| = 1.6029e-8 m
gain                  = 1.4124
```

The linear calculation therefore overestimated the gain, but the conclusion
did not change.

## Matching-layer thickness sweep

The direct-polymer focusing profile was held fixed while layer thickness was
swept in full FEM.  The target amplitude is maximized near `1.4 mm`:

```text
t = 0.0 mm  -> 1.2580e-8 m
t = 1.0 mm  -> 1.3567e-8 m
t = 1.4 mm  -> 1.4095e-8 m   (best sampled)
t = 1.5 mm  -> 1.3983e-8 m
t = 2.0 mm  -> 1.3764e-8 m
t = 3.0 mm  -> 1.1301e-8 m
t = 3.5 mm  -> 1.1093e-8 m
t = 4.0 mm  -> 1.0997e-8 m
t = 5.0 mm  -> 0.8500e-8 m
```

The best layer is only about `0.36` of the plane-P quarter wavelength.  This is
consistent with the full aperture containing oblique, guided and shear content
that a one-dimensional normal-incidence formula cannot match simultaneously.

## Recommended prototype

For maximum amplitude at the prescribed target, use the aluminium-to-aluminium
15-strip design in `selected_aluminium.csv`.  It uses a common `28 mm` length;
the central five strips and outermost pair are straight, while the intermediate
pairs use one- and two-period sinusoidal constrictions with `1.5--2.0 mm`
minimum gaps.  Treat the exact geometry as a numerical candidate until the
quadratic result is checked on one additional mesh.

If polymer printing is required, use the direct-polymer profile with an
approximately `1.4 mm` matching layer as the next continuation seed.  Its
profiles must then be re-optimized with the layer present; attaching the layer
after an isolated-cell optimization is not sufficient.

## Reproduction

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_sinusoidal_material_lens.jl
julia --startup-file=no --project=. src/code/metamaterial/run_material_lens_fem.jl \
  --solve aluminium:selected
julia --startup-file=no --project=. src/code/metamaterial/run_matching_layer_fem_sweep.jl \
  --thickness=1.4
```

Primary outputs are under:

- `tmp/sinusoidal_material_lens_242khz`;
- `tmp/sinusoidal_material_lens_fem_242khz`;
- `tmp/sinusoidal_material_lens_fem_q2_242khz`;
- `tmp/matching_layer_fem_sweep_242khz`.

## Five-cycle and defect continuation (P6)

The synthesis was repeated with the project-standard `5-cycle Hann @ 242 kHz`
pulse.  The direct-polymer and aluminium geometries were unchanged exactly.
Only the nominally matched polymer profile changed; its updated fine linear
FEM target amplitude is `8.086e-9 m`, versus `1.181e-8 m` for its straight
reference.  The matching layer therefore remains a phase-changing part of the
design rather than a universal post-processing improvement.

The aluminium winner was closed with an additional refined quadratic mesh:

```text
quadratic coarse:   selected 2.2639e-8 m, straight 1.6029e-8 m, gain 1.4124
quadratic refined:  selected 2.2713e-8 m, straight 1.6159e-8 m, gain 1.4056
```

The selected amplitude changes by only `+0.33%` and the gain by `-0.48%`.
At `242 kHz`, the bulk aluminium wavelengths are `lambdaP = 25.298 mm` and
`lambdaS = 12.743 mm`; the refined quadratic lens/output sizes are
`0.0471 lambdaS` and `0.1569 lambdaS` respectively.

For the first NDE pilot, a free cylindrical void of radius `2 mm` was placed
at `(88,0) mm`.  Relative to the matched straight aperture, the lens gives:

```text
defect-free local kinetic-energy density gain  1.9756
receiver-line RMS scattered ux gain            1.4406
equal-noise SNR proxy improvement               3.17 dB
```

The receiver line is `20 mm` before the defect.  Scattering is the complex
defect field minus the complex defect-free field on a common physical grid.
The SNR number assumes the same additive receiver noise; no stochastic noise
or strength model is included.

The permanent, reproducible result package is
`results/sinusoidal_material_transfer_p6_242khz/`.

## Broadband aluminium gate

The frozen Al→Al lens and matched straight reference were then solved at 17
full-FEM frequency points across `162.2--321.7 kHz`, with half spacing inside
`B_-6dB`.  The band contains `99.815%` of the five-cycle pulse energy.  Complex
interpolation to the FFT grid gives:

```text
selected peak       21.65 nm
straight peak       15.88 nm
G_peak               1.3634
carrier gain         1.6334
B_t                  0.8643
rho_pulse            0.9794
postcursor           0.0644
```

Frequency-grid refinement changes `G_peak` by only `-0.024%`.  Pulse fidelity
passes, but the amplitude gate `G_peak >= 2` fails.  An ideal scalar TTD mask
with the same aperture and focus gives `G_peak = 3.182`, so aperture size is
not the limiting mechanism; the real sinusoidal strip transfer functions are.
The next step is an eight-group physical response matrix and local geometry
Jacobian, not a dense blind sweep.

That diagnostic was completed at the carrier. Eight symmetric source-group
bases reconstruct the full lens and straight fields exactly. The lens has
coherence efficiency `0.5473`, maximum focal phase error `147.5 deg`, and a
phase-aligned same-amplitude ceiling of `2.984`. The pre-output 8x8 response
has off-diagonal-L2/diagonal `0.557`, so independent-cell synthesis is not
valid for this fused aluminium aperture.

A single coupled Jacobian step increased the gap of the worst group at
`|y|=49.2 mm` from `1.50` to `1.75 mm`. Gain fell `1.633 -> 1.541`, focus
fell `5.68%`, and its phase error worsened `-147.5 -> -155.0 deg`. The useful
gap direction is below the frozen minimum feature and has insufficient local
phase range. The gap regulator is therefore stopped; a new topology must
separate broadband phase/delay control from transmission.
