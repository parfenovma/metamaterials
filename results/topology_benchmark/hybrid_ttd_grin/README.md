# Hybrid TTD--GRIN topology benchmark

## Verdict

**`conditional`**, not a present winner.

The calibrated hybrid is physically compact, passes the isolated broadband
transfer and five-cycle pulse gates, and has a favourable reciprocal defect
proxy. It does not yet justify replacing a pure smooth Al TTD guide: on this
benchmark the calibrated pure path also fits inside `55 mm`, while the Al--Cu
insert increases pulse-weighted mode conversion from `0.30%` to `3.60%` and
has no measured loss/joint model. The strict proxy also places the hybrid echo
`7.32 dB` below its pessimistic device-ringing estimate, so G3 is not promoted
to a claimed monostatic PASS.

Continue only after a real graded Al--Cu coupon supplies measured complex
`S21/S11` and a coupled monostatic transient with the explicit `r=2 mm` void
confirms positive echo-to-ringing margin. Do not spend a full G4 sweep before
those two gates.

Validated code commit: `PENDING_FIRST_LOCAL_COMMIT`.

## G0 -- concrete geometry and manufacture

The 2D benchmark coupon is a stack of fifteen individually cut guide plates.
Each channel has a `3.2 mm` Al lead, a smooth zero-slope/zero-curvature
`sin^4` centreline, and a symmetric Al--Cu composition insert centred along
the route. All ports share the same `55 mm` input/output footprint.

| parameter | value |
|---|---:|
| aperture / channels | `71.6 mm / 15` |
| channel pitch | `4.885714 mm` |
| straight leads | `6 + 6 mm` |
| active axial route | `43 mm` |
| target maximum delay | `1.480305 us` |
| path / material split | `60 / 40%` |
| maximum extra path | `2.602573 mm` |
| maximum bend excursion | `6.199069 mm` |
| minimum inner bend radius | `5.955283 mm` |
| graded insert / Cu plateau | `10.200480 / 6.200480 mm` |
| minimum insert width | `1.384178 mm` |

The proposed lab route is:

1. tape-cast Al--Cu layers with increasing Cu volume fraction and consolidate
   them by hot pressing;
2. machine the impedance-compensating width profile and the Al carrier plates
   by wire EDM;
3. diffusion-bond or hot-press the insert into each carrier;
4. face-grind input/output planes in one setup and clamp the fifteen plates
   into a replaceable cassette.

The `0.4 mm` composition slices used by FEM are a laminate approximation of a
continuous tape-cast gradient, not a `0.4 mm` load-bearing ligament. The
minimum machined load path is `1.384 mm`, above the common `0.8 mm` gate.
Al--Cu intermetallic formation, porosity, residual stress and bond attenuation
are not represented and are major risks.

The mechanism is supported, but not validated for this coupon, by primary
work on [metal ultrasonic delay lines](https://doi.org/10.6028/jres.051.028),
[wide-passband FGM impedance transitions](https://doi.org/10.1016/j.physleta.2013.10.040),
and [tape-cast/hot-pressed Al--Cu graded material](https://doi.org/10.3390/ma15144834).
The use of cut and stacked metal plates as a manufacturable GRIN route also has
an experimental precedent in the [broadband pentamode lens](https://doi.org/10.1121/1.4985195),
although that device operated in water and is not evidence for solid-to-solid
bulk-wave performance here.

## Packing reduction

The first analytic seed used the old `4.793 km/s` Al-guide delay calibration
and over-delayed the channel. That uncalibrated raw run is intentionally kept
under `tmp/.../isolated_uncalibrated`: it measured `1.900 us` instead of
`1.480 us`. Exactly one measured calibration step was then applied.

On the calibrated common stand:

| comparison | pure smooth TTD | hybrid | reduction |
|---|---:|---:|---:|
| maximum extra path | `4.3376 mm` | `2.6026 mm` | `40.0%` |
| curvature-compliant excursion | `8.0548 mm` | `6.1991 mm` | `23.0%` |
| curvature-compliant depth | `53.9343 mm` | `55 mm` frozen footprint | none needed |

Relative to the historical wide-aperture all-Al seed (`23.03 mm` extra path),
the hybrid path is `88.7%` shorter, but this is context only: that seed used a
different `114.8 mm` aperture and must not be presented as the common-control
gain. The important result is that pure TTD is already feasible on the new
stand; the hybrid improves packing margin rather than unlocking feasibility.

## G1 -- isolated full-elastic channel

Nine full-elastic harmonic FEM anchors cover `162.2--321.7 kHz`. Rank-one
matched modal ports are used on the straight Al leads. Pulse metrics are
reconstructed with the common `5-cycle Hann @ 242 kHz` spectrum.

| variant | Tbar | Rbar | conversion | delay, us | ripple, us | Bt | rho | postcursor |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| straight | `1.00004` | `0.00010` | `<0.001%` | `0` | `0` | `1.0000` | `1.000000` | `0.05999` |
| path part only | `0.98523` | `0.01415` | `0.101%` | `0.92057` | `0.19310` | `1.00535` | `0.999959` | `0.06501` |
| hybrid | `0.95961` | `0.00977` | `3.602%` | `1.56771` | `0.31970` | `1.01070` | `0.999703` | `0.06064` |
| pure TTD | `0.96527` | `0.03213` | `0.299%` | `1.47926` | `0.15502` | `1.01070` | `0.999830` | `0.05920` |

All variants pass the declared isolated pulse gate. The hybrid is `5.9%`
over the delay target after the one allowed calibration step. No additional
fit is hidden in the aperture result.

## G2/G3 -- aperture and defect stand

A full coupled multi-material monostatic transient was not completed. G2/G3
therefore use a clearly labelled **reciprocal 2D scalar strict proxy**:

- the measured full-elastic G1 complex transfer is interpolated across each
  channel's delay fraction;
- `bare`, same-footprint straight, pure TTD and hybrid use equal incident
  energy;
- transmit and receive fields are multiplied by reciprocity;
- the defect is a point-scatterer approximation at `ka approximately 0.5`,
  not the explicit circular free boundary required by `void_r2_focus`;
- device ringing is reconstructed pessimistically from the finite-band FEM
  `S11`, set to zero outside the measured band. The resulting sinc tail is a
  screening bound, not a calibrated transducer baseline.

| topology | MF echo / bare | echo energy / bare | echo/device ringing | focus / bare | FWHM | sidelobe |
|---|---:|---:|---:|---:|---:|---:|
| bare | `1.000` | `1.000` | unavailable | `1.000` | `28.5 mm` | `0.500` |
| straight footprint | `0.635` | `0.403` | `+6.89 dB` | `0.798` | `29.5 mm` | `0.518` |
| pure TTD | `0.951` | `0.904` | `-15.14 dB` | `0.955` | `23.0 mm` | `0.322` |
| hybrid | `1.134` | `1.287` | `-7.32 dB` | `1.040` | `21.0 mm` | `0.216` |

The hybrid is the best focused/defect proxy of the manufactured candidates,
but it fails the echo-to-ringing screen. Because the proxy is not the explicit
common defect FEM, all unavailable primary NDT fields remain `NaN` in
`common_metrics.csv`; proxy values are kept separately in
`strict_proxy_metrics.csv`.

## Reproduction

From the repository root:

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_ttd_grin_isolated.jl --stage=mesh
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_ttd_grin_isolated.jl --stage=convert
METAMATERIALS_HYBRID_TTD_GRIN_WORKERS=2 \
  julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_ttd_grin_isolated.jl --stage=solve
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_ttd_grin_isolated.jl --stage=analyze
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_ttd_grin_aperture.jl
julia --startup-file=no --project=. test/topology_hybrid_ttd_grin.jl
julia --startup-file=no --project=. test/runtests.jl
```

Raw meshes, JSON models, per-frequency JLD2 points and plots remain under
`tmp/topology_benchmark/hybrid_ttd_grin/`. Only compact, reproducible summary
tables are committed here.

## Limitations

- The coupon and stand are 2D plane-strain surrogates, not axisymmetric rings
  or a full circular 3D aperture.
- The Al--Cu law is a Voigt mixture and has no measured damping, porosity,
  anisotropy, temperature or imperfect-bond model.
- The isolated FEM uses discrete constant-composition slices; a real FGM will
  differ.
- Common `none/void_r2_focus` monostatic transients, dead zone, explicit
  defect subtraction, CNR/POD and the robustness ensemble are not available.
- G4 was deliberately not run after the negative echo-to-ringing proxy.
