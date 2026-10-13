# Compact phase-wrapped Fresnel rings — `stop`

## Verdict

**`stop` for a phase-wrapped/Fresnel topology on the fixed benchmark.** The
negative result is geometric before it is numerical: the required edge-to-centre
phase correction is only `2.451 rad < pi`, and the first half-wave Fresnel
boundary is at `r = 40.962 mm`, outside the available `R = 35.8 mm` aperture.
There is therefore neither a binary Fresnel zone nor a `2pi` phase reset to
wrap. A compact kinoform/Fresnel plate degenerates into an ordinary unwrapped
annular delay plate in this stand.

The best manufacturable three-level plate was still evaluated as a useful
negative control. It preserves the five-cycle pulse and narrows the scalar
FWHM, but fails the local broadband transmission gate under even the optimistic
loss bracket (`min Tbar = 0.787 < 0.8`) and gives only `1.034x` carrier focus
amplitude relative to bare under nominal loss. G3 monostatic defect FEM was not
run after these G0/G1 no-go gates; all echo metrics are therefore `NaN`, not
substituted by a focal-field proxy.

This verdict does **not** reject all stepped annular delay plates. It says that
phase wrapping supplies no compactness advantage at `242 kHz`, `R=35.8 mm`,
`F=60 mm`; any continuation should be named and benchmarked as a conventional
quantised delay lens.

## Reproduction

Calculated implementation/results commit: `e777014` on
`codex/topology-fresnel-rings`.

```bash
cd /Users/mixx3/Desktop/betters/metamaterials/.worktrees/fresnel_rings
julia --startup-file=no --project=. test/topology_fresnel_rings.jl
julia --startup-file=no --project=. src/code/metamaterial/topology_fresnel_rings_run.jl
test -e tmp/gap_pilot_220khz || ln -s /Users/mixx3/Desktop/betters/metamaterials/tmp/gap_pilot_220khz tmp/gap_pilot_220khz
test -e tmp/period_length_pilot_220khz || ln -s /Users/mixx3/Desktop/betters/metamaterials/tmp/period_length_pilot_220khz tmp/period_length_pilot_220khz
julia --startup-file=no --project=. test/runtests.jl
```

The parameter screen is regenerated in
`tmp/topology_benchmark/fresnel_rings/manufacturing_screen.csv` and remains
ignored raw output. The small gate summaries in this directory are tracked.
The two symlinks only expose pre-existing ignored fixtures required by the
repository's unrelated legacy target-frequency tests; they are not inputs to
the Fresnel calculations. With them present, the full suite passes.

## G0 — geometry and manufacturing

For a planar aperture focused on-axis in Al 6061,

```math
\Delta \ell = \sqrt{F^2+R^2}-F=9.8687\ \mathrm{mm},
```

```math
\Delta \tau = \Delta \ell/c_P=1.61198\ \mathrm{us},\qquad
\Delta\phi=2\pi f_0\Delta\tau=2.45107\ \mathrm{rad}.
```

The carrier period is `4.13223 us`. Hence `floor(Delta tau/T0)=0`. More
strongly, the first phase-reversal boundary,

```math
r_1=\sqrt{(F+\lambda_P/2)^2-F^2}=40.9618\ \mathrm{mm},
```

lies outside the `35.8 mm` radius. The aperture contains only the central
Fresnel half-zone.

The screen covered 2, 3 and 4 levels, core depths `6.5--18 mm`, and matching
skins `1.5--3.5 mm`. The best lossless worst-ring transmissions were:

| levels | core | each matching skin | worst `Tbar`, lossless | reason |
|---:|---:|---:|---:|---|
| 2 | 7.0 mm | 2.375 mm | 0.874 | passes, but binary phase quantisation |
| 3 | 7.0 mm | 2.375 mm | 0.874 | selected; retains an intermediate delay |
| 4 | 9.0 mm | 2.625 mm | 0.676 | internal slab mismatch fails even lossless G1 |

### Fixed primary build

The selected all-solid coupon is a `71.6 mm` diameter, `7.0 mm` thick PEEK
core with concentric recesses turned from one face and vacuum-filled with cast
polyurethane. Two tungsten-loaded epoxy matching plates are bonded and lapped
flat. The assumed properties are:

| material | density | longitudinal speed | nominal attenuation at 242 kHz |
|---|---:|---:|---:|
| Al 6061 termination | 2700 kg/m3 | 6122.1 m/s | 0 dB/cm |
| PEEK | 1300 kg/m3 | 2530 m/s | 0.8 dB/cm |
| cast polyurethane | 1100 kg/m3 | 1500 m/s | 2.0 dB/cm |
| W-epoxy matching plate | 3000 kg/m3 | 2500 m/s | 0.5 dB/cm |

Attenuations are an explicit engineering bracket, not measured data. The
script also evaluates zero loss, `0.25x` (optimistic), and `3x` (pessimistic).
They must be replaced by coupon measurements before any fabrication decision.

| radial interval | PU depth | remaining PEEK | target relative delay |
|---|---:|---:|---:|
| `0--17.382 mm` | 5.939 mm | 1.061 mm | 1.612 us |
| `17.382--30.708 mm` | 2.970 mm | 4.030 mm | 0.806 us |
| `30.708--35.800 mm` | 0 mm | 7.000 mm | 0 us |

Each W-epoxy skin is `2.375 mm`; total depth is `11.75 mm`. The smallest
axial PEEK remainder is `1.061 mm`, and the narrowest radial ring is
`5.092 mm`, both above the `0.8 mm` manufacturing gate. A practical route is:

1. finish-turn the PEEK OD and two concentric recess depths in one setup;
2. plasma-treat, vacuum-fill with low-shrink PU, cure and face-grind;
3. mould/machine the two W-epoxy plates, bond in a parallel fixture, then lap;
4. clamp between the transducer shoe and Al specimen with the thinnest
   repeatable couplant film.

The calculation assumes perfect internal bonds and direct Al terminations; it
does not award any output renormalisation. Actual bond voids, shrinkage, and
contact films can only reduce confidence in the optimistic result.

## G1 — local broadband transfer

The model is an exact one-dimensional normal-incidence transfer-matrix stack
for scalar longitudinal pressure, evaluated over `162.2--321.7 kHz` and
weighted by the common five-cycle Hann spectrum. It includes every material
interface and attenuation but cannot represent elastic P/S conversion at ring
steps.

| loss bracket | ring 1 `Tbar` | ring 2 `Tbar` | ring 3 `Tbar` | worst `Bt` | worst `rho` | worst postcursor |
|---|---:|---:|---:|---:|---:|---:|
| lossless | 0.874 | 0.902 | 0.909 | 1.069 | 0.986 | 0.041 |
| optimistic (`0.25x`) | 0.787 | 0.830 | 0.864 | 1.059 | 0.986 | 0.038 |
| nominal | 0.578 | 0.650 | 0.744 | 1.059 | 0.985 | 0.031 |
| pessimistic (`3x`) | 0.261 | 0.345 | 0.500 | 1.040 | 0.982 | 0.019 |

The temporal gates (`Bt<=1.25`, `rho>=0.90`, postcursor `<=0.10`) pass. The
power gate does not: even the optimistic central ring is just below `0.8`, and
the nominal minimum is `0.578`. The nominal differential group delay between
centre and edge is `2.413 us`, rather than the geometric `1.612 us`, because
the multilayer Fabry response is dispersive, although the carrier phase itself
is close to the target.

## G2 — equal-input aperture proxy

G2 is an **axisymmetric scalar Rayleigh--Sommerfeld bulk-P proxy**, not the
shared 2D plane-strain elastic FEM and not a 3D solid model. The same uniform
incident aperture is used for every row. Local transmission remains in the
field integral. The `straight_footprint` cassette uses the same depth and
area-averaged PU inventory without radial phase variation. The
`unwrapped_continuous` control uses the same material stack with a continuous
unwrapped PU-depth profile.

| topology | weighted aperture transmission | focus amplitude / bare | focus energy / bare | FWHM | sidelobe | `Bt` | `rho` | postcursor |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| bare | 1.000 | 1.000 | 1.000 | 39 mm | n/a | 0.921 | 0.986 | ~0 |
| straight footprint | 0.642 | 0.842 | 0.710 | 39 mm | n/a | 1.040 | 0.990 | 0.021 |
| continuous unwrapped | 0.542 | 0.979 | 0.958 | 34 mm | 0.109 | 1.020 | 0.987 | 0.006 |
| selected 3-level plate | 0.658 | 1.034 | 1.069 | 33 mm | 0.125 | 1.030 | 0.990 | 0.016 |

The three-level plate is `1.227x` (`+1.78 dB`) in carrier amplitude relative
to its straight cassette, but only `1.034x` (`+0.29 dB`) relative to bare. Its
15% FWHM reduction is real within the scalar proxy, but is not enough to
override G0/G1 or establish defect detectability.

Chromatic observables at the fixed `60 mm` target are:

| frequency | focus amplitude / bare | spatial phase coherence | FWHM | sidelobe |
|---:|---:|---:|---:|---:|
| 162.2 kHz | 0.630 | 0.983 | 48 mm | 0.063 |
| 193.8 kHz | 0.760 | 0.975 | 42 mm | 0.094 |
| 242.0 kHz | 1.034 | 0.942 | 34 mm | 0.125 |
| 290.1 kHz | 1.036 | 0.900 | 28 mm | 0.156 |
| 321.7 kHz | 0.911 | 0.878 | 26 mm | 0.143 |

Thus the short pulse is not badly lengthened, but different spectral
components see materially different gain, spot width and coherence. The low
frequency half of the pulse is attenuated rather than usefully focused.

## G3/G4 and limitations

G3 was deliberately not run because:

- G0 proves there is no phase wrap or phase-reversal zone in the fixed aperture;
- G1 fails `Tbar>=0.8` under optimistic and nominal material loss;
- G2 offers only `+0.29 dB` carrier amplitude over bare in a scalar proxy.

Consequently `echo_mf_gain_bare`, `echo_energy_gain_bare`,
`echo_to_ringing_db`, dead-zone depth, CNR and POD remain `NaN`. No conclusion
is made about monostatic reciprocity, void scattering, P/S conversion, ring-edge
diffraction, or 3D finite-defect response. Those require coupled transient
elastic FEM, but are not justified for this stopped topology.

## Primary literature and transfer boundary

- Liu et al., [Design and experimental verification of phase-reversal Fresnel
  lens for contact stress characterization](https://doi.org/10.1016/j.ndteint.2024.103071),
  *NDT & E International* 144 (2024) 103071. This is the closest direct-contact
  solid NDT precedent: a five-ring copper lens was CNC-turned and experimentally
  tested. Its frequency, aperture, groove geometry and objective differ, so its
  gains are not imported into this benchmark.
- Dolmatov et al., [Application of Phase-Reversal Fresnel Zone Plates for
  High-Resolution Robotic Ultrasonic NDE](https://doi.org/10.3390/s21237792),
  *Sensors* 21 (2021) 7792. It experimentally supports PR-FZP use in NDE, but
  uses `5.55 MHz`, a localized coupling layer and a much smaller lens; it is not
  evidence for direct low-frequency Al-to-Al transmission here.
- Fuster et al., [Analysis of Fresnel Zone Plates Focusing Dependence on
  Operating Frequency](https://doi.org/10.3390/s17122809), *Sensors* 17 (2017)
  2809. Its measured underwater FZP establishes frequency-dependent focal
  position/distortion and motivates the chromatic G2 table.
- Pérez-López et al., [Transient Analysis of Fresnel Zone Plates for Ultrasound
  Focusing Applications](https://doi.org/10.3390/s20236824), *Sensors* 20
  (2020) 6824. It shows why transient duration must be checked for pulse-echo;
  this benchmark therefore uses `Bt`, correlation and postcursor rather than a
  carrier-only focus.

These papers establish that Fresnel lenses can work when multiple Fresnel zones
fit the aperture. They do not remove the present G0 condition. Increasing
frequency, increasing aperture, or reducing focus depth could create real
phase-reversal zones, but each would change the locked common stand and belongs
to a different research question.
