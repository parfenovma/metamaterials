# Hybrid TTD--GRIN: exact monostatic G3

## Verdict

**`conditional`**.

The explicit coupled plane-strain `none/void_r2_focus` calculation confirms
that the hybrid advantage in the earlier reciprocal scalar proxy was not a
sign error. At the refined decisive level, the hybrid matched-filter echo is
`0.970 dB` above the pure smooth TTD echo. The sign is stable from the coarse
level (`1.019 dB`), but the gain is below the predeclared `2.0 dB` CONTINUE
threshold.

The refined hybrid has positive echo-to-ringing (`+17.55 dB`, compared with
`+13.25 dB` for pure TTD) and essentially the same threshold-defined dead zone
(`54.738` versus `54.746 mm`). The hybrid shortens the *differential* geometric
delay path by `40%`. Including the `1.5 mm` common-mode path required only by
the non-intersecting 2D projection, the maximum excess path is `4.103 mm`
instead of `5.838 mm`, a still useful `29.7%` reduction. This satisfies the
predeclared `>=20%` packaging-benefit condition for a CONDITIONAL result.

Do not promote the topology to `continue` until a measured Al--Cu coupon or
joint supplies a complex loss model. A one-decibel numerical advantage is
small enough for unmodelled insert or bond loss to reverse it.

## What was actually solved

This is a full vector-elastic 2D plane-strain transient, not the scalar defect
proxy used in the prior hybrid report.

- Al 6061 specimen: `125 x 71.6 mm`, with `cP=6122.1024 m/s`,
  `cS=3083.8103 m/s`, and `rho=2700 kg/m3`.
- Explicit circular traction-free void: radius `2 mm`, centre at `x=60 mm`.
- Four topologies: `bare`, `straight_footprint`, `pure_ttd`, and `hybrid`.
- Paired `none` and `void_r2_focus` mesh and transient for every topology.
- `242 kHz`, five-cycle Hann traction burst.
- P/S characteristic absorbing traction on the source, lateral specimen, and
  downstream specimen boundaries; device and void surfaces are traction-free.
- Linear Newmark average acceleration (`beta=0.25`, `gamma=0.5`). The constant
  effective matrix is assembled/factorised once and reused for every step.
- The receiver is the aperture integral of normal velocity on the **same
  `Source` boundary that is driven**. No line near the defect enters a metric.
- Device sources consist of the fifteen `3.2 mm` guide entrances (`48 mm`
  active metal inside the `71.6 mm` aperture). Their pressure is multiplied by
  `sqrt(71.6/48)=1.221338`, so `p^2 A` and the nominal incident energy are the
  same as bare: `0.00419512 J/m`.

The source characteristic boundary can be interpreted as simultaneous
transmit/receive switching in a linear reciprocal model: the known drive is
removed by paired defect-free subtraction, while the aperture velocity is
retained throughout the record.

## Frozen analysis contract

The expected echo windows were chosen from source centroid plus geometric
two-way time of flight before comparing the topologies:

| topology | matched-filter / energy window |
|---|---:|
| bare | `24--40 us` |
| all `55 mm` devices | `45--68 us` |

For each topology,

```text
s0(t) = same-aperture normal velocity, defect-free
sd(t) = same-aperture normal velocity, explicit void
Delta s(t) = sd(t) - s0(t)
```

The matched filter is the normalised five-cycle source waveform. The ringing
denominator is the matched-filter peak of `s0` in exactly the same echo window;
echo energy and ringing energy use that same window. The dead zone is the
earliest point after the burst at which the analytic envelope of `s0` remains
below `-40 dB` of the driven maximum for at least one carrier cycle. Device
two-way delay is subtracted before converting this time to equivalent specimen
depth. Because no receiver-noise or material-clutter ensemble is present,
gain relative to bare is an equal-noise SNR proxy, not measured CNR or POD.

## Coarse common G3

The common control table uses `h_specimen=1.35 mm`, `h_guide=0.45 mm`,
`h_void=0.30 mm` and `18` samples per period.

| topology | MF echo / bare | echo energy / bare | echo/ringing | dead zone | MF peak time | G1 P/S conversion |
|---|---:|---:|---:|---:|---:|---:|
| bare | `1.000` | `1.000` | `+16.38 dB` | `63.245 mm` | `31.221 us` | `0%` |
| straight footprint | `1.089` | `1.194` | `+22.12 dB` | `57.644 mm` | `54.408 us` | `<0.001%` |
| pure TTD | `1.495` | `2.271` | `+12.68 dB` | `55.098 mm` | `55.785 us` | `0.299%` |
| hybrid | `1.681` | `2.830` | `+16.67 dB` | `55.968 mm` | `56.474 us` | `3.602%` |

Thus, on the common coarse stand the hybrid is `+4.51 dB` over bare and
`+1.02 dB` over pure TTD in matched-filter echo. `weighted_mode_conversion`
in `common_metrics.csv` remains the earlier full-modal isolated G1 value; the
new G3 table separately records the integrated tangential fraction on the
monostatic aperture and does not mislabel it as total converted power.

## Mesh/time-step convergence of the decisive comparison

The refined pair uses `h_specimen=1.00 mm`, `h_guide=0.32 mm`,
`h_void=0.22 mm` and `24` samples per period.

| topology | level | displacement DOF (`none`) | MF echo | echo/ringing | dead zone |
|---|---|---:|---:|---:|---:|
| pure TTD | coarse | `47,456` | `1.21205e-5` | `+12.68 dB` | `55.098 mm` |
| pure TTD | refined | `87,154` | `1.22441e-5` | `+13.25 dB` | `54.746 mm` |
| hybrid | coarse | `56,604` | `1.36293e-5` | `+16.67 dB` | `55.968 mm` |
| hybrid | refined | `96,690` | `1.36902e-5` | `+17.55 dB` | `54.738 mm` |

The pure echo changes by `+1.02%` and the hybrid echo by `+0.45%`. The hybrid
advantage changes only from `1.019` to `0.970 dB`; its sign is stable. Mesh and
time step are refined together, so this is a combined convergence gate rather
than separate error attribution.

## Geometry and the 2D packing surrogate

The previous single-channel G1 route can be stacked out of plane in separate
plates. A single 2D continuum cannot represent projected overlapping routes.
With zero common curvature, the fifteen variable-delay strips intersect near
the edge channels. Therefore both pure and hybrid device models receive the
same `1.5 mm` common-mode geometric excess path. It does not alter relative
delays, and straight remains straight; it only provides a non-intersecting
plane-strain projection. Minimum projected air clearance is `0.149 mm` for
pure and `0.721 mm` for hybrid.

Those small gaps are **not** proposed as the manufacturing route. The proposed
route remains individually wire-EDM-cut guide plates stacked out of plane.
Consequently this G3 tests coupled elastic pulse/void physics and relative
detectability, not final 3D plate packing or inter-plate cross-talk.

## Meshes and observed timings

| level / topology | nodes (`none`) | triangles (`none`) | Newmark steps | factorised solve time (`none`) |
|---|---:|---:|---:|---:|
| coarse bare | `6,459` | `12,537` | `349` | `3.23 s` |
| coarse straight | `23,031` | `41,588` | `349` | `7.29 s` |
| coarse pure | `23,728` | `42,899` | `349` | `4.51 s` |
| coarse hybrid | `28,302` | `47,901` | `349` | `5.83 s` |
| refined pure | `43,577` | `80,855` | `465` | `3.25 s` |
| refined hybrid | `48,345` | `86,152` | `465` | `3.22 s` |

The quoted solve time starts after spatial assembly and includes sparse LU plus
all time steps. Gmsh generation dominated this run: about `0.2--63 s` per
coarse mesh and `38--124 s` per refined device mesh on the test machine.

## Reproduction

From the repository root:

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_monostatic_run.jl \
  --stage=mesh --level=coarse --topology=all --case=all
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_monostatic_run.jl \
  --stage=convert --level=coarse --topology=all --case=all
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_monostatic_run.jl \
  --stage=solve --level=coarse --topology=all --case=all

julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_monostatic_run.jl \
  --stage=mesh --level=refined --topology=pure_ttd,hybrid --case=all
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_monostatic_run.jl \
  --stage=convert --level=refined --topology=pure_ttd,hybrid --case=all
julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_monostatic_run.jl \
  --stage=solve --level=refined --topology=pure_ttd,hybrid --case=all

julia --startup-file=no --project=. \
  src/code/metamaterial/topology_hybrid_monostatic_analysis.jl
julia --startup-file=no --project=. test/topology_hybrid_monostatic_g3.jl
julia --startup-file=no --project=. test/runtests.jl
```

Meshes, Gridap JSON models, JLD2 signals and manifests remain under
`tmp/topology_benchmark/hybrid_monostatic_g3/`. Compact committed outputs are:

- `common_metrics.csv`: common coarse controls in the benchmark schema;
- `g3_metrics.csv`: exact paired transient observables;
- `convergence_metrics.csv`: decisive coarse/refined comparison;
- `decision.csv`: frozen-gate verdict.

## Limitations and next physical gate

- This is 2D plane strain, not axisymmetric rings or a circular 3D aperture.
- The source is an ideal characteristic traction port; piezoelectric backing,
  electrical matching, contact layer and transmit/receive electronics are not
  represented.
- The Al--Cu insert uses the prior Voigt screening law. Porosity, anisotropy,
  intermetallics, attenuation and imperfect diffusion bonds are absent.
- Exact baseline subtraction suppresses deterministic geometry clutter; no
  stochastic grain/noise ensemble exists yet.
- The void is an infinite out-of-plane cylinder, not a spherical pore.
- Only the decisive pure/hybrid pair received combined mesh/time refinement.

The next justified experiment is therefore not another composition sweep. It
is a graded Al--Cu coupon with measured complex `S21/S11` across
`162--322 kHz`, including the actual bond. Those data should be inserted into
the refined G3 loss bracket. If the `~1 dB` hybrid advantage survives, proceed
to 3D plate packing/contact modelling; if it reverses, close this hybrid.
