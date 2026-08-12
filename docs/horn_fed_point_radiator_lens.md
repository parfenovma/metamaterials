# Horn-fed point-radiator lens at 242 kHz

## Objective

This experiment tests the research-plan B1 topology

```text
wide input collector → narrow guide → small radiator → common receiver
```

for a five-cycle elastic pulse.  The collector is intended to raise the local
velocity in a narrow guide while retaining useful radiated P-like energy.  It
is passive: all claims are therefore normalized either to an abrupt transition
with the same driven area or to a straight guide with the same unfolded path.

## Deterministic collector and guide

The photopolymer pilot uses `h_in = 3.2 mm` and `h_throat = 1.6 mm`.  Its
log-cosine height is

```math
q(x)=\frac{1-\cos(\pi x/L_h)}{2},\qquad
h(x)=h_{in}\exp\left(q(x)\ln\frac{h_{throat}}{h_{in}}\right).
```

The `B_-20dB` lower edge is `162.2 kHz`.  With `cP = 2340 m/s`, the WKB rule
in the plan gives `L_h = 25.0 mm` and `epsilon_ad = 0.049999`.  No horn-profile
sweep was used.

The first rounded guide packed a `30.42 mm` path into `20 mm` axially.  Its
`1.08 mm` inner radius retained too much P-to-S conversion.  One deterministic
correction increased the axial length to `26 mm`; the `3.93 mm` inner radius
then exceeded half the S wavelength at the lower pulse-band edge.  The path
length remained `30.42 mm`.

## Single-element transient gate

All cases use the same five-cycle Hann pulse, current Rayleigh losses, a common
homogeneous receiving block, and an on-axis probe `35 mm` from the small
radiator.  The equal-path log-cosine collector plus straight `30.42 mm` guide
is the gate reference.

| case | target peak / reference | target energy / reference | `B_t` | correlation | transverse / longitudinal energy |
|---|---:|---:|---:|---:|---:|
| abrupt `3.2 → 1.6 mm` | 0.626 | 0.550 | 1.500 | 0.896 | 0.0001 |
| collector + straight equal path | 1.000 | 1.000 | 1.000 | 1.000 | 0.0002 |
| collector + tight rounded guide | 0.769 | 0.596 | 0.941 | 0.995 | 0.157 |
| collector + gentle rounded guide | **0.940** | **0.901** | **1.020** | **1.000** | **0.047** |

The collector raises throat velocity by `1.625x` relative to the abrupt step
and the radiated on-axis peak by about `1.60x`.  The gentle delay guide passes
the plan's primary `peak >= 0.8` gate and the pulse-fidelity gates.  Local
throat amplification is reported only as a diagnostic; the gate is based on
the radiated pulse.

## Measured binary state and lens

For a common output plane at `x = 51 mm`, the zero state uses a straight
`26 mm` guide and the delay state uses the gentle `30.42 mm` path.  The complex
coefficient measured `2 mm` inside the receiver at `241.94 kHz` is

```text
gentle / device = 0.8592 ∠ -191.44 deg
group delay      = 2.94 us
```

A 15-radiator binary mask was selected at pitch `4.8 mm`, radiator width
`1.6 mm`, and focal distance `35 mm`.  The reduced measured-element model
predicted focus gain `2.226` and transverse FWHM `5.2 mm`.

The selected mask was then applied as complex traction to 15 independently
tagged small radiators on a common fully elastic photopolymer receiver.  The
comparison is a uniform point aperture with the same mesh, radiators, total
nominal load and absorbing boundary.

| FEM | lens `|ux|` at target | uniform `|ux|` | gain | transverse FWHM | first sidelobe / peak |
|---|---:|---:|---:|---:|---:|
| linear, fine mesh | 64.06 nm | 46.34 nm | 1.382 | 6.0 mm | 0.524 |
| quadratic, coarse geometry | 80.99 nm | 37.73 nm | 2.147 | 3.0 mm | — |
| quadratic, fine geometry | **81.16 nm** | **37.38 nm** | **2.171** | **3.0 mm** | **0.500** |

The quadratic target amplitude changed by only `0.21%` when the geometry mesh
was refined from `6852` to `14738` nodes.  The central peak is at `x = 34 mm`,
`y = 0`, close to the requested `(35, 0) mm`.  Its target-centred axial
`-3 dB` span is `11.5 mm`.

The `3 mm` central FWHM is accompanied by sidelobes around `0.50` of the main
amplitude.  It should be interpreted as a binary-aperture/superoscillatory
trade-off, not free sub-diffraction resolution.  The stated objective did not
penalize sidelobes.

## Scope and remaining gate

The common-receiver calculation is a full elastic aperture propagator, but it
is not yet an end-to-end 3D array of collectors and guides.  The measured
single-channel complex coefficients are imposed at the radiator patches.
Planar gentle paths have an approximately `8 mm` excursion and cannot be
stacked at `4.8 mm` pitch without overlap.  A physical lens therefore needs
out-of-plane routing, layered fabrication, or a new compact connection whose
neighbour coupling is explicitly simulated.

Absolute nanometre amplitudes must not be compared directly with the earlier
aluminium lens: receiver material, focal distance, aperture and driven areas
differ.  The robust comparison here is lens versus the matched uniform
point-aperture control.  The aperture calculation is harmonic; a full-array
five-cycle transient remains required even though both isolated binary states
passed their transient gates.

## Out-of-plane long-smooth continuation

The planar overlap was removed by routing bends in the `x-z` plane while the
aperture pitch remains along `y`. A short `sin^4` centerline removed the
curvature discontinuity at both ends but concentrated too much curvature in
the middle. Its five-cycle result (`peak=0.879`, transverse energy `0.225`)
was therefore rejected despite passing the first carrier neighbour screen.

The accepted profile was sized analytically, not swept. It preserves the
same `4.42 mm` extra path and uses the minimum axial length for which its inner
radius is no smaller than the successful cosine guide's `3.93 mm`:

```text
axial length       38.158 mm
unfolded path      42.578 mm
centerline height   7.797 mm
inner radius        3.930 mm
```

Against its straight equal-path control the long-smooth transient gives
`peak=0.967`, energy `0.956`, `B_t=0.990`, correlation `0.9988`, and transverse
energy ratio `0.0043`. The full local binary 3D carrier basis (`DDD`, `DSD`,
`DSS/SSD`, `SDD/DDS`, `SDS`, `SSS`) was then evaluated. In `D-S-D`, the
neighbour perturbation relative to isolated channels is only `-5.0%` and
`-4.6 deg`; the two straight neighbours change by less than `5.6%` and
`1.1 deg`.

The broadband `smooth/device` transient was anchored to these context-specific
3D carrier coefficients. Direct optimization of the five-cycle focus chose

```text
SDDSSDDDDDSSDDS
```

and produced `G_peak=2.059`, `B_t=1.091`, correlation `0.988`, postcursor
`0.065`, and impulse FWHM `6 mm` relative to the matched uniform device
aperture. A quadratic full-elastic common-receiver FEM using the same 15
context weights gave:

| metric | long-smooth context lens | previous cosine lens |
|---|---:|---:|
| `|ux(35,0)|` | **89.86 nm** | 81.16 nm |
| gain / uniform | **2.404** | 2.171 |
| transverse FWHM | 4.0 mm | 3.0 mm |
| sidelobe / peak | **0.317** | 0.500 |
| local axial peak | 33.5 mm | 34.0 mm |

This is a `10.7%` increase in target amplitude and gain. It deliberately
trades `1 mm` of transverse FWHM and `12.2 mm` of guide length for higher
amplitude, lower sidelobes, and much lower mode conversion. The remaining
qualification is broadband context dependence and a common-array transient;
the current context corrections are measured in 3D only at the carrier.

## Monotonic true-time-delay lens

The binary mask was a two-state phase-quantized proof, not the final geometry.
It has now been replaced by a symmetric monotonic delay law

```text
tau(y) = [sqrt(F^2 + y_edge^2) - sqrt(F^2 + y^2)] / c_receiver .
```

For `F=35 mm`, pitch `4.8 mm`, and 15 radiators this gives `0` delay at
`y=+-33.6 mm` and `5.777 us` at the centre. Delay fractions are converted to
physical extra path and each path length is inverted to a separate `sin^4`
bend amplitude. All outputs remain coplanar; the bend is routed in `x-z` and
the aperture pitch remains along `y`.

The amplitude-first continuation of the already qualified `4.42 mm` state
selected a maximum extra path of `9.724 mm`. The common section and its worst
geometry are:

```text
common axial length       50.664 mm
maximum unfolded length   60.388 mm
maximum bend height       13.746 mm
minimum inner radius       3.930 mm
```

The maximum state was then solved directly in the same full-elastic transient
model. Against its straight equal-path control it gives peak `0.977`, energy
`0.962`, `B_t=0.991`, correlation `0.9998`, postcursor `0.0928`, and transverse
energy ratio `0.0055`. Its measured group delay is `5.516 us`; this physical
endpoint, rather than the extrapolated reference state, calibrates all
intermediate monotonic states by log-transfer continuation.

The resulting five-cycle aperture gives `G_peak=2.904`, carrier scalar gain
`2.695`, `B_t=0.899`, correlation `0.99998`, postcursor `0.056`, impulse FWHM
`5 mm`, and sidelobe ratio `0.329`. The matched quadratic common-receiver FEM
comparison is:

| metric | monotonic TTD | best binary | uniform |
|---|---:|---:|---:|
| `|ux(35,0)|` | **93.17 nm** | 89.86 nm | 37.38 nm |
| gain / uniform | **2.493** | 2.404 | 1.000 |
| transverse FWHM | 5.0 mm | 4.0 mm | 40.0 mm |
| sidelobe / peak | **0.294** | 0.317 | -- |
| local axial peak | 34.0 mm | 33.5 mm | -- |

Thus the monotonic lens adds `3.68%` carrier amplitude over the best binary
lens and reduces its sidelobe, at the accepted cost of `1 mm` FWHM. The next
primary milestone is integration, not another single-channel qualification:
assemble one non-intersecting 15-channel 3D lens with a common receiver and
obtain one end-to-end carrier focus. The `80%` equal-path target, continuous
neighbour corrections and common-array transient remain iteration metrics
after this v0 lens exists.

That v0 array has now been assembled and solved. It contains all 15 physical
collectors and monotonic guides fused to one `80 x 24 x 55 mm` receiving block.
The coarse integration mesh has `59,751` nodes and `332,095` elements. All 15
input faces, the common elastic domain and the absorbing receiver boundary are
present in one mesh.

The end-to-end 3D carrier solve produces a central beam at the requested
`x=35 mm` plane: `|ux|=3.893 nm`, transverse FWHM `5 mm`, center/outer
amplitude ratio `5.28`, and sidelobe ratio `0.372`. The narrowest sampled waist
is at `x=30 mm` with FWHM `3 mm`; the output-phase fit itself corresponds to
`F=36.8 mm`, so the axial shift is not interpreted as a simple delay-scale
error on this coarse receiver mesh. Actual channel-output phases differ from
the reduced-model target by at most `12.7 deg`. This satisfies the v0
integration milestone. A full 3D uniform reference is the first comparison of
the next iteration, not a prerequisite for accepting the assembled lens. That
reference has now also been solved with 15 straight guides, identical sources,
axial thickness and receiver. It gives `|ux(35,0,0)|=2.645 nm` versus
`3.893 nm` for the monotonic lens, hence the matched end-to-end carrier gain is
`1.472`. The earlier `center/outer=5.28` remains a spatial contrast, not gain.

## Reproduction

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_point_radiator_pilot.jl

julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_binary_lens.jl --stage=design
julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_binary_lens.jl --stage=mesh
METAMATERIALS_HORN_LENS_ORDER=2 \
METAMATERIALS_HORN_LENS_FINE_MESH=true \
julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_binary_lens.jl --stage=solve
METAMATERIALS_HORN_LENS_ORDER=2 \
julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_binary_lens.jl --stage=plot
```

Primary outputs are in `tmp/horn_point_radiator_pilot` and
`tmp/horn_binary_lens_242khz`. The out-of-plane continuation is reproduced by
`run_horn_neighbor_3d_gate.jl`, `run_horn_long_smooth_pilot.jl`,
`run_horn_long_smooth_aperture.jl`, and `run_horn_long_smooth_fem.jl`; its
outputs are in `tmp/horn_neighbor_3d_242khz`, `tmp/horn_long_smooth_pilot`, and
`tmp/horn_long_smooth_aperture`.

The monotonic continuation is reproduced by
`run_horn_monotonic_aperture.jl`, `run_horn_monotonic_pilot.jl`, and
`run_horn_monotonic_fem.jl`. Primary outputs are in
`tmp/horn_monotonic_aperture` and `tmp/horn_monotonic_pilot`.

The full integrated v0 array is reproduced by
`run_horn_monotonic_array_3d.jl`; its mesh, carrier field, channel outputs and
focus map are in `tmp/horn_monotonic_array_3d`.
