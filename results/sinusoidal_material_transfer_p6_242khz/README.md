# P6 material transfer: sinusoidal lens in aluminium at 242 kHz

## Verdict

For the declared objective—maximum longitudinal displacement at a target
`60 mm` inside aluminium under the same `1 MPa` incident traction—the
`Al -> Al` sinusoidal lens is the numerical winner.  The selected aluminium
geometry is unchanged when the reduced synthesis is updated from a four-cycle
to the project-standard five-cycle Hann pulse.

On the same fine linear FEM prescription:

| configuration | selected `|ux|`, nm | matched straight, nm | gain |
|---|---:|---:|---:|
| photopolymer, direct bond | 12.580 | 8.296 | 1.516 |
| photopolymer + 3.5 mm nominal matching layer | 8.086 | 11.813 | 0.685 |
| aluminium elements into aluminium | 21.044 | 12.884 | 1.633 |

The fixed direct-polymer geometry peaks at `14.095 nm` when the matching layer
is swept to `1.4 mm`; this remains well below the aluminium candidate.  The
layer is not phase-neutral, so a polymer profile optimized without it cannot
simply receive a quarter-wave layer after the fact.

The leading aluminium design was checked with quadratic elements:

| level | lens mesh, mm | output mesh, mm | selected, nm | straight, nm | gain |
|---|---:|---:|---:|---:|---:|
| quadratic coarse | 0.75 | 2.60 | 22.639 | 16.029 | 1.4124 |
| quadratic refined | 0.60 | 2.00 | 22.713 | 16.159 | 1.4056 |

From coarse to refined quadratic FEM, the selected amplitude changes by
`+0.33%` and the gain by `-0.48%`.  The material ranking is therefore not a
linear-element or mesh artifact.

## Declared physical model

This package fixes the first aluminium application as a thick-body bulk-wave
problem, represented by a 2D plane-strain cross-section.  It does not claim to
represent a finite-thickness Lamb plate.

- isotropic 6061 aluminium: `rho = 2700 kg/m^3`, `cP = 6122.10 m/s`,
  `cS = 3083.81 m/s`;
- `f0 = 242 kHz`, hence `lambdaP = 25.298 mm` and `lambdaS = 12.743 mm`;
- 15 strips, each `7 mm` high, with `1.2 mm` slots; total aperture height
  `121.8 mm`;
- aluminium lens length `28 mm`; aluminium output length `105 mm`;
- target `(x,y) = (88,0) mm`, i.e. `60 mm` after the planar module/output
  interface;
- source proxy: matched longitudinal contact transducer represented by uniform
  normal traction `1 MPa` on all input strip faces;
- internal interfaces: perfect conforming bond; exterior profiled faces are
  traction-free; the source and remote output boundaries use elastic impedance
  terms;
- target measurement: complex longitudinal displacement `ux`.

The intended hardware interpretation is a separate profiled module attached
to the aluminium component at its planar exit, not slots machined into the
inspected region.  A bonded thickness-mode longitudinal PZT spanning the
active input faces is the selected source architecture for this numerical
gate.  A scanning in-plane vibrometer or a longitudinal receiver array can
measure `ux` for validation.

## Five-cycle synthesis audit

`design_5cycle/` contains a fresh reduced search for a `5-cycle Hann @ 242 kHz`
pulse.  The direct-polymer and aluminium CSV geometries are byte-identical to
the earlier four-cycle designs.  Only the nominally matched polymer profile
changes, and that profile was re-solved in full FEM before the material
comparison above.

The full continuum results in this package are carrier-frequency solves.  The
reduced transfer model performs the five-cycle pulse reconstruction.  The
dense broadband full-FEM continuation is documented below.

## Broadband full-FEM stop

The frozen aluminium lens and its matched straight reference were solved at
17 frequencies from `162.2` to `321.7 kHz`.  The spacing is approximately
`15.96 kHz` in the outer band and `7.98 kHz` throughout `B_-6dB`.  The same
fine linear meshes are reused at every frequency.  The solved band contains
`99.815%` of the five-cycle input spectral energy.

Complex amplitude and unwrapped phase were interpolated to the pulse FFT grid
without per-frequency phase adjustment.  The resulting physical focus signals
give:

```text
lens peak                         21.65 nm
matched straight peak             15.88 nm
impulse G_peak                     1.3634
carrier gain                       1.6334
broadening ratio B_t               0.8643
maximum-shift pulse correlation    0.9794
postcursor ratio                   0.0644
```

The pulse-fidelity metrics pass, but the amplitude target `G_peak >= 2` fails.
Doubling the frequency density inside `B_-6dB` changed `G_peak` from `1.36376`
to `1.36343` (`-0.024%`), so the verdict is not an interpolation-grid artifact.

An ideal nondispersive TTD mask with the same 15-strip aperture, pitch,
aluminium propagation and `60 mm` focus gives scalar `G_peak = 3.182`.  The
aperture is therefore capable of the target in principle; the shortfall is in
the complex transfer functions of the real sinusoidal strips.  A new blind
geometry sweep is not justified.  The next allowed diagnostic is an
eight-group physical response matrix followed by a local geometry Jacobian.

## Eight-group response matrix and local Jacobian

All 15 source faces were tagged separately, then combined into the centre
strip plus seven symmetric pairs. Eight independent carrier solves were run
for the lens and eight for the matched straight geometry. Summing the eight
complex focal contributions reconstructs the original full-drive fields to
floating-point precision:

```text
lens       21.044 nm
straight   12.884 nm
gain        1.6334
```

The carrier coherence efficiency `|sum C_i| / sum |C_i|` is only `0.5473`;
the largest focal contribution phase error is `147.5 deg`. If the measured
group amplitudes were retained but their phases aligned, carrier gain would be
`2.984`. Thus the target does not require more aperture or more group
amplitude first: it requires a regulator that can correct group phases.

At a plane `5 mm` after the strips, the response matrix has
off-diagonal-L2/diagonal `0.557` for the lens (`0.536` even for the straight
control). The independent-cell assumption used by the reduced design is not
valid for this fused aperture. Groups at `|y|=41.0` and `49.2 mm` make the
problem especially clear: they have identical two-period, `g=1.5 mm`
geometry, but very different device transmission and focal phase.

One predeclared local physical step tested the most erroneous group:

```text
group                      |y| = 49.2 mm (strips 2 and 14)
gap                         1.50 -> 1.75 mm
focus change                -5.68%
carrier gain                1.633 -> 1.541
group phase error          -147.5 -> -155.0 deg
d gain / d gap              -0.372 per mm
```

The coupled derivative outside the perturbed group is nonzero but small for
this focal observable (`1.09%` of squared derivative norm). More importantly,
the one-sided sign is wrong: increasing the gap worsens both focus and phase,
while the opposite direction would violate the frozen `1.5 mm` minimum and
does not offer enough linear phase range. Gap tuning is stopped after this
single step. The next topology needs an independent broadband phase/delay
regulator rather than a denser sweep of the same constriction.

## Defect pilot

A through cylindrical void of radius `2 mm` was placed at the future focus.
In this 2D model it is an out-of-plane cylindrical free boundary with
`diameter/lambdaS = 0.314`.  The refined quadratic lens and its matched
straight control were solved both with and without the void.

At the would-be defect centre in the defect-free fields, the time-averaged
kinetic-energy density

```text
<Ek> = rho * omega^2 * |u|^2 / 4
```

is `0.8051 J/m^3` for the lens and `0.4075 J/m^3` for the straight control:
an energy-density ratio of `1.9756`.

The complex defect-free field was subtracted from the complex defect field.
On the receiver line `x = 68 mm`, twenty millimetres in front of the void, the
RMS longitudinal scattered displacement is:

```text
lens       5.161 nm
straight   3.583 nm
ratio      1.4406
```

If additive receiver noise is unchanged between the two acquisitions, this
scattered-amplitude ratio is an SNR improvement of `3.17 dB`.  This is an
equal-noise SNR proxy, not a simulated electronics/noise model.  It must not be
interpreted as a material-strength limit or as a calibrated probability of
detection.

## Files

- `design_5cycle/`: reduced pulse synthesis, physical geometry CSVs and plots;
- `material_comparison/`: authoritative physical-unit material and convergence
  tables plus the combined plot;
- `full_fem_linear/`: compact fine-linear FEM fields used for the common-order
  material ranking;
- `convergence/`: actual coarse and refined quadratic aluminium meshes and
  fields;
- `matching_layer/`: full FEM layer sweep and the `1.4 mm` best sampled case;
- `defect/`: actual refined meshes, complex fields, receiver profiles, summary
  and plot.
- `broadband/`: 17-frequency complex FEM responses, both actual meshes,
  reconstructed physical waveforms, metrics and plot.
- `response_matrix/`: baseline and perturbed 8x8 pre-output matrices, focal
  response rows, actual tagged meshes, all basis fields and the local gap
  Jacobian.

## Reproduction

From the repository root:

```bash
METAMATERIALS_MATERIAL_LENS_OUTPUT=tmp/sinusoidal_material_lens_5cycle_242khz \
METAMATERIALS_MATERIAL_LENS_PULSE_CYCLES=5.0 \
julia --startup-file=no --project=. \
  src/code/metamaterial/run_sinusoidal_material_lens.jl

METAMATERIALS_MATERIAL_LENS_FEM_OUTPUT=tmp/sinusoidal_material_lens_fem_q2_refined_242khz \
METAMATERIALS_MATERIAL_LENS_FEM_ORDER=2 \
METAMATERIALS_MATERIAL_LENS_H_LENS_MM=0.60 \
METAMATERIALS_MATERIAL_LENS_H_OUTPUT_MM=2.00 \
julia --startup-file=no --project=. \
  src/code/metamaterial/run_material_lens_fem.jl --solve aluminium:selected

julia --startup-file=no --project=. \
  src/code/metamaterial/analyze_sinusoidal_material_transfer.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/analyze_material_lens_defect.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/run_material_lens_broadband.jl --all
julia --startup-file=no --project=. \
  src/code/metamaterial/run_material_lens_response_matrix.jl --all
julia --startup-file=no --project=. \
  src/code/metamaterial/analyze_material_lens_local_jacobian.jl
```

The defect run additionally sets:

```text
METAMATERIALS_MATERIAL_LENS_DEFECT_RADIUS_MM=2.0
METAMATERIALS_MATERIAL_LENS_DEFECT_X_MM=88.0
METAMATERIALS_MATERIAL_LENS_DEFECT_Y_MM=0.0
METAMATERIALS_MATERIAL_LENS_DEFECT_H_MM=0.35
```
