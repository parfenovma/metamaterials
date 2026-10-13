# Quantized round-hole Fresnel/zone plate at 305--321.7 kHz

## Bottom line

This branch found a real but modest **fine-mesh joint gain**, while failing the
campaign's focal-position gate.  The selected `f006` design at 321.7 kHz gives

| metric | f006 fine | straight reference |
|---|---:|---:|
| J1 nominal-focus gain | **1.0307** | 1.0000 |
| J2 equal-noise reduced CNR | **1.1512** (+1.223 dB) | 1.0000 |
| J2 clutter-dominated | **1.2167** | 1.0000 |
| exit-energy proxy | **0.8540** | 1.0000 |
| FWHM at x=115 mm | **29.52 mm** | -- |
| axial/lateral focus error | **9.0 mm** (`x=106`, `y=0`) | -- |
| nominal-plane side-lobe ratio | **0.0518** | -- |

The design passes the energy part of both manufacturing gates, but it passes
neither complete gate because `focus_error_mm > 5 mm`:

- prototype gate (`exit>=0.70`, error `<=5 mm`): **fail**;
- strict gate (`exit>=0.75`, error `<=3 mm`): **fail**.

Accordingly, `f006` is a useful topology direction, not a fabrication-ready
optimum.  No raw coarse gain in this folder is presented as a feasible result.

## Topology and manufacturing complexity

The aperture is split into nine mirror-symmetric transverse manufacturing
cells.  The continuous geometric-delay profile is quantized into 3--5 fill
levels and realised only with round through-drilled holes.  Three or four
longitudinal drill rows distribute each cell's void area.  A separate
outer-aperture porosity term provides an apodization DOF that can alter clutter
and side lobes without directly changing the central delay level.

The fine-selected `f006` has:

- 4 nominal phase/fill levels;
- 4 longitudinal rows and 36 holes;
- 5 actual drill diameters after the independent outer apodization is applied;
- a measured minimum hole/web feature of 2.613 mm, well above the 0.8 mm gate;
- a uniform equal-energy traction source and a same-frequency solid reference.

All variants are monolithic aluminium with through-cut circular voids.  No
resonant inserts, bonded layers or piezoelectric elements are assumed.

## DOE and Pareto result

The completed DOE contains exactly 24 coarse harmonic FEM cases: eight
geometries at 305, 315 and 321.7 kHz.  The axes are phase-level count, drill-row
count, porosity scale, independent outer apodization, delay-profile exponent
and two focal-depth retargets.

The strongest raw coarse points were:

| variant/frequency | J1 | J2 | exit | error | interpretation |
|---|---:|---:|---:|---:|---|
| f003 / 305 kHz | 1.2066 | 1.4661 | 0.7789 | 20.62 mm | raw Pareto, infeasible focus |
| f006 / 321.7 kHz | 1.1784 | 1.4755 | 0.8274 | 8.00 mm | selected for fine check |
| f008 / 321.7 kHz | 1.1040 | 1.2699 | 0.8643 | 9.00 mm | higher energy, lower gain |

None of the 24 points passes the prototype or strict gate.  Retargeting the
analytical design depth from 60 to 70 and 80 mm did not move the FEM peak to
115 mm.  The axial peak is therefore controlled mainly by the discrete elastic
scattering/interference pattern, not by the simple effective-delay initializer.

The `f006@321.7` coarse-to-fine change is material: J1 falls from 1.1784 to
1.0307 and J2 from 1.4755 to 1.1512, while exit energy rises from 0.8274 to
0.8540 and the focus error changes from 8 to 9 mm.  This is why only the fine
numbers above should be used outside this screening study.

## Honest comparison with ellipse c005

| fine result | J1 | J2 | exit | FWHM | error | side lobe | manufacture |
|---|---:|---:|---:|---:|---:|---:|---|
| Fresnel f006, 321.7 kHz | 1.0307 | 1.1512 | **0.8540** | **29.52 mm** | **9 mm** | **0.0518** | 36 circles / 5 diameters |
| Ellipse c005, 315 kHz | **1.2057** | **1.3859** | 0.7090 | 32.52 mm | 11 mm | 0.0821 | 27 ellipses / 10 shapes |

The round-hole plate is easier to specify, has a larger minimum feature, passes
more energy, narrows the nominal-plane lobe and lowers the same side-lobe
diagnostic.  It does **not** beat c005 in either primary objective, and neither
fine design passes the focal-position gate.  The comparison therefore supports
keeping this Fresnel family as the high-transmission/simple-manufacture branch,
not replacing c005 as the raw-gain baseline.

## Artifacts

- `data/campaign_metrics.csv`: common campaign columns first, followed by
  Fresnel DOFs, manufacturing complexity, side lobes and explicit gates.
- `data/doe_manifest.csv` and `data/f00*_zone_schedule.csv`: complete drill
  schedules and actual diameters.
- `data/fine_comparison_with_ellipse_c005.csv`: common fine comparison.
- `figures/doe_topologies.png`: all eight drill patterns.
- `figures/frequency_variant_heatmap.png`: J1/J2 frequency heatmaps.
- `figures/J1_J2_pareto.png`: raw two-objective Pareto plane.
- `figures/gate_diagnostics.png`: transmission versus focal error.
- `figures/best_field_f006_0321p7_fine.png`: fine field and transverse profile.
- `figures/mesh_f006_fine.png`: full and local conformal mesh views.

## Reproduction

```sh
julia --startup-file=no fresnel_zone_pipeline.jl --manifest
julia --startup-file=no fresnel_zone_pipeline.jl --run-coarse
julia --startup-file=no fresnel_zone_pipeline.jl --case f006 321.7 fine
julia --startup-file=no fresnel_zone_pipeline.jl --report
julia --startup-file=no fresnel_zone_pipeline.jl --plot-field f006 321.7 fine
julia --startup-file=no fresnel_zone_pipeline.jl --plot-mesh f006 fine
```

These are single-frequency, linear plane-strain carrier calculations.  They do
not predict transient ringing, group delay, nonlinear harmonics or real defect
echo CNR; those require the later multifrequency/transient campaign.
