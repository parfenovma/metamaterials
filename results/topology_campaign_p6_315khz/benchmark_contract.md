# Common benchmark contract: topology campaign at 300--322 kHz

Every candidate must be compared with the solid straight footprint at the same
frequency, mesh resolution and source-energy normalization.  Results that do
not satisfy this contract are diagnostic only and must not enter the common
Pareto set.

## Geometry and manufacturing gates

- Aluminium constants, aperture, device footprint and output domain are those
  in `ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl`.
- Nominal target: `(x,y)=(115,0) mm`.
- Minimum void dimension and minimum solid web: `0.8 mm`.
- The candidate must be realizable as a through-cut 2D pattern unless a branch
  is explicitly labelled an active upper-bound benchmark.

## Primary metrics

At each frequency `f`, with the matching straight solution denoted by `ref`:

- `J1 = |ux_lens(115,0;f)| / |ux_ref(115,0;f)|`.
- `J2_equal`: reduced reciprocal Born CNR gain from the spatial kernel
  `K=ux^2`, using `noise_variance = clutter_variance_ref`.
- `J2_clutter_dominated`: the same gain with
  `noise_variance = 0.25*clutter_variance_ref`.
- `exit_energy_proxy = sum(|ux_lens|^2) / sum(|ux_ref|^2)` on the common exit
  line.
- `focus_error_mm`: distance from `(115,0) mm` to the strongest point in the
  fixed local search window.
- `fwhm_mm`: lateral FWHM at the nominal target plane `x=115 mm`.

## Feasibility and ranking

Two gates are reported, not silently imposed:

- prototype gate: `exit_energy_proxy >= 0.70`, `focus_error_mm <= 5 mm`;
- strict gate: `exit_energy_proxy >= 0.75`, `focus_error_mm <= 3 mm`.

The principal joint score is `min(J1,J2_equal)`.  Pareto dominance is always
computed on the raw pair `(J1,J2_equal)`; the scalar score is used only to rank
otherwise nondominated designs.  Side-lobe, bandwidth and source dynamic-range
metrics remain explicit constraints/diagnostics, not hidden weights.

## Frequency and convergence

- Screening must include 305, 315 and 321.7 kHz or explain why a different
  bracket is necessary.
- Any claimed joint improvement needs a fine-mesh check against a fine straight
  reference at the same frequency.
- A single-frequency harmonic solution is a carrier response, not a transient
  echo prediction.  Group delay, ringing and pulse CNR require a later
  multifrequency/transient validation.

## Required campaign CSV columns

`topology,variant,resolution,frequency_khz,J1,J2_equal,J2_equal_db,J2_clutter_dominated,exit_energy_proxy,fwhm_mm,focus_error_mm,min_feature_mm,source_kind,manufacturing_class,artifact_path`

