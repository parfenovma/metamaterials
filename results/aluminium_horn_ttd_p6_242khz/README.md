# All-aluminium horn--TTD P6 package

This package contains the physical-unit meshes, fields, signals, response
matrices, tables and plots for the all-aluminium continuation around `242 kHz`.
The final topology is a `7.0 -> 3.5 mm` log-cosine collector, a smooth
sin-fourth delay guide, and a `3.5 -> 7.0 mm` output diffuser into an aluminium
receiver.

The isolated pilot in `isolated_pilot/` is a PASS. Its refined
`130 us / 40 samples-per-period` maximum-delay transient has envelope
peak/equal-path `0.9965`, energy/equal-path `0.9964`, `B_t=1.000`,
`rho=0.999972`, postcursor `0.09718`, and transverse energy ratio `0.00523`.
Measured delay calibration is `0.20862 us/mm` (`c_eff=4.79 km/s`).

The coupled aperture progression is preserved in `five_channel/`,
`nine_channel/`, and `fifteen_channel/`. The nine-channel carrier gate passes
with pressure/equal-power gains `1.902/1.649`. The 15-channel result is
`2.303 nm` versus `1.223 nm`, pressure gain `1.883`, equal-power gain `1.540`,
and FWHM `8 mm`. Its response-matrix phase-aligned ceiling is `1.988`; optimal
nonnegative equal-power weights give `1.834`. Therefore the abrupt-radiator
topology is stopped for the `2x` amplitude goal.

The continuation is stored in `output_diffuser_pilot/` and
`fifteen_channel_diffuser/`. The isolated output diffuser improves target
amplitude over the abrupt radiator by `1.2533x` at equal source work while
retaining `B_t=1.000`, `rho=0.999952`, postcursor `0.0924`, and transverse
energy `0.00618`. The unweighted full diffuser array reaches equal-power gain
`1.870` over the original abrupt straight reference.

The carrier response/power matrix gives a nonnegative equal-power ceiling
`2.0524`, reproduced by a direct 3D solve. Physical area mappings were tested
without sweeps in `mouth_jacobian/` and `output_jacobian/`; both fail because
the extreme transitions add `44--120 deg` phase shifts and do not reproduce
the requested amplitude contrast. They must not be presented as passive
solutions.

For the allowed lens-plus-segmented-drive system, five-cycle broadband weights
reach only `1.972`. A deterministic pulse-duration check finds the first
`2x` pass at a 30-cycle Hann pulse. The final five-frequency full-FEM
reconstruction over `226.04--257.94 kHz` contains `99.9485%` of pulse energy
and gives:

```text
G_peak                 2.1261
B_t                    1.1408
pulse correlation      0.9805
postcursor              0.0497
carrier equal-power     2.0503
carrier FWHM           16 mm
profile peak            y = 0 mm
```

Thus the primary amplitude target is achieved for the combined
double-horn TTD lens, segmented real nonnegative drive, and 30-cycle pulse.
The passive unweighted geometry remains at `1.870`.

The follow-up `passive_feed/` package tests whether the same pulse30 aperture
can be driven from one transducer. Absolute three-mode DtN ports replace the
invalid bulk-P power proxy. The best calibrated unequal Y has amplitude ratio
`0.41045` (target `0.41092`) and useful-mode transmission `0.96860`; the
required 50/50 Y retains `0.94205`. An optimistic full-tree cascade therefore
retains at most `0.87669` useful power, below the `0.88493` requirement. Its
power-only pulse ceiling is `1.99067` and routed carrier response-matrix gain
is `1.80498`. The passive one-input branch is a carrier-gate STOP; no
full-tree broadband solve was run.

Key figures:

- `isolated_pilot/aluminium_horn_ttd_geometry.png`
- `isolated_pilot/aluminium_horn_ttd_temporal_check.png`
- `nine_channel/nine_aperture_focus.png`
- `fifteen_channel/15_aperture_focus.png`
- `output_diffuser_pilot/diffuser_pilot_transient.png`
- `fifteen_channel_diffuser/15_weighted_pulse30_direct_focus.png`
- `fifteen_channel_diffuser/15_pulse30_local_band.png`
- `passive_feed/verdict/passive_feed_final_verdict.png`

Reproduction entry points:

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_aluminium_horn_ttd_pilot.jl --stage=mesh --case=maximum
julia --startup-file=no --project=. \
  src/code/metamaterial/run_aluminium_horn_ttd_pilot.jl --stage=solve --case=maximum
julia --startup-file=no --project=. \
  src/code/metamaterial/analyze_aluminium_horn_ttd_temporal_check.jl

julia --startup-file=no --project=. \
  src/code/metamaterial/run_aluminium_horn_ttd_small_aperture.jl --stage=compare
julia --startup-file=no --project=. \
  src/code/metamaterial/run_aluminium_horn_ttd_nine_aperture.jl \
  --stage=compare --channels=9
julia --startup-file=no --project=. \
  src/code/metamaterial/run_aluminium_horn_ttd_nine_aperture.jl \
  --stage=matrix_compare --channels=15
julia --startup-file=no --project=. \
  src/code/metamaterial/analyze_aluminium_horn_ttd_power_weights.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/optimize_aluminium_horn_ttd_broadband_weights.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/sweep_aluminium_horn_ttd_pulse_cycles.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/analyze_aluminium_horn_ttd_pulse30_local_band.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/design_aluminium_horn_ttd_passive_feed.jl
julia --startup-file=no --project=. \
  src/code/metamaterial/run_aluminium_horn_ttd_passive_splitter_modal.jl \
  --variant=calibrated_beat --frequency-hz=242000
julia --startup-file=no --project=. \
  src/code/metamaterial/analyze_aluminium_horn_ttd_passive_feed_verdict.jl
```
