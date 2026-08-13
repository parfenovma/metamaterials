# Metamaterials

Finite-element modelling and inverse design of solid elastic-wave
metamaterials and focusing meta-lenses.

The current design is a flat-ended all-aluminium true-time-delay lens built
from smooth `input collector -> delay guide -> output diffuser` channels.
The isolated channel preserves `99.65%` of the equal-path envelope peak. The
passive unweighted 15-channel aperture reaches `G=1.870` relative to the
original abrupt reference; a segmented nonnegative source and a 30-cycle Hann
pulse reach the verified system result `G_peak=2.126`, `B_t=1.141`, and pulse
correlation `0.980`.

The next physics question is source--topology matching, not another geometric
sweep. The planned concept uses one electrical waveform and a spatially
apodized piezoelectric layer to excite the measured time-reversal eigenchannel
of the passive TTD geometry. Independent electronic phase channels are outside
the present scope.

The current research roadmap and physical conventions are documented in:

- [`docs/research_plan.md`](docs/research_plan.md)
- [`docs/research_log.md`](docs/research_log.md) — archived experiment history
- [`docs/s_parameter_conventions.md`](docs/s_parameter_conventions.md)

The provisional impulse-risk gates are reproduced with:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_impulse_risk_pilot.jl
julia --startup-file=no --project=. src/code/metamaterial/run_port_band_risk.jl
julia --startup-file=no --project=. src/code/metamaterial/run_dispersive_delay_proxy.jl
julia --startup-file=no --project=. src/code/metamaterial/run_path_packing_bound.jl
```

They fix the five-cycle Hann spectrum, map the modal cut-on versus strip width,
compare the current and single-symmetric-mode ideal apertures, and test
dispersive geometric-delay proxies with delay quantization and Rayleigh loss,
and an optimistic packing bound for a literal in-plane meander.

The subsequent physical three-state gates are staged to keep Gmsh conversion
and quadratic modal FEM isolated:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_distributed_slow_wave_pilot.jl
julia --startup-file=no --project=. src/code/metamaterial/run_tapered_shunt_trim_pilot.jl --stage=mesh
julia --startup-file=no --project=. src/code/metamaterial/run_tapered_shunt_trim_pilot.jl --stage=convert
julia --startup-file=no --project=. src/code/metamaterial/run_tapered_shunt_trim_pilot.jl --stage=fine-lossless
julia --startup-file=no --project=. src/code/metamaterial/run_tapered_shunt_trim_pilot.jl --stage=fine-lossy
julia --startup-file=no --project=. src/code/metamaterial/run_tapered_shunt_trim_pilot.jl --stage=analyze-lossy
julia --startup-file=no --project=. src/code/metamaterial/run_measured_trim_cascade_proxy.jl
julia --startup-file=no --project=. src/code/metamaterial/run_perforated_channel_pilot.jl
```

The distributed slit topology fails the lossless transmission gate. The
tapered shunt state has high weighted transmission but fails pulse broadening;
even its optimistic uncoupled cascade reaches `G_peak = 2.01` only with
`B_t = 2.20` and `rho = 0.829`. Section 45 of the research log records the
no-go verdict; the six-state geometry library remains blocked.
The equal-area circular/elliptical through-hole diagnostic is also a no-go:
all three variants are almost fully reflecting at the carrier.

The current pipeline is split into:

1. `step1_mesher.jl` - parameterized geometry and Gmsh mesh;
2. `step1b_convert_models.jl` - Gmsh-to-Gridap conversion in an isolated process;
3. `step2_solver.jl` - transient elastic FEM and port observables;
4. `step3_analyzer.jl` - complex spectra and reference-calibrated response;
5. `step4_wavelet_analyzer.jl` - Morlet CWT diagnostics for an existing
   sample/reference signal pair. Its automatic frequency grid keeps at least
   half of the record outside the cone of influence and extends up to the
   smaller of `1 MHz` and `0.9` of Nyquist;
6. `plot_wavelet_analysis.jl` - scalograms with frequency on the horizontal
   axis, time on the vertical axis, ridges, COI limits and wavelet-delay plots;
7. `run_period_length_pilot.jl` - resumable 220 kHz pilot that separates the
   effects of the number of periods and corrugated-block length;
8. `run_gap_pilot.jl` - targeted minimum-gap sweep on selected period/length
   anchors; staggered-wall samples are recomputed and reuse only straight
   reference signals;
9. `run_lens_prototype.jl` - seven-element scalar diffraction prototype,
   element selection, focus plots and harmonic-field animation.
10. The same period/length pilot follows the rising `N=1` branch out to
    `L=46 mm`, while the gap pilot refines `L=33--35 mm` around its spectral
    transmission maximum.
11. `run_boundary_condition_study.jl` compares the standard absorbing right
    boundary with a free reflecting end, plots calibrated transmission for the
    seven pilot wall topologies, and generates their actual FEM-mesh gallery.
12. `run_working_band_boundary_study.jl` first confirms the common working band
    of the corrected staggered profiles with independent carrier runs, then
    generates absorbing/free-boundary scalogram pairs for every topology.

The gap-transmission study is summarized in
`docs/transmission_gap_study.md`.

The working-band scan followed by absorbing/free-boundary CWT for every pilot
topology is resumable with:

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_working_band_boundary_study.jl
```

The target-frequency complex library is rebuilt with:

```bash
julia --startup-file=no --project=. src/code/metamaterial/build_target_frequency_library.jl
```

The phase-efficiency envelope and lens-specific Pareto candidates are generated
with `src/code/metamaterial/analyze_phase_efficiency.jl`.

The focused FEM sweep for the missing `150--190 degree` phase sector is run
with `src/code/metamaterial/run_target_phase_search.jl`.

Independent localized constrictions and the coupled-mode pilot use
`src/code/metamaterial/harmonic_solver.jl`,
`src/code/metamaterial/run_bright_constriction_search.jl`, and
`src/code/metamaterial/run_coupled_pair_search.jl`.

The side-mass bright--dark cell and its two-cell neighbour check use
`src/code/metamaterial/side_mass_mesher.jl`,
`src/code/metamaterial/run_side_mass_pilot.jl`, and
`src/code/metamaterial/run_side_mass_neighbor_test.jl`.

The reproducible `R0 / B / BD` spectral ablation and harmonic-field animations
use `src/code/metamaterial/run_side_mass_ablation.jl` and
`src/code/metamaterial/render_side_mass_ablation.jl`.

The conservative `B / BD` eigenmode pilot, parity/localization classification,
clamped/free boundary sensitivity check, and standing-mode animation use
`src/code/metamaterial/modal_solver.jl` and
`src/code/metamaterial/run_side_mass_modal_pilot.jl`:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_side_mass_modal_pilot.jl
```

The fixed-interface component modes, dark-mass continuation, fine-mesh Fano
candidate near 242 kHz, and Rayleigh-loss check use
`src/code/metamaterial/side_mass_component_mesher.jl`,
`src/code/metamaterial/run_fano_dark_mass_tuning.jl`, and
`src/code/metamaterial/run_fano_loss_robustness.jl`. The verified lossless
candidate and its robustness data are stored in `tmp/fano_final_candidate` and
`tmp/fano_final_candidate_loss`; interpretation and limitations are recorded in
sections 3.3 of the research plan and 29 of the research log.

The ideal phase-only Huygens aperture search, used to fix the physical-cell
targets before a new FEM geometry sweep, runs with:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_ideal_huygens_aperture.jl
```

It compares aperture, focal distance and pitch choices and writes the selected
15-strip, 4.8 mm-pitch design to `tmp/ideal_huygens_aperture_242khz`.

The first physical broadband-cell ablation (two directly attached side-mass
pairs) uses `src/code/metamaterial/huygens_pair_mesher.jl` and
`src/code/metamaterial/run_huygens_pair_pilot.jl`. Its negative phase-coverage
result is stored in `tmp/huygens_pair_pilot_242khz` and section 31 of the
research log; it motivates the next balanced series-compliance/shunt-mass
topology.

The follow-up series-compliance/shunt-mass ablation uses
`src/code/metamaterial/balanced_huygens_mesher.jl` and
`src/code/metamaterial/run_balanced_huygens_pilot.jl`. Its component signs are
correct but its phase coverage remains insufficient; data and the reduced-lens
check are in `tmp/balanced_huygens_pilot_242khz` and section 32 of the log.

The 242 kHz lens prototype uses the rebuilt library with:

```bash
METAMATERIALS_LENS_FREQUENCY_HZ=242000 \
METAMATERIALS_LENS_LIBRARY=tmp/target_frequency_library/library_242khz.csv \
METAMATERIALS_LENS_OUTPUT=tmp/lens_prototype_242khz \
julia --startup-file=no --project=. src/code/metamaterial/run_lens_prototype.jl
```

The material-interface experiment for a 15-strip sinusoidal lens is screened
with the project-standard five-cycle pulse and a fixed incident traction using:

```bash
METAMATERIALS_MATERIAL_LENS_PULSE_CYCLES=5.0 \
julia --startup-file=no --project=. src/code/metamaterial/run_sinusoidal_material_lens.jl
```

It compares direct photopolymer-to-6061 aluminium coupling, a geometric-mean
impedance layer, and aluminium elements bonded to aluminium.  Selected and
uniform control apertures can then be meshed and checked with the heterogeneous
full-lens harmonic FEM, one resumable case per process:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_material_lens_fem.jl \
  --solve polymer_matched:selected
```

The reduced model optimizes the time-domain pulse peak; the FEM check is a
carrier-frequency validation and is deliberately reported separately. The P6
material transfer is now closed for the bulk-longitudinal plane-strain pilot:
Al→Al reaches `22.713 nm` on the refined quadratic mesh with matched gain
`1.4056`; a `2 mm`-radius void at the focus gives a `+3.17 dB` equal-noise
scattered-field SNR proxy. Results, actual meshes and limitations are in
[`results/sinusoidal_material_transfer_p6_242khz/`](results/sinusoidal_material_transfer_p6_242khz/README.md).

The subsequent all-aluminium horn--TTD branch is archived in
[`results/aluminium_horn_ttd_p6_242khz/`](results/aluminium_horn_ttd_p6_242khz/README.md).
Its isolated `7.0 -> 3.5 mm` collector plus smooth delay guide preserves
`99.65%` of the equal-path envelope peak with `rho=0.999972`. A deterministic
9-channel aperture reaches equal-pressure/equal-power gains `1.902/1.649`,
but the abrupt 15-channel version stops at `1.883/1.540`. A smooth
`3.5 -> 7.0 mm` output diffuser raises the unweighted equal-power gain over
the original straight reference to `1.870`. Simple passive input/output-area
apodizers fail their complex physical Jacobian gates. With a segmented
nonnegative drive and a 30-cycle Hann pulse, the final five-frequency full-FEM
gate reaches `G_peak=2.126`, `B_t=1.141`, `rho=0.980` and postcursor `0.050`;
the local band contains `99.9485%` of pulse energy. The result is explicitly a
lens-plus-drive system referenced to the original abrupt control, not a purely
passive apodized lens. At the carrier its matched-diffuser gain is `1.9276`.

The one-input passive-feed branch uses absolute multimode ports. A straight
7 mm control gives `T_fund=1.00009`; the calibrated unequal Y gives the
requested ratio `0.41045` and `T_fund=0.96860`, while the 50/50 Y retains
`0.94205` in the useful longitudinal modes. The routed balanced tree fails
with efficiency `0.87669`, power-only pulse ceiling `1.99067`, and carrier
gain `1.80498`. This is not a global no-go: topology-optimal but unbuilt
surrogates give only a marginal `2.006--2.008` power-only ceiling. The tested
feed implementations are stopped; the class remains physically unresolved.
Details are in
[`results/aluminium_horn_ttd_p6_242khz/`](results/aluminium_horn_ttd_p6_242khz/README.md).

The horn-fed point-radiator branch implements research-plan B1 with a
deterministic `3.2 -> 1.6 mm` log-cosine collector, straight equal-path
reference, rounded delay guide and homogeneous receiving block:

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_point_radiator_pilot.jl
```

The wavelength-scale gentle bend passes the isolated five-cycle gate with
peak `0.940`, energy `0.901`, `B_t=1.020` and correlation `1.000`.  Its
measured binary coefficient is assembled into a 15-point aperture and checked
with a fully elastic common-receiver propagator using:

```bash
julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_binary_lens.jl --stage=design
METAMATERIALS_HORN_LENS_ORDER=2 \
METAMATERIALS_HORN_LENS_FINE_MESH=true \
julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_binary_lens.jl --stage=mesh
METAMATERIALS_HORN_LENS_ORDER=2 \
julia --startup-file=no --project=. \
  src/code/metamaterial/run_horn_binary_lens.jl --stage=solve
```

The planar-overlap continuation routes the guide in `x-z` and uses a
curvature-limited `sin^4` path (`38.158 mm` axial, `42.578 mm` unfolded). Its
isolated transient gives peak/equal-path `0.967`, `B_t=0.990`, and transverse
energy `0.0043`. A complete 3D binary triplet library calibrates the
context-aware mask `SDDSSDDDDDSSDDS`. The reconstructed five-cycle aperture
gives peak gain `2.059`; the fine quadratic receiver FEM gives `89.86 nm`,
gain `2.404`, FWHM `4.0 mm`, and sidelobe ratio `0.317` versus the matched
uniform point aperture.

The current winner replaces that binary mask with a monotonic unwrapped
true-time-delay law. The 15 channels have continuously increasing extra path
from `0` at the aperture edges to `9.724 mm` at the centre. The common delay
section is `50.664 mm` long, the largest `sin^4` excursion is `13.746 mm`, and
the minimum inner radius remains `3.93 mm`. The maximum-delay transient gives
peak/equal-path `0.977`, `B_t=0.991`, correlation `0.9998`, and transverse
energy `0.0055`. Reconstructing the full five-cycle aperture from this measured
endpoint gives peak gain `2.904`; the same quadratic receiver FEM gives
`93.17 nm`, gain `2.493`, FWHM `5.0 mm`, and sidelobe ratio `0.294`. This is a
`3.68%` carrier-amplitude increase over the best binary lens while retaining
the amplitude-first trade-off.

Broadband neighbour dependence and one common-array transient remain open.
The full assumptions and limitations are in
[`docs/horn_fed_point_radiator_lens.md`](docs/horn_fed_point_radiator_lens.md).

The first integrated 3D monotonic array is also complete. It fuses all 15
physical horn/guide channels to one receiving block (`59,751` nodes,
`332,095` elements) and is driven only through the 15 common-phase input
faces. The end-to-end `242 kHz` solve produces a central beam at `x=35 mm`
with `|ux|=3.893 nm`, FWHM `5 mm`, and center/outer amplitude `5.28`; actual
output phases remain within `12.7 deg` of the reduced design. Reproduce it
with `run_horn_monotonic_array_3d.jl`; outputs are in
`tmp/horn_monotonic_array_3d`. A matched full-3D straight-array reference gives
`2.645 nm` at the same target versus `3.893 nm` for the lens: the first honest
end-to-end carrier gain is therefore `1.472` (`+47.2%` amplitude).

`METAMATERIALS_LENS_ELEMENT_COUNT`, `METAMATERIALS_LENS_FOCAL_DISTANCE_MM`, and
`METAMATERIALS_LENS_MINIMUM_AMPLITUDE` optionally control the aperture and
candidate filter.

The resumable frequency sweep runs all four stages for seven representative
profiles at `160:20:280 kHz` and writes
`tmp/pilot_220khz/library_frequency_sweep.csv`:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_pilot_sweep.jl
```

Use `--list` to inspect its profiles without starting FEM. Re-running the
pipeline keeps completed meshes, models and signal files. FEM profiles run in
four isolated Julia processes by default; set `METAMATERIALS_PILOT_JOBS` to
change that number.

The period/length sensitivity pilot uses a sinusoidal profile at `220 kHz`,
screens `N=1:4` at `17 mm`, refines the length grid for `N=1,2`, and writes its
table and plots to `tmp/period_length_pilot_220khz`. It can also be resumed:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_period_length_pilot.jl
```

Set `METAMATERIALS_GEOMETRY_JOBS` to change its number of parallel FEM
processes.

The targeted gap pilot uses the staggered-wall convention and runs every sample
geometry; only the matching straight reference signals are reused:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_gap_pilot.jl
```

Build the reduced seven-element lens prediction and animation from the gap
library:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_lens_prototype.jl
```

Use the pinned Julia environment for every command:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. src/code/metamaterial/step1_mesher.jl
julia --project=. src/code/metamaterial/step1b_convert_models.jl
julia --project=. src/code/metamaterial/step2_solver.jl
```
