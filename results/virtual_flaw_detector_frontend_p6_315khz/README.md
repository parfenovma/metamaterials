# Virtual flaw-detector receive chain — 315 kHz

This branch converts the completed reflective-backwall pulse-echo FEM traces into an instrument-like A-scan and a stochastic detection experiment. It does not run new FEM.

## Signal chain

1. Existing reflective-backwall PZT stress traces at 315 kHz, five-cycle burst.
2. Open-circuit voltage surrogate `V = g33 * t_pzt * coupling * stress`, with `g33 = 24.8e-3 V m/N`, `t_pzt = 6 mm`, and an explicit uncalibrated receive-coupling factor of `0.02`.
3. Smooth 120–600 kHz FFT band-pass.
4. Twenty-two randomly located weak reflectors, normalized by their Gate-A RMS, plus band-limited electronic noise at 2% of the clean 2-mm echo peak.
5. A 1 V full-scale ADC. The sample screen uses 21.20 dB receive gain, chosen so the active-aperture Gate-A indication is about 80% FSH; the front pulse and backwall can clip.
6. Detection statistic: peak display envelope in Gate A, 90–135 mm.
7. Threshold: the 99th percentile of 600 no-defect realizations, giving empirical `Pfa = 1%`.

The same absolute post-receive clutter/noise realization is added to uniform and active traces. This isolates the irreducible receiver/output-clutter ceiling, but it intentionally does **not** credit the active aperture with spatial suppression of individual microstructural scatterers.

## Detection result

| Clutter RMS / clean echo peak | Uniform Pd | Active Pd |
|---:|---:|---:|
| 0% | 100% | 100% |
| 10% | 100% | 100% |
| 20% | 97.2% | 95.0% |
| 35% | 40.0% | 40.0% |
| 50% | 13.7% | 11.2% |
| 100% | 0.8% | 0.8% |

At `Pfa = 1%`, the current 2-mm reflector crosses a provisional `Pd >= 90%` acceptance criterion only while Gate-A structural clutter stays at roughly 20% of the clean echo peak or below. At 35%, increasing display gain cannot recover detection because the threshold rises with the clutter.

The active controller does not improve this particular ROC. That is an important result: after clutter has become spatially indistinguishable at the receiver output, aperture optimization has no remaining degree of freedom to remove it. A future microstructure-resolved FEM or measured multi-channel dataset is required to test whether spatial filtering changes that conclusion.

## Economic go/no-go interpretation

Before machining another passive topology, measure a real coupon with the intended probe and a known 2-mm reference reflector:

- if `RMS_clutter / peak_echo <= 0.20`, this model predicts a credible route to `Pd >= 90%` at `Pfa = 1%`, so a physical controller/lens comparison is justified;
- if the ratio is near `0.35` or above, the present single-gate detector is not economically rescued by more lens optimization alone; the next investment should be spatially resolved acquisition, coded excitation, SAFT/TFM, or a different frequency;
- the cutoff is a simulation-derived experiment-design criterion, not a certified POD claim.

## Reproduction and artifacts

Run:

```bash
julia --startup-file=no results/virtual_flaw_detector_frontend_p6_315khz/frontend_pipeline.jl
julia --startup-file=no results/virtual_flaw_detector_frontend_p6_315khz/build_visualization.jl
```

Artifacts:

- `data/roc_by_clutter.csv` — threshold, empirical Pfa, and Pd versus clutter;
- `data/pd_by_echo_scale.csv` — Pd versus an echo-amplitude surrogate;
- `data/frontend_metrics.csv` — electrical scaling and default-screen metrics;
- `data/sample_adc_traces.csv` — clean voltage traces and one reproducible noise/clutter realization;
- `data/frontend_ensemble.jld2` — compact machine-readable bundle;
- `frontend_pipeline.jl` — reproducible receive-chain and ensemble calculation;
- `build_visualization.jl` — generates the interactive A-scan fragment.

## Limitations

- The PZT conversion omits electrical impedance, cable, backing, matching layers, pulser waveform, and a calibrated transimpedance/load model.
- The reflector is a 2-D circular void; echo-amplitude scaling is not equivalent to changing a real flaw radius.
- Random scatterers are a phenomenological clutter ensemble rather than grains or inclusions explicitly meshed in FEM.
- The display envelope is a causal rectified-decay surrogate; ROC scores use the analytic envelope in Julia.
- No experimental probability-of-detection claim should be made from these numbers.
