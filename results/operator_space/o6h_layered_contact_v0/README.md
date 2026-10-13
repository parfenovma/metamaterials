# O6h-B exact layered-contact motion replay

The frozen P12/P17 detector is driven by an exact normal-incidence thin-layer transfer matrix with all layer multiple reflections. The nominal 61.16 um couplant response is calibrated once; exact per-stop, per-frequency residual transfer multiplies signal and nuisance without an affine approximation.

| Layout | Max passing scale | Thickness RMS | Z log-RMS | c log-RMS | q95 non-affine residual | First failed gate |
|---|---:|---:|---:|---:|---:|---|
| `motion_compact_center_p12_native` | `1.24463` | `38.06 um` | `12.45%` | `6.22%` | `6.184%` | `pfa` |
| `motion_uniform_full_p17_a8` | `0.98279` | `30.05 um` | `9.83%` | `4.91%` | `4.739%` | `pfa` |

The fitted equivalent log-amplitude and delay RMS are reported in the threshold CSV to compare this physical layer with O6g. The calculation is exact for normal-incidence longitudinal plane waves, but omits oblique P/S conversion, lateral squeeze-out and finite-footprint lift-off; those define the justified 2D spot check.
