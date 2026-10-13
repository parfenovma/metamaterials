# O6h-C exact 2D angular contact spot check

**FAIL: the normal-incidence P12 contact qualification does not survive 2D angular P/S physics.**

The exact solid--fluid--solid boundary solve conserves lossless power to `3.997e-15` and reproduces the 1D normal-incidence transmitted power to `2.776e-15`.

| Scenario | Min layer/direct power, 0--45 deg | Max Delta S fraction | 55-deg min power | Pass |
|---|---:|---:|---:|---|
| `minus_sigma` | `0.7245` | `0.2327` | `0.6661` | `false` |
| `nominal` | `0.6541` | `0.2327` | `0.5977` | `false` |
| `plus_sigma` | `0.5970` | `0.2327` | `0.5432` | `false` |

This is an exact 2D plane-wave spot check over the operator band and ROI-relevant angles. It still omits lateral thickness texture, squeeze-out edges and a finite contact footprint; a failure is therefore a physical blocker, while a pass is only permission for a finite-footprint prototype model.
