# O6f synthesis and O6g handoff

O6f measured the missing physical-risk coordinate on both nondominated O6e
branches. It did not introduce a scalar cost or reoptimize either operator.

## Evidence obtained

The q8 branch now has a concrete reciprocal factorization: three input ports,
three near-unit singular attenuators and a 45-coupler adjacent-Givens output
lattice. The zero-parasitic reconstruction error is `7.0e-16`. Its optimistic
component ceilings are `0.047376 dB` uniform insertion loss and `15.252 ns`
uniform delay per output coupler. With independent errors on all 45 couplers,
the largest common 99%-Wilson scale is `4.0732`, corresponding to
`0.4073 deg` angle RMS and `8.1465 ns` delay RMS.

These values are distributed-network bounds. They are not isolated-component
datasheet requirements and do not include the two input transforms, dump-port
loss, winding resistance, dielectric loss or piezo impedance.

The motion branch was replayed with the saved P12/P17 layouts, templates,
training nuisance precision and thresholds. The replay exactly reproduced
the original signal, localization and H0 counts. Under paired encoder and
contact errors, both layouts pass every preregistered scale through `2.0`:

| Layout | Pfa upper 95% | alpha90 / nominal | Pdl at reference |
|---|---:|---:|---:|
| P12 | `0.00627` | `1.0611` | `0.8777` |
| P17 | `0.00927` | `1.0834` | `0.9140` |

Scale `2.0` means `0.10 lambda` position RMS, `10%` contact log-amplitude RMS
and `0.04` carrier-cycle contact-delay RMS. The nearest gate is relative
alpha90 for both layouts. P17 retains better absolute detection but has the
smaller relative robustness margin.

## Decision

Both branches remain physically credible and nondominated. O6f therefore
does not authorize a prototype selection. It closes the preregistered audit
and exposes the next calculable boundaries:

1. allocate a joint q8 loss/error budget across all physical blocks, rather
   than treating the output lattice alone;
2. locate and attribute the motion failure boundary above scale `2`, rather
   than reporting an unbounded pass.

The next stage is O6g. After it, further reduced-surrogate sweeps are stopped;
the surviving risk must be tested at a higher physical evidence level.

## Reproduction

```text
julia --startup-file=no --project=. src/code/metamaterial/run_operator_o6f_beamformer_network.jl
julia --startup-file=no --project=. src/code/metamaterial/run_operator_o6f_motion_risk.jl
```

Results:

- `results/operator_space/o6f_beamformer_network_v0/`
- `results/operator_space/o6f_motion_risk_v0/`
