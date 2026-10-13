# O8 PZT coupon v1 — BUILD

This package authorizes fabrication of **one stationary four-state de-risking
coupon**.  It is not a production flaw detector and it does not authorize a
POD claim.  Its single purpose is to determine whether two passive RL-shunted
PZT pairs create four repeatable elastic measurement operators in aluminium.

## Frozen article

- Al 6061-T6 block: `50.60 x 202.38 x 8.00 mm` (`x, y, z`);
- four identical CTS 3195HD plates: `20.00 x 8.00 x 1.00 mm`;
- inner pair A at `y = +/-20.07 mm`, wired with equal polarity;
- outer pair B at `y = +/-91.19 mm`, one plate electrically reversed;
- state A target: `L_total=70.15 uH`, `R_total=17.89 ohm`;
- state B target: `L_total=73.27 uH`, `R_total=27.57 ohm`;
- open/shunted relay states: `q00`, `q10`, `q01`, `q11`.

The total load values include inductor ESR, relay contact resistance, wiring
and the external resistor.  They are recalculated from the **measured** paired
capacitance before the final resistor is fitted.

## Frozen numerical manufacturing gate

The final rounded geometry was rerun with the CTS 3195HD public material
values `epsilon_r=1900`, `kt=0.49`, `Qm=80` over `193.6...290.4 kHz`.

| Metric | Result | Gate |
|---|---:|---:|
| minimum output/reference amplitude | 0.7367 | >= 0.70 |
| maximum output/reference amplitude | 0.8054 | <= 1.50 |
| minimum four-state distance at 242 kHz | 0.0769 | >= 0.05 |
| maximum raw mirror error | 0.01083 | <= 0.02 |
| maximum reciprocity error | 6.98e-16 | <= 1e-6 |

**Verdict: numerical manufacturing gate passed.**  No geometry tuning is
allowed after this freeze.

## Package contents

- `coupon_block.step` — machinable aluminium solid;
- `coupon_assembly.step` — block plus four PZT solids;
- `dimensioned_drawing.svg` — dimensions and placement;
- `pzt_placement.csv` — exact patch coordinates and parity;
- `bom.csv` — procurement and build quantities;
- `pzt_rfq.md` — text for the ceramic vendor;
- `wiring.svg` and `shunt_targets.csv` — two passive switching channels;
- `assembly_and_tuning.md` — controlled build sequence;
- `bench_acceptance.md` — preregistered build/kill experiment;
- `source_points.csv` and `receiver_scan_points.csv` — bench coordinates;
- `coupon_validation.msh`, `model_validation_band.csv`, and
  `model_validation.jld2` — reproducible model evidence.

Regenerate the CAD and numerical evidence with:

```bash
julia --startup-file=no --project=. src/code/metamaterial/export_o8_pzt_coupon.jl
```

## Terminal decision after the bench test

- **BUILD** if every electrical and acoustic gate in `bench_acceptance.md`
  passes.  The next hardware article must preserve this geometry and parity.
- **KILL** if the complete test still fails after the one permitted
  capacitance/parasitic retune.  Do not rescue the coupon by moving patches or
  lowering the state-distance gate.

## Vendor basis

CTS lists 3195HD at `epsilon_r=1900`, `Qm=80`, and `kt=0.49`, and offers
custom bulk plates and electrode configurations:

- <https://www.ctscorp.com/Products/Piezoelectric/Bulk-Piezoceramics/Ceramic-Material-Properties>
- <https://www.ctscorp.com/Files/DataSheets/Piezoelectric/Bulk-Piezoceramics/CTS-Piezoelectric-Bulk-Piezoceramics-Plates-Datasheet.pdf>

The selected bond adhesive is EPO-TEK 301-2.  Its current data sheet gives
low viscosity and room-temperature or `80 C / 3 h` cure options:

- <https://www.epotek.com/docs/en/Datasheet/301-2.pdf>

