# O8 handoff: plane-strain stationary PZT boundary

> Completed 15 August 2026.  Preregistered contract:
> [`operator_o8_plane_strain_boundary_contract.md`](operator_o8_plane_strain_boundary_contract.md).
> Numerical bundle:
> [`results/operator_space/o8_plane_strain_v0/`](../results/operator_space/o8_plane_strain_v0/README.md).

## Verdict

The frozen two-bit PZT branch passes O8.  The frozen mechanical-stiffness
branch fails before detector replay because its four physical operators are
not sufficiently distinct.

This is a conditional numerical pass, not a prototype recommendation.  The
PZT detector margin is extremely small at the active detection gate:
`Pd=0.900245` against `0.90`.

## Physical plane-strain result

The model is a `202.38 mm x 50.60 mm` aluminium plane-strain block, obtained
directly from the frozen `8 lambdaP x 2 lambdaP` O7 normalization at `242 kHz`.
It has 17 output bins and three fixed DCT traction ports.  The base meshes have
about 3,000 nodes; the `lambdaS/8` checks have about 5,300 nodes.

The frozen section is shown in
[`pzt_head_schematic.png`](../results/operator_space/o8_plane_strain_v0/pzt_head_schematic.png).

| Branch | amplitude/reference | centre min state distance | max raw mirror error | max reciprocity error | max refinement projector error | Physical gate |
|---|---:|---:|---:|---:|---:|---|
| PZT/RL | 0.721–0.779 | 0.0766 | 0.01135 | 7.47e-16 | 0.0498 | pass |
| mechanical stiffness | 0.978–0.992 | 0.0104 | 0.01116 | 5.71e-16 | 0.0575 | **fail** |

The unstructured mesh reflection error falls on refinement and is removed by
the fixed `C_s` group projector declared in the contract addendum.  The
projected error is zero to working precision.  No candidate parameter depends
on this projection.

Across the full band, the PZT pairwise projector distance spans
`0.0479...0.2152`; the smallest high-edge pair is below the centre-only `0.05`
gate but remains part of a detector-passing broadband operator.  The
mechanical range is only `0.0080...0.0217`: preserving almost all amplitude is
not useful if switching barely changes the measurement subspace.

The mechanical O7 pass was therefore an artefact of giving its reduced
stiffness perturbations more modal leverage than the explicit local Robin
boundary provides.  That frozen mechanism is stopped.  Revisiting it would
require a new resource allocation or a new latch topology, not parameter
tuning inside O8.

## Fresh detector replay

Only the physically passing PZT branch entered a new O8 split.  The O7 final
test was not reused.

| Architecture | channels | states | events | Pfa95 | alpha90 | Pd@FMC ref | Pdl@FMC ref | Final gate |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| same FEM head, q00 only | 3 | 1 | 3 | 0.229% | 421.65 | 0.810 | 0.778 | fail |
| FEM PZT two-bit | 3 | 4 | 12 | 0.129% | 316.14 | 0.900 | 0.884 | **pass** |
| ideal K3Q4 | 3 | 4 | 12 | 0.155% | 277.94 | 0.940 | 0.937 | pass |
| FMC | 17 | 1 | 17 | 0.075% | 219.19 | 0.968 | 0.966 | pass |

Four physical states improve `alpha90` by `25.02%` relative to the identical
head using frequency samples in `q00` only.  The FEM candidate remains
`13.74%` worse than ideal K3Q4 and `44.23%` worse than FMC in alpha90.

The uniform A-scan had `Pfa95=6.80%` on this fresh background shift and is not
a valid quantitative gap anchor.

## What is now physically specified

- one stationary footprint;
- two even and one odd electronic source/receive sectors;
- inner equal-polarity PZT pair controlling the even sector;
- outer opposite-polarity PZT pair controlling the odd sector;
- four `1 mm` PZT boundary patches, `2 nF` each;
- two series-RL shunts shared within mirrored pairs;
- `94.41 uH + 24.07 ohm` and `98.61 uH + 37.10 ohm`;
- open/shunted states `q00`, `q10`, `q01`, `q11`;
- no negative capacitance, cell-wise control or probe motion.

## Remaining evidence gap

O8 still imposes the three DCT traction ports rather than realizing their
transducers.  Each PZT patch is a finite boundary interval, but its electric
field is condensed to a one-port constitutive Schur complement.  Patch
geometry through thickness, actual electrode area, adhesive shear lag,
dielectric loss, BVD motional branches, switch parasitics, terminal wiring,
thermal drift and the specimen contact are absent.

The fixed `2 nF` capacitance is not yet reconciled with the actual out-of-plane
width of a domain-resolved patch.  The `202 mm` aperture is also a real size
and cost issue, not merely a plot scale.

## Next authorized stage: O9

Only the PZT branch advances.  O9 must preserve the O8 geometry, parity wiring
and loads and add:

1. explicit PZT and adhesive domains in plane strain;
2. consistent electrode area/capacitance and piezoelectric tensor;
3. series-BVD dielectric/motional loss plus switch `R/L/C` parasitics;
4. finite aluminium-head/specimen contact;
5. a preregistered component/adhesive tolerance ensemble;
6. one fresh detector replay after the physical model is frozen.

Because the O8 `Pd` margin is only `2.45e-4`, O9 failure is plausible.  O8 by
itself did not authorize fabrication, a POD claim, or geometry optimization on
the O8 test split.

## Manufacturing addendum

The later material-realism freeze `O8M` removed one narrow part of that stop.
The frozen geometry was rerun with public CTS 3195HD values and then with four
identical `20 x 8 x 1 mm` plates.  The rounded article retained a carrier
minimum state distance of `0.07689` and an amplitude/reference range of
`0.7367...0.8054`, so fabrication of **one mechanism-only coupon** is now
authorized.  See
[`manufacturing/o8_pzt_coupon_v1/`](../manufacturing/o8_pzt_coupon_v1/README.md).
Product fabrication and a POD claim remain unauthorized until the missing O9
physics and detector evidence are supplied.
