# O8 frozen plane-strain boundary realization

Evidence: plane-strain elastic FEM with finite control zones. PZT electrodes/load are represented by an exact one-port constitutive Schur complement; the patch electric field is not domain-meshed.

| Branch | min amplitude | max amplitude | carrier min state distance | max raw mirror error | max reciprocity error | Physical gate |
|---|---:|---:|---:|---:|---:|---|
| pzt | 0.7212 | 0.7791 | 0.0766 | 0.01135 | 7.465835747902474e-16 | true |
| mechanical | 0.9778 | 0.9921 | 0.0104 | 0.01116 | 5.71094215202699e-16 | false |

- `pzt`: failed gates = none.
- `mechanical`: failed gates = carrier_state_distance.

Only a physically passing branch may enter the fresh O8 detector replay. The O7 final-test split was not used.

## Fresh detector replay

| Architecture | alpha90 | Pfa95 | Pd@ref | Pdl@ref | final gate |
|---|---:|---:|---:|---:|---|
| o8_a_scan | 569.455 | 6.804% | 0.714 | 0.201 | false |
| o8_pzt_k3q1 | 421.651 | 0.229% | 0.81 | 0.778 | false |
| o8_pzt_k3q4 | 316.136 | 0.129% | 0.9 | 0.884 | true |
| o8_ideal_k3q4 | 277.94 | 0.155% | 0.94 | 0.937 | true |
| o8_fmc | 219.187 | 0.075% | 0.968 | 0.966 | true |

PZT Q4 alpha90 gain over the same FEM head in q00 only: 25.02%.

Verdict: the frozen PZT boundary passes O8 and may enter an explicit patch/adhesive piezoelectric FEM.
