# O7 stationary global-state gate

Evidence level: coupled-modal surrogate; not coupled piezo FEM.

| Frozen candidate | alpha90 | Pfa95 | Pd@ref | Pdl@ref | Q4 gain over matched Q1 | Final gate |
|---|---:|---:|---:|---:|---:|---|
| pzt_k3q4_frozen | 285.147 | 1.134% | 0.987 | 0.851 | 17.84% | true |
| mechanical_k3q4_frozen | 270.672 | 0.796% | 0.996 | 0.95 | 34.16% | true |

Positive controls: PZT `sigma=0.9799999999999999`, reflection error `2.4824555417333764e-15`, minimum state distance `0.34023607565617514`; mechanical `sigma=0.9800000000000002`, reflection error `2.3977649908846795e-15`, minimum state distance `0.45115401875123007`.

The RLC values in `component_translation.csv` assume two 2 nF mirrored patches in parallel per bit. They translate the reduced resonance only; they are not a matched BVD terminal design.

Verdict: at least one stationary two-bit mechanism passes the reduced O7 gate. The next authorized calculation is a frozen coupled plane-strain physical model.
