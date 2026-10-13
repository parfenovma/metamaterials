# O6e passive-beamformer versus motion resource fork

**The passive-beamformer branch clears its algebraic pre-geometry gates; neither branch dominates the full resource vector.**

## Passive 50-tap branch

The nominal q8 map has useful amplitude `0.999053`, leaving at most `0.913 dB` one-way common insertion loss before the O4 floor. Its frozen F1d replay passes: `Pfa95=0.089%`, `Pd=0.982`, `Pdl=0.976`, `alpha90=2.4298`, gap closure `0.909`; alpha90 changes by `0.137%` from ideal DCT3.

The largest common manufacturing-error scale with Wilson lower-95 pass probability >=99% is `0.7516`: log-amplitude RMS `2.255%`, delay RMS `62.081 ns`. At that point observed pass is `99.076%`, lower bound `99.003%`, q01 useful amplitude `0.930920`, leaving `0.299 dB` additional common-loss headroom.

At the full O4 reference perturbation (`3%`, `82.6 ns`) only `20.694%` of
draws pass jointly, so the q8 network does not inherit the entire coefficient
tolerance budget: fabrication and network errors together must stay below the
reported scale.

Single-open-tap audit: `2/50` faults still pass O4 after optimistic fixed port recalibration. This is diagnostic only; no fault tolerance is credited to the architecture.

## Unsummed resource frontier

| Branch | Ch | Events | Positions / moves | T electronics / settle / travel / slow | Columns | Passive taps | Gap closure* |
|---|---:|---:|---:|---:|---:|---:|---:|
| `beamformer_q8` | 3 | 3 | 1 / 0 | 3.0 / 3.0 / 3.0 / 3.0 | 17 | 50 | 0.909 |
| `motion_compact_center_p12_native` | 1 | 12 | 12 / 11 | 12.0 / 23.0 / 28.5 / 67.0 | 1 | 0 | 0.868 |
| `motion_uniform_full_p17_a8` | 1 | 17 | 17 / 16 | 17.0 / 33.0 / 41.0 / 97.0 | 1 | 0 | 0.931 |

`*` Gap closures are normalized inside different evidence campaigns: full-wave frozen F1d for beamformer q8, reduced fresh Born O3 for motion. They are not a shared ROC number.

Relative to P12/P17, the static beamformer spends two extra electronic channels and 50 calibrated passive taps, but saves `9/14` transmit events and `11/16` moves. Under the existing time models it saves `9--64` normalized units versus P12 and `14--94` versus P17. P12/P17 avoid the distributed network but retain unresolved encoder, couplant-reset and probe-footprint risks.

## Decision

Both branches remain open. The next paired physical-risk stage is: one realizable passive-network topology with measured/estimated parasitics for q8, and one positioning/contact perturbation campaign for P12/P17. A scalar winner remains forbidden until those missing coordinates exist.
