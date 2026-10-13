# O4 transfer target: static K4,Q2

> This is a realizable transfer-function contract, not a geometry and not a 17-channel electronics prescription.

The selected O3 operator is `trained_k4_q2_unconstrained`: four electronic ports, two sequential passive states and eight transmit events. Each state is a real orthonormal `17 x 4` map from four device ports to a discretized continuous contact aperture.

The target is

```text
H_q(omega) = g_q(omega) exp(-i omega tau_q) [S_q + E_q(omega)], q=1,2
reverse transfer = transpose(H_q)             (reciprocal)
sigma_max(H_q) <= 1                           (passive)
```

`S_q` is stored in `state_coefficients.csv`; the common scalar gain and delay are calibrated gauges. Only differential errors change the operator shape. The qualified band is `193.6--290.4 kHz` (seven anchors, five-cycle 242 kHz pulse).

## Hard acceptance numbers

- coefficient amplitude error: <= 3% RMS; differential group-delay error: <= 82.6 ns RMS (7.2 deg at carrier);
- normalized transfer fidelity >= 0.98; maximum port coherence <= 0.15; column imbalance <= 0.75 dB;
- useful one-way amplitude >= 0.9631, hence one-way insertion loss <= 0.327 dB and round-trip loss <= 0.654 dB;
- reflected plus parasitic one-way power <= 7.25%; reciprocity residual <= 3%;
- calibrated defect-dictionary fidelity >= 0.90 and the frozen F1c acceptance replay must still pass.

The loss bound is derived from the observed O3 margin `alpha90=2.7824` at `reference alpha=3.0`. Because pulse echo traverses the network twice, the one-way amplitude floor is the square root of the round-trip margin.

The 10,000-draw tolerance envelope is deterministic (`seed=0x04a42026`). It is a consistency envelope for the already-tested 3%/0.02-cycle perturbation law, not additional evidence or a replacement for the F1c replay.

## Manufacturability diagnostic

The exact unconstrained target is not a low-order lens field. Its four state-to-state principal angles are `32.5, 50.1, 68.7, 79.6 deg`; individual ports need `11--17` DCT aperture modes for 95% energy and contain `5--9` sign changes. O5 must therefore test a geometry-aware smooth projection before any dense geometry optimization.
