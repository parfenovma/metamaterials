# O6f-A concrete reciprocal q8 lattice

**PASS as an optimistic topology bound; electrical load and input-hybrid parasitics remain.**

The q8 `17 x 3` contraction is reconstructed to relative error `7.008e-16` as three singular-value attenuators, two ideal calibrated `3 x 3` hybrid blocks and a `45`-coupler adjacent-Givens output lattice. The small hybrids, dump-port loss and piezo impedance are counted but not yet assigned parasitics.

With every other error zero, the O4 ceilings are `0.047376 dB` identical insertion loss per output coupler (active gate: amplitude RMS) and `15.252 ns` identical delay per output coupler (active gate: differential delay). These are network-distributed ceilings, not datasheet limits for an isolated part.

For independent component errors, the largest common scale with Wilson lower-95 pass >=99% is `4.0732`: angle RMS `0.4073 deg`, delay RMS `8.1465 ns`; observed pass `99.120%`, lower bound `99.005%`, q01 useful `0.995964`.

This is a concrete reciprocal factorization, but not an electrical prototype. The remaining common loss budget must still pay for two input hybrids, three attenuator/dump paths, winding resistance, dielectric loss and mismatch to 17 capacitive piezo loads.
