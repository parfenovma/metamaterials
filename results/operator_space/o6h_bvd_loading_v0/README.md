# O6h-A reciprocal q8 network with loaded BVD terminals

**PASS for a finite broadband-matched BVD region.**

The q8 contraction is embedded in an exact `23 x 23` reciprocal lossless scattering matrix with three matched dump ports; isometry error is `3.687e-16`.

With all 17 series-BVD loads matched at f0, the maximum effective Q across the 0.8--1.2 band is `1.48148`; active gate `public_return_loss`. The mismatch audit is run at half this value, Q=`0.74074`.

At audit Q, the 99%-Wilson load-error threshold is scale `2.2544`: resonance RMS `2.254%`, resistance log-RMS `4.509%`, Q log-RMS `11.272%`; Wilson lower-95 `99.005%`.

The frozen intrinsic-Q=35 load is `FAIL`: worst-band public return loss `0.069 dB`, useful amplitude `0.1261`. It therefore requires broadband/multiresonant matching to reach the effective-Q ceiling; this is not a global passive-network no-go.

The model includes reciprocal multiple reflections and matched dump ports, but remains normalized: a concrete matching circuit, transformer leakage and electrical-noise model are still required before fabrication.
