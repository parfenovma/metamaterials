# O6j solid-contact attachment qualification

> Frozen after O6i selected one robust interior solid-contact property point.
> O6j is the final pre-prototype numerical qualification. P17 and liquid
> contact remain stopped; q8 remains blocked on a specified lossless or
> multi-resonant matching topology.

## Frozen P12 contact handoff

```text
ZP / sqrt(Zprobe Zsample)       0.75
cP / cP_Al                     1.00
cS / cP                        0.25
rho                            2053.61 kg/m3
physical thickness             505.96 um = 0.02 lambdaP_layer
```

No contact property or P12 layout is reoptimized.

## A. Angular attachment corners

The exact 2D three-solid P/S layer solve is repeated for all `2^4=16`
one-sigma sign corners:

```text
physical thickness log variation     +/-10%
P-wave impedance log variation       +/-5%
P-wave speed log variation           +/-2%
S-wave speed log variation           +/-5%
```

Every corner is checked at seven frequencies and angles
`0,15,30,45,55 deg`. Qualification over `0--45 deg` retains the O6h gates:
minimum layer/direct power `>=0.80`, change in transmitted S fraction
`<=0.10`; energy balance error must remain `<=1e-8`.

## B. Frozen detector replay with attachment plus encoder error

Only compact P12 is replayed. Encoder position error is fixed at its original
qualification RMS `0.05 lambda` using the frozen O6f standard draws. Solid
attachment thickness and impedance use the other frozen standard draws;
P-speed uses a fresh deterministic O6j draw with seed `0x06f0e117`:

```text
thickness log RMS       10% * attachment scale
P impedance log RMS      5% * attachment scale
P speed log RMS          2% * attachment scale
```

For each stop and frequency, an exact normal-incidence three-solid transmitted
P velocity is divided by the frozen nominal-contact value and squared for
pulse-echo. This complex response is applied to signal and nuisance at the
actual encoder-perturbed stop position; no affine amp/delay reduction is used.

Attachment scales `0,0.5,1,1.5,2,3,4` are extended and bisected through first
failure. The original P12 dictionary, train nuisance precision and threshold,
60 signal contexts, fresh O6f noise streams and 12,000 H0 trials per point are
frozen. Gates are unchanged:

```text
Pfa Wilson upper 95%               <= 1.25%
alpha90 / zero-error O6f replay    <= 1.10
Pdl at frozen reference amplitude  >= 0.80
```

## Decision

A conditional motion-prototype recommendation requires all 16 angular corners
and attachment scale `1.0` to pass. The finite detector attachment ceiling is
reported as a build tolerance. A pass does not certify wear, lateral
squeeze-out or a commercial material; it authorizes material identification,
finite-footprint modelling and a benchtop coupon. A failure leaves the project
without a prototype-ready architecture at the current evidence level.
