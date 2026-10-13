# O6h loaded-network and physical-contact contract

> Frozen after O6g and before O6h calculation. O6h raises the evidence level;
> it may not optimize layouts, q8 coefficients or detector thresholds.

## O6h-A: reciprocal q8 scattering completion with BVD loads

The frozen q8 transfer `B` is embedded into an exactly isometric forward map
by adding three matched dump ports:

```math
F = [B; (I-B^T B)^{1/2}], \qquad F^T F=I.
```

The reciprocal lossless scattering completion is

```math
S = [0, F^T; F, I-FF^T].
```

Only the first 17 output ports carry piezo loads; the three dump ports remain
matched. For diagonal load reflection `Gamma`, the loaded output and public
reflection are calculated with the full multiple-reflection resolvent

```math
x=(I-(I-FF^T)Gamma)^{-1}Fu,
R_public=F^T Gamma x.
```

Every piezo terminal uses a normalized series-BVD impedance, already matched
at `f0`:

```math
z_i(r)=R_i [1 + i Q_i (r/r_i-r_i/r)],
Gamma_i=(z_i-1)/(z_i+1),
t_i=2/(1+z_i).
```

The mechanical drive observable is `diag(t) x` on the 17 physical ports. The
deterministic effective-Q sweep is `0, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 35`,
with bisection at first failure. Gates are the frozen O4 transfer gates plus
worst-band public return loss of at least `10 dB` (`sigma_max(R)<=0.31623`).

At half the deterministic Q ceiling, `20,000` common draws use independent
load errors:

```text
resonance-frequency RMS       1% * scale
matched-resistance log RMS    2% * scale
effective-Q log RMS           5% * scale
```

The largest scale with one-sided Wilson lower-95 `>=99%` joint pass is
reported. No per-load retuning, post-load rebalance or global contraction is
allowed. This is a normalized BVD/matching-bandwidth gate, not a component
datasheet claim.

## O6h-B: exact thin-layer contact transfer for motion

The contact is an exact normal-incidence three-medium transfer matrix,
including all multiple reflections in the layer. Frozen nominal values are:

```text
probe-face impedance          17.0 MRayl
Al sample impedance           rho_Al * cP_Al = 16.5297 MRayl
couplant impedance            1.5 MRayl
couplant wave speed           1480 m/s
carrier                       242 kHz
nominal layer thickness       0.01 lambda_c = 61.16 um
```

For layer matrix `M`, the pressure-wave transmission is converted to a
power-wave amplitude and divided by the zero-thickness direct-interface
value. Pulse-echo receives the squared one-way transfer. A nominal-layer
frequency response is calibrated once and divided out; only stop-to-stop
departures remain. The exact complex transfer, not an affine amplitude/delay
fit, multiplies signal and all nuisance atoms at every frequency.

The physical variation scale uses the frozen O6f standard draws:

```text
layer-thickness RMS           0.005 lambda_c * scale
couplant-impedance log RMS    10% * scale
couplant-speed log RMS        5% * scale
```

Thickness is reflected at zero rather than allowed to become negative. The
scale grid `0, 0.5, 1, 2, 3, 4` is extended and bisected through the first
failure for frozen P12/P17. The same 60 signal contexts, fresh O6f noise
streams and 12,000 H0 trials per point are reused. Gates remain Pfa upper
`<=1.25%`, alpha90 ratio `<=1.10`, and reference Pdl `>=0.80`.

The transfer residual relative to its best per-stop log-amplitude plus linear
delay fit is also reported. A large residual demonstrates physics omitted by
O6g; a small residual validates the earlier abstraction in this regime.

## O6h decision

- q8 survives only if a finite effective-Q/mismatch region passes both O4 and
  `10 dB` return loss. A high-Q BVD failure is a requirement for broader or
  multi-resonant matching, not a no-go for every passive network.
- Motion survives only if the nominal physical contact specification passes.
- Only after these gates may one branch receive a 2D angular/full-wave spot
  check. A failed branch is not rescued by another reduced parameter sweep.

## O6h-C addendum: preregistered 2D angular spot check

O6h-B passed the nominal layer specification for P12 and failed it for P17;
therefore only P12 is eligible for this spot check. Before calculation, the
normal-incidence layer is promoted to an exact 2D solid--fluid--solid
plane-wave boundary solve. It includes reflected/transmitted P and SV waves,
forward/backward fluid waves, tangential-wavenumber conservation and all
layer reflections.

The probe face uses the frozen `17 MRayl` impedance with Al-like wave speeds;
the sample uses `rho=2700 kg/m3`, `cP=6122.1024 m/s`, `cS=3083.8103 m/s`.
The fluid properties and nominal thickness remain those of O6h-B. Frequencies
are the seven operator anchors and incident P angles are
`0, 15, 30, 45, 55 deg`. A zero-thickness direct solid--solid interface is
the paired reference at every frequency and angle.

Numerical positive controls:

```text
lossless energy balance |R+T-1|       <= 1e-8
normal-incidence transmission          agrees with 1D layer <= 1e-8
```

Qualification gates over `0--45 deg` (55 deg is diagnostic):

```text
minimum layer/direct transmitted power >= 0.80
maximum change in transmitted S fraction <= 0.10
```

The same observables are reported at the P12 nominal one-sigma adverse
anchors (`d0 +/- 30.58 um`, `Zc * exp(+/-10%)`, `cc * exp(+/-5%)`) without
forming an unregistered worst-sign optimizer. Failure blocks a motion
prototype claim until contact geometry or calibration explicitly includes
angular P/S physics.
