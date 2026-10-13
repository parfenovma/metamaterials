# O6 physical twin contract: direct DCT3 modal plate

> Frozen before evaluating the dense-band transfer. The purpose of O6 is to
> test one physical realization of the O4 target, not to restart a geometry
> search.

## 1. Candidate

The first `translation_kz` realization is a homogeneous plane-strain Al 6061
transition plate:

```text
three fixed distributed reciprocal piezoelectric couplers
    -> 71.6 mm continuous aperture x 3.0 mm Al transition plate
    -> perfectly bonded Al specimen.
```

The electrical ports impose the first three continuous DCT displacement
profiles on the input face: uniform even, first odd and second even. Each is a
single global electrode/polarization pattern. There are three pulser/receiver
channels, one passive state and three transmit events; there are no
element-wise phases or independently wired aperture cells.

The plate is deliberately uniform. Its only job is to transport the modal
profiles into the specimen. Any spatial synthesis is assigned to the three
distributed couplers, so a successful result does not authorize a taper,
resonator library or topology optimization.

## 2. Why the length is fixed at 3.0 mm

For a scalar longitudinal lateral harmonic `m` in a strip of width `W`,

```math
k_{x,m}=\sqrt{(\omega/c_P)^2-(m\pi/W)^2},\qquad
\tau_{g,m}/L=\frac{1}{c_P\sqrt{1-(f_{c,m}/f)^2}},
```

with `f_c,m=m c_P/(2W)`. At `W=71.6 mm`, `c_P=6122.10 m/s` and the lower
O4 band edge `193.6 kHz`, the worst required harmonic `m=2` gives a maximum
length of about `4.4 mm` at the O4 differential-delay gate `82.6 ns`.

The candidate length is frozen at `3.0 mm`, below `0.7` of that optimistic
upper bound. It is not selected from a sweep. The exact elastic gate below
uses all P/S-mixed free-strip modes and may still reject it.

## 3. Modal model

At every frequency the existing quadratic elastic port pencil is solved on
the full `71.6 mm` cross-section. Free upper/lower faces are natural boundary
conditions. All right-going propagating modes in the matching reflection
sector are retained; the prescribed DCT input is represented by the
minimum-norm outgoing modal radiation solution. Propagation over `3.0 mm`
uses each computed complex wavenumber.

This is an optimistic model of the distributed piezo coupler but a complete
model of propagation inside the uniform transition plate. It therefore has
two separate gates:

1. `coupler/source fidelity >= 0.98` throughout the band;
2. the propagated, passivity-contracted `H(omega)` must satisfy all O4
   transfer gates.

Failure of the first gate rejects the assumed three global electrode
profiles. Failure only after propagation rejects the uniform plate.

## 4. Frozen evaluation

- band: `193.6--290.4 kHz`;
- dense grid: 49 uniformly spaced frequencies, including the seven O4
  anchors;
- port cross-section: 160 linear elements;
- output sampling: the frozen 17 O4 aperture bins;
- fixed carrier drive calibration: equal incident modal power at `242 kHz`;
- one global contraction applied to the whole band so
  `max_omega sigma_max(H)<=1`;
- no per-frequency or per-cell correction.

Common complex gain and common delay are calibrated gauges. Differential
amplitude, delay, column coherence and imbalance are not removed.

## 5. Gates and stopping rule

The physical twin passes only if, on the dense grid:

```text
source-profile fidelity             >= 0.98
normalized propagated fidelity      >= 0.98
maximum port coherence              <= 0.15
column imbalance                    <= 0.75 dB
differential amplitude RMS          <= 3%
differential delay RMS              <= 82.6 ns
one-way useful amplitude            >= 0.8994
maximum singular value              <= 1
reciprocity residual                <= 3%
```

The modal result authorizes exactly one 2D FEM confirmation of the same
rectangle. It does not authorize tuning its length. If the modal gate fails,
no geometry is built and the next allowed action is to revise the coupler
mechanism or use the frozen `K4,Q1` fallback.

After a successful FEM transfer check, the measured `H(omega)` is inserted
into the frozen expanded F1d replay. Prototype machining is forbidden until
that replay passes.
