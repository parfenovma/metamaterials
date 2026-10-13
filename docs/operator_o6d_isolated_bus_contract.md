# O6d mechanically isolated three-bus bound

> Frozen before the O6d calculation. O6c showed that reducing the effective
> lateral stiffness of one continuous layer improves O4 but does not remove
> differential spatial propagation. O6d separates the mechanical-isolation
> question from the electrical-coding burden before any kerf geometry exists.

## 1. Frozen object

The load face is represented by the same `17` spatial bins and the same three
orthonormal DCT0--2 target columns. The hypothetical transducer has exactly
three public reciprocal electrical terminals, one passive state and three
transmit events. It has no per-pixel phase, delay, switching, pulser or ADC.

Mechanical columns are treated as ideally load-local for the electrical
coding bound. Conversely, the mechanical-isolation bound uses the measured
O6c transfer while granting an ideal exact DCT bus map. These two optimistic
bounds deliberately do not compound unrelated penalties.

## 2. Mechanical residual bound

At every O6c frequency, split the calibrated transfer into its common target
part and differential residual:

```text
H_o6c(omega) = g(omega) S + D(omega)
H_kappa(omega) = g(omega) S + kappa D(omega),  0 <= kappa <= 1.
```

`kappa=1` is the measured continuous O6c layer; `kappa=0` is perfect local
isolation at the same common useful amplitude. The largest `kappa` satisfying
all O4 gates is the required residual-coupling ceiling. This interpolation is
an operator requirement, not a kerf-width or material law.

## 3. Electrical connection classes

An ideal lossless three-port rotation is allowed before the internal buses;
this is generous to sparse wiring and costs no additional electronic channel.
For each isolated column, retain the best continuous signed coupling to:

1. at most one internal bus (`<=17` active taps);
2. at most two internal buses (`<=34` active taps);
3. all three buses (`<=51` taps; exact distributed-beamformer upper bound).

For a fixed rotation and sparsity, orthogonal projection supplies the optimal
least-squares taps. A deterministic `200,000`-rotation SO(3) campaign with
local coordinate refinement searches the first two classes. Because this is
not a certified global optimizer, a close miss cannot establish a no-go; a
large miss may stop only the tested low-complexity family. The full three-bus
map is analytic and exact.

Every candidate is granted fixed column balancing and a final contractive
scale. Reverse transfer is its transpose, network loss is zero and
differential electrical delay is zero. These are optimistic upper-bound
assumptions: failure cannot be repaired by a real passive network.

## 4. Weight-resolution audit

The exact full map is additionally quantized to `1, 2, 4, 8, 16` nonzero
magnitude levels per DCT port using deterministic one-dimensional Lloyd
quantization. Signs and symmetry zeros are retained. This records how much
analog electrode/transformer resolution is required; it is not a geometry
sweep.

## 5. Gates and decision

All O4 transfer gates remain unchanged: fidelity `>=0.98`, coherence
`<=0.15`, imbalance `<=0.75 dB`, differential amplitude RMS `<=3%`, useful
amplitude `>=0.8994` and `sigma_max<=1`. Electrical reciprocity and delay pass
by construction in this ideal bound.

The selected low-complexity O6d class passes only if a map with at most two
bus connections per isolated column passes all gates. A full three-bus map is
reported as a passive `3-to-17` beamformer fallback, not as a pass of the
monolithic low-complexity class already selected in O5.

- Low-complexity failure blocks kerf FEM and prototype geometry.
- A low-complexity pass authorizes exactly one macro-column twin at the
  calculated isolation requirement.
- A full-map-only pass reopens an explicit product-level choice between the
  passive beamformer burden and another point on the architecture frontier; it
  does not silently rename the beamformer a lens.
