# O8: frozen plane-strain boundary realization — preregistered contract

> Frozen before O8 mesh generation and solves on 15 August 2026.

## Question and evidence level

Can the two frozen O7 stationary mechanisms survive an explicit
two-dimensional plane-strain elastic block when their boundary physics and
resources are represented rather than fitted as modal perturbations?

O8 uses a full plane-strain finite-element host and finite boundary zones.  A
PZT zone is condensed to its exact one-electrical-port constitutive Schur
complement, so electrodes and the common passive load are explicit but the
electric field through the patch is not spatially meshed.  This is a higher
evidence level than O7, but still one level below a domain-resolved
piezoelectric FEM.

## Frozen common head

- aluminium 6061 plane-strain block;
- carrier `242 kHz`, seven anchors `0.8...1.2 f0`;
- aperture `8 lambdaP` with 17 output bins;
- block length `2 lambdaP`;
- three fixed source tractions in reflection sectors: DCT0 and DCT2 even,
  DCT1 odd;
- matched absorbing backing on the source face and matched normal/shear
  radiation on the output face; free lateral edges;
- one fixed contact footprint and no motion.

The O7 normalized positions and widths become finite source-face intervals.
`position * (aperture/2)` is the centre and `width * (aperture/2)` is the
interval half-width.  Clipping at the aperture edge is recorded, not repaired
by optimization.

## O8-A: PZT/RL boundary

The O7 PZT geometry is frozen:

```text
pair A centre |y|/(H/2)=0.1983573, half-width=0.0992270
pair B centre |y|/(H/2)=0.9026094, half-width=0.1044452
```

Each patch has `C0=2 nF`, thickness `1.0 mm`, density `7750 kg/m3`, short
surface stiffness from `c33E=117 GPa`, thickness coupling `kt^2=0.45` and
mechanical loss factor `1/55`.  Both mirrored patches in a pair share one
electrical voltage and one passive series-`RL` shunt across their combined
`4 nF` capacitance.

The loads are frozen from O7:

```text
pair A: L=94.4099 uH, R=24.0716 ohm
pair B: L=98.6120 uH, R=37.0997 ohm
```

`bit=0` is open circuit and `bit=1` connects the passive `RL` shunt.  To cover
both reflection sectors without extra channels, pair A is wired with equal
polarity (even electrical functional) and pair B with opposite polarity (odd
functional).  This polarity choice is frozen analytically before the solve.
There remain exactly two switch bits and two loads.

For mechanical displacement `u`, patch coupling vector `h`, pair capacitance
`C` and load admittance `Y`, the electrical contribution is

```math
K_e(omega)=\frac{i\omega}{Y(omega)+i\omega C}\,h h^T,
\qquad Y=0\;\text{or}\;(R+i\omega L)^{-1}.
```

No negative capacitance, fitted complex gain or state-wise normalization is
allowed.

## O8-B: mechanical boundary

The O7 mechanical intervals, strengths and loss factors are frozen.  The
dimensionless on-state strength maps to surface stiffness through

```math
k_s=s\,\rho_{Al}c_P^2/\lambda_P.
```

The off state retains `2%` of the on stiffness and the actuator carries a
fixed `3.9 kg/m2` surface mass.  A switched pair applies the same local Robin
stiffness independently at `+y` and `-y`; it therefore preserves reflection
while acting in both parity sectors.  No actuator geometry or latch is tuned.

## Normalization and comparisons

Every source profile has unit `L2` traction norm.  Raw transfer, accepted
source work and output displacement are saved.  For operator replay, one
common dimensional scale per mechanism maps the maximum singular value across
all four states and seven frequencies to `0.98`.  State-wise and
frequency-wise rescaling is forbidden.  This scale is only the unresolved
common transduction calibration; it may not hide differential loss.

A same-mesh unmodified aluminium reference supplies the physical amplitude
ratio.  For every state/frequency, output Frobenius amplitude relative to the
reference must remain in `[0.70, 1.50]`.

## Mesh and gates

The frozen base mesh target is `lambdaS/6`; a carrier/edge spot check uses
`lambdaS/8`.  No geometric quantity is changed between meshes.

### Numerical symmetry addendum (before the band campaign)

A carrier-only implementation check, performed before any band ranking or O8
detector split was opened, showed `~0.9%` reflection error from Gmsh's
unstructured triangular mesh while port-work reciprocity was already below
`7e-16`.  Because the continuum geometry and boundary coefficients commute
exactly with reflection, every extracted transfer is henceforth projected by
the fixed group projector

```math
H_{Cs}=\tfrac12(H+R_{out}H R_{port}),
\qquad R_{port}=\operatorname{diag}(1,1,-1).
```

Both raw and projected errors are saved.  The raw base-mesh residual must be
`<=2%` and decrease or remain below `2%` on refinement; the preregistered
`1e-8` gate applies to the group-projected operator.  This projector is common
to both branches and all states/frequencies and contains no fitted parameter.

Mandatory full-wave controls:

```text
port-work reciprocity residual            <= 1e-6
raw mesh reflection residual               <= 2e-2
projected reflection equivariance residual <= 1e-8
all R, L, C, density, stiffness, loss      positive
minimum centre-frequency state distance   >= 0.05
zero switched coupling collapses to Q1    <= 1e-10
raw amplitude/reference                    0.70...1.50
refined/base projector RMS                 <= 0.10
refined/base amplitude relative error      <= 10%
```

If a branch passes the physical gates, its frozen `17 x 3 x 4 x 7` transfer
is replayed with new O8 nuisance, signal and H0 seeds.  Detector gates remain:

```text
Wilson upper 95% Pfa       <= 1.2%
Pd at fresh FMC reference  >= 0.90
Pdl at fresh FMC reference >= 0.85
finite alpha90
alpha90 gain over matched Q1 >= 5%
```

The old O7 final-test split is never opened.  A branch failing a physical
gate is not detector-tested and is stopped before a domain-resolved PZT or
actuator FEM.  Passing O8 authorizes O9 with explicit patch/adhesive geometry
and tolerance ensemble, not fabrication.
