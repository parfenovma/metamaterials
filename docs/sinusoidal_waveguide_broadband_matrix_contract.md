# Sinusoidal Al waveguides: full broadband response-matrix contract

> Frozen after the carrier reconstruction pilot and before solving the other
> sixteen frequencies on 4 October 2026.

## Question

Across the already frozen `162.2--321.7 kHz` five-cycle band, does the selected
15-strip sinusoidal aluminium lens provide more power-normalized spatial modes
or a better signed time-reversal focus than the matched straight-strip
aperture?

## Fixed physical model

- the selected and straight geometries, materials, `1 MPa` traction convention,
  linear element order, mesh-size prescription, absorbing boundaries, and 17
  frequencies are inherited unchanged from the permanent P6 package;
- the tagged carrier response-matrix meshes are reused at every frequency;
- all 15 strips are excited separately, producing 15 right-hand sides after
  one matrix assembly and LU factorization per geometry/frequency;
- longitudinal displacement is sampled at all 15 signed strip-centre
  coordinates on `x=33 mm`, five millimetres after the strips;
- no geometry or frequency is selected from the new result.

The campaign therefore requires `17 x 2 = 34` factorizations and 15 backsolves
per factorization, rather than 510 independent assemblies/factorizations.

## Carrier control and numerical symmetry

At `242 kHz`, summing the individual strip columns into the old centre-plus-
seven-pair basis must reconstruct the permanent `8 x 8` carrier matrix with
relative error `<=1e-10`. The observed control errors are `2.74e-15` for the
sinusoidal lens and `1.69e-15` for the straight aperture.

The unstructured meshes introduce a raw reflection residual of about one
percent even though the continuum geometry is exactly symmetric and input
work reciprocity is near machine precision. The raw output and port-work
matrices are retained and must satisfy

```text
raw reflection residual       <= 2e-2
port-work reciprocity         <= 1e-8
minimum diagonal accepted work > 0
```

For modal and time-reversal analysis, both matrices use the fixed reflection
projector

```math
H_s=\tfrac12(H+RHR),
\qquad
W_s=\tfrac12(W+RWR).
```

No data-dependent smoothing or per-column normalization is allowed.

## Physical input-power normalization

Let `W=B^T A^{-1}B` be the complex source work/compliance matrix. With the
solver's harmonic convention, the accepted-power quadratic form is

```math
Q(f)=-\frac{\omega}{2}\frac{W_s-W_s^\dagger}{2i}.
```

`Q` must be positive definite within numerical tolerance. The transfer used
for modal comparison is

```math
\widetilde H(f)=H_s(f)Q(f)^{-1/2}.
```

Thus a unit Euclidean input coefficient has the same accepted power even when
strip coupling is non-negligible. The eigenvalue range and any regularization
required by `Q` are reported. A negative eigenvalue beyond `1e-8` of the
largest eigenvalue stops the analysis.

## Frozen diagnostics

The five-cycle Hann amplitude weights include trapezoidal quadrature on the
nonuniform 17-point grid,

```math
w_f^2\propto |S(f)|^2\Delta f_{\mathrm{trap}},
\qquad \sum_f w_f^2=1.
```

For each geometry,

```math
A=[w_f\widetilde H(f)]_f,
\qquad K=AA^\dagger.
```

The audit reports:

1. broadband singular spectrum, stable rank, entropy rank, and observable mode
   counts at 20, 30, 40, and 60 dB;
2. the same spectra separately in the eight-dimensional even and
   seven-dimensional odd parity sectors;
3. signed focal-bin validity, sampled FWHM, focus contrast, and sidelobes of
   `K` over all 15 positions;
4. frequency-wise mode counts and conditioning;
5. all numerical gates, pulse weights, and the complete complex transfer,
   work, whitened-transfer, and kernel arrays.

This is a deterministic information/focusing audit. It does not introduce a
new stochastic noise draw or update probability of detection.

## Post-failure refinement addendum

The blind band run exceeded the raw `2%` reflection gate at the upper
sinusoidal anchors (first observed before analysis). The gate is not relaxed.
Before opening modal/time-reversal results, both geometries at the fixed worst-
case upper edge `321.7 kHz` will be rebuilt with linear mesh targets reduced
from `0.42/1.55 mm` to `0.32/1.20 mm` in the lens/output.

The projected base-grid campaign is authorized for analysis only if all hold:

```text
refined raw reflection residual                <= 2e-2
refined port-work reciprocity                   <= 1e-8
refined accepted-power metric                   positive definite
base/refined normalized power-whitened kernel   <= 10% relative Frobenius
base/refined normalized singular spectrum       <= 10% relative L2
```

This addendum was written after seeing the raw symmetry failure but before
running the refinement or any broadband modal/time-reversal analysis. If it
fails, the result remains provisional and no book-guided conclusion is made.
