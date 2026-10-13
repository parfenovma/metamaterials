# Sinusoidal aluminium waveguides: carrier time-reversal audit

## Result

At `242 kHz`, the selected sinusoidal lens does **not** add spatial channels relative to the matched straight strips. Both matrices retain all eight reflection-even modes within 20 dB, while the sinusoidal spectrum is less balanced.

Under the nominal source-area normalization:

| geometry | stable rank | entropy rank | condition number | modes @20 dB | median contrast / dB | median sidelobe outside FWHM / dB |
|---|---:|---:|---:|---:|---:|---:|
| sinusoidal | 2.245 | 5.035 | 5.197 | 8 | 9.32 | -4.7 |
| straight | 2.844 | 5.667 | 3.201 | 8 | 10.07 | -5.97 |

The sinusoidal stable rank is 78.9% of the straight value, its condition number is `1.62x` larger, and its median time-reversal contrast changes by -0.75 dB. The sign of all three comparisons is the same under the original as-solved traction convention.

This does not contradict the existing `1.633x` carrier focus gain under the fixed uniform drive. That design optimized one coherent focal combination; it did not optimize the full singular spectrum or general-purpose channel conditioning. In book terminology, it redirects energy into one useful wavefront but does not increase the carrier information capacity of this even-sector basis.

## What the kernel means

`K = H H^dagger` is computed on the plane five millimetres after the strips. Its columns are virtual phase-conjugate refocusing profiles over `|y|`. Because every source is a mirror-symmetric strip pair, an off-axis target represents two equal physical foci at `+y` and `-y`; signed lateral localization has not been demonstrated.

The sampled half-aperture spacing is `8.2 mm = 0.324 lambdaP`. Every sampled half-maximum width is one bin, so the present matrix cannot compare sub-bin resolution. Contrast and sidelobe metrics remain resolved.

## Scope and next physical calculation

This is a one-frequency result. A genuine broadband time-reversal test needs the same eight source-group responses and eight output samples on the frozen 17-frequency grid, plus complex source work for exact power normalization. Odd source groups are also required before claiming signed lateral localization. Until those data exist, no conclusion about temporal compression, Q, ringing, or broadband super-resolution is made.
