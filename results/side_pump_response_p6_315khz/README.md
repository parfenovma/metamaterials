# Symmetric side-pump screen at 315 kHz

## Question

Can one inexpensive symmetric side-patch channel partly decouple harmonic
focusing from pulse-echo detectability for the c005 ellipse plate?

The tested hardware grammar has one central front PZT channel and one common
channel feeding equal 12 mm patches on the top and bottom plate edges.  The
sweep varies side mechanical-energy fraction, side/front carrier phase,
traction angle, longitudinal patch position and carrier frequency.  Every
harmonic candidate is normalized to the same mechanical input power as the
straight front-driven reference.

## Numerical campaign

- 2 topologies, 3 frequencies and 13 harmonic RHS per case;
- 2,976 coarse algebraic candidates;
- 992 fine harmonic candidates at 315 kHz;
- 26 fine time-domain sin/cos response histories per paired defect/free mesh;
- 496 pulse-echo candidates reconstructed from the transient basis;
- coarse c005 mesh: 29,127 nodes / 56,944 triangles;
- fine c005 mesh: 60,404 nodes / 118,881 triangles;
- transient mesh: 8,711 nodes defect-free and 8,836 nodes with the 4 mm void.

## Harmonic result

The target-locked coarse winner survives the independent fine mesh without
retuning.

| metric | coarse | fine |
|---|---:|---:|
| side effort fraction | 0.25 | 0.25 |
| phase | 225 deg | 225 deg |
| traction angle | 30 deg | 30 deg |
| patch center | 47 mm | 47 mm |
| J1 focus amplitude gain | 1.3939 | 1.4550 |
| harmonic bistatic J2 proxy | 1.6119 | 1.7422 |
| exit-energy proxy | 1.2391 | 1.3019 |
| FWHM | 17.48 mm | 16.81 mm |
| focus error | 0 mm | 1 mm |
| direct source-to-PZT coupling proxy | 0.708 | 0.700 |

This is a real fine-mesh harmonic pass.  It is not yet a defectoscope pass.

## Paired pulse-echo result

The pulse model uses a five-cycle 315 kHz burst, a mechanical PZT-epoxy-Al
stack, a 4 mm through-void at x=115 mm and paired defect/free subtraction.
Signals are compared at equal integrated mechanical input work.

The best target-locked absolute system point in the tested two-channel grammar
is eta=0.25, phase=225 deg, traction angle=45 deg and x=47 mm:

| metric | versus straight | versus c005 front-only |
|---|---:|---:|
| harmonic J1 | 1.4496 | 1.3128 |
| harmonic J2 proxy | 1.7685 | 1.3306 |
| matched-filter echo | 0.4849 | 0.8044 |
| pulse-echo CNR amplitude ratio | 0.00727 | 1.3771 |
| focus error | 1 mm | — |

Thus the side pump improves both harmonic focus and actual CNR relative to the
same c005 plate, but it still loses matched-filter echo amplitude and remains
far below the monolithic straight pulse-echo reference.  No tested point
crosses the absolute pulse-echo acceptance line.

The gentler eta=0.10, phase=270 deg, normal-traction point preserves 95.5% of
the c005 matched-filter echo while improving c005 CNR by 13.6% and harmonic J1
by 25.7%; it is the least-worst incremental compromise.

## Verdict

Side pumping is a useful second degree of freedom and the harmonic effect is
not a mesh artifact.  It does not rescue the present perforated c005 plate as
a defectoscope.  The dominant blocker is long-lived lens ringing, not lack of
focal amplitude.  The honest next topology should reduce stored resonant
energy first (weak isolated perturbations, lossy/matched boundary treatment or
a non-perforated phase transformer), then reuse the same side-pump response
basis as a small correction channel.

The 225 degree carrier relation is equivalent to a -135 degree phase.  At
315 kHz this is about 1.19 microseconds, or roughly 7.3 mm of longitudinal path
in aluminium.  It can be made with a narrowband passive delay or an analog
phase network, but broadband operation must be revalidated.

## Limits

- two-dimensional plane-strain model;
- mechanical isotropic PZT surrogate; no electromechanical/electrical matching;
- ideal symmetric prescribed side traction and ideal coupling;
- deterministic defect subtraction, not a calibrated probability of detection;
- one defect size/location and one pulse length;
- harmonic J2 is a bistatic field proxy and is demonstrably optimistic relative
  to the paired time-domain CNR.

## Reproduction

```bash
julia --startup-file=no side_pump_response.jl --run-coarse
julia --startup-file=no side_pump_response.jl --report coarse
julia --startup-file=no side_pump_response.jl --basis straight 315 fine
julia --startup-file=no side_pump_response.jl --basis c005 315 fine
julia --startup-file=no side_pump_response.jl --fine-report-315
julia --startup-file=no side_pump_pulse_echo.jl --report
```

Primary data are in `data/side_pump_sweep_coarse.csv`,
`data/side_pump_sweep_fine_315khz.csv` and
`data/side_pump_transient_sweep_315khz.csv`.  Mesh, field, Pareto and virtual
defectoscope figures are in `figures/`.
