# O6f-B frozen motion physical-risk replay

The saved P12/P17 GLRT thresholds, nominal templates and training nuisance model were reconstructed without retraining. Signal, localization and H0 false-alarm counts reproduce the saved baseline exactly.

Each row uses 60 frozen signal contexts and 12,000 fresh aggregate H0 trials. Encoder, contact-amplitude and contact-delay errors are paired across architectures and nested across scales.

- `motion_compact_center_p12_native`: maximum tested joint-pass scale `2.00` = position RMS `0.1000 lambda`, contact log-amplitude RMS `10.00%`, contact delay RMS `0.0400 cycles`; all tested scales pass: `true`.
- `motion_uniform_full_p17_a8`: maximum tested joint-pass scale `2.00` = position RMS `0.1000 lambda`, contact log-amplitude RMS `10.00%`, contact delay RMS `0.0400 cycles`; all tested scales pass: `true`.

This remains a reduced Born replay: it includes encoder and multiplicative/delay contact errors, but not wear, lift-off-dependent mode conversion, surface waves or scan-time economics.
