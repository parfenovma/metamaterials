# Assembly and shunt tuning

## 1. Machine and baseline the block

1. Machine the Al 6061-T6 block from `coupon_block.step`.
2. Verify `50.60 x 202.38 x 8.00 mm`, overall tolerance `+/-0.05 mm`.
3. Verify source-face flatness `<=0.03 mm` and `Ra <=1.6 um`.
4. Mark the geometric centre and the four patch intervals from
   `pzt_placement.csv`; do not scribe inside a bond area.
5. Before bonding PZT, acquire the bare-block transfer baseline using the
   complete source fixture and the 17 receiver points in the acceptance
   protocol.  This is the amplitude reference and cannot be recreated later.

## 2. Incoming PZT inspection

1. Photograph every serialized plate and record the vendor polarity mark.
2. Measure dimensions, capacitance and loss.  Reject cracks, chips inside the
   active area, electrode delamination, or ambiguous polarity.
3. From the six plates choose two matched pairs.  Within each pair,
   `abs(C1-C2)/mean(C1,C2) <= 5%`.
4. Record the measured combined capacitance `C_pair=C1+C2`.  Do not force the
   design to the nominal `5.383 nF`.

## 3. Bond the plates

1. Mask the four `20.00 x 8.00 mm` bond areas.  Lightly abrade aluminium with
   a fresh fine abrasive, vacuum debris, then clean with reagent-grade IPA.
2. Clean the ceramic only by the vendor-approved method.  Do not bend or
   ultrasonically clean it.
3. Mix EPO-TEK 301-2 by the current manufacturer data sheet.  Apply enough to
   produce a continuous cured line of `25...50 um`; no dry islands or edge
   voids are allowed.
4. Locate all four plates to `+/-0.10 mm`.  The two B plates are flush with the
   top and bottom block edges.  Apply uniform dead-weight pressure only; do
   not point-clamp the ceramic.
5. Cure either `48 h at 23 C` or `3 h at 80 C`, following the adhesive data
   sheet.  Use the same cure for all four plates.
6. Inspect the perimeter under magnification.  Reject visible cracks or a
   bond void extending more than `1 mm` into an active edge.

## 4. Attach leads and freeze parity

1. Attach 36 AWG twisted pairs to the wrap-around electrode tabs with
   silver-filled conductive epoxy.  Add nonconductive strain relief after the
   electrical joint cures.
2. Wire A- and A+ in parallel with equal electromechanical polarity.
3. Wire B- and B+ in parallel with one plate reversed.  Capacitances still add;
   the mechanical coupling changes sign.
4. Keep each PZT-to-relay path below `50 mm`.  Photograph and continuity-check
   the finished polarity before potting or strain relief hides it.

## 5. Tune the two passive loads

For each measured pair capacitance, use

```text
L_total = 1 / ((2*pi*f_target)^2 * C_pair)
R_total = 2*pi*f_target*L_total / Q_target
```

with:

| pair | f_target | Q_target | nominal L_total | nominal R_total |
|---|---:|---:|---:|---:|
| A | 258.989 kHz | 6.382 | 70.15 uH | 17.89 ohm |
| B | 253.411 kHz | 4.232 | 73.27 uH | 27.57 ohm |

1. Start A with `68 uH + 2.2 uH` and B with `68 uH + 4.7 uH`.
2. Measure the assembled load at the PZT terminal plane near the target
   frequency.  Bin or replace the trim inductor until `L_total` is within 2%.
3. Measure inductor ESR, relay contact resistance and wire resistance together.
4. Fit thin-film resistors so their sum plus all measured parasitic resistance
   equals `R_total` within 5%.
5. Confirm that the relay-open state adds less than `20 pF` at the PZT port.
6. Label the two relay inputs A and B and verify the state truth table in
   `wiring.svg`.

Only one post-assembly retune of L/R is permitted by the acceptance protocol.
Patch positions, ceramic parity and statistical gates remain frozen.

