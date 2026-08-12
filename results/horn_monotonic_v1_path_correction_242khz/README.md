# Monotonic horn lens v1 — physical path-correction verdict

## Verdict

Первая и единственная разрешённая path-length correction проверена сначала на
изолированном максимальном `sin^4`-канале, затем в полной 15-канальной 3D
линзе. Изолированный pulse-fidelity gate проходит, но carrier gate полной
линзы проваливается с большим запасом. V1 отвергнута; общий 3D transient не
разрешён и геометрический sweep не запускается.

## Maximum-delay transient

Исправленный maximum-state имеет осевую длину `50.664 мм`, развёрнутый путь
`57.620 мм`, extra path `6.957 мм` и bend amplitude `11.374 мм`. Три
полноупругих 2D transient-модели использовали одинаковый `5-cycle Hann @
242 кГц`, материал, источник, temporal step и receiver:

```text
peak / straight equal path      0.9691    >= 0.80
energy / straight equal path    0.9619
broadening ratio                1.0000    <= 1.25
pulse correlation               0.99988   >= 0.90
postcursor ratio                0.09246   <= 0.10
transverse energy ratio         0.00227   <= 0.10
```

Pulse-fidelity gate проходит. Однако measured delay относительно прямого
device равна только `3.987 мкс` против target `5.777 мкс`. Это было первым
признаком, что центральная монотонная экстраполяция broadband slope не
переносится в физический maximum-state.

## 15-channel carrier gate

V1 была построена на symmetry-constrained half-domain. Phase4 mesh содержит
`184 716` узлов, но два независимых direct solve воспроизводимо упали в
системном macOS `_xzm` allocator после LU assembly. Оба неполных JLD2 были
проверены как unreadable и не используются.

Carrier verdict получен на `phase3`:

```text
path h       0.387 мм   (h/lambda_S = 0.08)
focus h      0.774 мм
receiver h   1.209 мм
nodes        142 904
elements     817 287
```

Matched straight reference геометрически идентичен сертифицированному v0
reference: path correction не меняет mouths, source faces, axial thickness,
receiver или mesh prescription. Равенство `common_axial_length_mm` проверено,
поэтому использовано существующее phase3 straight field того же уровня.

| Carrier metric | V1 | Gate | Result |
|---|---:|---:|---|
| Lens focus | 1.809 нм | — | −24.89% к v0 |
| Matched gain | 1.604 | ≥2.0 | fail |
| Central/outer channel | 0.711 | ≥0.85 | fail |
| Maximum output phase error | 165.1° | ≤10° | fail |
| FWHM | 10 мм | ≤6 мм | fail |
| Sidelobe/peak | 0.552 | ≤0.40 | fail |
| Mirror amplitude mismatch | 0 | ≤3% | pass |
| Mirror phase mismatch | 0° | ≤3° | pass |

V0 на том же уровне имеет focus `2.409 нм`, matched gain `2.135` и maximum
phase error `23.88°`. Уменьшение paths `9.724 -> 6.957 мм` изменило carrier
phase почти на противоположную и расширило фокус. Разница настолько велика,
что неизвестность phase3→phase4 порядка процента не может изменить stop.

## Outlet/receiver coupling audit

Полные `8x8` response matrices на `193.8 / 242 / 290.1 кГц` были проверены
двумя observables: diagonal preout transfer и вклад каждого источника в
фокус. Для семи надёжных групп их group-delay estimates отличаются RMS на
`0.607 мкс`, то есть broadly подтверждают один и тот же broadband slope.

Carrier cross-talk мал:

```text
maximum lens off-diagonal L2 / diagonal       0.0643
maximum straight off-diagonal L2 / diagonal   0.0622
dominance threshold                            0.20
```

Следовательно, failure вызван не сильной соседней связью. Причина —
несовместимость broadband phase slope и carrier phase intercept для
path-only регулятора. Центральная группа нелинейна по частоте, а коррекция
остальных задержек разрушает carrier phasing. По правилу P5.1D кандидат
возвращается в полосовую модель как dispersively incorrect; независимый
wrapped `2pi` trim не добавляется.

## Содержимое

- `horn_monotonic_maximum_transient.png` и `maximum_delay_summary.csv` —
  corrected maximum-state transient gate;
- `signals/*.jld2` — три transient signal sets;
- `transient_meshes/*.msh`, `transient_models/*.json` — воспроизводимые 2D
  transient inputs;
- `horn_monotonic_v1_carrier_gate.png` — focal profile, channel amplitudes и
  phase error;
- `v1_carrier_gate_summary.csv`, `v1_carrier_channels.csv` — carrier gates;
- `horn_monotonic_broadband_coupling.png` и `broadband_coupling_*.csv` —
  outlet/focus/cross-talk audit;
- `field_half_lens_phase3.jld2` — сохранённое комплексное v1 carrier field;
- `mesh_half_lens_phase3.msh`, `mesh_half_lens_phase3.jld2` — фактически
  использованная phase3 mesh и её размерные/безразмерные параметры;
- `corrected_monotonic_lens_v1.jld2` — входная геометрия.

Полный 15-канальный transient намеренно отсутствует: carrier gate равен
`false`, следовательно его запуск нарушал бы заранее записанный порядок.
