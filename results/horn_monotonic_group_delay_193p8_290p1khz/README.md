# Monotonic horn lens — broadband delay correction v1

## Verdict

Для замороженной 15-канальной геометрии v0 рассчитаны независимые комплексные
отклики на `193.8`, `242.0` и `290.1 кГц`. Новые anchors потребовали 32
отдельных FEM-решения: восемь симметричных групп источников, lens/straight и
две дополнительные частоты. Использован calibration-уровень `phase2`:

```text
path h       0.483 мм
focus h      0.774 мм
receiver h   1.209 мм
```

Из восьми неотрицательных групп семь дают однозначный монотонный excess group
delay. Центральный канал не проходит phase-linearity gate: его три точки дают
residual `86.46 deg`, тогда как максимум остальных каналов равен `19.89 deg`.
Центральная оценка `-0.240 мкс` поэтому отвергнута, а не используется как
физическая отрицательная задержка.

Для валидных каналов измерена эффективная скорость группы `1147.4 м/с`.
Исходная линза пере-задерживает большую часть апертуры:

| y, мм | target, мкс | measured, мкс | path v0, мм | path v1, мм |
|---:|---:|---:|---:|---:|
| 0.0 | 5.777 | rejected | 9.724 | 6.957 |
| 4.8 | 5.637 | 7.983 | 9.488 | 6.796 |
| 9.6 | 5.224 | 7.794 | 8.794 | 5.846 |
| 14.4 | 4.560 | 6.980 | 7.676 | 4.900 |
| 19.2 | 3.674 | 5.271 | 6.184 | 4.352 |
| 24.0 | 2.598 | 3.799 | 4.373 | 2.995 |
| 28.8 | 1.364 | 2.060 | 2.296 | 1.497 |
| 33.6 | 0.000 | 0.000 | 0.000 | 0.000 |

Для центра path v1 продолжен от ближайшего валидного канала по той же
измеренной скорости группы. Это детерминированная монотонная экстраполяция,
а не подмена отвергнутой фазовой производной.

Первая рассчитанная коррекция:

```text
maximum extra path       9.724 -> 6.957 мм
maximum bend amplitude  13.746 -> 11.374 мм
common axial length     50.664 -> 50.664 мм
minimum inner radius     3.930 -> 4.917 мм
monotonic gate           passed
curvature gate           passed
```

Predicted delay RMS после correction равен нулю только внутри локальной
линейной модели `delta L = c_g (tau_target - tau_measured)`. Это не FEM-
верификация новой линзы и не заявка на фактический gain. Следующий gate —
пятицикловый transient максимального исправленного delay-state против прямого
equal-path контроля, затем одна 15-канальная v1 и matched straight reference —
[выполнен здесь](../horn_monotonic_v1_path_correction_242khz/README.md). Pulse
fidelity прошла, но полная v1 провалила carrier gate и была остановлена.

## Передача и ограничение полосы

Диагональное displacement-отношение `|lens/straight|` не проваливается в ноль:
минимум по всем каналам и anchors равен `0.511`, по семи валидным delay-
каналам — `0.599`. На несущей соответствующие минимумы равны `0.684` и
`0.714`. Центральная проблема — сильная нелинейность фазы по полосе, а не
нулевой carrier response. Перед физическим принятием central path полезны
более близкие к несущей anchors либо transient arrival-time check.

## Содержимое

- `horn_monotonic_group_delay.png` — delay, path, phase и transmission в
  физических единицах;
- `corrected_monotonic_lens_v1_geometry.png` — 3D-вид 15 исправленных
  `sin^4`-траекторий;
- `group_delay_path_correction.csv` — измеренные задержки и первая поправка;
- `frequency_anchor_transfer.csv` — wrapped/unwrapped phase всех anchors;
- `group_delay_summary.csv` — основные физические и безразмерные метрики;
- `corrected_monotonic_lens_v1_geometry.csv` — полная 15-канальная геометрия;
- `corrected_monotonic_lens_v1.jld2` — mesher-compatible v1 design;
- `response_matrix_*` — шесть комплексных `8x8` матриц;
- `basis_*_193p8khz.jld2`, `basis_*_290p1khz.jld2` — 32 независимых FEM-
  столбца новых частот.

Воспроизводимые драйверы:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_response_matrix.jl --stage=solve --variant=lens --source=1 --frequency-khz=193.8
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_response_matrix.jl --stage=aggregate --variant=lens --frequency-khz=193.8
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_group_delay.jl
```
