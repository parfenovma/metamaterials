# Многопараметрический screening эллиптической пластины, 242 кГц

## Задача

Проверить, можно ли развязать две цели конструкции:

- `J1` — harmonic gain амплитуды `|u_x|` в номинальном фокусе относительно
  монолитного straight footprint при одинаковой интегральной энергии возбуждения;
- `J2` — reduced defect-CNR gain относительно того же straight reference.

Рассчитано 59 полноупругих 2D plane-strain coarse FEM-конфигураций и три
выбранных кандидата на fine-сетке.

## Степени свободы

| Параметр | Смысл | Диапазон screening |
|---|---|---:|
| `q` | `a_x/a_y` наружных согласующих эллипсов | 1.15–3.50 |
| `s_c` | масштаб площади отверстий центрального ряда | 0.08–0.35 |
| `q_c` | форма центральных эллипсов при постоянной площади | 0.50–2.00 |
| `d/L` | расстояние наружного ряда от центра / длина пластины | 0.16–0.42 |
| `gamma` | показатель поперечного delay/porosity-профиля при сохранении суммарной площади | 0.55–2.00 |
| `alpha` | нормированная амплитудная аподизация возбуждения | 0–0.8 |

Первые 13 точек — контролируемые однофакторные гипотезы, 19 — Halton screening,
27 — регулярный refinement `q × s_c × d/L`. Во всех точках refinement сохранены
все 27 отверстий и manufacturing gate `0.8 mm`.

## Определение reduced J2

Для моностатического point-scatterer surrogate использована взаимность. Двухходовой
Born-kernel определён как

```math
K(x,y)=u_x^2(x,y).
```

В окне вокруг ожидаемой глубины дефекта

```math
\mathrm{CNR}_{red}=
\frac{|K(x_d,y_d)|}
{\sqrt{\sum_{\Omega_c}|K(x,y)|^2+\sigma_n^2}},
\qquad
J_2=\frac{\mathrm{CNR}_{red,lens}}
{\mathrm{CNR}_{red,straight}}.
```

Default-окно: `x_d ± 15 mm`, `|y| ≤ 30 mm`, исключение радиуса `4 mm` вокруг
дефекта. Default noise variance равна clutter variance straight reference.
Проверены 18 сочетаний размера окна и noise/clutter ratio.

Это физически мотивированный harmonic screening CNR, но не окончательный
transient CNR: модель пока не содержит пьезоэлектрический ringing, контактный
слой, реальные зерновые реализации clutter и электронику.

## Screening effects

Стандартизованные линейные эффекты глобального DOE:

| Параметр | `J1` | `J2` | Вывод |
|---|---:|---:|---|
| `d/L` | +0.582 | +0.586 | самый сильный параметр; выраженный интерференционный максимум |
| `alpha` | -0.321 | -0.357 | taper увеличивает exit-energy, но ухудшает обе цели |
| `s_c` | -0.264 | -0.250 | в среднем ухудшает амплитуду; создаёт отдельную clutter-ветвь |
| `q` | +0.150 | +0.117 | слабый положительный эффект только совместно с правильным offset |
| `gamma` | -0.069 | -0.052 | слабый эффект |
| `q_c` | +0.014 | +0.024 | практически нейтрален в этом диапазоне |

## Регулярная heatmap и fine-проверка

Coarse-refinement формально дал максимум
`r025 = (q=1.8, s_c=0.12, d/L=0.38)`:

```text
J1 = 1.01087
J2 = 1.00302
```

Однако fine-сетка сняла этот небольшой выигрыш:

| Кандидат, fine | `J1` | `J2`, equal noise/clutter | Exit-energy | FWHM |
|---|---:|---:|---:|---:|
| `r025` | 0.99991 | 0.98295 | 0.74137 | 37.54 mm |
| `r016` | 0.99871 | 0.98448 | 0.71873 | 37.30 mm |

Следовательно, coarse-выигрыш порядка 1% был сеточной погрешностью. Геометрия
практически восстанавливает амплитуду straight plate, но default CNR остаётся
ниже на `0.14–0.15 dB`.

При clutter-dominated режиме появляется реальная, хотя маленькая развязка.
Кандидат `c005` с увеличенным центральным рядом (`s_c=0.20`) на fine-сетке даёт:

```text
J1 = 0.93931
J2_clutter = 1.00825  (+0.071 dB)
```

То есть около 6% амплитуды обмениваются лишь на 0.8% CNR. Для практической
системы это слишком слабый Pareto-компромисс, но он подтверждает сам механизм:
усиленный центральный ряд снижает интегральный clutter kernel сильнее, чем
дефектный отклик.

## Verdict

`stop as standalone lens / retain as matching corrector`.

Проверенные геометрические и пьезоэлектрические степени свободы не дают
значимого одновременного роста `J1` и `J2`. Следующая топология должна вводить
настоящую независимую фазовую координату — TTD/GRIN-каналы — и отдельную
приёмную координату: сегментированный электрод/beamforming. Эллиптические ряды
можно сохранить только как согласующие элементы около входа и выхода.

## Артефакты

- `data/doe_manifest.csv` — все 59 геометрий;
- `data/multidof_metrics.csv` — coarse-метрики;
- `data/fine_validation_metrics.csv` — три fine-проверки;
- `data/screening_effects.csv` — стандартизованные эффекты;
- `data/cnr_sensitivity.csv` — устойчивость к noise/clutter и окну;
- `figures/regular_q_sc_offset_heatmap.png` — регулярная heatmap `J1/J2`;
- `figures/regular_decoupling_heatmap.png` — `J2/J1^2`;
- `figures/J1_J2_pareto.png` — coarse Pareto plane;
- `ellipse_multidof_pipeline.jl` — воспроизводимый pipeline.

## Воспроизведение

```bash
julia --startup-file=no --project=. \
  results/ellipse_multidof_p6_242khz/ellipse_multidof_pipeline.jl --manifest

julia --startup-file=no --project=. \
  results/ellipse_multidof_p6_242khz/ellipse_multidof_pipeline.jl --case c001

julia --startup-file=no --project=. \
  results/ellipse_multidof_p6_242khz/ellipse_multidof_pipeline.jl --report
```
