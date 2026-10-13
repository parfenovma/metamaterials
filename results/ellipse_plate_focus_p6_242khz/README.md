# Первый FEM-прототип монолитной пластины с эллиптическими пустотами

## Краткий результат

Проверен изготовимый аналитический инициализатор `9 поперечных зон × 3
продольных ряда` сквозных эллиптических пустот. Полный 2D plane-strain FEM
показывает, что эта конкретная трёхрядная топология **не является пригодной
фокусирующей линзой** на `242 kHz`.

- малые пустоты (`scale=0.25`) сохраняют амплитуду в номинальной точке
  (`J1=1.037`), но дают широкое поле `FWHM=36.0 mm`, максимум на `x=93 mm`
  вместо `115 mm` и не формируют локализованный фокус;
- исходный TTD-инициализатор (`scale=1.00`) создаёт более узкий осевой лепесток
  `FWHM=23.8 mm`, но на `x=137 mm`, а амплитуда в номинальном фокусе падает до
  `J1=0.213`; proxy выходной продольной энергии равен только `0.255` от
  сплошного straight reference;
- результаты `scale=0.25` по `J1` и FWHM хорошо сходятся между двумя сетками;
  для неоднородного поля `scale=1.00` положение лепестка и energy proxy
  сходятся, но FWHM меняется с `17.2` до `23.8 mm`. Это ещё один повод не
  интерпретировать узкий coarse-лепесток как доказанный результат.

Verdict: **stop для этой трёхрядной effective-medium топологии**. Она полезна
как диагностический no-go и стартовая сетка, но не как «оптимальная линза».

## Постановка

| Параметр | Значение |
|---|---:|
| Модель | 2D plane strain, полноупругая `u_x,u_y` |
| Материал | Al-6061 surrogate, однородный и без потерь |
| `rho` | `2700 kg/m^3` |
| `cP` | `6122.1024 m/s` |
| `cS` | `3083.8103 m/s` |
| Частота | `242 kHz` |
| Апертура | `71.6 mm` |
| Длина пластины | `55 mm` |
| Длина образца после пластины | `125 mm` |
| Номинальный фокус | `60 mm` от выходной грани, `x=115 mm` |
| Источник | равномерная нормальная тяга `1 MPa` по входной грани |
| Reference | сплошной Al-блок тех же внешних размеров |
| Боковые/выходная границы образца | first-order absorbing impedance |
| Поверхности пустот | traction-free |

Piezo, контактный слой и вода/другая связующая среда не моделировались.
`J1` здесь — отношение гармонических `|u_x|` при одинаковой заданной входной
тяге, а не transient peak пятициклового импульса.

## Как построена геометрия

Для центров девяти зон вычислена требуемая невёрнутая задержка

```math
\Delta t(y)=
\frac{\sqrt{F^2+(D/2)^2}-\sqrt{F^2+y^2}}{c_P}.
```

Она отображена в начальную долю пустот ограниченной нерезонансной формулой

```math
\frac{c_{eff}}{c_P}=\sqrt{\frac{1-\phi}{1+\phi}}.
```

Половина площади пустот помещена в центральный продольный ряд, по четверти —
во входной и выходной согласующие ряды. Полуоси имеют отношение
`b_y/a_x=1.10`. Это формула **только для инициализации**; все численные метрики
получены из полноупругой FEM, а не из effective-medium модели.

Для `scale=1.00`:

- 27 сквозных эллипсов;
- пористость зон от `3.49%` до `16.36%`;
- минимальная полуось `1.051 mm`;
- максимальная полуось `3.540 mm`;
- минимальный зазор между эллипсами `0.963 mm`;
- наружная перемычка не менее `0.8 mm`.

Для диагностического `scale=0.25` минимальная полуось `0.525 mm`, а
минимальный зазор `4.459 mm`. Обе геометрии проходят заданные manufacturing
ограничения `semi-axis >= 0.4 mm`, `gap/ligament >= 0.8 mm` и не содержат
закрытых полостей.

## Сетка и сходимость

Сетка конформна эллиптическим свободным поверхностям. Gmsh локально сгущает её
до `0.42/0.28 mm` у пустот на coarse/fine уровнях. В выходном образце размеры
элементов равны примерно `1.45/1.00 mm`, то есть около `lambdaP/17` и
`lambdaP/25`.

| Вариант | Сетка | Узлы линзы | Треугольники линзы | `J1` nominal | `x_peak`, mm | FWHM, mm | exit-energy proxy |
|---|---|---:|---:|---:|---:|---:|---:|
| `scale=0.25` | coarse | 28 865 | 56 655 | 1.0347 | 93 | 35.41 | 0.6117 |
| `scale=0.25` | fine | 59 719 | 117 867 | 1.0369 | 93 | 36.00 | 0.5876 |
| `scale=1.00` | coarse | 28 918 | 56 371 | 0.2436 | 137 | 17.16 | 0.2534 |
| `scale=1.00` | fine | 60 043 | 117 912 | 0.2128 | 137 | 23.78 | 0.2547 |

`exit-energy proxy = sum_y |u_x(56 mm,y)|^2` относительно straight reference.
Он не является полным `S21`: в нём нет явной проекции на P-моду, потока
мощности и вклада `u_y`. `DOF_-3dB` не замкнулся внутри расчётного окна и
оставлен `NaN`. Sidelobe ratio также не считается сходимой метрикой для этой
неудачной многолепестковой конфигурации.

Полная таблица: [`data/focusing_metrics.csv`](data/focusing_metrics.csv).

## Главные изображения

- [`figures/ellipse_plate_topology.png`](figures/ellipse_plate_topology.png) —
  параметризованная геометрия `scale=1.00`;
- [`figures/mesh_full_fine.png`](figures/mesh_full_fine.png) — полная fine-сетка;
- [`figures/mesh_ellipses_closeup_fine.png`](figures/mesh_ellipses_closeup_fine.png) —
  локальное сгущение у эллипсов;
- [`figures/focusing_field_normalized_scale_1.00_fine.png`](figures/focusing_field_normalized_scale_1.00_fine.png) —
  форма поля сильной пластины, нормированная на собственный максимум;
- [`figures/focusing_field_map_scale_1.00_fine.png`](figures/focusing_field_map_scale_1.00_fine.png) —
  то же поле в абсолютной нормировке на straight nominal focus;
- [`figures/transverse_profile_scale_1.00_fine.png`](figures/transverse_profile_scale_1.00_fine.png) и
  [`figures/axial_profile_scale_1.00_fine.png`](figures/axial_profile_scale_1.00_fine.png) — профили;
- соответствующие файлы `scale_0.25` показывают слабо перфорированную ветку.

## Воспроизведение

Расчёт выполнен на commit `775244a2424508d890fea7012d50fa655712258f`.
Из корня репозитория:

```bash
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case straight coarse 0.25
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case lens coarse 0.25
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case straight coarse 0.50
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case lens coarse 0.50
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case straight coarse 0.75
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case lens coarse 0.75
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case straight coarse 1.00
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case lens coarse 1.00
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case straight fine 0.25
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case lens fine 0.25
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case straight fine 1.00
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --case lens fine 1.00
julia --startup-file=no --project=. results/ellipse_plate_focus_p6_242khz/ellipse_plate_pipeline.jl --report
```

Скрипт не перезаписывает существующие JLD2 без явного вызова Julia API с
`force=true`; сетки при повторном расчёте строятся детерминированно.

## Что не рассчитано

- transient `5-cycle Hann`, pulse fidelity и postcursor;
- `S11/S21` с модовой проекцией;
- P–S mode conversion;
- дефект, pulse-echo, CNR и ringing;
- 3D конечная толщина, piezo, backing, bond layer;
- manufacturing perturbations.

Эти стадии сознательно не запускались: кандидат не прошёл дешёвый harmonic G2
gate. Следующая содержательная ветка должна сначала получить локальные
комплексные `S21/group delay` для одной-двух эллиптических ячеек и использовать
больше продольных рядов или геометрическую TTD-трассу. Простое масштабирование
трёх рядов уже показало монотонный компромисс «передача против фазового
воздействия» и не заслуживает дальнейшего brute-force sweep.
