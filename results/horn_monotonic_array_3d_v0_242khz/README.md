# Horn + monotonic waveguide lens v0 — 242 kHz

> **Обновление convergence:** `1.472×` ниже воспроизводится как результат
> исходной v0-сетки, но больше не считается сеточно сертифицированным gain.
> Full/half-domain refinement сохранил положительную фокусировку, однако не
> прошёл количественный gain/FWHM gate. См. [пакет convergence](../horn_monotonic_array_3d_convergence_242khz/README.md).

## Verdict

Связка **плавная воронка → узкий плавно изогнутый волновод → малый
излучатель** принята как рабочая базовая топология линзы. Она прошла три
последовательных уровня проверки:

1. Максимальный delay-state в полноупругом transient сохранил
   `0.9775` амплитуды и `0.9620` энергии относительно прямого волновода той же
   развёрнутой длины; `B_t=0.9906`, `rho=0.99983`, поперечная энергия
   `0.00550`.
2. Редуцированная 15-элементная импульсная апертура дала `G_peak=2.904`,
   FWHM `5 мм` и postcursor `0.056`.
3. Полная 15-канальная 3D-модель v0 дала в целевой точке `|ux|=3.893 нм`
   против `2.645 нм` у согласованной прямой решётки: baseline-выигрыш
   **`1.472×` (`+47.2%`)**. Последующая проверка показала, что его точное
   значение ещё зависит от сетки.

`5.280×` в отчёте полной линзы — это отношение амплитуд выходов центрального
и внешнего каналов, а не gain относительно исходного импульса или прямой
решётки.

## Постановка v0

- материал каналов и общего принимающего блока: фотополимер;
- частота: `242 кГц`;
- 15 одинаковых синфазных источников по `1 МПа`;
- шаг каналов по апертуре: `4.8 мм`;
- log-cosine воронка: `3.2 → 1.6 мм`;
- монотонный out-of-plane `sin^4`-изгиб без складывания пути;
- extra path монотонно растёт от `0` на краях до `9.724 мм` в центре;
- фокусная плоскость: `35 мм` после плоскости излучателей;
- сетка линзы: `59 751` узел, `332 095` элементов.

На целевой плоскости FWHM равен `5 мм`, sampled waist находится около
`x=30 мм` с FWHM `3 мм`, боковой лепесток равен `0.372` от главного максимума,
а максимальная ошибка выходной фазы относительно редуцированной целевой
модели — `12.72°`.

## Содержимое

- `full_array_focus.png`, `full_array_summary.csv` — поле и метрики полной
  монотонной линзы;
- `matched_full_3d_gain.png`, `matched_gain_summary.csv` — прямое сравнение с
  согласованной uniform-решёткой;
- `channel_outputs.csv` — комплексные выходы всех 15 каналов;
- `lens_carrier_field.jld2`, `uniform_carrier_field.jld2` — компактные поля
  обеих 3D-моделей на несущей;
- `mesh_summary_lens.csv`, `mesh_summary_uniform.csv` — размеры сеток;
- `monotonic_lens_geometry.*` — закон задержек и геометрия каналов;
- `maximum_delay_*` — transient максимального delay-state;
- `reduced_aperture_impulse*` — полосовая импульсная модель апертуры;
- `receiver_fem_*` — предыдущая редуцированная полноупругая модель приёмника.

Полные `.msh` не продублированы в этом пакете: они занимают около `26 МБ` и
детерминированно пересобираются исходным скриптом. Рабочие копии остаются в
`tmp/horn_monotonic_array_3d/`.

## Воспроизведение полного 3D-сравнения

Из корня проекта после построения монотонного aperture-design:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d.jl --stage=mesh --variant=lens
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d.jl --stage=solve --variant=lens
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d.jl --stage=report --variant=lens
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d.jl --stage=mesh --variant=uniform
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d.jl --stage=solve --variant=uniform
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d.jl --stage=compare --variant=lens
```

Основной код: [`run_horn_monotonic_array_3d.jl`](../../src/code/metamaterial/run_horn_monotonic_array_3d.jl),
[`horn_monotonic_array_3d_mesher.jl`](../../src/code/metamaterial/horn_monotonic_array_3d_mesher.jl),
[`horn_monotonic_array_3d_harmonic_solver.jl`](../../src/code/metamaterial/horn_monotonic_array_3d_harmonic_solver.jl).

## Ограничения результата

Это carrier-level v0 на сравнительно грубой 3D-сетке. Перед печатью ещё нужны
полный transient линза/reference, mesh convergence, потери и разброс свойств,
контактный/согласующий слой и перенос фокусировки в алюминий. Эти проверки
должны улучшать или опровергать текущую реализацию, но поиск базовой топологии
на данном этапе закрыт: следующая итерация развивает воронку с волноводом.
