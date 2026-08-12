# Monotonic horn-waveguide lens — 3D convergence at 242 kHz

## Verdict

Связка **воронка → плавный монотонный волновод → малый излучатель** прошла
phase-path refinement в symmetry-constrained half-domain. На двух последних
уровнях сетка тракта уменьшена с `0.387` до `0.338 мм`; matched target gain
изменился только с `2.135` до `2.114` (`0.994%`), а FWHM остался `4 мм`.
Предварительно заданный gate `5% / 1 scan step` пройден.

Рабочий carrier-результат v0 при `242 кГц`: `2.226 нм` у линзы против
`1.053 нм` у прямого согласованного контроля, то есть gain `2.114×`.
Это сертифицирует относительное усиление для следующей итерации. Абсолютная
амплитуда линзы между двумя последними уровнями меняется сильнее (`7.59%`),
поэтому абсолютную калибровку в нанометрах всё ещё следует считать отдельной
задачей receiver/boundary refinement.

## Единицы и нормализация

Модель теперь задаёт сетку внутри через `h/λS` и `kS h`, а геометрические
масштабы — через `L/λP` или `L/λS`. При `242 кГц` для текущего фотополимера
`λP=9.669 мм`, `λS=4.835 мм`.

Пользовательские результаты намеренно остаются физическими: координаты и
FWHM в миллиметрах, амплитуды в нанометрах, частоты в килогерцах. В CSV
сначала идут физические колонки, после них — normalized diagnostics. Полный
перевод размеров приведён в `physical_and_dimensionless_scales.csv`.

## Full-domain screening

| Уровень | Path mesh, мм | Lens, нм | Straight, нм | Target gain | Lens FWHM, мм | Max mirror mismatch |
|---|---:|---:|---:|---:|---:|---:|
| coarse | 1.10 | 4.858 | 1.795 | 2.706 | 3 | 17.37% / 6.25° |
| reference | 0.90 | 3.893 | 2.645 | 1.472 | 5 | 11.36% / 6.13° |
| refined | 0.80 | 3.928 | 2.466 | 1.593 | 6 | 12.81% / 3.79° |

Переход `0.90 → 0.80 мм` меняет target gain на `7.60%`; разрешено не более
`5%`. Full-domain тетраэдры также не проходят зеркальный gate `3% / 3°`.

## Symmetry-constrained half-domain

Половина `y ≥ 0` решена с точным условием `u_y=0` на центральной плоскости.
Это исключает случайную лево-правую асимметрию сетки.

| Уровень | Path mesh, мм | Lens, нм | Straight, нм | Target gain | Plane-peak gain | Window-amplitude gain | Lens FWHM, мм |
|---|---:|---:|---:|---:|---:|---:|---:|
| reference | 0.900 | 3.456 | 2.130 | 1.623 | 1.322 | 1.475 | 8 |
| refined | 0.800 | 3.444 | 2.382 | 1.446 | 1.624 | 1.733 | 8 |
| fine | 0.700 | 4.220 | 2.268 | 1.861 | 1.662 | 1.422 | 4 |
| phase1 | 0.580 | 2.518 | 1.646 | 1.530 | 1.541 | 1.548 | 6 |
| phase2 | 0.483 | 2.794 | 1.230 | 2.271 | 1.941 | 1.796 | 4 |
| phase3 | 0.387 | 2.409 | 1.128 | 2.135 | 1.912 | 1.692 | 4 |
| phase4 | 0.338 | 2.226 | 1.053 | 2.114 | 1.907 | 1.653 | 4 |

`Window-amplitude gain` — квадратный корень из отношения энергий `|ux|²` в
окне `|y|≤3 мм` целевой плоскости. Уровни `phase3` и `phase4` меняют только
сетку фазочувствительного тракта; фокусная и остальная приёмная зоны остаются
на `0.774` и `1.209 мм`. Переход `0.387 → 0.338 мм` меняет target gain на
`0.994%`, plane-peak gain на `0.283%`, window-amplitude gain на `2.30%` и не
меняет FWHM. Half-domain gate пройден.

## Содержимое

- `full_domain_convergence.*`, `full_domain_gate.csv`,
  `full_domain_mirror_pairs.csv` — full-domain screening;
- `half_domain_convergence.*`, `half_domain_gate.csv` — симметричная проверка;
- `half_domain_mesh_statistics.csv` — физические размеры и фактическое число
  узлов/элементов новых сеток;
- `physical_and_dimensionless_scales.csv` — физические значения и их
  внутренние нормированные аналоги;
- `full_*.jld2`, `half_*.jld2` — компактные фокальные профили и выходы
  каналов для всех решённых уровней, включая `phase1`–`phase4`.

Полные `.msh` не дублируются. Они остаются в `tmp/` и воспроизводятся
драйверами:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d_convergence.jl --stage=mesh --level=all --variant=all
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d_convergence.jl --stage=solve --level=all --variant=all
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d_convergence.jl --stage=analyze

julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d_half_symmetry.jl --stage=mesh --level=all --variant=all
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d_half_symmetry.jl --stage=solve --level=phase4 --variant=lens
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d_half_symmetry.jl --stage=solve --level=phase4 --variant=uniform
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_array_3d_half_symmetry.jl --stage=analyze
```

Тяжёлые lens/reference solve запускаются раздельно, чтобы прямой решатель не
сохранял фрагментированную память между двумя миллионно-элементными задачами.
