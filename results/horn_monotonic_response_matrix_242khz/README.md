# Monotonic horn lens — response matrix at 242 kHz

## Verdict

Восемь независимых симметричных групп источников рассчитаны для линзы и
matched straight reference. Сумма столбцов воспроизводит исходные равномерные
поля с относительной ошибкой `2.7e-15` для линзы и `1.4e-14` для контроля.

Глобальный неотрицательный амплитудный optimum при фиксированном номинальном
энергетическом бюджете

```text
sum(m_i s_i^2) = 15,  m = [1, 2, 2, 2, 2, 2, 2, 2]
```

поднимает фокус линзы на response-matrix сетке на `4.45%`, сохраняет FWHM
`4 мм` и даёт sidelobe/peak `0.260`. Но matched gain падает
`2.271 -> 1.761`, потому что те же веса ещё сильнее увеличивают прямой
контроль.

Независимая проверка одного рассчитанного набора весов на финальной сетке
подтвердила направление и величину прогноза:

| Метрика | Equal drive | Weighted | Изменение |
|---|---:|---:|---:|
| Фокус линзы, нм | 2.226 | 2.326 | +4.48% |
| Фокус straight reference, нм | 1.053 | 1.409 | +33.84% |
| Matched gain | 2.114 | 1.650 | -21.93% |
| FWHM линзы, мм | 4 | 4 | 0 |

Proxy gate `>=15%` улучшения matched gain и `gain>=2.2` не пройден.
Переносить эти веса в физические horn-mouth areas нельзя. Следующая ветка —
полосовые anchors `193.8 / 242 / 290.1 кГц`, group-delay error и одна
геометрическая коррекция невёрнутой длины пути — завершена в
[`horn_monotonic_group_delay_193p8_290p1khz`](../horn_monotonic_group_delay_193p8_290p1khz/README.md).

## Сетки и ограничение памяти

Response matrices рассчитаны на неравномерном уровне `phase2`:

```text
path h       0.483 mm
focus h      0.774 mm
receiver h   1.209 mm
```

Он используется только для направления линейной коррекции. Финальная
проверка выполнена на `phase4`, где path `h=0.338 мм`; её результат является
приёмочным.

LU-факторизация `phase3/phase4` почти полностью занимает доступную память
macOS, поэтому повторные RHS в одном процессе вызывают системный allocator
crash. Каждый столбец `phase2` решён в отдельном свежем Julia-процессе.
Точный boundary-work integral после LU по той же причине не вычислялся;
текущий proxy сохраняет номинальный pressure-energy budget через
`sum(m_i s_i^2)`. Комплексная введённая мощность остаётся отдельным
boundary-аудитом и явно отмечена как unavailable в summary CSV.

## Содержимое

- `response_matrix.csv` — амплитуды в нанометрах и фазы всех `8x8` откликов;
- `source_weights.csv` — рассчитанные множители давления;
- `precompensation_proxy_summary.csv` — response-matrix prediction;
- `phase4_verification_summary.csv` — независимая финальная проверка;
- `response_matrix_*.jld2` — комплексные матрицы, фокальные отклики и профили;
- `basis_*.jld2` — шестнадцать независимых столбцов;
- `weighted_*_phase4.jld2` — финальные weighted поля линзы и контроля;
- `horn_monotonic_response_matrix.png` — график в физических единицах.

Воспроизводимый драйвер:

```bash
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_response_matrix.jl --stage=solve --variant=lens --source=1
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_response_matrix.jl --stage=aggregate --variant=lens
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_response_matrix.jl --stage=analyze
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_response_matrix.jl --stage=verify --variant=lens
julia --startup-file=no --project=. src/code/metamaterial/run_horn_monotonic_response_matrix.jl --stage=verify --variant=uniform
```
