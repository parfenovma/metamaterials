# Контракт O3: paired pulse-echo и структурный clutter

> Зафиксировано 15 августа 2026 года до расчёта O3. Этот документ задаёт
> следующий stage gate после положительного O2c. Он не является заявлением
> POD и не разрешает синтез геометрии.
>
> Статус исполнения: reduced campaign завершён и не прошёл equal-time
> mechanical non-inferiority. Результат находится в
> [`paired_clutter_born_v0/`](../results/operator_space/paired_clutter_born_v0/README.md);
> исходные gates ниже не изменялись после test-run.

## 1. Решение, которое должен принять O3

O3 отвечает на один вопрос:

> Сохраняет ли выбранный малоканальный класс преимущество над A-scan и
> механическим сканированием, когда решение принимается по сырым
> defect/no-defect данным с ringing, backwall, contact variation и material
> clutter, а не по идеальному defect-minus-reference словарю?

В O3 входит один основной кандидат:

```text
four_channel_two_static_states
electronic channels     4
sequential states       2
state matrices          real, reciprocal, contractive
symmetry                2 even + 2 odd ports in each state
energy normalization    1/sqrt(2) across states
2D -> 3D status         group_sector + translation_kz
```

Точные матрицы состояний берутся из
`o2c_multistate_codebooks.jld2`; их нельзя дооптимизировать на test ensemble.
Второй прошедший O2c codebook хранится как frozen резерв, но не образует ещё
один архитектурный класс.

Комплексный true-time-delay класс не проходит O2c по calibration fidelity и в
O3 не продвигается. Односостоянийный real-класс сохраняется как отрицательный
control.

## 2. Обязательные baselines

Все операторы собираются из одной элементарной Tx/Rx-библиотеки:

1. `a_scan_uniform` — один неподвижный преобразователь;
2. `mechanical_scan_equal_time` — один канал и столько положений, сколько
   допускает тот же acquisition-time budget;
3. `four_channel_two_static_states` — основной кандидат;
4. `paut_fmc_oracle` — верхняя граница полной виртуальной апертуры.

Дополнительно сохраняется `four_channel_one_static_state`, чтобы измерить
цену и пользу второго состояния. Нельзя сравнивать кандидат с A-scan без
одновременного FMC oracle: иначе не измеряется доля реально закрытого разрыва.

## 3. Общая нормировка ресурсов

Для всех архитектур фиксируются:

```text
physical aperture
center frequency and usable bandwidth
five-cycle incident pulse
total incident energy per inspected point
common ADC bandwidth and electronics-noise PSD
maximum inspection time
calibration data budget
```

Четырёхпортовый кандидат получает два последовательных состояния, но не
двойную энергию: каждый state dictionary умножается на `1/sqrt(2)`.
Mechanical scan получает явный штраф по времени/позиционированию. FMC остаётся
oracle и не получает fictitious BOM parity.

В metadata раздельно хранятся:

```text
K = electronic channels
Q = transmit/state events
P = mechanical positions
T = normalized inspection time
C = calibration samples
```

## 4. Paired data model

Для одного background realization `eta`, гипотезы дефекта `xi` и амплитуды
рассеяния `alpha` генерируются связанные наблюдения

```math
y_0(\eta)=b(\eta)+r(\eta)+n_0,
```

```math
y_1(\xi,\alpha,\eta)
=b(\eta)+r(\eta)+\alpha a(\xi;\eta)+n_1.
```

Здесь `b` — backwall/material/contact background, `r` — front-end ringing,
`n_0,n_1` — независимый electronics noise. Общий `eta` делает пару
физически честной, но алгоритм на test не получает `y_0` для идеального
вычитания, если такая reference-процедура не включена в ресурсный контракт.

Training/calibration, threshold-validation и final-test используют разные
random seeds и разные background realizations. Codebook и порог замораживаются
до final-test.

## 5. Последовательные nuisance levels

| Level | Состав | Назначение |
|---|---|---|
| `C0` | общая electronics noise PSD | проверка decision pipeline |
| `C1` | `C0` + direct ringing + backwall variability | coherent nuisance |
| `C2` | `C1` + distributed elastic scatterers | structural clutter |
| `C3` | `C2` + O2c material/contact/calibration perturbations | model mismatch |

В каждом уровне сохраняются результаты предыдущего. Провал между соседними
уровнями локализует механизм потери качества.

Первый reduced/Born campaign использует nuisance atoms, сформированные в той
же полной виртуальной апертуре до сжатия архитектурами. Это обеспечивает общий
физический background. Независимо сгенерированный удобный шум отдельно для
каждой архитектуры запрещён.

Минимальный `C3` ensemble содержит:

```text
at least 100 material/contact realizations
at least 10^4 no-defect trials for Pfa = 1%
all 136 frozen ROI locations
at least 5 scattering-amplitude levels
isotropic and oriented P/S scattering tensors
```

## 6. Decision pipeline

Первый общий estimator — covariance-regularized matched-filter/GLRT bank на
замороженном defect dictionary. Для каждого trial он выдаёт:

```text
maximum detection score
estimated defect hypothesis
localization error
out-of-distribution / calibration residual
```

Порог выбирается только на no-defect validation set. Primary operating point:

```text
Pfa = 1%
```

Обязательно сохраняется верхняя односторонняя доверительная граница Pfa, а не
только наблюдаемая доля. Более сложный learned estimator разрешён лишь как
второй анализ с одинаковыми train/test данными для всех архитектур.

## 7. Primary metrics

Для каждого nuisance level и архитектуры сохраняются:

1. `Pd(alpha)` при фиксированном `Pfa=1%`;
2. `alpha90` — минимальная амплитуда для `Pd >= 0.90`;
3. `Pdl(alpha)` — вероятность обнаружить и локализовать не хуже
   `0.5 lambda_P0`;
4. median и Q90 ошибки локализации среди detected trials;
5. calibration false-reject rate;
6. inspection time и полный resource vector.

Доля закрытого разрыва до PAUT определяется для метрики, где меньше лучше:

```math
g_{PAUT}=
\frac{\alpha_{90,A}-\alpha_{90,candidate}}
{\alpha_{90,A}-\alpha_{90,FMC}}.
```

Если знаменатель статистически неотличим от нуля, эта метрика не заявляется.

## 8. Предварительно зафиксированные O3 gates

Основной кандидат проходит reduced O3 только при одновременном выполнении:

```text
test Pfa upper 95% bound                  <= 1.25%
Pd at preregistered reference alpha       >= 0.90
Pdl within 0.5 lambda_P0                  >= 0.80
Q10 PAUT-gap closure over C3 backgrounds  >= 0.50
Q90 localization error                    <= 1.0 lambda_P0
no full recalibration per contact
resource vector                           K=4, Q=2
```

Порог `reference alpha` выбирается до final-test как точка, где FMC достигает
`Pd=0.95` на validation ensemble. Так абсолютная шкала Born-амплитуды не
маскируется произвольной ручной настройкой.

Отдельно проверяется non-inferiority к `mechanical_scan_equal_time`. Если
кандидат выигрывает у A-scan, но проигрывает дешёвому механическому скану при
том же времени контроля, продуктовая ветка не считается выбранной.

## 9. Переход reduced model -> full-wave

Положительный reduced/Born O3 разрешает ровно один следующий дорогой шаг:

```text
paired 2D/2.5D transient full-elastic confirmation
    with one no-defect and one defect mesh family
    for the frozen entrant and all baselines.
```

Full-wave проверка обязана включать конечный дефект, mode conversion,
front/backwall echoes и одинаковый receive chain. До неё результат называется
`architectural surrogate`, а не POD.

Только успешный full-wave O3 открывает O4 и формирование целевого
многопортового `H(omega)` с допусками. Геометрия начинается после O4, а не
после красивой ROC-кривой reduced model.

## 10. Правила остановки

- Не менять codebook после просмотра final-test.
- Не вычитать defect-free reference, если его получение не оплачено временем
  и калибровкой.
- Не давать кандидату больше импульсной энергии из-за двух состояний.
- Не считать частотные bins независимыми observations.
- Не подбирать clutter отдельно под каждую архитектуру.
- Не переходить к геометрии при провале C2/C3; сначала пересмотреть состояние,
  ROI, task prior или calibration burden.
- Не называть reduced ROC промышленным POD.
