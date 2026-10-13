# План исследования малоканального упругого дефектоскопа

> Обновлено 15 августа 2026 года. Это основной управляющий документ проекта.
> Полное пространство архитектур и определения операторных метрик находятся в
> [`measurement_operator_design_space.md`](measurement_operator_design_space.md),
> безразмерный контракт первого pilot — в
> [`operator_task_contract.md`](operator_task_contract.md), история физических
> расчётов — в [`research_log.md`](research_log.md), текущий paired benchmark —
> в [`operator_o3_benchmark_contract.md`](operator_o3_benchmark_contract.md).
> Новые геометрии не синтезируются до gates `O0--O3` и явного handoff
> `O4 -> O5`.

# Часть I. Основная амбиция проекта

## 1. Продуктовая и физическая цель

Разработать ультразвуковой дефектоскоп, который в ограниченной области
контроля сопоставим с PAUT по обнаружению и локализации заданного класса
дефектов, но использует существенно меньше электронных каналов и более
дешёвую архитектуру:

```text
1--4 электронных канала
    + пассивное или малоканальное упругое кодирование
    + калиброванный алгоритм решения обратной задачи
    -> полезная доля качества PAUT внутри ограниченной ROI.
```

Проект не ставит целью построить ещё одну фазированную решётку. Дорогие
независимые element-wise phase channels должны быть заменены пространственным,
временным, частотным, модовым или механическим кодированием. Цена системы
включает не только ADC/pulser, но также PZT, изготовление, калибровку, контакт,
время измерения и повторяемость посадки.

Линза, метаматериал, пассивная сеть, несколько symmetry-портов или
механически перестраиваемая поверхность являются возможными реализациями
измерительного оператора, но ни одна геометрия не выбирается заранее.

## 2. Центральная операторная постановка

Для дефекта `xi`, архитектуры `D`, возбуждения `u_k` и приёмного функционала
`r_l` измерение записывается как

```math
y_{kl}(t;\xi,D)=
\mathcal R_{D,l}\,\mathcal G(\xi)\,\mathcal T_{D,k}u_k(t)
+b_{D,kl}(t)+n_{kl}(t).
```

После одинакового feature extraction:

```math
z=A_Dx+\varepsilon.
```

Столбец `a_j` матрицы `A_D` — весь измеряемый отклик гипотезы дефекта
`xi_j`. Проектируется тройка

```math
(D,\{u_k\},\mathcal A_{decision}),
```

а не только механическая структура. Норма столбца характеризует
обнаружимость, сходство столбцов — неоднозначность локализации, singular
spectrum — число устойчиво наблюдаемых степеней свободы.

Фокус, FWHM и полевая амплитуда являются диагностикой механизма. Они не
ранжируют дефектоскопы без paired pulse-echo, noise/clutter и общей
нормировки ресурсов.

## 3. Базовые инварианты

1. Сравнение ведётся одновременно с `A-scan baseline`, механически
   сканируемым одиночным преобразователем и `PAUT/FMC oracle` при общей
   апертуре, полосе и падающей энергии.
2. Основные расчёты живут в 2D/2.5D. Каждый результат обязан иметь статус
   `translation_kz`, `axisymmetric_m`, `group_sector`,
   `plane_strain_surrogate` или `full_3d_required`.
3. P- и S-волны, их разные скорости, затухание и mode conversion входят в
   оператор; S-составляющая не объявляется заранее только потерей или только
   полезным каналом.
4. Денежная цена не назначается до построения физического Pareto-фронта.
   Ресурсный вектор хранит каналы, transmit events, положения, состояния,
   fabrication и calibration burden раздельно.
5. Положительный Born/reduced результат разрешает следующий операторный gate,
   но не fabrication и не продуктовый POD claim.

## 4. Stage gates проекта

| Gate | Вопрос | Текущий статус |
|---|---|---|
| `O0 — task contract` | Что обнаруживаем, где и относительно каких baselines? | **pass** |
| `O1 — information ceiling` | Разрешима ли задача полной виртуальной апертурой? | **pass** |
| `O2 — compression frontier` | Сколько кодов/портов сохраняют полезную информацию? | **pass** |
| `O2c — causal robustness` | Сохраняется ли результат после causal `H(omega)`, pulse covariance и perturbations? | **pass: 4 channels x 2 states** |
| `O3 — architecture frontier` | Выигрывает ли класс в paired pulse-echo и clutter ensemble? | **pass: expanded F1d, smooth K3,Q1** |
| `O4 — transfer target` | Какой многопортовый `H(omega)` и с какими допусками нужно изготовить? | **complete: DCT3 modal H(omega)** |
| `O5 — geometry handoff` | Какая физическая топология реализует выбранный оператор? | **complete: 3-port monolithic mode transformer** |
| `O6 — physical twin` | Может ли один физический механизм реализовать O4 без скрытой решётки? | **O6j: conditional compact-P12 solid-contact coupon; q8 matching blocked** |
| `O7 — stationary global states` | Можно ли заменить запрещённое движение четырьмя глобальными состояниями одной головки? | **pass: PZT-2bit и mechanical-2bit на coupled-modal gate** |
| `O8 — plane-strain boundary` | Выживают ли frozen O7-состояния в конечном упругом блоке с физическими границами и нагрузками? | **conditional pass: PZT; mechanical stopped** |
| `O8M — manufacturing freeze` | Сохраняется ли механизм у закупаемой PZT и округлённой геометрии? | **pass: BUILD one CTS-3195HD coupon** |

## 5. Доказанное состояние O1/O2

Первый geometry-free pilot использует `lambda_P0=c_P=f0=1`,
`c_S/c_P=0.5`, апертуру `8 lambda_P0`, ROI
`x/lambda_P0=2.5--6.0`, `y/lambda_P0=-2.0--2.0` и семь частотных anchors.
Полный отчёт и словари находятся в
[`results/operator_space/dimensionless_elastic_born_v0/`](../results/operator_space/dimensionless_elastic_born_v0/README.md).

Физические positive controls пройдены:

```text
reciprocity error             9.54e-17
mirror equivariance error     0.0
FMC numerical rank            136
FMC stable rank               29.804
```

Основные результаты:

| Оператор | Каналы | Min detectability | Nonlocal coherence | Stable rank |
|---|---:|---:|---:|---:|
| uniform A-scan | 1 | 0.489 | 1.000 | 3.070 |
| fixed reciprocal focus | 1 | 0.069 | 1.000 | 2.713 |
| asymmetric dispersive port | 1 | 0.261 | 0.975 | 3.423 |
| cosine two-port | 2 | 0.397 | 0.842 | 5.462 |
| cosine four-port | 4 | 0.820 | 0.732 | 7.140 |
| FMC oracle | 17 | 1.123 | 0.159 | 29.804 |

Fixed focus повышает detectability в своей точке в `8.79x` относительно
uniform A-scan, но имеет точную неоднозначность `y <-> -y` и плохой минимум
по ROI. Это сильный `known-zone detection` оператор, а не imaging architecture.

P/S-ablation подтверждает, что поперечная мода несёт информацию: stable rank
FMC растёт с `25.130` для P-only до `29.804` для P+S. Для двух портов рост
составляет `4.989 -> 5.462`, для четырёх — `6.674 -> 7.140`.

Идеальная SVD-компрессия сохраняет `44.5%`, `72.9%` и `94.2%` энергии
словаря при `16`, `32` и `64` абстрактных измерениях. Это aggregate-information
ceiling, а не требование такого же числа электронных каналов.

## 6. Текущее архитектурное решение O2b

На замороженной библиотеке проверены reciprocal real/complex codebooks с
сохранением и нарушением reflection sectors. По `300` детерминированных
кандидатов на семейство дали:

| Порты | Codebook | Min detectability | Nonlocal coherence | Stable rank |
|---:|---|---:|---:|---:|
| 2 | cosine baseline | 0.397 | 0.842 | 5.462 |
| 2 | best complex sector | 0.280 | 0.615 | 8.566 |
| 4 | cosine baseline | **0.820** | 0.732 | 7.140 |
| 4 | best real sector | 0.587 | 0.362 | **13.337** |
| 4 | best complex sector | 0.563 | **0.343** | 12.259 |

Текущий лидер `L2` — четырёхпортовый codebook, сохраняющий even/odd сектора
группы отражения. Он примерно вдвое уменьшает дальнюю неоднозначность и почти
вдвое увеличивает stable rank относительно простых cosine-портов.

Два порта остаются условной дешёвой веткой. Ни один из `1200` испытанных
двухпортовых вариантов не сохранил worst-case detectability A-scan; это не
математический no-go, но достаточная причина не выбирать два порта по
умолчанию.

## 7. Результат O2c

Цель — проверить, остаётся ли четырёхпортовый результат после снятия трёх
оптимистических допущений `O2b`.

### 7.1. Causal transfer

Сравниваются:

1. `carrier_oracle` — замороженная матрица портов;
2. `static_real_sector` — real amplitudes/signs во всей полосе;
3. `delay_lift` — carrier phases заменены положительными true-time delays;
4. `contractive_delay_lift` — frequency-dependent матрица дополнительно
   масштабируется так, чтобы `sigma_max(W(omega)) <= 1`.

Комплексная carrier-матрица является только верхней границей. В `O4` может
передаваться лишь causal и contractive target.

### 7.2. Time-bandwidth contract

Семь harmonic anchors не считаются семью независимыми приборами. Проверяются:

```text
1 central anchor
3 band anchors
7 oversampled anchors
five-cycle pulse weighting
correlated frequency-feature covariance.
```

Заявляется только размерность, устойчивая к изменению частотной сетки и
covariance model.

### 7.3. Perturbation ensemble

Первый безразмерный ensemble включает:

```text
cP and cS                  +/-1%
attenuation                +/-20%
port amplitudes            3% RMS
port delay/phase           0.02 carrier cycles RMS
lateral registration       +/-0.10 lambda_P0
```

Для каждой реализации сохраняются operator metrics и calibration fidelity
между nominal и perturbed defect dictionaries.

### 7.4. Provisional O2c gates

Для четырёхпортового кандидата:

```text
nominal minimum detectability       >= 1.0 x A-scan
Q10 minimum detectability           >= 0.8 x A-scan
nominal nonlocal coherence          <= 0.50
Q90 nonlocal coherence              <= 0.60
nominal stable rank                 >= 10
Q10 stable rank                     >= 8
Q10 calibration fidelity            >= 0.90
max sigma(W(omega))                 <= 1 + numerical tolerance
```

Провал gate не отправляет проект к новой геометрии. Сначала пересматривается
число портов, codebook, ROI, разрешение или calibration model.

### 7.5. Одно состояние не прошло; два состояния прошли

Сначала были проверены выбранные на O2b одно-состоянийные targets. Затем ещё
`1000` четырёхпортовых symmetry-codebooks были оптимизированы непосредственно
по причинному pulse operator. Лучший номинальный gate score остался равным
`0.952`; следовательно, провал не списывается на один случайный carrier
codebook.

Минимальное успешное расширение сохраняет четыре электронных канала, но
использует два последовательных статических вещественных состояния. Энергия
между состояниями делится множителем `1/sqrt(2)`. Оба состояния reciprocal,
contractive и сохраняют even/odd сектора отражения.

| Кандидат | Channels | States | Q10 min / A-scan | Q90 coherence | Q10 stable rank | Q10 calibration fidelity | O2c |
|---|---:|---:|---:|---:|---:|---:|---|
| real one-state O2b | 4 | 1 | 0.843 | 0.568 | 7.364 | 0.909 | fail |
| complex delay one-state | 4 | 1 | 0.915 | 0.603 | 7.475 | 0.884 | fail |
| real two-state A | 4 | 2 | **1.120** | 0.493 | **10.074** | 0.908 | **pass** |
| real two-state B | 4 | 2 | 0.995 | **0.489** | 9.846 | **0.910** | **pass** |

Точный текущий архитектурный лидер:

```text
4 electronic channels
    x 2 sequential passive static real states
    x 2 even + 2 odd symmetry ports per state
    -> O3 entrant.
```

Complex true-time-delay ветка не продвигается: её лучший двухсостоянийный
вариант сохраняет высокую обнаружимость и ранг, но имеет Q10 calibration
fidelity `0.884 < 0.90`.

## 8. Текущий этап O3 и условие перехода к геометрии

O2c закрыт; первый четырёхканальный класс был передан в O3, после чего
пространство расширено на static и motion ветви. Все классы обязаны пройти
общий paired transient pulse-echo benchmark с defect/no-defect, ringing,
backwall, common receive chain и stochastic clutter. Полный preregistered
benchmark задан в
[`operator_o3_benchmark_contract.md`](operator_o3_benchmark_contract.md).

Обязательные comparators:

```text
uniform A-scan
equal-time mechanical scan
four-channel one-state negative control
four-channel two-state entrant
PAUT/FMC oracle.
```

Основные O3 observables — `Pd(alpha)` при фиксированном `Pfa=1%`, амплитуда
`alpha90`, вероятность совместного detection/localization, ошибка локализации
и доля закрытого разрыва A-scan -> FMC на общем clutter ensemble. Codebook,
threshold-validation и final-test разделены; идеальное вычитание defect-free
reference не допускается без оплаты соответствующего измерения.

### 8.1. Результат reduced O3

Frozen campaign выполнен на независимых train/validation/test nuisance
libraries, `20 000` validation и `20 000` test no-defect trials и `100`
material/contact/calibration contexts. Полный отчёт находится в
[`paired_clutter_born_v0/`](../results/operator_space/paired_clutter_born_v0/README.md).

| Архитектура | Pfa upper 95% | Pd at reference | Pdl at reference | alpha90 |
|---|---:|---:|---:|---:|
| uniform A-scan | 0.022% | 0.539 | 0.180 | 644.1 |
| equal-time mechanical scan | 0.926% | **0.990** | **0.983** | **177.8** |
| four-channel one-state | 2.408% | 0.953 | 0.936 | 233.4 |
| four-channel two-state | **0.745%** | 0.942 | 0.937 | 235.3 |
| FMC oracle | 6.246% | 0.979 | 0.972 | 138.4 |

Двухсостоянийный кандидат проходит собственные detection, localization, Pfa
и PAUT-gap gates. Он закрывает `80.8%` aggregate разрыва A-scan -> FMC; Q10
по всем `100/100` C3 contexts равен `0.775` при gate `0.50`.

Reduced O3 тем не менее **не пройден**: при тех же восьми transmit events
одноканальный механический скан имеет `alpha90=177.8` против `235.3` у
статического кандидата. Второй passive state в основном улучшает rejection
неизвестного background: Pfa upper падает с `2.41%` для one-state до `0.75%`,
но detectability почти не меняется.

### 8.2. Параллельная разметка обеих ветвей

Вместо немедленного продуктового выбора обе ветви пересчитаны независимо на
общем **новом** контракте с base seed `0x04f0a117`. Старый O3 final-test не
использовался для выбора codebooks или scan layouts. Обе кампании используют
одинаковые nuisance powers, `10 000/10 000` H0 trials, `24/60`
validation/test signal contexts и одну амплитудную сетку.

Static campaign охватывает `2 304` codebooks, `192` nuisance-aware finalists
и шестнадцать заранее замороженных `(K=1..4,Q=1..4)` winners. Motion campaign
охватывает `32` layouts при `P=1,2,3,4,6,8,12,17`: uniform, compact,
Chebyshev, aperture sweep и train-only D-optimal.

Полная объединённая карта находится в
[`combined_frontier_v0/`](../results/operator_space/combined_frontier_v0/README.md).

### 8.3. Static frontier

Три статические архитектуры проходят все intrinsic gates:

| Static operator | Channels `K` | States `Q` | Tx events | alpha90 | Pfa upper 95% |
|---|---:|---:|---:|---:|---:|
| unconstrained | 4 | 2 | 8 | 185.7 | 1.167% |
| unconstrained | 3 | 3 | 9 | 183.9 | 0.906% |
| balanced even/odd sectors | 3 | 4 | 12 | **168.6** | **0.862%** |

`K3,Q4` является лучшим проверенным static endpoint по `alpha90`. `K4,Q2`
использует минимальное число событий, `K3,Q3` уменьшает channels и states.
Кажущийся лидер `K4,Q4` с `alpha90=148.0` отвергнут из-за
`Pfa upper=3.21%`: это background over-specialization, а не Pareto-точка.

Balanced reflection sectors выигрывают полезный high-information endpoint;
pure even или pure odd family часто сохраняет точную неоднозначность. Для
аналитического 2D -> 3D lift это аргумент в пользу совместного использования
irreducible sectors, а не одного симметричного поля.

### 8.4. Motion frontier

| Motion operator | Channels | Positions/events | alpha90 | Reference gate |
|---|---:|---:|---:|---|
| train-D-optimal | 1 | 8 | 193.6 | narrow fail (`Pd=0.895`) |
| compact span `5.5 lambda` | 1 | 12 | 183.4 | pass |
| full monostatic aperture | 1 | 17 | **159.7** | pass |

При optimistic `T=events` первый aggregate crossover относительно прежнего
frozen `K4,Q2` возникает около `P=12`. Но новый static `K3,Q4` при тех же
двенадцати событиях лучше: `168.6` против `183.4`. Только полный `P17` scan
обгоняет лучший проверенный static endpoint примерно на `5.3%`, используя
семнадцать событий и шестнадцать перемещений.

При ненулевой цене settle/travel static branch быстро выходит вперёд. При
жёстком ограничении одного электронного канала motion остаётся сильным классом.

### 8.5. Объединённый ресурсный вывод

Одного победителя нет. Главные точки взаимно недоминируемы:

```text
static K4,Q2   4 channels,  8 events, 1 switch,  alpha90=185.7
static K3,Q3   3 channels,  9 events, 2 switches, alpha90=183.9
static K3,Q4   3 channels, 12 events, 3 switches, alpha90=168.6
motion P12     1 channel,  12 events, 11 moves,   alpha90=183.4
motion P17     1 channel,  17 events, 16 moves,   alpha90=159.7
FMC            17 channels, 17 events,            alpha90≈134.
```

Ресурсный вектор поэтому не сворачивается в одну назначенную цену. В
[`combined_resource_frontiers.csv`](../results/operator_space/combined_frontier_v0/combined_resource_frontiers.csv)
сохранены четыре sensitivity conventions: `events_only`,
`balanced_overhead`, `motion_expensive`, `state_expensive`.

Обе ветви остаются открытыми. Если не выбирать между ними, следующий
разрешённый physics step — paired full-wave сравнение не более четырёх frozen
representatives:

```text
static K4,Q2  — minimum-event valid static
static K3,Q4  — best valid static alpha90
motion P12    — first valid motion class
motion P17    — best tested low-channel alpha90.
```

### 8.6. Full-wave carrier F0

Четыре frozen representatives перенесены в одну 2D plane-strain
full-elastic постановку. Полный `17 x 17` элементный отклик решается для
общей пары «заполненный образец / свободная круглая полость `r=2 mm`»;
static, motion, A-scan и FMC затем являются только линейными проекциями одной
физической матрицы. Контракт переноса находится в
[`operator_o3_fullwave_contract.md`](operator_o3_fullwave_contract.md),
а прошедший refinement-run — в
[`fullwave_carrier_commonmesh_h1p0_v0/`](../results/operator_space/fullwave_carrier_commonmesh_h1p0_v0/README.md).

Номинальные независимые сетки выявили численный blocker: paired mirror
residual был `25.8%`. Refinement независимых сеток дал `10.2%`, общая
conforming-сетка `h=1.5 mm` — `4.70%`, и общая сетка `h=1.0 mm` —
`1.805% < 3%`. Исходный gate не менялся. На финальной сетке reciprocity
residual равен `5.85e-14`; tangential/normal norm для scattering равен `1.25`,
то есть mode conversion не исчез из модели.

Одночастотный центральный круглый дефект почти одномодов и потому особенно
благоприятен для uniform A-scan. Его carrier contrast нельзя превращать в
архитектурный verdict: здесь нет ROI-словаря, локализации, pulse separation
или clutter ensemble. Содержательный вывод F0 — физический оператор и
нормировки перенесены, а не то, что A-scan выиграл O3.

F0 разрешает следующий уровень `F1`: broadband пятицикловый paired
pulse-echo с конечным defect ensemble, backwall/ringing/contact/material
variations и исходными `Pfa/Pd/Pdl/alpha90` gates. До положительного F1 O4 и
синтез геометрии остаются заблокированы.

### 8.7. Direct transient F1a и stochastic F1b

Прямой multi-RHS Newmark solver строит `17 x 17 x 510` response за одну
факторизацию на mesh. Пятицикловый `242 kHz` burst рассчитан до `70 us` с
шагом `0.138 us`; свободный backwall и P/S impedance boundaries сохранены.
Полный `F1a`-отчёт находится в
[`fullwave_transient_v0/`](../results/operator_space/fullwave_transient_v0/README.md).

`F1a` прошёл physics gates:

```text
reciprocity residual              3.73e-14
paired mirror residual            0.549%       (gate <= 3%)
precursor energy fraction         2.63e-11     (gate <= 1%)
paired tangential/normal norm     0.662
```

Затем рассчитан finite full-wave словарь из шести классов: два размера void,
три глубины, lateral offset и наклонная щель. Frozen stochastic protocol
использует train-only nuisance subspace, общий FMC-derived electronics noise,
независимые `10 000/10 000` H0 validation/test trials и contact/material/
ringing variations. Отчёт находится в
[`fullwave_f1b_pilot_v0/`](../results/operator_space/fullwave_f1b_pilot_v0/README.md).

Ни один static/motion кандидат не прошёл F1b pilot:

| Architecture | Pfa upper 95% | Pd ref | Pdl ref | Min class Pd | alpha90 |
|---|---:|---:|---:|---:|---:|
| static `K4,Q2` | 1.532% | 0.858 | 0.851 | 0.183 | 2.702 |
| static `K3,Q4` | 1.756% | 0.849 | 0.838 | 0.238 | 2.713 |
| motion `P12` | **1.113%** | 0.581 | 0.562 | 0.058 | 5.046 |
| motion `P17` | 1.393% | 0.603 | 0.584 | 0.071 | 5.067 |
| FMC | 1.575% | **0.978** | **0.976** | **0.892** | 1.523 |

Главный blocker ещё фундаментальнее: uniform A-scan имеет `alpha90=0.635`,
то есть лучше FMC на данном узком prior. Изолированные дефекты различаются по
времени, а coherent uniform transmit получает SNR advantage; в словаре нет
structural clutter или нескольких одновременных scatterers. Поэтому
`A-scan -> FMC` gap не определён, и pilot не является валидным пространством
для выбора приближения к PAUT.

Это было no-go для frozen F1b и требование пересобрать task ensemble: плотная
lateral ROI, full-wave multi-scatterer background и явный PAUT positive
control. Возвращение было сделано к operator prior, а не к линзе.

### 8.8. Structural clutter F1c и expanded calibration F1d

F1c использует пять одинаковых targets `r=1 mm` на глубине `60 mm` и
lateral-позициях `-28,-14,0,14,28 mm`; каждый physical context содержит пять
настоящих traction-free background voids. H0 и пять H1 строятся из одной
conforming master mesh, поэтому multiple scattering решается direct transient
без independent mesh jitter. Контракт находится в
[`operator_o3_f1c_structural_contract.md`](operator_o3_f1c_structural_contract.md).

Первый four-context pilot создал требуемое пространственное преимущество FMC
и выбрал unconstrained `K4,Q2`, но следующий fresh context сорвал Pfa даже у
FMC. Это показало недостаток числа физических calibration contexts, а не
исчезновение spatial information.

В revised F1d старые пять наблюдавшихся contexts явно переведены в
development/train. Threshold заморожен на двух новых validation, final metrics
получены на двух других test contexts. Общий electronics noise сохранён.
Полный отчёт находится в
[`fullwave_f1d_expanded_v0`](../results/operator_space/fullwave_f1d_expanded_v0/README.md).

| Architecture | Pfa upper 95% | Pd@3 | exact Pdl@3 | Min Pd | alpha90 | Gap closure | Gate |
|---|---:|---:|---:|---:|---:|---:|---|
| A-scan | **0.027%** | 0.537 | 0.221 | 0.000 | 9.697 | 0.000 | fail |
| smooth `K3,Q1` | 0.089% | 0.982 | 0.976 | 0.957 | 2.427 | 0.910 | **pass** |
| smooth `K4,Q1` | 0.312% | **0.993** | **0.983** | **0.988** | 1.911 | 0.974 | **pass** |
| unconstrained `K4,Q2` | 1.361% | 0.980 | 0.979 | 0.947 | 2.120 | 0.948 | fail |
| FMC | 1.026% | 0.989 | 0.983 | 0.980 | **1.706** | 1.000 | oracle |

PAUT positive control прошёл: FMC template stable rank `4.300`, coherence
`0.146`, exact Pdl `0.983`, Pfa upper `1.026%`. По заранее заданному resource
priority выбран `K3,Q1`: три электронных канала, одно состояние и три events.
`K4,Q1` остаётся robust fallback.

### 8.9. O4 transfer target и O5 handoff

Положительный O3 сформировал O4:

```text
target multiport H(omega)
amplitude / phase / group-delay tolerances
allowed insertion loss and reflection
port isolation / reciprocity
minimum useful calibration fidelity.
```

Финальный target находится в
[`o4_k3q1_transfer_target_v0`](../results/operator_space/o4_k3q1_transfer_target_v0/README.md).
Это первые три гладких ортонормированных DCT aperture modes (`even/odd/even`)
в полосе `193.6--290.4 kHz`. Hard gates: amplitude error `<=3% RMS`,
differential delay `<=82.6 ns RMS`, transfer fidelity `>=0.98`, port coherence
`<=0.15`, reciprocity residual `<=3%`. Без drive compensation loss ограничен
`0.921 dB` one way.

O5 handoff зафиксирован в
[`operator_o5_geometry_handoff.md`](operator_o5_geometry_handoff.md): один
трёхпортовый монолитный нерезонансный mode transformer с непрерывной
апертурой, без element-wise phases и без state switch. Первый lift —
`translation_kz`; objective будущей единственной геометрии — совпадение
измеренного `H(omega)`, не focal peak.

### 8.10. O6a: короткая однородная modal plate

Первый physical twin был заморожен до расчёта: три распределённых DCT-порта,
однородная полоса Al 6061 шириной `71.6 mm` и длиной `3.0 mm`, затем
согласованное продолжение в алюминий. Скалярная P-оценка для второй чётной
гармоники на нижней границе полосы даёт предельную длину `4.413 mm` по
`82.6 ns`; поэтому длина `3.0 mm` выбрана с запасом, а не из sweep.
Контракт находится в
[`operator_o6_physical_twin_contract.md`](operator_o6_physical_twin_contract.md),
полный результат — в
[`o6a_modal_gate_v0`](../results/operator_space/o6a_modal_gate_v0/README.md).

Полноупругий quadratic port pencil выявил `14--21` открытых правых P/S-mixed
ветвей в полосе. Энергия и низкая межпортовая coherence сами по себе не стали
проблемой, но требуемые пространственные профили не являются тремя чистыми
модами такого широкого свободного guide:

| Метрика dense band | Результат | Gate |
|---|---:|---:|
| source-profile fidelity | `0.96099` | `>=0.98` |
| propagated fidelity | `0.98190` | `>=0.98` |
| maximum port coherence | `0.06896` | `<=0.15` |
| column imbalance | `0.366 dB` | `<=0.75 dB` |
| differential amplitude RMS | `16.29%` | `<=3%` |
| differential delay RMS | `>1 us` около cut-on | `<=82.6 ns` |
| useful one-way amplitude | `0.93783` | `>=0.8994` |

Targeted discretization check `120/160/200` элементов по сечению даёт
source fidelity `0.96103/0.96109/0.96111` и amplitude RMS
`15.72/15.26/15.27%`. Фазовая производная около открытия ветви сходится
хуже, но verdict от неё не зависит: два независимых shape-gate уже устойчиво
провалены.

По заранее заданному stopping rule 2D FEM прямоугольника и as-built F1d replay
не запускались. Отвергнут не оператор `K3,Q1`, а конкретный механизм
`DCT profiles -> laterally guiding uniform Al plate`. Переход к `K4,Q1` не
лечит этот механизм: он сохраняет те же первые три профиля и добавляет ещё
одну более высокую гармонику.

Следующий разрешённый `O6b` — load-local modal transducer: три фиксированных
распределённых пьезоэлектрических coupling-профиля действуют непосредственно
на контактной плоскости и не проходят через латеральный Al-guide. Сначала
нужен coupled electromechanical/radiation upper bound с общей нагрузкой и
reciprocity; геометрия слоя разрешается только если этот bound проходит O4.

### 8.11. O6b: direct load-local PZT coupler

O6b сохранил `K3,Q1` и перенёс три DCT-профиля непосредственно к нагрузке.
Контракт находится в
[`operator_o6b_direct_coupler_contract.md`](operator_o6b_direct_coupler_contract.md).

Сначала рассчитан необходимый power-wave upper bound: точный outgoing
plane-strain P/S surface impedance Al для `k_y=q*pi/W`, продольный mechanical
port PZT-5A и заранее замороженный `100 um` epoxy bond. Dense `49`-точечный
bound прошёл весь O4:

| Upper-bound metric | Result | Gate |
|---|---:|---:|
| fidelity | `0.999993` | `>=0.98` |
| imbalance | `0.072 dB` | `<=0.75 dB` |
| differential amplitude RMS | `0.353%` | `<=3%` |
| differential delay RMS | `0.705 ns` | `<=82.6 ns` |
| useful amplitude | `0.90457` | `>=0.8994` |

Это авторизовало одну finite geometry без sweep: `6.0 mm` изотропный
mechanical PZT-5A surrogate, `0.10 mm` epoxy и `35 mm` matched Al load.
Три непрерывные DCT-тяги приложены к общей задней поверхности PZT; выходной
`H(omega)` измерен в 17 bins на bonded Al interface на `13` частотах.
Полный пакет находится в
[`o6b_stack_twin_v0`](../results/operator_space/o6b_stack_twin_v0/README.md).

Finite twin не прошёл O4:

| Finite-twin metric | Result | Gate |
|---|---:|---:|
| fidelity | `0.97331` | `>=0.98` |
| port coherence | `0.04499` | `<=0.15` |
| imbalance | `0.230 dB` | `<=0.75 dB` |
| differential amplitude RMS | `19.59%` | `<=3%` |
| differential delay RMS | `846 ns` | `<=82.6 ns` |
| useful amplitude | `0.88250` | `>=0.8994` |
| reciprocity residual | `1.31e-15` | `<=3%` |

Refinement `PZT/Al h=0.60/1.00 -> 0.45/0.75 mm` на трёх band anchors
сохранил worst fidelity `0.97331 -> 0.97262` и amplitude RMS
`19.59 -> 19.84%`. Следовательно, failure не является mesh artifact.

Физический вывод: load-local coupling как предел возможен, но сплошная
изотропная PZT-пластина сама становится коротким латеральным multimode guide и
возвращает механизм O6a. Согласно stopping rule electrical terminal model и
F1d replay не запускались.

Следующий разрешённый `O6c` — один laterally-decoupled `1-3` composite /
effective-anisotropic layer, который физически приближает прошедший O6b upper
bound. Это изменение механизма, а не sweep толщины: поперечная effective
stiffness должна быть подавлена при сохранении thickness coupling. Сначала
замораживается один homogenized tensor и проверяется тот же finite `H(omega)`;
явная pillar/kerf geometry разрешается только после такого pass.

### 8.12. O6c: непрерывный homogenized 1-3 composite

До solve был заморожен один `65%` ceramic random-fibre effective medium:
`rho=5461 kg/m3`, измерительный axial impedance `17 MRayl`,
`(C11,C22,C12,C66)=(52.921,17.465,8.674,4.391) GPa`, `Q=35`. Axial stiffness
получена из импеданса; lateral, cross и shear stiffness — harmonic mixtures
проектных PZT и photopolymer. Ни одно число не подбиралось по O4. Контракт —
[`operator_o6c_composite_contract.md`](operator_o6c_composite_contract.md),
результат —
[`o6c_composite_twin_v0`](../results/operator_space/o6c_composite_twin_v0/README.md).

На неизменной геометрии `6.0 mm composite -> 0.10 mm epoxy -> 35 mm Al` и
той же 13-точечной полосе получено:

| Finite-twin metric | O6c | Gate |
|---|---:|---:|
| fidelity | `0.99384` | `>=0.98` |
| port coherence | `0.01748` | `<=0.15` |
| imbalance | `0.173 dB` | `<=0.75 dB` |
| differential amplitude RMS | `9.02%` | `<=3%` |
| differential delay RMS | `641 ns` | `<=82.6 ns` |
| useful amplitude | `0.92970` | `>=0.8994` |
| maximum singular value | `0.97390` | `<=1` |
| reciprocity residual | `5.94e-16` | `<=3%` |

Lateral decoupling дала существенное, но недостаточное улучшение относительно
O6b: fidelity `0.97331 -> 0.99384`, amplitude error `19.59% -> 9.02%`, delay
`846 -> 641 ns`. Худший отдельный профиль — DCT2 (`10.65%`, `784 ns`).
Common raw compliance монотонно уменьшается по sampled band, поэтому данные не
позволяют объяснить failure только номинальным half-wave resonance слоя около
`259.4 kHz`: провал остаётся дифференциальным пространственным эффектом
непрерывной среды.

По preregistered stopping rule mesh refinement, explicit homogenization cell,
electrical terminal model и F1d replay не запускались. Закрыт именно класс
`continuous 6 mm effective layer`, а не идеальный load-local bound.

Следующий допустимый механизм `O6d` обязан разорвать латеральный механический
путь топологически, а не ещё раз уменьшать effective stiffness. Разрешён
только low-rank isolation bound для механически разделённых колонок с тремя
общими reciprocal electrical buses: без individual pixel phases, delay laws
и 17 независимых pulser/ADC. До доказательства пассивной реализации трёх
DCT-профилей и оценки fabrication/calibration burden явная kerf geometry и
FEM запрещены. Это отдельная проверка на «скрытую решётку», а не автоматический
переход к устройству.

### 8.13. O6d: mechanically isolated three-bus bound

O6d разделил две причины сложности до какой-либо kerf geometry. Контракт —
[`operator_o6d_isolated_bus_contract.md`](operator_o6d_isolated_bus_contract.md),
результат —
[`o6d_isolated_bus_bound_v0`](../results/operator_space/o6d_isolated_bus_bound_v0/README.md).

Сначала measured O6c transfer был разложен как
`H=g(omega)S+D(omega)`. В optimistic isolated model оставлена доля
`kappa*D`. Все O4 gates проходят только при
`kappa<=0.136698`: требуется подавить не менее `86.33%` differential
оператора O6c (`-17.28 dB` residual amplitude). Активным ограничением на
границе является delay `82.6 ns`; amplitude RMS там уже `1.236%`. Это
требование к оператору изоляции, не предсказание ширины kerf.

Независимый electrical upper bound разрешил идеальный lossless `3x3` port
rotation, continuous signed weights, fixed balancing, нулевые loss/delay и
точную reciprocity. После `200000` детерминированных SO(3) rotations и local
refinement получен frontier:

| Внутренняя карта | Активные связи | Fidelity | Amplitude RMS | Useful | O4 |
|---|---:|---:|---:|---:|---|
| одна bus на колонку | `17` | `0.95519` | `25.55%` | `0.91228` | fail |
| до двух bus на колонку | `34` | `0.99919` | `4.04%` | `0.99727` | fail |
| полная три-bus карта | `50` | `1.00000` | `0%` | `1.00000` | pass upper bound |

У two-bus карты проходит всё, кроме amplitude RMS `4.04% > 3%`. Поиски на
gate score, fidelity и minimum amplitude error сошлись к одному результату.
Это сильный numerical result, но не certified global optimum; поэтому no-go
относится к протестированному sparse family, а не ко всем возможным пассивным
сетям.

Аналитически точная full map требует `50` nonzero bus-to-column taps и `18`
различных continuous magnitude ratios. Quantization проходит O4 впервые при
`8` уровнях на DCT-port (`17` различных ratios по полной карте); четыре уровня
дают amplitude RMS `7.37%`. Такая реализация имеет три внешних канала, но
является именно пассивным `3-to-17` beamformer, который O5 заранее пометил как
fallback из-за fabrication/calibration burden.

По stopping rule kerf FEM и F1d replay не запущены. Следующий этап обязан
сохранить обе ветви открытыми и сравнимыми:

1. **beamformer branch:** сначала electrical loss/tolerance/cost bound полной
   `50-tap`, 8-level карты; geometry разрешается только после него;
2. **architecture-return branch:** вернуть на общий resource frontier уже
   проверенные mechanically scanned `P12/P17` варианты, не кредитуя им цену
   как трёхканальному static operator.

Выбирать между ними до общей таблицы fabrication, calibration, acquisition
time и detection quality запрещено.

### 8.14. O6e: passive-beamformer versus motion resource fork

O6e сохранил обе O6d-ветви и запретил искусственный scalar cost. Контракт —
[`operator_o6e_resource_fork_contract.md`](operator_o6e_resource_fork_contract.md),
результат —
[`o6e_resource_fork_v0`](../results/operator_space/o6e_resource_fork_v0/README.md).

#### Полная q8 passive map

Первая прошедшая O6d-квантизация (`50` active taps, `8` magnitude levels на
DCT-port, `17` различных ratios) имеет useful amplitude `0.999053`. Поэтому
номинальный common one-way loss ceiling до O4 равен `0.913 dB`.

Frozen F1d replay на тех же test contexts и random draws прошёл:

| Pfa upper | Pd | Pdl | Min Pd | alpha90 | Gap closure |
|---:|---:|---:|---:|---:|---:|
| `0.089%` | `0.982` | `0.976` | `0.957` | `2.4298` | `0.909` |

Относительно ideal DCT3 alpha90 ухудшился лишь на `0.137%`. Это проверяет
quantized operator, но ещё не parasitic physical network.

Fresh O6e tolerance ensemble использовал `50000` общих draws на каждой точке,
13 частот, единый worst-band contractive scale и joint O4 gate. Максимальный
общий multiplier с one-sided Wilson lower-95 `>=99%` равен `0.7516`, то есть:

```text
tap log-amplitude RMS       <= 2.255%
tap differential delay RMS <= 62.081 ns
observed joint pass         99.076%
Wilson lower 95             99.003%
q01 useful amplitude        0.930920
remaining common loss       0.299 dB
```

При полном исходном O4 perturbation budget `3% / 82.6 ns` jointly проходят
только `20.694%` draws. Следовательно, coefficient tolerance нельзя целиком
отдать отдельным taps: суммарные fabrication и network errors должны быть
ниже. Single-open audit проходит только для `2/50` самых слабых DCT2 taps;
fault tolerance архитектуре не засчитывается.

#### Несжатый resource frontier

| Architecture | Ch | Events | Positions / moves | T electronic / settle / travel / slow | Columns | Passive taps | Gap closure* |
|---|---:|---:|---:|---:|---:|---:|---:|
| q8 beamformer | `3` | `3` | `1 / 0` | `3 / 3 / 3 / 3` | `17` | `50` | `0.909` |
| compact P12 | `1` | `12` | `12 / 11` | `12 / 23 / 28.5 / 67` | `1` | `0` | `0.868` |
| uniform P17 | `1` | `17` | `17 / 16` | `17 / 33 / 41 / 97` | `1` | `0` | `0.931` |

`*` Beamformer использует full-wave frozen F1d, motion — reduced fresh Born
O3. Gap closures нормированы внутри своих evidence campaigns и не являются
одним ROC-экспериментом.

Точный обмен: beamformer тратит на два electronic channels больше и добавляет
`50` калибруемых passive taps, но экономит `9/14` transmit events и `11/16`
moves относительно P12/P17. Motion избегает distributed network, но сохраняет
непроверенные encoder error, couplant reset и probe-footprint риски. Поэтому
ни одна ветвь не доминирует полный resource vector.

Следующий разрешённый `O6f` снова парный:

1. **O6f-A:** одна реализуемая reciprocal q8 network topology с parasitic
   loss, bandwidth, impedance и tolerance accounting; без geometry sweep;
2. **O6f-B:** positioning/contact perturbation campaign для frozen P12/P17,
   без повторной оптимизации layout на test.

Scalar winner и prototype fabrication запрещены, пока обе отсутствующие
координаты не измерены на сопоставимых gates.

### 8.15. O6f: paired physical-risk audit

Контракт был заморожен до расчёта обеих ветвей:
[`operator_o6f_physical_risk_contract.md`](operator_o6f_physical_risk_contract.md).
Сводный вывод и следующий gate находятся в
[`operator_o6f_synthesis.md`](operator_o6f_synthesis.md).

#### O6f-A: reciprocal q8 lattice

Квантованный `17 x 3` оператор точно разложен на три почти единичных
singular-value attenuators, два малых `3 x 3` input transforms и
`45` adjacent-Givens output couplers. Численная ошибка реконструкции равна
`7.0e-16`. В optimistic модели, где input blocks и нагрузки пока идеальны:

```text
uniform loss ceiling / output coupler       0.047376 dB
active loss gate                            amplitude RMS = 3%
uniform delay ceiling / output coupler      15.252 ns
active delay gate                           differential delay = 82.6 ns
99%-Wilson angle RMS ceiling                0.4073 deg
99%-Wilson component delay RMS ceiling      8.1465 ns
```

Задержка первоначально была подписана как `passivity` из-за численной ничьей
`sigma_max=1+7e-16`; audit показал, что bisection в действительности
останавливает differential-delay gate. Результат исправлен без изменения
границы или Monte Carlo.

Это topology-bound PASS, но ещё не electrical prototype: общий loss budget
должен оплатить два входных hybrid blocks, dump paths, winding/dielectric
loss и mismatch к `17` пьезоэлектрическим нагрузкам.

#### O6f-B: frozen P12/P17 motion replay

Номинальные dictionaries, train nuisance precision и сохранённые thresholds
восстановлены без переобучения. Исходные `alpha90`, reference `Pd/Pdl` и
целые числа H0 false alarms воспроизведены точно. После этого выполнен paired
прогон на тех же `60` signal contexts и `12000` fresh H0 trials на точку.

| Layout | Scale | Position RMS | Contact log-amp RMS | Delay RMS | Pfa upper | alpha90 ratio | Pdl ref |
|---|---:|---:|---:|---:|---:|---:|---:|
| compact P12 | `2.0` | `0.10 lambda` | `10%` | `0.04 cycles` | `0.627%` | `1.0611` | `0.8777` |
| full P17 | `2.0` | `0.10 lambda` | `10%` | `0.04 cycles` | `0.927%` | `1.0834` | `0.9140` |

Оба варианта проходят все три frozen gates до верхней протестированной точки.
Ближайшее ограничение у обоих — относительный `alpha90 <= 1.10`; у P17 запас
меньше, несмотря на лучшую абсолютную чувствительность.

O6f не даёт scalar winner. Следующий `O6g` обязан закрыть две конкретные
неопределённости без переоптимизации операторов:

1. **O6g-A:** joint end-to-end budget полной q8 сети, включая оба input
   transforms, attenuator/dump paths и одновременные loss/tolerance errors;
2. **O6g-B:** найти coupled failure boundary P12/P17 выше `scale=2` и
   one-at-a-time attribution encoder/contact-amplitude/contact-delay.

После O6g дешёвый reduced surrogate для выигравшей координаты считается
исчерпанным: следующий положительный шаг должен повышать evidence level до
loaded RLC/piezo network или full-wave contact/lift-off, а не добавлять ещё
один Born sweep.

### 8.16. O6g: joint budgets and finite failure boundaries

O6g был последним разрешённым reduced-model этапом. Контракт —
[`operator_o6g_joint_boundary_contract.md`](operator_o6g_joint_boundary_contract.md),
сводный handoff —
[`operator_o6g_synthesis.md`](operator_o6g_synthesis.md).

#### O6g-A: complete q8 network

Обе signed input rotations разложены ещё на `3+3` Givens couplers; вместе с
output lattice полная сеть содержит `51` rotation couplers и три
attenuator/dump paths. Ошибка реконструкции q8 равна `6.70e-16`.

| Budget | Zero-error ceiling | Active gate |
|---|---:|---|
| insertion loss каждого из 51 coupler | `0.045074 dB` | amplitude RMS |
| excess loss каждого attenuator path | `0.912715 dB` | useful amplitude |
| multiplier луча через половину обоих ceilings | `1.18759` | useful amplitude |

При zero loss совместный `99%`-Wilson manufacturing scale равен `2.1572`:
`0.2157 deg` coupler-angle RMS, `4.314 ns` coupler/attenuator delay RMS и
`2.157%` attenuator log-amplitude RMS. На preregistered allocation
`0.022537 dB/coupler + 0.456357 dB/attenuator path` остаётся scale `1.1472`:

```text
coupler angle RMS             0.1147 deg
coupler delay RMS             2.294 ns
attenuator log-amplitude RMS  1.147%
attenuator delay RMS          2.294 ns
Wilson lower 95               99.005%
```

Это end-to-end ideal-network PASS, но не loaded electrical PASS.

#### O6g-B: P12/P17 motion thresholds

Те же contexts, error draws и noise streams O6f были продолжены через первый
отказ и 12-step bisection. Никакого reoptimization не выполнялось.

| Layout | Scenario | Max passing scale | Physical ceiling | First failed gate |
|---|---|---:|---|---|
| P12 | coupled | `2.35522` | `0.1178 lambda`, `11.78%`, `0.0471 cycles` | alpha90 |
| P17 | coupled | `2.23883` | `0.1119 lambda`, `11.19%`, `0.0448 cycles` | alpha90 |
| P12 | position only | `2.76808` | `0.1384 lambda` | alpha90 |
| P17 | position only | `2.57225` | `0.1286 lambda` | alpha90 |
| P12 | delay only | `3.14356` | `0.0629 cycles` | Pfa |
| P17 | delay only | `2.69906` | `0.0540 cycles` | Pfa |
| P12 | amplitude only | `5.52696` | `27.63%` | Pfa |
| P17 | amplitude only | `4.53560` | `22.68%` | Pfa |

P17 сохраняет лучшую absolute detection, но имеет меньший relative physical
margin по каждому сценарию. Для обеих архитектур coupled предел ниже любого
one-at-a-time: комбинацию encoder/contact phase/contact amplitude нельзя
сертифицировать независимыми scalar tolerances.

Reduced evidence теперь исчерпан. Следующий `O6h` не имеет права продолжать
идеальные Born/coupler sweeps:

1. q8 переносится в reciprocal scattering completion с реально нагруженными
   piezo/BVD terminals и broadband matching;
2. motion переносится в физическую layered-contact/lift-off transmission
   model с multiple reflections, а затем только прошедший вариант — в 2D
   full-wave spot check.

Prototype selection остаётся запрещённым до O6h.

### 8.17. O6h: loaded network and physical contact

O6h впервые повысил evidence level за пределы unilateral network и заданной
вручную пары contact amplitude/delay. Контракт —
[`operator_o6h_loaded_physics_contract.md`](operator_o6h_loaded_physics_contract.md).

#### Loaded q8

q8 embedded в точную reciprocal lossless `23 x 23` scattering matrix:
`3` public, `17` piezo и `3` matched dump ports. Isometry error `3.69e-16`.
Полный multiple-reflection resolvent с series-BVD loads дал:

```text
maximum passing effective Q       1.48148
active gate                       10 dB public return loss
mismatch audit Q                  0.74074
99%-Wilson resonance RMS          2.254%
99%-Wilson resistance log RMS     4.509%
99%-Wilson Q log RMS              11.272%
```

Intrinsic `Q=35` не проходит: worst-band return loss `0.069 dB`, useful
amplitude `0.1261`. Следовательно, q8 требует реального broadband /
multiresonant matching, а не только точных couplers.

#### Layered motion contact

Exact normal-incidence transfer matrix для жидкого слоя с multiple reflections
дал физические границы:

| Layout | Max scale | Thickness RMS | Z log RMS | c log RMS | First gate |
|---|---:|---:|---:|---:|---|
| P12 | `1.24463` | `38.06 um` | `12.45%` | `6.22%` | Pfa |
| P17 | `0.98279` | `30.05 um` | `9.83%` | `4.91%` | Pfa |

P17 уже не проходит nominal scale `1.0`. P12 проходит normal gate, но exact
2D solid--fluid--solid P/S spot с energy error `4.0e-15` его останавливает:
nominal minimum layer/direct power до `45 deg` только `0.654`, изменение
S-fraction `0.233`. Даже favourable one-sigma anchor даёт `0.725`.
Жидкостный moving contact для текущей ROI закрыт.

### 8.18. O6i: mitigation-class map

Контракт — [`operator_o6i_mitigation_contract.md`](operator_o6i_mitigation_contract.md).

Простое resistive broadening q8 от `Q=35` до `1.48148` требует `23.625x`
total resistance и оставляет только `4.233%` полезной мощности, amplitude
`0.2057 < 0.8994`. Этот mitigation STOP; lossless/multiresonant matching
остаётся отдельной незаданной topology-задачей.

Для motion exact three-solid 2D map проверил `168` density-valid property
points. `122` прошли исходный angular gate, `84` — строгий interior gate.
Frozen handoff по preregistered rule:

```text
ZP / sqrt(Zprobe Zsample)      0.75
cP / cP_Al                    1.00
cS / cP                       0.25
rho                           2053.6 kg/m3
thickness                     0.02 lambdaP = 505.96 um
minimum power, 0--45 deg      0.9441
maximum Delta S fraction      0.0317
```

Это property target для shear-transmitting dry/solid contact, не выбранный
коммерческий материал.

### 8.19. O6j: final numerical attachment qualification

Контракт —
[`operator_o6j_attachment_qualification_contract.md`](operator_o6j_attachment_qualification_contract.md),
сводный handoff —
[`operator_o6j_final_handoff.md`](operator_o6j_final_handoff.md).

Все `16/16` exact 2D attachment corners прошли: worst minimum power `0.9139`,
worst `Delta S=0.04985`, maximum energy-balance error `1.33e-15`.

Frozen P12 detector replay одновременно включал fixed encoder RMS
`0.05 lambda` и exact solid-contact response. На nominal attachment scale
`1.0`:

```text
Pfa Wilson upper 95%       0.441%
alpha90 / O6f zero error   1.00874
reference Pdl              0.9055
```

Finite attachment boundary равен `9.06836`: `90.68%` thickness log-RMS,
`45.34%` impedance log-RMS и `18.14%` cP log-RMS; первым ломается Pfa.
Это большой numerical margin, а не предлагаемый manufacturing tolerance.

**Текущий архитектурный verdict:** разрешён только conditional
compact-P12 solid-contact coupon/prototype: `1` electronic channel, `12`
events/positions, `11` moves, span `5.5 lambda`. q8 остаётся research fallback,
заблокированным до конкретного broadband/multiresonant matching circuit; P17
и жидкий moving contact остановлены. Следующие обязательные работы уже не
являются ещё одним geometry-free sweep: identification реального contact
material, finite-footprint/edge model, coupon bench, wear/debonding и
scan-time economics. Product POD claim по текущим расчётам запрещён.

### 8.20. O7: stationary two-bit global boundary

После запрета движения головки относительно образца прежний compact-P12
перестаёт быть продуктовым кандидатом. Новый контракт
[`operator_o7_global_state_contract.md`](operator_o7_global_state_contract.md)
проверил одну неподвижную `K3,Q4` головку с двумя глобальными бинарными
граничными управлениями. Подробный handoff —
[`operator_o7_stationary_state_handoff.md`](operator_o7_stationary_state_handoff.md).

На новом seed/split были независимо заморожены две ветви: две зеркальные пары
PZT-патчей с двумя пассивными RLC-нагрузками и две зеркальные пары
переключаемой механической жёсткости. Частотное кодирование осталось честным
matched `K3,Q1` comparator. Final test был открыт только после freeze.

| Frozen architecture | Pfa95 | alpha90 | Pd ref | Pdl ref | alpha90 gain над matched Q1 |
|---|---:|---:|---:|---:|---:|
| PZT `K3,Q4` | `1.134%` | `285.15` | `0.987` | `0.851` | `17.84%` |
| mechanical `K3,Q4` | `0.796%` | `270.67` | `0.996` | `0.950` | `34.16%` |
| ideal frozen `K3,Q4` | `0.229%` | `259.44` | `0.979` | `0.976` | — |
| FMC | `0.116%` | `220.19` | `0.984` | `0.982` | — |

Обе физические ветви проходят preregistered reduced gate и показывают, что
четыре состояния дают информацию сверх семи частот одного фиксированного
состояния. PZT-проход близок к границе по `Pdl` и `Pfa`; mechanical имеет
больший запас. Positive controls дали `max sigma(W)=0.98`, reflection error
`<2.5e-15`, точную взаимность и точное схлопывание `Q4 -> Q1` при нулевой
управляющей связи.

Это coupled-modal surrogate, не coupled piezo FEM. Он не содержит явную
геометрию блока, пьезотензор, клей, BVD terminals, switch parasitics или
конечную механику clamp. Поэтому следующий разрешённый `O8` сохраняет обе
замороженные ветви и повышает evidence level до plane-strain coupled physical
model. Старый O7 test повторно для выбора геометрии не используется.

### 8.21. O8: frozen plane-strain boundary realization

Контракт —
[`operator_o8_plane_strain_boundary_contract.md`](operator_o8_plane_strain_boundary_contract.md),
сводный handoff —
[`operator_o8_plane_strain_handoff.md`](operator_o8_plane_strain_handoff.md).

O8 перенёс обе O7-ветви в конечный `202.38 x 50.60 мм` aluminium
plane-strain block. Четыре управляющие зоны стали конечными интервалами,
source и output остались неподвижны. PZT-пары получили одну общую электрическую
степень свободы на пару и точный passive RL Schur complement; mechanical
ветвь — конечные Robin stiffness и actuator mass.

| Frozen branch | Amplitude/reference | Carrier min state distance | Max refinement projector error | Physical verdict |
|---|---:|---:|---:|---|
| PZT/RL | `0.721--0.779` | `0.0766` | `0.0498` | **pass** |
| mechanical stiffness | `0.978--0.992` | `0.0104` | `0.0575` | **stop** |

Mechanical O7 pass не пережил явную локальную границу: состояния почти не
меняют измерительное подпространство, хотя почти не теряют амплитуду. Эта
замороженная ветвь остановлена до явного нового resource/topology contract.

Только PZT открыл новый O8 detector split:

| Architecture | Pfa95 | alpha90 | Pd ref | Pdl ref | Verdict |
|---|---:|---:|---:|---:|---|
| same FEM head `K3,Q1` | `0.229%` | `421.65` | `0.810` | `0.778` | fail |
| FEM PZT `K3,Q4` | `0.129%` | `316.14` | `0.90025` | `0.884` | **pass** |
| ideal `K3,Q4` | `0.155%` | `277.94` | `0.940` | `0.937` | pass |
| FMC | `0.075%` | `219.19` | `0.968` | `0.966` | pass |

Четыре PZT-состояния улучшают `alpha90` на `25.02%` относительно той же
FEM-головки только в `q00`. Но активный `Pd>=0.90` пройден с запасом лишь
`0.000245`; результат считается conditional, а не robust.

Изначально после O8 fabrication была запрещена до O9 из-за несогласованной
ёмкости и слишком оптимистического `kt`. Этот запрет заменён только для одного
ограниченного de-risking coupon решением O8M ниже. Product prototype и POD
claim по-прежнему требуют domain-resolved PZT/adhesive, BVD/switch parasitics,
finite specimen contact и component tolerance ensemble.

### 8.22. O8M: material-realistic manufacturing freeze

Перед закупкой обнаружено, что O8 использовал `kt^2=0.45`, то есть
`kt=0.671`, тогда как публичные данные реальных bulk-PZT дают примерно
`0.49--0.55`. Замороженная O8-геометрия была повторена без оптимизации с
паспортными `epsilon_r`, `kt` и `Qm`; RL пересчитаны только для согласования с
фактической площадью электрода и ёмкостью.

| Material anchor | kt | C одной 20.25 x 8 x 1 мм пластины | Carrier min state distance | Amplitude/reference | Gate |
|---|---:|---:|---:|---:|---|
| CTS 3195HD bulk | `0.49` | `2.726 nF` | `0.0785` | `0.732--0.801` | pass |
| CTS 3203HD bulk | `0.55` | `5.451 nF` | `0.0817` | `0.729--0.794` | pass |
| PI PIC255 1-3 lower bound | `0.60` | `1.782 nF` | `0.1108` | `0.783--0.869` | pass |

Для изготовления выбран не composite, а CTS 3195HD: это обычная bulk-плита,
её меньшая ёмкость даёт более высокие сопротивления шунта и меньшую
чувствительность к relay/contact parasitics. Все четыре патча округлены до
одинаковых `20.00 x 8.00 x 1.00 мм`; outer pair сдвинута внутрь на доли
миллиметра и сделана flush с краями блока. Повторный full-band FEM уже этой
manufacturing geometry дал:

```text
amplitude/reference                 0.7367--0.8054
carrier minimum state distance      0.07689 >= 0.05
maximum raw mirror error            0.01083
maximum reciprocity error           6.98e-16
manufacturing gate                  PASS
```

**Терминальное решение: BUILD one coupon.** Пакет
[`manufacturing/o8_pzt_coupon_v1/`](../manufacturing/o8_pzt_coupon_v1/README.md)
содержит STEP, размерный чертёж, BOM, RFQ на шесть PZT-пластин, две RL-цепи,
сборку и preregistered bench acceptance. Номинальные нагрузки после
согласования площади равны `70.15 uH + 17.89 ohm` и
`73.27 uH + 27.57 ohm`; финальные значения пересчитываются один раз по
измеренной pair capacitance и полным паразитикам.

Coupon имеет ровно один разрешённый repair: электрическую ретуновку по
измеренным `C/L/R`. Если после неё minimum state distance остаётся ниже
`0.05`, топология получает **KILL** без движения патчей, снижения gate или
нового geometry sweep. Если coupon проходит, геометрия и parity сохраняются
для стационарной головки; coupon сам по себе ещё не является POD-валидированным
дефектоскопом.

# Часть II. Доказательная физическая TTD-ветка

Ниже сохранён прежний план линзовой ветки. Он фиксирует уже доказанные
ограничения пассивных задержек, полосы, модовой конверсии и изготовления, но не
является текущим порядком выбора продуктовой архитектуры.

## 1. Цель физической ветки

Исходная цель ветки — разработать плоскую с обоих торцов упругую линзовую
прослойку, которая из
одного электрического импульса формирует внутри твёрдого тела пространственно
и временно локализованный фокус:

```text
один электрический вход + пространственно согласованный преобразователь
    + пассивная TTD-топология
    -> согласованные пространственные задержки
    -> сфокусированный импульс с G_peak >= 2
    -> минимальное временное и спектральное размытие.
```

Фотополимер использовался для первого воспроизводимого прототипа. Текущий
физический кандидат — алюминиевая система
`input horn -> smooth delay guide -> output diffuser`, сопряжённая с
пространственно аподизированным пьезопреобразователем. Целевое прикладное
направление — диагностика алюминиевых деталей и другие задачи неразрушающего
контроля.

Механическая часть является пассивным аналогом `true-time-delay`-маски и имеет
заранее заданный фокус. Пространственный профиль источника разрешено задавать
геометрией, локальной пьезоэлектрической связью или долей активного материала,
но не независимыми фазовыми каналами. Электронное сканирование и перестройка
фокуса не входят в текущую задачу.

Под «высокой добротностью фокуса» понимаются высокий локальный gain, малые
боковые лепестки и сохранение формы импульса. Высокий резонансный `Q` не
является целью: узкий резонанс может дать звон, большую дисперсию и растянуть
импульс.

## 2. Физическая постановка

Пусть входной импульс имеет спектр `S_in(omega)`, а комплексное прохождение
`i`-й полосы линзы в целевой направляемой моде равно

```math
t_i(\omega)=A_i(\omega)e^{i\phi_i(\omega)}.
```

С учётом распространения от полосы до фокуса

```math
U_f(\omega)=S_{in}(\omega)
\sum_i a_i t_i(\omega)G_i(\omega),
```

где `G_i` — модовый или полноупругий propagator, а `a_i` учитывает освещение
апертуры.

Для когерентного импульсного фокуса во всей полезной полосе требуется

```math
\Theta_i(\omega)=\arg t_i(\omega)+\arg G_i(\omega)
\simeq \Theta_{common}(\omega)+2\pi n_i,
```

причём общая фаза должна быть почти линейной:

```math
\Theta_{common}(\omega)\simeq \Theta_0-\omega t_f.
```

Эквивалентные условия:

```math
\tau_{cell,i}(\omega)+\tau_{prop,i}(\omega)\simeq t_f,
```

```math
\frac{d^2\Theta_i}{d\omega^2}\simeq 0
```

в области, где энергия входного спектра существенна. Одновременно нужны
высокое и достаточно ровное `A_i(omega)`, малое отражение и малая конверсия в
нежелательные моды.

Отсюда следует центральное изменение постановки: восемь фаз на `242 кГц` —
только монохроматический ориентир. Для импульса проектируются **невёрнутые
задержки и полосовые передаточные функции**, а не фазы по модулю `2pi` в одной
точке.

## 3. Обязательная спецификация входа и референсов

### 3.0. Геометрическая конвенция стеночных профилей

Для `SinusoidalProfile`, `ExponentialProfile` и `PowerProfile` верхняя
стенка сдвинута на полпериода относительно нижней. На общем
профилированном участке обе границы являются вертикальными копиями,
поэтому толщина канала постоянна:

```math
g = H-A.
```

Прямые входной и выходной порты соединяются с профилем гладкими
полупериодами. `LegacySinusoidalProfile` не изменяется и служит
только для воспроизводимости исходного кода.

Сетки новой геометрии имеют суффикс `_STG`. Все полученные ранее
сигналы, графики прохождения, задержки и CWT для профилей без
`_STG` считаются архивными и не используются в новых выводах.

Для проверки гладкости экспоненциального семейства добавляется мягкий вариант
`κ=0.5`. Предел `κ->0` совпадает с синусоидальным профилем, поэтому отдельного
«ещё более мягкого экспоненциального» предела внутри этой формулы нет.
Сравнение выполняется не только визуально, но и по `max|dy/dx|`, кривизне у
прямого стыка и `max|d2y/dx2|`, после чего `κ=0.5` проходит тот же working-band
и CWT pipeline.

Пилот завершён: относительно синусоидального предела `κ=0.5` уменьшает
стыковую кривизну на `23%`, но увеличивает глобальные максимумы наклона на
`2%` и кривизны на `27%`. На `140 кГц` получено `|Hcal|=1.870` против
`1.882`, практически та же полоса и тот же envelope shift `-5.0 мкс`.
Дальнейший sweep по малым `κ` не выполняется; для реального уменьшения всех
производных нужно увеличивать длину перехода или уменьшать число периодов.

Ограниченный тест `κ=0.5`, `N=1` завершён при прежних `L=17 мм`, `A=2.5 мм`,
`g=4.5 мм`. Один период растянул полупереход примерно в `1.67` раза
относительно `N=2`: максимальный наклон уменьшился до `0.60`, а стыковая и
максимальная кривизны — до `0.36` от двухпериодных значений. Волновой результат
оказался противоположным геометрическому: спектральный максимум сместился с
`140` к `130 кГц`, пиковый `|Hcal|` упал с `1.870` до `0.625`, а на
`140 кГц` — до `0.569`. Отражающая граница увеличивает интеграл выходной
скорости в `9.35` раза против `3.37` для `N=2`, что указывает на более сильное
накопление/переотражение. Поэтому дальнейшее растягивание этого семейства не
выполняется: глобальная гладкость стенки сама по себе не гарантирует хорошего
модового и импедансного согласования.

Отдельным диагностическим семейством вводится `RoundedNotchProfile`: четыре
односторонних U-образных углубления чередуются между нижней и верхней стенками,
их входы почти вертикальны, а внутренний конец имеет радиус
`notch_width_mm/2`. В отличие от гладких `_STG`-профилей этот канал не имеет
постоянной толщины и проверяет именно резкий zig-zag путь. Он допускается в
сравнение только через тот же working-band gate; большая задержка при малом
прохождении не считается положительным результатом.

Первый вариант `A=2.5 мм`, `N=4`, `w=1.4 мм`, `g=4.5 мм` проверен на
`120--160 кГц`. При `140 кГц` он дал `|Hcal|=1.063`, но только `0.40` от
прямого референса по пику выходной скорости и `0.21` по интегралу `v^2`;
envelope shift равен `-0.71 мкс`. Поэтому эта конкретная alternating
U-геометрия не считается delay-кандидатом. Её продолжение разрешено только
после постановки полного multimode `S(omega)`, поскольку в порте `7 мм` она
нарушает зеркальную симметрию и может конвертировать падающую моду в нечётные.

До новой оптимизации геометрии фиксируется `PulseConfig`:

```text
waveform
f0
число периодов либо длительность
амплитуда и полная падающая энергия
окно/огибающая
полосы B_-6dB и B_-20dB
целевая измеряемая компонента: displacement/velocity/stress/energy
частота дискретизации и временное окно
```

Если экспериментальный импульс ещё не задан, численный baseline —
пятипериодный Hann-windowed tone burst с `f0 = 242 кГц`. Трёх- и
семипериодные импульсы используются как проверка чувствительности, а не как
три независимые задачи оптимизации.

Используются два референса:

1. **Физический референс** — однородный образец без линзы при той же падающей
   энергии, той же моде и том же положении измерения.
2. **Идеальный референс** — апертура с единичным прохождением и точными
   true-time delays. Она задаёт достижимый для выбранной апертуры потолок.

Старый узкий `R0` сохраняется только для абляции исторических резонаторов и не
используется как нормировка прохождения или усиления.

## 4. Метрики и рабочие критерии

### 4.1. Пространственный фокус

До сравнения фиксируется наблюдаемая `q`: для первого демонстратора это
амплитуда смещения целевой моды; для NDE она может быть заменена на скорость,
напряжение или плотность энергии, но не после просмотра результата. Основная
метрика:

```math
G_{peak}=
\frac{\max_t |q_{focus}^{lens}(t)|}
     {\max_t |q_{focus}^{uniform}(t)|}.
```

Она вычисляется при одинаковой падающей энергии. Дополнительно сохраняются:

- энергия в заданном пространственно-временном окне фокуса;
- поперечная FWHM;
- осевая `DOF_-3dB`;
- максимальный боковой лепесток;
- положение максимума и его смещение с частотой;
- Strehl ratio относительно идеальной импульсной апертуры.

Цель первой линзы:

```text
G_peak >= 2
DOF_-3dB <= 25--30 мм
без доминирующих боковых фокусов.
```

### 4.2. Сохранение импульса

Длительность определяется минимальным интервалом, содержащим `90%` энергии
аналитической огибающей. Вводятся

```math
B_t=\Delta t_{90}^{lens}/\Delta t_{90}^{uniform},
```

```math
\rho_{pulse}=\max_\tau
\operatorname{corr}
\left[u_{focus}^{lens}(t),u_{focus}^{ideal}(t-\tau)\right],
```

и отношение энергии хвоста после главного импульса к энергии главного окна.

Начальные инженерные gates, которые окончательно фиксируются после расчёта
референсного импульса:

```text
B_t <= 1.25
rho_pulse >= 0.90
postcursor energy / main-pulse energy <= 0.10.
```

Если они физически недостижимы одновременно с `G_peak >= 2`, строится явный
Pareto-фронт gain--broadening; пороги не ослабляются молча.

### 4.3. Модовые и спектральные характеристики

Для входного спектра определяются взвешенные характеристики, например

```math
\bar T_{00}=
\frac{\int |S_{in}(\omega)|^2T_{00}(\omega)d\omega}
     {\int |S_{in}(\omega)|^2d\omega},
```

и аналогично `R`, диссипация и суммарная модовая конверсия `C`.

Рабочий gate одиночного элемента:

```text
Tbar_00 >= 0.8
Rbar_00 <= 0.1
Cbar_unwanted <= 0.1
замкнутый энергетический баланс
нет узкого длинного ringing внутри B_-20dB.
```

Условия проверяются также поточечно на несущей и около спектральных максимумов.
Большая групповая задержка в области малого `T_00` не считается полезным
результатом.

### 4.4. Скейлограммы

CWT Морле строится для входа, однородного референса, идеального фокуса и
реальной линзы. Сохраняются:

- time-frequency ridge и его кривизна;
- временная ширина энергии на каждой существенной частоте;
- интегральный временно-частотный разброс;
- относительная энергия, корреляция огибающих и cone-of-influence mask.

CWT является независимой диагностикой формы и прихода импульса. Она не
заменяет комплексные `S(omega)` и не интерпретируется в областях слабой
энергии.

## 5. Доказанное текущее состояние

| Блок | Результат | Решение |
|---|---|---|
| Портовая модель | В полосе `4.2 мм` на `242 кГц` распространяются три моды; целевая чётная ветвь имеет `k = 1152.1 1/м`, `lambda_g = 5.45 мм` | Использовать направляемые моды, не bulk-P-приближение |
| Provisional-импульс | Для `5-cycle Hann @ 242 кГц`: `B_-6dB = 193.8--290.1 кГц`, `B_-20dB = 162.2--321.7 кГц` | Проектировать по фактической широкой полосе, не по `222--262 кГц` |
| Модовый риск | Для полосы `4.2 мм` extra-symmetric cut-on равен `252.9 кГц`, выше него лежит `35.2%` спектральной энергии импульса; при `3.2 мм` cut-on поднят до `332.0 кГц` | Provisional-апертура переводится на полосы `3.2 мм`; запас до верхней границы `B_-20dB` равен `10.3 кГц` |
| Абсолютный port-gate | Quadratic FEM на однородной полосе дал `T00 = 0.9999`, `R00 ~ 1e-7` | Пройден для симметричной rank-one-задачи |
| Идеальная фазовая апертура | `15` полос, `D = 71.4 мм`, `F = 35 мм`, `G_ideal = 2.998` | Только монохроматический baseline; пересчитать как true-time-delay aperture |
| Идеальная импульсная апертура | Оптимистический scalar-bulk baseline для `19 x 3.2 мм`, `D = 71.6 мм`, `F = 35 мм`: `G_peak = 3.52`, `DOF_-3dB = 23.5 мм`, FWHM `5.0 мм`, диапазон задержек `5.96 мкс` | Системный aperture/delay gate пройден только как верхняя граница; требуется полноупругий propagator |
| Дисперсионный delay-proxy | Шесть phase-matched path states длиной до `10.42 мм` при текущих Rayleigh-потерях дали `G_peak = 2.22`, `B_t = 1.08`, `rho = 0.971`; при `1.5x` потерях gain падает ниже `2` | Нерезонансная ветка физически перспективна, но почти без transmission-margin |
| Delay + phase trim proxy | Шесть состояний с путём до `6.20 мм` и независимым carrier phase-trim до `3.67 рад` дали `G_peak = 2.66` при текущих и `2.28` при удвоенных потерях | Основной робастный кандидат для следующей физической ячейки |
| In-plane packing | В шаге `3.8 мм` для полосы `3.2 мм` остаётся только `0.3 мм` centerline amplitude; даже оптимистический синусоидальный меандр требует длину `34.7 мм` для extra path `6.2 мм` и `57.1 мм` для `10.42 мм` | Буквальный плоский меандр закрыт как компактная ячейка; использовать slow-wave лигаменты, распределённую комплаенсность или третий размер |
| Global folded path | При одинаковом развёрнутом пути `30.42 мм` sharp fold дал относительный пик `0.284`, rounded fold — `0.565`; у обоих `E_y/E_x ~= 1.3` на выходе | Широкий изгиб закрыт: основная потеря полезного сигнала вызвана `P -> S`-конверсией, а не материальным затуханием |
| Horn-fed point radiator | Монотонный TTD-profile из 15 `sin^4`-каналов имеет extra path `0--9.724 мм`, общую осевую длину `50.664 мм`, максимальный вынос `13.746 мм` и `R_inner=3.93 мм`. Максимальный канал дал peak/equal-path `0.977`, энергию `0.962`, `B_t=0.991`, `rho=0.9998`, `E_y/E_x=0.0055`; измеренная group delay `5.516 мкс` против требуемых `5.777 мкс` | Центральный transient-gate пройден. Endpoint-calibrated 15-точечный импульс дал `G_peak=2.904`; quadratic receiver FEM — `93.17 нм`, gain `2.493`, FWHM `5 мм`, sidelobe `0.294`, то есть `+3.68%` амплитуды к лучшей бинарной линзе |
| All-Al double-horn TTD | Изолированный канал сохраняет `0.9965` envelope peak и `0.9964` энергии. Полный неаподизированный массив даёт `1.870` относительно исходного abrupt straight reference и `1.758` относительно matched-diffuser carrier control | Delay/throughput primitive доказан; главный оставшийся дефицит — пространственное распределение амплитуды источника |
| Согласованный сегментированный drive | Общие неотрицательные веса и `30-cycle Hann` дали five-frequency `G_peak=2.1261`, `B_t=1.1408`, `rho=0.9805`, postcursor `0.0497`; carrier gain равен `2.0503` относительно abrupt control и `1.9276` относительно matched-diffuser control | Системный `2x` gate пройден, но это результат `источник + линза`, а не чисто пассивной геометрии |
| Одновходовый passive feed | Реализованный balanced-tree surrogate даёт efficiency `0.87669`, power-only ceiling `1.99067` и routed carrier `1.80498`. Однако topology-optimal суррогаты дают `0.89056 / 2.00635` для contiguous tree и `0.89160 / 2.00752` без ограничения порядка | Рассчитанная feed-топология остановлена; глобальный no-go не доказан, но остаточный запас меньше `0.4%` и не оправдывает инженерный тюнинг без нового механизма |
| Distributed slow-wave FEM | Четыре симметричные прорезные секции дали `Tbar = 0.451`, `Rbar = 0.549` и delay `8.10 мкс`; предельно слабый печатный вариант дал `Tbar = 0.551`, delay `3.53 мкс` | Связность закрыта: delay получается только вместе со стоп-полосами |
| Tapered shunt trim FEM | Лучшее одиночное состояние на fine-grid: `Tbar = 0.977`, `Rbar = 0.023`, phase `+0.968 рад`, но `B_t = 1.38` lossless и `1.41` с текущими потерями; stop на `318 кГц` устойчив при lead `12 -> 24 мм` | Годится как диагностический phase actuator, но не проходит импульсный gate |
| Измеренный trim-cascade bound | Оптимистический каскад без взаимных отражений дал lossless `G_peak = 1.85`, `B_t = 1.60`, `rho = 0.845`; при Rayleigh x1 `G_peak = 2.01`, но `B_t = 2.20`, `rho = 0.829` | Каскад закрыт: достижение peak-gain сопровождается недопустимым разрушением импульса; новая библиотека запрещена |
| Равноплощадные through-holes | Круг и эллипсы двух ориентаций в прямом канале имеют на `242 кГц` `T00 = 0.0005--0.0104`, а во всей coarse-полосе `max R00 = 0.998--0.9995` | Диагностический механизм закрыт: ориентация меняет фазу, но не снимает почти полный stop |
| Fano/bright--dark | Узкая lossless-линия доказана, но исчезает уже при умеренных потерях | Не использовать как основную импульсную линзу |
| Две shunt-пары | Полезное фазовое покрытие `0.646 рад`, `G = 1.47` | Закрыто |
| Короткое series-сужение + shunt | Покрытие `0.472 рад`, `G = 1.375` | Закрыто |
| Компонент `M` | Чётная продольная мода `244.55 кГц`, mass fraction `0.996`, longitudinal fraction `0.966` | Компонент физически подтверждён |
| Компонент `K/2` | Положительная динамическая жёсткость в `222--262 кГц`, ближайшая чётная мода `380.5 кГц` | Компонент прошёл локальный gate |
| Открытая `K/2--M--K/2` | На `242 кГц`: `T00 = 0.085`, `R00 = 0.913`; даже на `246 кГц` сильное отражение | Текущая сборка не является фазовым элементом; dense sweep запрещён |
| Разнесённые стеночные `_STG`-профили | Общий максимум `|H_cal|` на `140 кГц`; оценка общей `-3 dB`-полосы `129.4--145.8 кГц` | Использовать `140 кГц` как несущую для сравнения границ; не переносить на неё старые выводы при `220 кГц` |

Дополнительные ограничения:

- около `253 кГц` открываются дополнительные чётные портовые ветви;
- rank-one DtN допустим только пока полоса и симметрия исключают их возбуждение;
- в полном массиве, при асимметрии и для широкого импульса нужен multimode DtN
  либо проверенный PML;
- переход от отдельных полос к общей выходной детали является отдельной
  модовой задачей.

## 6. Кандидатные топологии

### A. Mass--stiffness `K/2--M--K/2`

Это вторичный диагностируемый кандидат, а не выбранная универсальная
конструкция. Его локальный якобиан разрешён только как короткий контроль после
первого нерезонансного path-delay FEM; он не должен задерживать основную ветку.

Ветка продолжается только если существует хорошо обусловленное направление,
которое одновременно:

- уменьшает спектрально взвешенное отражение;
- сохраняет высокое прохождение;
- изменяет невёрнутую задержку в нужную сторону;
- не вводит сильную фазовую кривизну и ringing.

Если такого направления нет, меняется связность, а не плотность sweep.

### B. Нерезонансная true-time-delay-топология

Это текущая основная ветка. Оптимистический дисперсионный proxy показал, что
шесть состояний могут пройти импульсные gates при текущих потерях. Теперь нужно
реализовать добавочный путь и, при необходимости, отдельный слабый phase-trim
с абсолютным контролем отражения:

- несколько слабых согласованных секций вместо одного сильного резонанса;
- печатаемые slow-wave прорези или лигаменты с сохранённым прямым несущим
  путём; буквальный in-plane меандр внутри шага `3.8 мм` не используется;
- изменение пути либо эффективной фазовой скорости без полюсов внутри
  `B_-20dB`;
- зеркальная монолитная конструкция без закрытых полостей и поддержек.

Эта ветка потенциально длиннее, но естественнее для сохранения короткого
импульса. Она сравнивается с веткой A по delay--bandwidth, потерям, размерам и
печатным допускам, а не только по фазе на `f0`.

### B1. Входной collector, узкая delay-линия и точечный излучатель

Новая разрешённая связность использует каждый элемент линзы не как широкий
сквозной порт, а как пассивный канал

```text
плоская входная площадка
    -> нерезонансный horn / modal taper
    -> узкий стержнеподобный delay-guide
    -> короткий прямой выход
    -> малая площадка-излучатель в общей выходной детали.
```

Обратный horn до исходной ширины не требуется: малая выходная площадка является
одним из пространственно распределённых излучателей линзы. Один общий выход не
заменяет апертуру; будущая линза должна сохранить набор выходных фазовых центров
по всему диаметру. Входные площади и коэффициенты деления мощности разрешено
выбирать неравными, чтобы пассивно компенсировать различное затухание длинных и
коротких delay-ветвей.

Воронка не создаёт энергию. Рост локальной скорости в горле примерно как
`sqrt(S_in/S_throat)` используется только как диагностический признак сбора
мощности. Полезность элемента определяется P-импульсом, который малая выходная
площадка реально излучает в однородный принимающий материал. Обязательный
референс имеет те же входной horn, горло, выходную площадку и развёрнутую длину,
но прямой узкий guide. Кандидат отличается только delay-трактом.

Первый входной horn задаётся **без геометрического sweep**:

1. Для planar-пилота принять `h_in = 3.2 мм`, `h_throat = 1.6 мм`. Это даёт
   отношение площадей `R_A = 2` и идеальный локальный amplitude transformer
   `sqrt(R_A) ~= 1.414`. В axisymmetric-геометрии та же процедура применяется
   к площади `S`, а радиус восстанавливается из `S = pi r^2`.
2. Использовать log-cosine профиль площади

   ```math
   q(x)=\frac{1-\cos(\pi x/L_h)}{2},
   \qquad
   S(x)=S_{in}\exp\!\left[q(x)\ln\frac{S_{throat}}{S_{in}}\right],
   ```

   у которого `dS/dx = 0` на обоих стыках. Для 2D-модели постоянной глубины
   достаточно заменить `S(x)` на высоту `h(x)`.
3. Длину определить из локального adiabatic/WKB-критерия, а не перебирать:

   ```math
   \epsilon_{ad}=
   \frac{\max_x |d\ln Z_m(x)/dx|}{2k_{min}} \le 0.05,
   ```

   где `Z_m` — force/velocity impedance целевой продольной моды, а `k_min` —
   минимальный действительный волновой номер этой моды вдоль taper во всём
   `B_-20dB`. В первом одномерном недисперсионном приближении
   `Z_m ~= rho c_L S`; для реального guide используется вычисленный
   force/velocity modal impedance. Поэтому для указанного log-cosine профиля

   ```math
   L_h \ge
   \frac{\pi |\ln R_A|}{4\epsilon_{ad}k_{min}}
   \simeq 1.73\lambda_{max}
   \quad (R_A=2,\ \epsilon_{ad}=0.05).
   ```

4. Horn является travelling-wave переходом и не подгоняется к `lambda/2` или
   собственной частоте. После синтеза выполняется один полноупругий FEM с
   modal/PML-нагрузкой горла. Ориентиры: `T_P >= 0.9`, `R <= 0.1`, не менее
   `90%` прошедшей мощности в целевой моде и отсутствие заметного ringing.
5. Если вычисленная адиабатическая длина физически не помещается, единственный
   разрешённый fallback — синтезировать Klopfenstein impedance taper для
   фиксированных `Z_in`, `Z_throat`, нижней частоты полосы и заданного
   `R_max`, затем один раз проверить его FEM. Перебор linear/exponential/
   catenoidal/Bezier-профилей не выполняется.

Источники и границы применимости:

- [Guiman & Rosca, *A New Approach on Vibrating Horns Design*](https://doi.org/10.1155/2017/8532021) — Webster equation для продольной волны в solid horn переменного сечения; из работы берётся одномерная механика, но не резонансная настройка;
- [Rosca et al., *Experimental and numerical study on an ultrasonic horn...*](https://doi.org/10.1016/j.apacoust.2015.02.009) — переход от Webster/transfer-matrix модели к обязательной 3D FEM-проверке чистоты продольной моды;
- [Dahmani et al., *Piezoelectric Transduction of a Wavelength-Scale Mechanical Waveguide*](https://doi.org/10.1103/PhysRevApplied.13.024069) — пример modal decomposition и широкополосной передачи через elastic horn в узкий single-mode guide; материал и SH-поляризация отличаются от нашей задачи;
- [Klopfenstein, *A Transmission Line Taper of Improved Design*](https://doi.org/10.1109/JRPROC.1956.274847) — аналитический minimum-length impedance taper при заданном passband reflection; перенос допустим только после проверки одномодовости и отображения `Z_m -> S`.

Первый системный пилот включает принимающий однородный блок. Основные метрики:

- пик и энергия первого P-пакета в заданной точке на оси/в будущем фокусе;
- отношение `peak_P(delay)/peak_P(straight equal-path)` как диагностическая
  метрика; ориентир `0.8` проверяется после сборки первой полной линзы и не
  блокирует геометрическую интеграцию v0;
- невёрнутая задержка, `B_t`, корреляция и postcursor;
- отдельно энергия и время прихода S-компоненты;
- локальная амплитуда в горле только как вспомогательная диагностика.

Только если узкий delay-guide проходит этот end-to-end gate, разрешаются один
симметричный `1 -> 2` делитель и затем малый набор распределённых выходов.

### C. Поверхностная или присоединённая линза для алюминия

После выбора прикладной моды проверяется отдельный печатный модуль,
приклеенный или прижатый к алюминиевой детали. Возможные реализации —
поверхностная GRIN-структура, столбики или распределённая задерживающая маска.

Эта ветка не требует прорезать исследуемую деталь и поэтому предпочтительна
для дефектоскопии. Контакт, согласующий слой и модовая конверсия входят в
модель явно.

### D. Глобальный или суперъячеечный дизайн

Topology optimization и совместная настройка всего массива разрешаются только
если локальная библиотека нарушается из-за сильной связи соседей. Гладкие
сквозные отверстия и печатные ограничения встраиваются в параметризацию.
Чёрный ящик не используется до появления корректной полосовой прямой модели.

## 7. Универсальность

Универсальной должна быть топология и процедура синтеза, а не размеры.

```text
MaterialConfig:
    rho, c_p, c_s, loss_model(f)

ManufacturingConfig:
    min_wall, min_gap, min_radius, thickness, orientation

WaveConfig:
    pulse, target_mode, dispersion, aperture, focal_point
```

Неизменными остаются связность, физические роли регуляторов, правила
изготовления и последовательность gates. Пересчитываются размеры,
компонентные частоты, задержки и число необходимых состояний.

Основной масштаб задаётся дисперсией целевой направляемой моды:

```math
\lambda_g(\omega)=2\pi/\operatorname{Re}k_g(\omega),
```

а не только bulk-длиной `c_p/f0`. Простое геометрическое масштабирование не
сохраняет задачу при изменении `c_s/c_p`, коэффициента Пуассона, потерь,
толщины и печатного разрешения.

После первого успешного прототипа один и тот же synthesis-pipeline применяется
минимум к трём контрастным изотропным наборам свойств. Разрешён автоматический
пересчёт размеров; запрещена ручная смена связности под каждый материал.
Результатом может быть как печатаемая конструкция, так и честный `infeasible`
с указанием границы применимости.

## 8. Порядок работ и gates

### P0. Зафиксировать импульс и прикладной канал

- [x] Создать воспроизводимый provisional `PulseConfig` для
      `5-cycle Hann @ 242 кГц`.
- [x] Зафиксировать `B_-6dB = 193.8--290.1 кГц` и
      `B_-20dB = 162.2--321.7 кГц`.
- [ ] Зафиксировать модовую падающую энергию и полное временное окно.
- [x] Выбрать provisional-ширину полосы фотополимерного демонстратора
      `3.2 мм`: первая дополнительная чётная мода открывается только около
      `332.0 кГц`.
- [x] Для первого алюминиевого gate выбрана bulk longitudinal P-постановка
      толстого тела в 2D plane strain, а не конечная Lamb-пластина. При
      `242 кГц` получено `lambda_P=25.298 мм`, `lambda_S=12.743 мм`;
      finite-thickness plate требует отдельного повторного выбора моды.
- [x] Для NDE зафиксирован отдельный профилированный модуль с плоским
      bonded-контактом к алюминиевой детали. Источник — контактный
      thickness-mode longitudinal PZT, в модели matched port с равномерной
      нормальной тягой `1 МПа`; целевая компонента — `u_x`.

**Gate:** источник, референс, целевая мода и полоса однозначно определены.

### P1. Полосовой solver и восстановление импульса

- [ ] Продолжить все распространяющиеся и близкие эванесцентные портовые моды
      по `B_-20dB`; отметить cut-on/cut-off.
- [ ] Реализовать multimode DtN или верифицированный PML, если полоса пересекает
      открытие дополнительных чётных мод либо геометрия нарушает симметрию.
- [ ] Вычислять абсолютные комплексные `S(omega)` на адаптивной сетке частот.
- [ ] Реализовать общий конвейер
      `S(omega) * S_in(omega) -> IFFT -> envelope -> CWT -> metrics`.
- [ ] На однородной полосе и одной тестовой ячейке сравнить восстановленный
      импульс с прямым transient-FEM.
- [ ] Проверить причинность, достаточную ширину частотного окна, zero-padding и
      отсутствие FFT wrap-around.

**Gate:** одинаковые arrival time, waveform и энергия в frequency-domain и
transient-расчётах в пределах заранее записанного численного допуска.

### P2. Идеальная импульсная апертура

- [x] Пересчитать геометрические невёрнутые задержки всех полос; не применять
      phase wrapping по `2pi`.
- [x] Для `D = 71.6 мм`, `F = 35 мм`, шага `3.8 мм` построить оптимистический
      scalar-bulk идеальный сфокусированный импульс и пространственные метрики.
- [ ] Заменить scalar kernel на полноупругий homogeneous propagator и построить
      CWT в фокусе.
- [ ] Сравнить несколько апертур, шагов и фокусов по `G_peak`, FWHM, DOF,
      боковым лепесткам и требуемому диапазону задержек.
- [ ] Определить минимальное число уникальных delay states и допустимые ошибки
      амплитуды, задержки и дисперсии.
- [ ] Оценить пассивный потолок gain при реальном освещении апертуры.

**Gate:** существует идеальная печатно-разрешимая апертура с достаточным
запасом относительно `G_peak = 2` и временных критериев.

### P3. Выбор топологии одиночного элемента

- [x] Проверить оптимистические дисперсионные proxies для phase-matched path и
      group-matched path с отдельным carrier phase-trim.
- [x] Выполнить ограниченный диагностический пилот прямого канала со
      сквозными отверстиями: сплошной референс, круг и два равноплощадных
      эллипса разной ориентации. Сравнивать только абсолютные `T00`, `R00`,
      относительную фазу и excess group delay во всём `B_-20dB`; пилот не
      превращать в плотный sweep до положительного transmission--delay verdict.
      Для восьми отверстий площадью `4.02 мм^2` получена глубокая stop-band
      около `242 кГц`: `T00=0.0005--0.0104`, поэтому эта конфигурация не является
      широкополосной delay-ячейкой. Разрешён только отдельный ограниченный
      radius/pitch-detuning, если ветку отверстий решено продолжить.
- [x] Построить одну печатаемую нерезонансную path-delay-ячейку с прямым либо
      плавно изогнутым несущим путём в полосе `3.2 мм`.
- [x] Проверить состояния `0`, среднее и максимальное по абсолютным
      `T00(omega)`, `R00(omega)`, delay, GDD и ringing во всём `B_-20dB`.
- [ ] Только если path-delay не проходит gate, рассчитать ограниченный
      physical Jacobian `K/2--M--K/2` и один predicted-vs-FEM step.
- [x] До ветки Люнеберга выполнить бинарный global-fold gate: прямой образец
      толщины устройства, прямой guide той же развёрнутой длины, sharp
      mirror-like fold и rounded guided fold. Использовать один и тот же
      пятицикловый импульс; основной относительный критерий — пик полезной
      выходной компоненты fold/unfolded не ниже `0.8`, затем `B_t`, waveform
      correlation, postcursor и поперечная компонента. Этот пилот является
      полноупругим transient-screening, а не заменой будущего multimode
      `S(omega)`. При одинаковом пути `30.42 мм` sharp proxy дал peak ratio
      `0.284`, rounded guide — `0.565`; у обоих отношение выходной поперечной
      энергии к продольной около `1.3`. Наивный planar fold не проходит gate.
- [x] Проверить один horn-fed point-radiator без profile sweep: синтезировать
      входной log-cosine taper `3.2 -> 1.6 мм` по adiabatic-критерию; сравнить
      `horn -> straight narrow guide -> point radiator` и
      `horn -> rounded narrow delay guide -> point radiator` при одинаковой
      развёрнутой длине и одинаковом принимающем блоке. Gate ставится по
      излучённому P-импульсу, а не по локальной амплитуде горла. Tight bend
      дал peak `0.769` и не прошёл; один физически мотивированный gentle-step
      с `R_inner=3.93 мм` дал peak `0.940`, energy `0.901`, `B_t=1.020`,
      `rho=1.000` и прошёл gate.
- [x] Перенести delay-guide в плоскость `x-z`, сохранив шаг апертуры по `y`,
      и закрыть полный бинарный локальный triplet-базис. Curvature-limited
      `sin^4` state с осевой длиной `38.158 мм` и путём `42.578 мм` дал в
      transient peak/equal-path `0.967`, energy `0.956`, `B_t=0.990`,
      `rho=0.9988`, transverse energy `0.0043`; в representative `D-S-D`
      carrier neighbour penalty составил `-5.0%` и `-4.6 deg`.
- [x] Сравнить кандидатов по единому полосовому и печатному gate и выбрать
      рабочую baseline-топологию. Выбрана связка
      `log-cosine horn -> smooth sin^4 waveguide -> point radiator`: она
      сохранила импульс максимального delay-state и дала положительный gain в
      полной 15-канальной 3D-сборке.

**Стратегический verdict:** связка
`log-cosine horn -> smooth sin^4 waveguide -> point radiator` является
нормально работающей базовой стратегией проекта. Воронка собирает импульс в
узкий канал, плавный волновод создаёт монотонную true-time delay, а малый
выход формирует апертуру. Максимальный delay-state сохранил
peak/equal-path `0.977`, `B_t=0.991`, `rho=0.9998` и transverse energy
`0.0055`; первая полная 15-канальная 3D-сетка сфокусировала `|ux|=3.893 нм`
против `2.645 нм` у согласованной прямой решётки, то есть дала baseline
carrier gain `1.472` (`+47.2%`). Последующее phase-controlled уточнение
symmetry-constrained half-domain довело сетку тракта до `0.338 мм`
(`h/lambda_S=0.07`) и установило matched carrier gain `2.114`: переход
`0.387 -> 0.338 мм` изменил gain на `0.994%`, FWHM остался `4 мм`. Relative
gain gate пройден; абсолютная амплитуда в нанометрах остаётся отдельной
receiver/boundary-калибровкой.
[Пакет результатов v0](../results/horn_monotonic_array_3d_v0_242khz/README.md)
содержит карты полей, CSV-метрики, выходы 15 каналов и компактные JLD2-поля
линзы и согласованного референса.

Distributed slow-wave, tapered shunt trim и простой planar global fold не
проходят одновременно transmission--delay--pulse gates. После бинарного
proof-of-concept монотонная невёрнутая TTD-линза дала extra path от `0` на
краях до `9.724 мм` в центре. Endpoint-calibrated пятицикловая апертура дала
`G_peak=2.904`, `B_t=0.899`, `rho=0.99998`, postcursor `0.056`, FWHM `5 мм`.
Fine-mesh quadratic receiver FEM дал `|ux(35,0)|=93.17 нм`, gain `2.493`,
FWHM `5 мм` и боковой лепесток `0.294`: это `+3.68%` амплитуды к лучшей
бинарной линзе. Следующие риски — полный transient общего массива,
continuous-state neighbour correction, mesh convergence и перенос через
согласующий слой в алюминий.

**Gate для линзы v0:** выбранная топология даёт монотонный набор требуемых
невёрнутых задержек, собирается в единую непересекающуюся 15-канальную
геометрию и даёт положительную aperture-to-focus фокусировку. `Tbar`, `Rbar`,
P/S-конверсия и ориентир `peak_P(delay)/peak_P(straight equal-path) >= 0.8`
сохраняются как метрики итеративной доводки, но не блокируют сборку первого
полного массива.

### P4. Библиотека и взаимодействие соседей

- [x] Получить delay states continuation по физическим параметрам, без
      декартова геометрического sweep.
- [x] Построить 15-канальную монотонную геометрию по невёрнутому
      true-time-delay law и проверить максимальный delay-state отдельным
      full-elastic transient против straight equal-path.
- [ ] Для каждого состояния сохранить `S(omega)`, impulse response, CWT,
      максимальные напряжения и производственные размеры.
- [x] Рассчитать полный бинарный локальный triplet-базис в 3D на несущей:
      `DDD`, `DSD`, `DSS/SSD`, `SDD/DDS`, `SDS`, `SSS`; edge contexts
      откалиброваны отдельно.
- [ ] Измерить изменение `S(omega)`, delay, ringing и cross-talk.
- [ ] При сильной связи перейти к суперячейке или глобальному малому массиву.

**Gate:** библиотека сохраняет характеристики в окружении реальных соседей
либо заменена явно коллективной моделью.

### P5. Полный фотополимерный демонстратор

- [x] **P0:** собрать единую 3D-геометрию всех 15 монотонных horn-fed каналов
      с общим выходным фронтом и принимающим блоком; проверить отсутствие
      пересечений и построить сетку. V0 содержит `59 751` узел и `332 095`
      элементов; все 15 источников входят в единый `Domain`.
- [x] **P0:** выполнить один end-to-end carrier solve полной линзы и убедиться,
      что фокус существует. Это первичный интеграционный milestone; gain,
      transmission и sidelobes улучшаются следующими итерациями. На плоскости
      `x=35 мм` получены центральный максимум, FWHM `5 мм`, center/outer
      `5.28`; waist находится около `x=30 мм`, FWHM `3 мм`.
- [x] Рассчитать matched full-3D uniform reference с теми же 15 источниками,
      осевой длиной и принимающим блоком: `|ux|=2.645 нм` против `3.893 нм`
      у монотонной линзы, carrier baseline `1.472`. Convergence P5.1A
      впоследствии показал, что это число ещё не сеточно сертифицировано.
- [x] Собрать редуцированную полосовую модель апертуры с измеренным transient
      `smooth/device` и carrier context coefficients из 3D triplets.
- [x] Восстановить пятицикловый импульс в фокусе и проверить quadratic
      full-elastic receiver FEM на несущей.
- [ ] Построить CWT контекстно-калиброванного фокуса.
- [ ] Выполнить полный transient-FEM только для выбранного массива и
      однородного референса.
- [ ] Проверить потери, `+-0.05/0.10 мм`, свойства материала, контакт и
      ориентацию печати.
- [ ] Откалибровать напечатанный материал на прямых образцах нескольких длин.
- [ ] Печатать линзу только после прохождения gain и pulse-fidelity gates.

### P5.1. Детерминированная предкомпенсация апертуры

**Цель:** поднять честный end-to-end gain полной 15-канальной линзы с
carrier-baseline `1.472` к системному `G_peak >= 2`, не меняя выбранную
связность `log-cosine horn -> smooth sin^4 guide -> point radiator`. Коррекция
должна компенсировать фактическое падение выходной амплитуды к центру
(`0.646` от среднего края) и ошибку true-time delay, а не подгонять только
фазу на `242 кГц`.

**Соглашение о масштабах:** внутри synthesis и convergence все размеры сетки
задаются через `h/lambda_S`, фазовое разрешение — через `k_S h`, геометрия —
через `L/lambda_P` либо `L/lambda_S`, частота — через `f/f0`, задержка — через
`f0 tau` и `2 pi f0 tau`. Первичные графики и таблицы для пользователя
остаются в физических единицах: `мм`, `мкс`, `кГц`, `нм`, `МПа`; gain и
отношения энергий по определению безразмерны. Каждый CSV сначала содержит
физические столбцы, затем служебные normalized columns.

Плотный геометрический sweep запрещён. Разрешены измерение линейного отклика,
один рассчитанный physical step и не более одной Jacobian-коррекции, если
первый step качественно следует предсказанию.

#### A. Сначала отделить физику от ошибки сетки

- [x] Пересчитать v0 lens и matched straight reference минимум на двух
      согласованных сетках. Предпочтительный дополнительный check —
      half-domain `y >= 0` с условием зеркальной симметрии `u_y=0` на
      центральной плоскости. Выполнены три full-domain уровня и три
      symmetry-constrained half-domain уровня.
- [x] Для пар каналов `i` и `16-i` сохранить разность комплексных выходов.
      До синтеза требуются `<= 3%` по амплитуде и `<= 3 deg` по фазе либо
      явная оценка численной неопределённости, если эти пороги дороже разумной
      3D-сетки. Full refined дал `12.81% / 3.79 deg` и не прошёл; half-domain
      навязывает точную зеркальность и использован как диагностический control.
- [x] Принять carrier gain v0 численно устойчивым, только если refinement
      меняет matched gain не более чем на `5%`, а FWHM — не более чем на один
      шаг пространственного scan. Старые уровни gate не проходили, поэтому
      добавлены четыре phase-controlled уровня. Финальный переход тракта
      `0.387 -> 0.338 мм` (`h/lambda_S: 0.08 -> 0.07`) меняет matched gain
      `2.135 -> 2.114`, то есть на `0.994%`; FWHM сохраняется `4 мм`.

**Stop:** геометрическая предкомпенсация не строится из асимметричных v0
выходов. Если half-domain и full-domain не сходятся, сначала исправляется
mesh/boundary/evaluation pipeline.

**Результат P5.1A:** relative-gain stop снят для symmetry-constrained v0.
Финальный уровень даёт `2.226 нм` у линзы против `1.053 нм` у matched straight
reference: target gain `2.114`, plane-peak gain `1.907`, amplitude-equivalent
gain окна `1.653`, FWHM `4 мм`. Между двумя последними уровнями matched gain
меняется на `0.994%`, plane-peak gain — на `0.283%`, window-amplitude gain —
на `2.30%`. Можно переходить к response matrix замороженной симметричной v0.
Абсолютное `|ux|` меняется на `7.59%`, поэтому оно не используется как
абсолютно откалиброванная передаточная функция до отдельного
receiver/boundary refinement. Все новые уровни по-прежнему задаются через
`h/lambda_min` и `k h/p`, а графики остаются в физических единицах.

#### B. Измерить response matrix замороженной v0

- [x] Для одной замороженной постановки решить независимые правые части для
      восьми уникальных по симметрии источников. Сохранить
      `C_ji(omega)` — выход канала `j` от источника `i`, и `g_i(omega)` —
      вклад источника `i` в целевую точку. Выполнено на phase2 отдельно для
      lens и matched straight reference; из-за предела памяти каждый RHS
      решён в свежем процессе. Суперпозиция закрылась с ошибкой
      `2.7e-15 / 1.4e-14`.
- [x] Первый дешёвый этап выполнить на `f0=242 кГц`. Для delay-коррекции
      добавить как минимум `B_-6dB` anchors `193.8` и `290.1 кГц`, развернуть
      фазу и оценить group delay; одна carrier phase не используется как
      true-time-delay observable. Рассчитаны 32 дополнительных FEM-RHS;
      семь из восьми уникальных групп прошли delay gate. Центральная группа
      отвергнута по phase-linearity residual `86.46 deg`, максимум остальных
      равен `19.89 deg`.
- [x] По `C` отделить собственное падение transmission длинных каналов от
      neighbour/cross-talk и общего receiver response. Проверить зеркальность
      response matrix в пределах численного допуска. Reciprocity проверяется
      только для отдельной reverse-задачи с work-conjugate port observables и
      не приписывается автоматически матрице `input face -> preout probe`.
      Carrier off-diagonal L2/diagonal не превышает `0.0643` у lens и
      `0.0622` у matched straight, поэтому cross-talk не доминирует.
      Diagonal-preout и source-to-focus delay observables на семи надёжных
      группах отличаются RMS на `0.607 мкс` и подтверждают один и тот же
      dispersive slope.
- [x] На неизменной геометрии найти симметричные неотрицательные множители
      входа `s_i` constrained least-squares-задачей. Цель — приблизить
      выбранную выходную амплитудную огибающую и максимизировать фокус при
      ограничениях на FWHM и sidelobes, а не обязательно сделать все выходы
      строго равными. Глобальный carrier optimum при `sum(m_i s_i^2)=15`
      найден и независимо проверен на phase4.

Основной proxy сохраняет входной бюджет:

```text
s_i = s_(16-i) >= 0
sum(A_i) = sum(A_i_v0)
same total incident energy for lens and reference
```

Равная сумма площадей является только геометрическим ограничением. При
pressure drive отдельно вычисляется комплексная входная мощность
`1/2 Re integral(t dot conj(v)) dA`; carrier и transient сравнения
перенормируются к одинаковой полной введённой энергии. Результат при
одинаковых `1 MPa` сохраняется как дополнительная source-condition метрика.

Отдельный `available-piezo-area` сценарий разрешён только как системная
абляция. Для него увеличение собранной площади и gain фокусировки
отчитываются раздельно; matched reference получает те же входные площади и
ту же полную падающую энергию.

**Proxy gate:** рассчитанная коррекция должна давать не менее `15%` к
устойчивому v0 gain без ухудшения FWHM выше `6 мм` и sidelobe/peak выше
`0.40`. Целевой запас перед physical step — carrier gain `>= 2.2`; если даже
линейный source-weight proxy не достигает `G=2`, изменение horn mouths не
считается достаточным механизмом.

**Результат amplitude proxy:** gate не пройден. Response matrix предсказала
`+4.45%` к фокусу линзы, FWHM `4 мм`, sidelobe `0.260`; phase4 подтвердила
`2.226 -> 2.326 нм` (`+4.48%`). Но matched straight вырос
`1.053 -> 1.409 нм`, поэтому gain упал `2.114 -> 1.650`. Перенос этих весов
в horn mouths остановлен без physical sweep. Точный boundary-work integral
после LU не поместился в доступную память; текущий proxy сохраняет
номинальный pressure-energy budget, а отдельный power audit остаётся
обязательным до любого изменения площадей.

**Результат broadband delay anchors:** диагональный `lens/straight` transfer
даёт валидную монотонную excess group delay для `y=4.8--33.6 мм` и
эффективную скорость группы `1147.4 м/с`. Минимальный amplitude ratio по
валидным каналам и трём anchors равен `0.599`, поэтому их phase slope не
вычисляется в transmission null. Центр имеет ratio не ниже `0.511`, но
сильно нелинейную фазу и не используется напрямую. Полный полосовой
cross-talk audit матрицы остаётся отдельным подпунктом.

#### C. Перенести амплитудные веса в физические horn mouths

- [ ] Менять входную высоту horn в `z`, сохраняя channel depth `3.2 мм`, шаг
      по `y=4.8 мм`, throat `1.6 мм` и плоский входной торец. Для одинакового
      давления первый прогноз использует
      `A_i/A_0 = h_i/h_0 ~= s_i`; это только начальный Jacobian, а не
      утверждение о точной линейности horn efficiency.
- [ ] Проверить этот Jacobian двумя физическими anchor-моделями: edge/straight
      и maximum-delay channel с рассчитанной центральной mouth height. Не
      строить sweep промежуточных высот.
- [ ] Назначить всем каналам общую horn length, рассчитанную по наибольшему
      expansion ratio из условия `epsilon_ad <= 0.05`. Это сохраняет плоские
      торцы; увеличение общей толщины допустимо, уменьшение guide axial length
      с нарушением curvature gate — нет.
- [ ] Для основного кандидата сохранить `sum(A_i)=sum(A_i_v0)`. Если
      требуемая высота нарушает throat/min-wall/общую толщину, решить
      constrained projection один раз и повторно оценить response-matrix
      proxy.

#### D. Скорректировать невёрнутую задержку

- [x] Из полосового phase slope получить фактические `tau_i`, сравнить их с
      геометрическим `tau_i_target` и вычислить первую поправку
      `delta L_i = c_g * (tau_i_target - tau_i)`. Знак проверить отдельным
      тестом увеличенного прямого пути. Знак согласуется с ранее измеренным
      positive delay максимального transient-state и закрыт exact-delay unit
      test. Для отвергнутого центра выполнено монотонное продолжение от
      ближайшей валидной группы, а не использована нефизическая оценка
      `-0.240 мкс`.
- [x] Повторно обратить arc-length equation для каждой `sin^4` centerline,
      сохранив зеркальность, монотонность задержки, общий выходной фронт и
      `R_inner >= 3.93 мм`. Первая v1 уменьшает maximum extra path
      `9.724 -> 6.957 мм`, maximum bend `13.746 -> 11.374 мм`, сохраняет
      axial length `50.664 мм` и поднимает minimum inner radius до
      `4.917 мм`; monotonic и curvature gates пройдены.
- [x] Не добавлять независимый `2pi` carrier trim. Если phase intercept и
      group delay требуют несовместимых поправок, кандидат возвращается в
      полосовую модель как дисперсионно некорректный. V1 дала maximum output
      phase error `165.1 deg` и `-24.89%` к v0 focus; carrier phase intercept
      и broadband group delay несовместимы для path-only регулятора.
- [x] Максимальный скорректированный delay-state повторно проверить против
      straight equal-path пятицикловым transient. Сохраняются gates
      peak ratio `>= 0.8`, `B_t <= 1.25`, `rho >= 0.90`, postcursor `<= 0.10`
      и transverse energy ratio `<= 0.10`. Получено peak `0.969`, energy
      `0.962`, `B_t=1.000`, `rho=0.99988`, postcursor `0.0925`, transverse
      energy `0.00227`; pulse-fidelity gate проходит. Изолированная delay
      `3.987 мкс` против target `5.777 мкс` подтверждает риск central
      extrapolation.

#### E. Одна физическая линза v1 и matched references

- [x] Собрать одну 15-канальную precompensated v1 с рассчитанными mouths и
      path lengths. До solve автоматически проверить плоские торцы,
      отсутствие пересечений, min wall/gap, `R_inner` и зеркальную симметрию.
      Amplitude mouths не перенесены после failed proxy, поэтому сохранены
      исходные равные mouths; изменены только path lengths. Symmetry,
      monotonicity, mesher и `R_inner=4.917 мм` gates пройдены.
- [x] Построить matched straight v1 reference с теми же пятнадцатью mouths,
      source faces, общей осевой толщиной, receiver и mesh prescription. Для
      сравнения с v0 дополнительно сохранить gain при одинаковой полной
      падающей энергии, а не при одинаковом номинальном давлении и разной
      площади. V1 сохраняет identical mouths и axial thickness, поэтому
      сертифицированный phase3 straight reference совпадает точно и повторно
      использован после config equality check.
- [x] Выполнить carrier solve линзы и reference. Целевые метрики v1:

```text
matched carrier gain >= 2.0
central channel amplitude / mean edge >= 0.85
maximum output phase error <= 10 deg
transverse FWHM <= 6 mm
sidelobe / peak <= 0.40
mirror-pair mismatch <= 3% amplitude and <= 3 deg phase
```

  Phase4 mesh (`184 716` узлов) построена, но direct solve дважды
  воспроизводимо упал в системном macOS `_xzm` allocator после assembly.
  Carrier stop получен на phase3 (`142 904` узла, `h/lambda_S=0.08`): сразу
  пять метрик провалены с большим запасом, поэтому отсутствие phase4 числа не
  может изменить знак verdict и не оправдывает дальнейший memory sweep.

- [x] Только после carrier gate выполнить общий `5-cycle Hann @ 242 кГц`
      transient линзы и matched reference при одинаковой падающей энергии.
      Финальный системный gate остаётся `G_peak >= 2`, `B_t <= 1.25`,
      `rho_pulse >= 0.90`, postcursor `<= 0.10`, без доминирующего бокового
      фокуса. Carrier gate равен `false`, поэтому full transient по правилу
      не авторизован и не запускался.
- [x] Если physical v1 улучшает gain в правильную сторону, но промахивается
      по рассчитанной амплитуде, разрешена одна correction step из измеренного
      physical Jacobian. Если знак или механизм улучшения не совпадает с
      response-matrix prediction, ветка mouth precompensation останавливается
      без sweep и пересматривается outlet/receiver coupling. V1 phase3 дала
      gain `1.604`, FWHM `10 мм`, sidelobe `0.552`, centre/edge `0.711` и
      phase error `165.1 deg`; branch остановлена. Coupling audit показал,
      что причина — phase/group-delay incompatibility, а не cross-talk.

Минимальный воспроизводимый пакет:

```text
run_horn_monotonic_array_3d_convergence.jl
run_horn_monotonic_response_matrix.jl
design_horn_aperture_precompensation.jl
run_horn_precompensated_array_3d.jl
tmp/horn_monotonic_precomp_v1/
```

### P6. Перенос на алюминиевую деталь

- [x] Построить дисперсионную и модовую модель реальной алюминиевой геометрии.
      Для первого thick-body gate выбраны nondispersive bulk P/S ветви
      изотропного Al 6061 в plane strain: `c_P=6122.10 м/с`,
      `c_S=3083.81 м/с`; апертура `121.8 мм`, модуль `28 мм`, фокус
      `60 мм` за интерфейсом. Это не модель Lamb-мод конечной пластины.
- [x] Выбрать преобразователь, контакт и целевую компоненту измерения.
      Выбран longitudinal contact PZT с равномерной `1 МПа` port-тягой,
      perfect bonded planar contact и измерение комплексного `u_x`.
- [x] Тем же synthesis-pipeline пересчитать выбранную топологию. Повторный
      reduced synthesis для `5-cycle Hann @ 242 кГц` оставил Al→Al-профиль
      неизменным. Fine linear FEM даёт `21.044 нм` против `12.884 нм`
      straight; refined quadratic — `22.713` против `16.159 нм`, gain
      `1.4056`. Coarse→refined quadratic: `+0.33%` по амплитуде линзы и
      `-0.48%` по gain.
- [x] Сравнить defect-free и defect-сценарии по локальной энергии, рассеянию и
      улучшению signal-to-noise. Для свободного цилиндрического дефекта
      `r=2 мм` в фокусе local kinetic-energy ratio равен `1.9756`, RMS
      scattered `u_x` на линии за `20 мм` до дефекта — `5.161` против
      `3.583 нм`, equal-noise SNR proxy `+3.17 dB`.
- [x] Не интерпретировать усиление дефектного отклика как прямое измерение
      предела прочности без отдельной калибровочной модели. Текущий результат
      явно помечен как линейный carrier SNR proxy без electronics noise,
      probability-of-detection и strength calibration.

Материальный verdict: Al→Al устойчиво выигрывает. На общей fine-linear
постановке direct polymer даёт `12.580 нм`, обновлённый five-cycle polymer с
номинальным слоем `3.5 мм` — `8.086 нм`, Al→Al — `21.044 нм`. Лучший
фиксированный polymer-профиль со слоем `1.4 мм` достигает `14.095 нм` и не
меняет ранжирование. Полный пакет:
`results/sinusoidal_material_transfer_p6_242khz/`.

Последующий 17-точечный full-FEM broadband gate замороженной Al→Al линзы
закрыт: диапазон `162.2--321.7 кГц` содержит `99.815%` энергии входного
5-cycle pulse, удвоение плотности внутри `B_-6dB` меняет `G_peak` только на
`-0.024%`. Получено `21.65` против `15.88 нм`, `G_peak=1.3634`,
`B_t=0.8643`, `rho=0.9794`, postcursor `0.0644`. Форма импульса проходит,
но амплитудный gate `G_peak>=2` — нет. Идеальный TTD-предел той же апертуры
равен `3.182`, поэтому следующий шаг — восьмигрупповая physical response
matrix и локальный geometry Jacobian, а не изменение апертуры или blind sweep.

Carrier response-matrix/Jacobian gate также завершён. Восемь симметричных
source-базисов точно восстанавливают full drive. Coherence efficiency линзы
равна `0.547`, maximum focal phase error `147.5 deg`; при выравнивании тех же
измеренных амплитуд carrier ceiling равен `2.984`. Pre-output offdiag/diag
`0.557` доказывает нарушение independent-cell approximation. Единственный
gap-step для худшей группы `g_7:1.50->1.75 мм` ухудшил focus на `5.68%`, gain
`1.633->1.541` и phase error `-147.5->-155.0 deg`. Gap-регулятор закрыт:
полезный знак требует выйти ниже min gap и всё равно не даёт нужного phase
range. Следующий topology step должен ввести отдельный broadband phase/delay
регулятор при сохранении прохождения.

#### P6.1. Трёхсостоянийный Al horn--TTD pilot

Текущая sinusoidal-связность больше не продолжается. Следующий физический
pilot переносит уже проверенную стратегию
`collector -> smooth guide -> radiator` в постановку Al→Al и отделяет
throughput от delay:

```text
Al log-cosine horn 7.0 -> 3.5 mm
    -> Al smooth sin^4 guide
    -> малая Al площадка
    -> bonded Al receiver
```

- [x] Зафиксировать воронку по `epsilon_ad <= 0.05` на нижней границе
      `162.2 кГц`: `L_horn=65.4 мм`, area ratio `2.0`.
- [x] Для `F=60 мм`, `y_max=57.4 мм` и `c_P=6122.10 м/с` зафиксировать
      дополнительные пути состояний `0 / 14.65 / 23.03 мм` и maximum delay
      `3.76 мкс`.
- [x] Выбрать первый осевой seed направляющей `132.5 мм` без sweep — как
      безразмерное продолжение polymer pilot по отношению скоростей волн.
- [x] Построить состояния zero/middle/maximum и прямой equal-path control;
      проверить кривизну, minimum inner radius и общую сетку.
- [x] Прогнать `5-cycle Hann @ 242 кГц` при одинаковой падающей энергии.
      Gate: `peak/equal >= 0.8`, `B_t <= 1.25`, `rho >= 0.90`,
      `postcursor <= 0.10`, transverse energy `<= 0.10`.
- [x] Только после прохождения pilot собирать малую Al-апертуру; полный
      15-канальный массив и плотный sweep до этого запрещены.

Уточнённый `130 мкс / 40 samples-per-period` gate максимального состояния
прошёл: peak/equal `1.004` (envelope `0.9965`), energy `0.9964`, `B_t=1.000`,
`rho=0.999972`, postcursor `0.09718`, transverse energy `0.00523`.
Фактическая пакетная скорость узкого Al-направляющего равна примерно
`4.79 км/с`; поэтому aperture paths калибруются по `0.20862 мкс/мм`, а не по
bulk `c_P`.

#### P6.2. Coupled Al-апертура и stop текущего radiator

- [x] Пятиканальный spanning pilot `y=0 / 28.7 / 57.4 мм`: начальный
      equal-pressure gain `1.040`, equal-power `0.948`. Трёхгрупповая
      response matrix дала один path step; он поднял pressure gain до `1.226`
      и реализовал `93.7%` matrix ceiling, но equal-power gain остался `1.069`.
- [x] Увеличить sampling density без sweep: девять каналов, pitch `14.35 мм`.
      Получено pressure gain `1.902`, equal-power gain `1.649`, ideal
      efficiency `0.539`, FWHM `8 мм`, осевой peak — промежуточный gate PASS.
- [x] Проверить полный 15-канальный массив с pitch `8.2 мм`. Получено
      `2.303` против `1.223 нм`, pressure gain `1.883`, equal-power gain
      `1.540`, ideal efficiency `0.447`, FWHM `8 мм` — целевой gate FAIL.
- [x] Восьмигрупповая matrix точно восстановила full field: coherence `0.947`,
      phase-aligned pressure ceiling `1.988 < 2`. Оптимальные неотрицательные
      веса при равной активной мощности дают upper bound `1.834 < 2`.

Verdict: `collector -> narrow sin^4 delay guide -> abrupt small radiator`
доказан как высококачественный broadband delay channel, но остановлен как
амплитудная линза с целью `2x`. Следующий разрешённый pilot сохраняет delay law
и добавляет отдельный плавный output diffuser/radiator `3.5 -> 7.0 мм`.
Сначала снова isolated equal-path/max transient; aperture и sweep до его
throughput gate запрещены. Полный пакет:
`results/aluminium_horn_ttd_p6_242khz/`.

#### P6.3. Double-horn output diffuser

- [x] Изолированный `3.5 -> 7.0 мм`, `65.41 мм` output diffuser прошёл
      refined transient: maximum/equal peak `1.0032`, energy `0.9963`,
      `B_t=1.000`, `rho=0.999952`, postcursor `0.0924`, transverse energy
      `0.00618`. По сравнению с abrupt radiator амплитуда при
      равной source work выросла в `1.2533` раза.
- [x] Полный 15-channel diffuser дал `2.973` против `1.633 нм`
      matched-diffuser control. По отношению к исходному abrupt straight
      equal-power gain равен `1.870`; FWHM `12 мм`, peak `y=0`.
- [x] Восьмигрупповая power matrix восстановила мощность с ошибкой
      `4.6e-13`. Неотрицательные carrier weights дают equal-power upper
      bound `2.0524`; прямой 3D solve воспроизвёл его с точностью
      `1e-13`, осевым peak и FWHM `12 мм`.

#### P6.4. Физическая аподизация и узкополосный импульс

- [x] Перенос carrier weights в площади входных mouths проверен двумя
      крайними physical anchors без sweep. Входы `1.408 / 11.484 мм`
      дали амплитуды `0.461 / 0.883` вместо `0.201 / 1.641` и
      фазовые сдвиги `48.3 / 120.4 deg`. Input-area branch STOP.
- [x] Аналогичные output-height anchors дали `0.580 / 1.037` вместо
      `0.201 / 1.641` и `44.2 / 102.8 deg`. Простой пассивный
      area regulator также STOP: он меняет комплексный импеданс, а не
      только амплитуду.
- [x] Постоянные вещественные веса сегментированного излучателя
      проверены на `193.8 / 242.0 / 290.1 кГц`. Для five-cycle pulse
      broadband optimum только `1.972`: цель `2x` не проходит.
- [x] Поскольку временное/пространственное разрешение вторично,
      выполнен deterministic cycle sweep. Минимальный проход в anchor model
      получен на `30 cycles`; `20 cycles` дают `1.989`.
- [x] Финальный five-frequency full-FEM gate `226.04--257.94 кГц`
      покрывает `99.9485%` энергии 30-cycle Hann и даёт
      `G_peak=2.1261`, `B_t=1.1408`, `rho=0.9805`, postcursor `0.0497`.
      Carrier direct solve: equal-power gain `2.0503`, FWHM `16 мм`, peak
      `y=0`, transverse energy `1.78e-4`.

**Verdict:** цель `2x` достигнута для системы
`double-horn TTD lens + segmented nonnegative drive + 30-cycle Hann`. Это
не чисто пассивная аподизация: без сегментированного drive текущая
геометрия даёт `1.870`. Если следующая цель требует именно
пассивную линзу, нужен отдельный phase-neutral broadband attenuator;
простое изменение площади входа/выхода больше не разрешено.

#### P6.5. Одновходовый пассивный feed

- [x] Синтезировано lossless corporate tree для 15 pulse30-весов:
      `14` делителей, maximum depth `4`, minimum child fraction `0.14446`.
      Для сохранения `G_peak>=2` feed должен передать не менее `0.88493`
      полезной мощности.
- [x] Bulk-P absorber заменён на абсолютный многомодовый DtN-порт. Прямой
      `7 мм` Al-контроль дал `R=9.99e-8`, `T_fund=1.00009`, power balance
      `1.00009`; прежние абсолютные оценки около `0.4` признаны port-artifact.
- [x] Лучший симметричный `50/50` Y (`131 мм`) почти не отражает, но из-за
      возбуждения двух нечётных мод сохраняет только `T_fund=0.94205`.
      Даже идеальное lossless dyadic tree с одним листом на канал имеет
      full-FEM carrier ceiling `1.898 < 2`.
- [x] Неравный узел откалиброван двумя physical anchors и секущей:
      geometric small fraction `0.19558`, junction width `1.369 мм` дают
      amplitude ratio `0.41045` при цели `0.41092`. Beat-компенсированная
      длина `144.478 мм` повышает `T_fund` до `0.96860`, но оставляет
      carrier phase offset `79.81 deg`.
- [x] Альтернативы проверены до разрешённого gate: общая воронка
      `7 -> 121.8 мм` перед 15 каналами передаёт `0.99276` всей мощности,
      но только `0.05280` в полезные fundamental modes, с coherence к
      целевому профилю `0.131`; reduced MMI width sweep `15.2--40 мм`
      даёт ceiling coherence только `0.683`.
- [x] Рассчитанный balanced tree использует measured `0.94205` для двух
      equal nodes, `0.96860` для остальных двенадцати nodes и идеальную phase
      compensation. Получено total useful efficiency `0.87669 < 0.88493`,
      power-only pulse ceiling `1.99067` и routed carrier response-matrix gain
      `1.80498`. Broadband solve именно этого дерева не запускался.
- [x] Проведён topology audit того же optimistic cascade surrogate. Лучшее
      дерево с непрерывными группами каналов имеет efficiency `0.89056` и
      power-only ceiling `2.00635`; unrestricted binary tree — `0.89160` и
      `2.00752`. Это отменяет глобальный no-go, но не является physical PASS:
      разные split ratios не проверены FEM, фазы и маршрутизация идеализированы,
      а амплитудный запас меньше `0.4%`.

**Verdict:** реализованный balanced feed и проверенные Y/manifold/MMI
топологии остановлены, но весь класс одновходовых пассивных feed не закрыт.
Стратегия `воронка -> волновод -> диффузор` остаётся рабочей delay/throughput
ячейкой; риск возникает при многократном распределении мощности в присутствии
изгибных мод. Blind feed sweep запрещён. Возврат к дереву разрешён только при
новом физическом механизме подавления mode conversion или после появления
достаточного intrinsic gain margin.

#### P6.6. Time-reversal-matched piezo source

Следующая ветка меняет объект исследования с механического делителя мощности
на согласованный собственный канал всей системы:

```text
один общий электрический waveform
    -> пространственно аподизированный пьезослой
    -> одинаковые input collectors
    -> smooth TTD guides
    -> одинаковые output diffusers
    -> фокус в алюминии
```

Амплитудный профиль задаётся локальной эффективностью электромеханического
преобразования — например, долей PZT в `1-3`-композите, локальной
поляризацией или активной площадью общего электрода. Скорость и задержка не
должны кодироваться в пьезослое: их по-прежнему задаёт пассивная TTD-геометрия.

1. [ ] Поставить виртуальный продольный источник в фокус и по взаимности
   back-propagate поле до плоскости пьезослоя во всей полосе импульса.
2. [ ] Посчитать два upper bounds при одинаковой энергии: произвольный
   комплексный matched source и его проекцию на общий waveform с
   вещественным неотрицательным пространственным профилем.
3. [ ] Проверить без геометрических sweep текущую линзу для `5-cycle` и
   `30-cycle Hann`. Если constrained projection не даёт `G_peak >= 2`,
   topology/source concept останавливается.
4. [ ] Только после прохода заменить prescribed traction на гомогенизированный
   coupled piezoelectric layer с общим напряжением и пространственным
   `e_33(y)` или fill fraction. Проверить электрическую энергию, фазу,
   lateral coupling и чувствительность к дискретизации профиля.
5. [ ] Gate физического концепта: один электрический режим возбуждает нужную
   механическую аподизацию, `G_peak >= 2` при равной электрической энергии,
   а pulse-fidelity остаётся в текущих пределах.

## 9. Ближайший расчётный пакет

Текущий этап завершён. Доказаны три разных результата, которые нельзя
смешивать:

1. Пассивный all-Al double-horn TTD-массив без аподизации даёт `G=1.870`
   относительно исходного abrupt reference (`1.758` относительно
   matched-diffuser carrier reference).
2. Та же линза с сегментированным неотрицательным drive и `30-cycle Hann`
   проходит системный gate: `G_peak=2.1261`; matched-diffuser carrier gain
   равен `1.9276`.
3. Проверенный balanced passive feed не проходит, однако topology audit
   оставляет лишь теоретический power-only запас `2.006--2.008`; поэтому это
   открытый физический риск, а не доказанный глобальный no-go.

Следующий разрешённый расчёт — только P6.6: reciprocity/time-reversal upper
bound для общего waveform и пространственно аподизированного источника.
Новые механические sweeps, ещё один вариант Y-tree и ручная подстройка
длин/площадей до этого gate запрещены. Полная история старых пакетов остаётся
в [`research_log.md`](research_log.md), а воспроизводимый итог P6 — в
[`results/aluminium_horn_ttd_p6_242khz/README.md`](../results/aluminium_horn_ttd_p6_242khz/README.md).

## 10. Критерии остановки и запреты

- не оптимизировать только фазу на `f0`;
- не считать фазы, различающиеся на `2pi`, эквивалентными для импульса;
- не максимизировать резонансный `Q` или групповую задержку без ограничения на
  transmission, bandwidth и ringing;
- не интерпретировать фазовую производную при малом прохождении как время
  прихода пакета;
- не запускать плотный sweep до физического якобиана и проверки его ранга;
- не считать optimistic path-delay proxy доказательством прохождения без FEM
  изгибов/переходов и абсолютного `T/R`;
- не использовать rank-one port вне доказанной полосы и симметрии;
- не считать изолированные элементы независимыми без димерного теста;
- не называть одну миллиметровую геометрию универсальной;
- не прорезать диагностическую мета-линзу в исследуемой алюминиевой детали;
- не переходить к нелинейности, surrogate или topology optimization до
  верификации линейной полосовой прямой задачи.

## 11. Формат каждого нового результата

Каждая серия должна сохранять:

- `MaterialConfig`, `ManufacturingConfig`, `WaveConfig` и версию топологии;
- фактическую сетку, порядок FEM и проверку сходимости;
- полный комплексный `S(omega)` в существенной полосе;
- энергетический баланс и модовую конверсию;
- входной, референсный и выходной импульсы при одинаковой падающей энергии;
- огибающие, CWT и доверительные маски;
- `G_peak`, энергию фокального окна, FWHM, DOF и боковые лепестки;
- `B_t`, `rho_pulse`, postcursor ratio, delay и GDD;
- максимальные напряжения и минимальные печатные размеры;
- короткий verdict по текущему gate и следующий разрешённый шаг;
- ссылки на код и каталог воспроизводимых артефактов.

Подробные промежуточные кривые и история не переносятся обратно в этот файл:
они добавляются в [`research_log.md`](research_log.md).
