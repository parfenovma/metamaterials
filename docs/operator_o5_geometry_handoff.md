# O5 geometry handoff: трёхпортовый гладкий модовый преобразователь

> Финальное решение зафиксировано 15 августа 2026 года после expanded F1d:
> пять development/train, два fresh validation и два fresh final-test
> full-wave structural contexts.

## 1. Переданный оператор

В геометрию передаётся не исходный случайный `K4,Q2`, а прошедший
geometry-aware оператор:

```text
architecture       K3,Q1
electronic ports   3
passive states     1
transmit events    3
aperture modes     DCT 0, 1, 2
symmetry sectors   even, odd, even
band               193.6--290.4 kHz
lift               translation_kz, затем finite-width 3D check
```

Полный `H(omega)`, коэффициенты и допуски находятся в
[`o4_k3q1_transfer_target_v0`](../results/operator_space/o4_k3q1_transfer_target_v0/README.md).
Семнадцать строк файла коэффициентов — численная дискретизация непрерывной
апертуры, а не 17 независимых pulser/ADC.

## 2. Почему handoff изменился после первого O4

Первый F1c pilot выбрал unconstrained `K4,Q2`, но его карты требовали 11--17
пространственных мод, 5--9 sign changes и 116 Stiefel-параметров. Smooth
projection показала, что информация находится не в этой случайной сложности:
уже первые три DCT-моды сохраняют нужную detection/localization.

Fresh O5a test выявил недостаток двух train-контекстов по Pfa. После
расширения calibration исходный `K4,Q2` не прошёл `Pfa_upper=1.361%`, тогда
как гладкие архитектуры прошли:

| Architecture | Pfa upper | Pd@3 | exact Pdl@3 | Min Pd | alpha90 | Gap closure |
|---|---:|---:|---:|---:|---:|---:|
| smooth `K3,Q1` | **0.089%** | 0.982 | 0.976 | 0.957 | 2.427 | 0.910 |
| smooth `K4,Q1` | 0.312% | **0.993** | **0.983** | **0.988** | **1.911** | **0.974** |
| unconstrained `K4,Q2` | 1.361% | 0.980 | 0.979 | 0.947 | 2.120 | 0.948 |
| FMC oracle | 1.026% | 0.989 | 0.983 | 0.980 | 1.706 | 1.000 |

Выбран `K3,Q1` по заранее заданному resource priority. `K4,Q1` остаётся
robust fallback, если реальная трёхпортовая реализация не выдержит O4.

## 3. Выбранный topology class

Основной класс proof-of-physics:

```text
три фиксированных reciprocal modal couplers
    -> один монолитный нерезонансный elastic manifold / wedge
    -> непрерывная плоская контактная апертура.
```

Порты возбуждают:

1. uniform even mode;
2. first odd cosine mode;
3. second even cosine mode.

Это могут быть три распределённых piezoelectric electrode/polarization
pattern или три механических входных модовых порта. Между элементами нет
управляемых фазовращателей, состояния не переключаются, а signed weights
фиксированы симметрией мод. Поэтому система не является малоканальной PAUT
решёткой, замаскированной под пассивную сеть.

Первая геометрия живёт в 2D plane strain и экструдируется по `kz=0`.
Reflection sectors `even/odd/even` дают аналитические symmetry checks.
Finite-width 3D нужен только после совпадения измеренного 2D `H(omega)` с O4.

## 4. Отвергнутые реализации

| Класс | Решение | Причина |
|---|---|---|
| Одиночная focusing lens | reject | rank-1 порт не сохраняет lateral information advantage |
| 17-channel array или активная phase network | reject | нарушает исходный cost/architecture invariant |
| Дискретная 3-to-17 пассивная beamformer network | fallback only | реализуема, но прячет fabrication/calibration complexity в couplers и wiring |
| Резонансный metamaterial / длинный labyrinth | reject first | лишняя dispersion и loss при пятицикловой полосе |
| Монолитный low-order mode transformer | selected | напрямую реализует три гладких group modes с одной геометрией |

## 5. Геометрические степени свободы первого прототипа

Разрешены только верхнеуровневые параметры, влияющие на измеренный transfer:

```text
общая длина и толщина manifold
три входные port widths/positions
два symmetry-preserving taper envelopes
распределение electrode/coupling strength для трёх mode shapes
один общий contact-layer impedance.
```

Не разрешены per-cell phase variables, 17 независимых weights, resonator
library или blind topology optimization. Objective — расстояние между
измеренным и target `H(omega)`, включая loss/reciprocity, а не focal peak.

## 6. O5 acceptance gates

До defect simulation геометрия должна пройти empty-load multiport test:

```text
normalized H fidelity              >= 0.98
maximum port coherence             <= 0.15
column imbalance                   <= 0.75 dB
coefficient amplitude error RMS    <= 3%
differential delay RMS             <= 82.6 ns
reciprocity residual               <= 3%
sigma_max(H)                       <= 1
```

Без drive/noise compensation one-way useful amplitude должна быть `>=0.8994`
(`<=0.921 dB` loss). Большая потеря допустима только при явно оплаченном
pulser headroom и сохранении output-referred noise budget.

После transfer gate as-built calibrated `H(omega)` подставляется в frozen F1d
replay. Требуются прежние `Pfa/Pd/Pdl/min-class/gap` gates. Лишь затем
разрешаются prototype machining и finite-width 3D.

## 7. Состояние handoff

O5 архитектурно завершён: выбран один topology class, одна symmetry group,
минимальный набор параметров и hard acceptance contract. Следующий физический
шаг — **один** mode-transformer prototype/digital twin по этому контракту, а
не новый обзор геометрий.
