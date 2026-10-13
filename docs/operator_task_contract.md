# O0: безразмерный task contract операторного pilot

> Этот контракт не задаёт стенд, частоту или материал прототипа. Он фиксирует
> первую сравнительную задачу в единицах длины P-волны и предназначен только
> для gates `O1--O2` из
> [`measurement_operator_design_space.md`](measurement_operator_design_space.md).

## Цель pilot

Проверить, сколько независимых измерительных мод требуется для:

1. обнаружения слабого дефекта в заранее известной малой зоне;
2. локализации слабого дефекта в ограниченной двумерной ROI;
3. различения зеркально расположенных дефектов;
4. оценки выигрыша от собственной P/S-физики относительно scalar-P модели.

Sizing, classification и сертификационный POD в первый pilot не входят.

## Безразмерная физика

Единицы выбираются как

```math
\lambda_{P,0}=1,\qquad c_P=1,\qquad f_0=1.
```

Первый representative elastic ratio равен

```math
c_S/c_P=0.5.
```

Он близок к металлам с коэффициентом Пуассона около одной трети, но не
обозначает конкретный сплав. Частоты задаются через `f/f0`; все выводы pilot
относятся к dimensionless operator, а не к килогерцам или миллиметрам.

## Виртуальная апертура и ROI

Используется линейная апертура длиной `8 lambda_P0` из нечётного числа
виртуальных элементов. Нечётность сохраняет точный центральный элемент и
действие группы отражения `C_s`.

Первая ROI:

```text
x/lambda_P0 = 2.5 ... 6.0
y/lambda_P0 = -2.0 ... 2.0
```

Это ограниченная область, а не попытка заменить PAUT во всём полупространстве.
Геометрия апертуры и ROI может масштабироваться вместе; первичны отношения
`D/lambda`, `x/lambda`, `y/lambda` и `B/f0`.

## Дефектный словарь

На `O1--O2` используется single-scattering/Born surrogate. Дефект задаётся
симметричным tensor `Q`:

```math
H_{rt}(\omega;\xi)
=e_x^T G(r_r,\xi,\omega)\,Q_\xi\,
G(\xi,r_t,\omega)e_x.
```

В первой localization-задаче `Q=I`. Затем добавляются ориентированные
анизотропные `Q(theta)`. Это словарь слабых point-scatterer hypotheses, а не
точная модель свободной полости или контактной трещины.

Green tensor раскладывается на P/S-проекторы. Поэтому mode conversion и
угловые диаграммы присутствуют на уровне линейного surrogate, но near-field,
finite-defect resonance и нелинейный crack contact отсутствуют.

## Частота и noise level

Используется несколько общих anchors внутри

```text
f/f0 = 0.8 ... 1.2.
```

На `N0` каждый сохранённый complex feature имеет единичную white-noise
variance. Это algebra/debug ceiling. Частотные anchors нельзя интерпретировать
как независимые физические измерения до введения общей pulse covariance на
`N1--N3`.

## Нормировка ресурсов

- Пространственный вес каждого порта имеет `sum(abs2(w)) = 1`.
- Для `Q` последовательных transmit-кодов полная падающая энергия фиксируется;
  каждый код получает множитель `1/sqrt(Q)`.
- Receive outputs не перенормируются по собственной выходной энергии.
- Число ADC, transmit events, frequency anchors и механических состояний
  сохраняется раздельно.

`Ideal-K` является линейным compression bound полного словаря. Он не считается
физически реализуемым Tx/Rx-кодом и не участвует в BOM-сравнении.

## Обязательные операторы

```text
a_scan_uniform
fixed_focus_reciprocal
coded_even_single_port
coded_asymmetric_single_port
few_port_2
few_port_4
ideal_k_code
paut_fmc_oracle
```

`coded_even` обязан показать точную зеркальную неоднозначность. Это positive
control для coherence metric. `coded_asymmetric` проверяет, сколько информации
добавляет нарушение `C_s`, сохраняя трансляционный 2D→3D lift. `few_port_2`
использует even/odd symmetry sectors и проверяет альтернативу нарушению
симметрии.

## Gates

### `O1 pass`

- reciprocity error полного response tensor не превышает численный допуск;
- mirror equivariance выполняется для `Q=I`;
- PAUT/FMC dictionary имеет numerical rank больше единицы;
- mirrored hypotheses не обязаны быть различимы в symmetric one-port control.

### `O2 result`

Строятся:

- singular-value spectrum полного оператора;
- retained dictionary energy для `ideal K`;
- minimum/median detectability;
- maximum pair coherence и худшая пара дефектов;
- stable rank и information log-det;
- ресурсная таблица каналов и transmit events.

`O2` не имеет заранее назначенного победителя. Его результат — compression и
architecture frontiers, после которых определяется состав более дорогого
pulse-echo/clutter campaign `O3`.

## Ограничение вывода

Положительный результат означает только, что архитектурный механизм способен
создать информативный оператор в reciprocal elastic-wave surrogate. Он не
разрешает синтез геометрии до проверки finite bandwidth, contact variation,
ringing и paired-defect pulse-echo на `O3`.
