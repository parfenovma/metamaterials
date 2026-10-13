# Контракт переноса O3 в full-wave

> Зафиксировано 15 августа 2026 года до первого расчёта нового FEM-класса.
> Этот контракт не открывает синтез геометрии и не заменяет paired transient
> gate из [`operator_o3_benchmark_contract.md`](operator_o3_benchmark_contract.md).

## 1. Зачем нужен промежуточный уровень F0

Reduced O3 сравнивал архитектуры на общем Born-операторе. Следующий дорогой
уровень должен проверить, что вывод не был создан точечным дефектом,
асимптотическим Green tensor или идеальной границей. При этом запускать полный
transient ensemble до проверки самой FEM-постановки нельзя: ошибка в разметке
апертуры или backwall сделает весь ensemble бесполезным.

Поэтому full-wave O3 делится на два вложенных уровня:

```text
F0  paired carrier operator
    -> проверяет физику, границы, reciprocity и перенос frozen operators;

F1  paired broadband/transient ensemble
    -> проверяет pulse echo, ringing, clutter и исходные detection gates.
```

Положительный F0 разрешает F1, но сам по себе не считается успешным O3.

## 2. Общая физическая модель F0

```text
dimension                    2D plane strain, full isotropic elasticity
material                     nominal 6061 aluminium
carrier                      242 kHz
specimen                     125.0 x 71.6 mm
aperture                     full left face, 17 equal boundary elements
top/bottom                   first-order P/S absorbing boundary
input face                   same P/S absorbing boundary + applied traction
backwall                     free reflecting boundary
defect                       free circular void, radius 2.0 mm
defect centre                x=60.0 mm, y=0.0 mm
receive observable           integrated normal displacement per element
secondary observable         integrated tangential displacement per element
mesh                         h_bulk=1.5 mm, h_defect=0.35 mm
```

Одна комплексная FEM-матрица факторизуется один раз на mesh/frequency. Затем
решаются 17 правых частей — по одной равномерной нормальной traction-нагрузке
на каждый элемент. Это даёт физическую матрицу `R[receive, transmit]`. Все
архитектуры собираются из неё линейно; ни одна не получает отдельную сетку,
материал или receive chain.

Дальний торец намеренно не входит в absorbing boundary. Иначе из постановки
исчезает обязательный backwall echo.

## 3. Frozen representatives

F0 не оптимизирует codebook по FEM-результату. Из свежих reduced campaigns
переносятся ровно четыре заранее выбранных представителя:

| ID | Ресурс | Роль |
|---|---:|---|
| `trained_k4_q2_unconstrained` | `K=4, Q=2, E=8` | minimum-event static pass |
| `trained_k3_q4_reflection_balanced` | `K=3, Q=4, E=12` | best static pass |
| `motion_compact_center_p12_native` | `P=12, E=12` | first valid motion crossover |
| `motion_uniform_full_p17_a8` | `P=17, E=17` | best low-channel motion point |

Дополнительные comparators: uniform A-scan и `17×17` FMC. Статические
матрицы загружаются без изменения из `static_frontier_v0`; motion indices
равны соответственно `3:14` и `1:17`.

## 4. Нормировка операторов

Для статического состояния с вещественной матрицей портов `W_q`:

```math
y_q = W_q^T R W_q / \sqrt{KQ}.
```

Для motion берётся диагональ `R` на выбранных позициях и делится на
`sqrt(P)`. FMC делится на `sqrt(17)`. Uniform A-scan использует единичный
нормированный пространственный порт. Таким образом во всех классах полная
incident-energy convention совпадает с reduced O3.

Paired scattering operator определяется как

```math
Delta R = R_void - R_free.
```

Это вычитание используется только как physics diagnostic F0. Оно не является
доступным detector input в F1 и не отменяет запрет на бесплатный reference.

## 5. F0 gates

До сравнения архитектур одновременно должны выполняться:

```text
all 17 source/receiver tags present                 yes
free backwall is distinct from absorbing boundary  yes
maximum relative reciprocity residual               <= 1e-6
maximum centred-case mirror residual                <= 3%
paired normal scattering norm                       finite and > 0
tangential receive response                         finite and reported
same mesh/material/receiver for all architectures   yes
```

Architecture metrics F0 — captured paired-scattering energy, raw
defect/free contrast и их отношения к A-scan/FMC. Они описывают перенос
оператора, но не `Pd`, `Pfa`, POD или промышленную чувствительность.

## 6. Условие перехода F0 -> F1

При пройденных physics gates выполняется широкополосная проверка на общей
частотной сетке, после чего из неё строится пятицикловый pulse-echo response.
Минимальный F1 обязан добавить mesh-convergence, contact/material variations,
ringing и независимые no-defect trials. Только F1 возвращается к исходным
O3 gates (`Pfa`, `Pd`, `Pdl`, `alpha90`, localization и resource vector).

Если F0 меняет порядок static/motion представителей, они всё равно не
переоптимизируются: F1 получает те же frozen четыре точки, а изменение порядка
фиксируется как результат смены физического уровня.

## 7. Результат numerical remediation

Первая независимая номинальная пара сеток провалила только paired mirror
gate: при вычитании двух отдельно триангулированных carrier-матриц усилился
mesh jitter. Gate не ослаблялся. Была введена общая conforming-триангуляция:
в free case круглая подобласть заполнена тем же алюминием, а в void case
удаляются только её элементы. Convergence-run при `h_bulk=1.0 mm`,
`h_defect=0.25 mm` дал paired mirror residual `1.805%` и прошёл исходный
порог `3%`. Полный audit trail сохранён в
`fullwave_carrier_convergence_v0`.
