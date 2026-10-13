# Контракт O3 F1c: structural clutter и PAUT positive control

> Зафиксировано 15 августа 2026 года до построения первой F1c scene. Этот
> уровень исправляет task-prior failure F1b; он не меняет frozen operators и
> не разрешает геометрию до повторного O3 gate.

## 1. Почему F1c обязателен

В F1b шесть изолированных дефектов хорошо разделялись по времени, а uniform
A-scan получил больший coherent SNR, чем sequential FMC. Поэтому не
существовал положительный разрыв `A-scan -> FMC`, относительно которого можно
оценивать дешёвое приближение к PAUT.

F1c сначала строит задачу с обязательной пространственной неоднозначностью:

```text
одинаковый target radius       1.0 mm
одинаковая target depth        60.0 mm
lateral target positions       -28, -14, 0, 14, 28 mm
background per context         5 actual traction-free circular voids
background radii               0.6--1.1 mm
background ROI                 x=25--105 mm, y=-27--27 mm
```

Target classes имеют одинаковую форму и почти одинаковое время прихода.
Отличить их можно прежде всего по пространственной структуре aperture
response, а не по TOF или размеру.

## 2. Actual full-wave scene family

Каждый structural context строится как одна conforming master-триангуляция,
содержащая пять target subdomains и пять background subdomains. Из неё без
ремешинга записываются:

```text
H0      background voids removed; all target subdomains filled with aluminium
H1_j    same background voids + target j removed; other targets filled
```

Таким образом H0 и пять H1 имеют идентичную внешнюю и background mesh,
а paired difference не содержит independent mesh jitter. Background holes
являются настоящими свободными границами FEM, поэтому multiple scattering и
backwall interaction решаются direct transient, а не складываются из Born
atoms.

Physics/time discretization сохраняется из прошедшего F1a:

```text
2D plane strain isotropic elasticity
common h_bulk/h_local = 1.0/0.25 mm
five-cycle Hann 242 kHz burst
30 samples per cycle, record 70 us
17 x 17 multi-RHS response
free backwall, P/S impedance input/lateral boundaries.
```

## 3. Frozen context split

```text
base seed              0x051c2026
train contexts         2
validation contexts    1
test contexts          1
target classes         5 per context
```

Background positions/radii выводятся только из context seed с rejection по
границам, target disks и минимальному separation. Test context не используется
для выбора projector, templates, threshold или noise scale.

Четыре context — physics pilot, а не POD sample size. Они достаточны для
проверки того, появилась ли в задаче пространственная information advantage;
после positive control stochastic counts расширяются без изменения сцен.

## 4. Positive-control gate до сравнения кандидатов

Сначала сравниваются только uniform A-scan и FMC на одинаковой полной
incident energy, одном electronics-noise-per-observation и train-only
background calibration.

F1c task prior считается валидным только при одновременном выполнении:

```text
FMC template stable rank                         >= 3.0
FMC maximum off-diagonal template coherence      <= 0.80
FMC alpha90                                      <= 0.80 x A-scan alpha90
FMC Pdl at reference alpha                       >= 0.90
test Pfa upper 95% for FMC                       <= 1.25%
```

Reference alpha — первая validation point, где FMC `Pd>=0.95`. Localization
в positive-control pilot требует правильного lateral class; tolerance по
координате дополнительно репортируется, но не заменяет exact-class gate.

Если positive control не проходит, static/motion frontier не интерпретируется:
сначала усиливается structural task prior. Если проходит, те же frozen
`K4,Q2`, `K3,Q4`, `P12`, `P17` оцениваются без переоптимизации.

## 5. Условие выхода из O3

После positive control кандидат должен сохранить прежние gates:

```text
Pfa upper 95%                     <= 1.25%
Pd(reference alpha)               >= 0.90
exact-class Pdl                   >= 0.80
minimum target-class Pd           >= 0.75
A-scan -> FMC alpha90 closure      >= 0.50
```

Только такой результат открывает O4. Иначе проект возвращается к operator
class/task prior; новая линза или метаматериал по-прежнему запрещены.

