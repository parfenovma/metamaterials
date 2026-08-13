# Concentric composite rings: topology benchmark

## Verdict

`conditional`

Семь концентрических колец из W--эпоксидного композита дают требуемый
профиль задержки в изготовимой геометрии и устойчиво выигрывают у
`straight_footprint`. Однако выигрыш относительно `bare` сохраняется только
при затухании core не выше примерно **0.62 dB/cm на 242 kHz** при контактном
слое `0.10 mm`. Прямая линейная экстраполяция опубликованных измерений при
`1 MHz` даёт около `2.2 dB/cm`; в этом source-informed сценарии
point-scatterer echo proxy равен только `0.352 x bare`.

Следующий расходующий вычисления gate -- не более подробный FEM, а измерение
купонов при `162--322 kHz`. Если для рецептур `10--20 vol% W` получено
`alpha_242 <= 0.60 dB/cm`, ветку имеет смысл продолжать в общем paired-defect
FEM. При `alpha_242 >= 0.8 dB/cm` bulk DP190/W-кольца следует остановить и
перейти к низкопотерянному термопластическому backup.

Кодовый commit, на котором получены приведённые числа:
`fa32a3a5a2f8daf4b01a19d14ea97e8ab43efaef`. Общий benchmark взят из
`775244a2424508d890fea7012d50fa655712258f`.

## Общий стенд и честный уровень модели

Контракт `docs/defect_topology_benchmark.md` не менялся:

- Al 6061 surrogate, `cP=6122.1024 m/s`, `rho=2700 kg/m3`;
- `242 kHz`, `5-cycle Hann`, aperture `71.6 mm`;
- focus/defect `60 mm` от интерфейса;
- обязательные контроли `bare` и `straight_footprint`;
- максимальная длина устройства `55 mm`.

Реально рассчитаны:

1. normal-incidence complex transfer matrix каждого кольца, включая оба
   face sheet, контактные слои, `S21`, `S11`, потери и group delay;
2. 2-D scalar Huygens aperture surrogate с абсолютными потерями относительно
   того же unit incident aperture у `bare`;
3. broadband `5-cycle Hann` pulse;
4. строго помеченный reciprocal point-scatterer proxy `H(f)^2`, используемый
   только как дешёвый G3 gate.

Не рассчитаны full-elastic P/S conversion, рассеяние на радиальных швах и
paired monostatic `none`/`void_r2_focus` FEM. Поэтому common echo/CNR/POD
метрики оставлены `NaN`; point-scatterer proxy не переименован в CNR или POD.

## G0 -- выбранная конструкция

### Основная архитектура

Плоский диск состоит из:

- core `18.0 mm`: семь равных по радиальной ширине колец;
- два общих face sheet по `1.65 mm` из измеренного `20.22 vol% W` состава;
- два явно моделируемых контактных слоя DP190 по `0.10 mm`;
- суммарная осевая длина `21.5 mm`;
- минимальная ширина кольца `5.114 mm`.

Профиль получается из

```math
\tau(r)-\tau(r_o)=
\frac{\sqrt{F^2+r_o^2}-\sqrt{F^2+r^2}}{c_{Al}},
\qquad
c_i=\frac{L}{L/c_o+\Delta\tau_i}.
```

Внешнее кольцо использует `9.92 vol% W`; к центру доля W возрастает до
`19.48 vol%`. Полученный диапазон:

| Ring, center to edge | `cP`, m/s | `Z`, MRayl | W, vol% | delay, us |
|---:|---:|---:|---:|---:|
| 1 | 1718 | 8.265 | 19.480 | 1.393 |
| 2 | 1730 | 8.184 | 19.061 | 1.323 |
| 3 | 1754 | 8.011 | 18.204 | 1.183 |
| 4 | 1790 | 7.733 | 16.898 | 0.976 |
| 5 | 1839 | 7.324 | 15.113 | 0.707 |
| 6 | 1902 | 6.747 | 12.807 | 0.380 |
| 7 | 1982 | 5.946 | 9.920 | 0.000 |

Al 6061 имеет `Z=16.53 MRayl`. Без matching одиночная граница пропускает по
нормали примерно `0.78--0.89` мощности. Общие `20.22 vol% W` face sheet
(`cP=1698 m/s`, `Z=8.41 MRayl`) выбраны как reproducible zero-porosity anchor,
а не более импедансный `31.52 vol%` образец: в опубликованном последнем
образце уже измерено `10.8%` пористости.

### Маршрут изготовления

1. Использовать 3M DP190 и crushed W grit `~250 um`; для точного соответствия
   измеренным autoclave anchors добавить `1 vol%` стеклянных шариков `300 um`.
2. Отмерить семь рецептур по объёму. Это соответствует приблизительно
   `65.6--80.7 wt% W`, поэтому обычное ручное смешивание без pressure cure
   неприемлемо.
3. Отлить семь annular blanks в разъёмных силиконовых/металлических формах,
   провести vacuum degassing и cure `4 bar, 29 degC, 4 h`, как у primary
   property source.
4. Проточить OD/ID колец, вложить их в общий fixture; радиальные швы заполнить
   тем же соседним композитом с контролем `<=0.10 mm`, повторно pressure-cure.
5. Отлить/обработать два общих face sheet `20.22 vol% W`, соединить и совместно
   притереть торцы. Контролировать параллельность, пустоты (X-ray/UT) и
   фактическую толщину.
6. До сборки полной линзы изготовить цилиндрические coupons каждой рецептуры
   для измерения `cP(f)`, `alpha(f)` и `Z` на `162--322 kHz`.

Оценка диска: core `72.5 cm3`, face sheets `13.3 cm3`, около `341 g` общей
массы и `~0.25 kg` W. BOM-класс -- `medium`: материалы доступны, но нужны
точные кольцевые формы, токарная/притирочная операция и pressure pot/autoclave.

### Единственный backup

Низкопотерянный nested thermoplastic stack: PEI/PEEK в центре и
коротковолокнистый CFR-PEEK с поперечной ориентацией к краю, core около
`51 mm` плюс тонкие W--эпоксидные matching skins. Для PEEK опубликованы
`cL=2536 m/s`, `Z=3.26 MRayl`, `BUA=1.33 dB/(cm MHz)`; для transverse CFR-PEEK
`cL=2693 m/s`, `Z=3.59 MRayl`, `BUA=1.45 dB/(cm MHz)`. Это даёт намного лучше
подтверждённый low-frequency loss, но требуемая скорость центрального кольца
около `2508 m/s` лежит чуть ниже neat PEEK, а непрерывная серия рецептур пока
не подтверждена. Поэтому backup не выдаётся за готовый дизайн: сначала нужны
coupons PEI/PEEK/CFR-PEEK и проверка bond-line scattering.

## Материальный отбор

- **DP190/W, primary:** Powell et al. измерили density, `cL`, `Z`, attenuation,
  porosity и bond strength для `0--35 vol% W`; pressure cure устранил поры у
  `9.92` и `20.22 vol%` и дал именно нужные `cL=1982--1698 m/s`,
  `Z=5.99--8.41 MRayl`. DOI
  [10.1177/07316844261436540](https://doi.org/10.1177/07316844261436540).
- **Независимый W/epoxy loss anchor:** Franco et al. при `1 MHz` и
  `10--20 vol% W` показывают по Figure 5 примерно `8--10 dB/cm`; скорость,
  плотность и impedance также измерены. [Открытый PDF](https://www.abcm.org.br/anais/cobem/2005/PDF/COBEM2005-2186.pdf).
- **Process sensitivity:** ball-milled W/epoxy показал изменение
  `c=2.2--1.4 km/s`, `Z=8.31--2.30 MRayl` и рост damping вместе с загрузкой из-за
  porosity/agglomeration; это причина не переносить свойства между рецептами.
  DOI [10.1016/j.rinma.2023.100503](https://doi.org/10.1016/j.rinma.2023.100503).
- **PEEK backup:** Fitch, Hoffmeister & de Ana дали скорости, impedance и BUA
  neat/CFR-PEEK. DOI
  [10.1007/s10853-010-4428-1](https://doi.org/10.1007/s10853-010-4428-1).
- **Alumina/silica epoxy:** доступные диапазоны `c~2.65--2.8 km/s` и
  `Z~3--5.5 MRayl` дают слишком малый delay span при `L<=55 mm` и худшее
  согласование с Al. Измеренный alumina/EPO-TEK пример: DOI
  [10.1109/TUFFC.2011.1825](https://doi.org/10.1109/TUFFC.2011.1825).
- **Hollow microspheres, foams, rubber/polyurethane:** скорость понижается
  легко, но одновременно падает `Z`, растут scattering и attenuation. Они
  полезны как backing, а не как power-preserving transmission lens.
- **Сплошные стекло/керамика:** слишком быстрые; отдельные кольца дают резкие
  `Z`-ступени и не решают центральную положительную задержку без второго
  медленного материала.

## G1 -- local transfer и pulse

Потери заданы не константой, а

```math
\alpha(f)=\alpha_0\left(f/242\,\mathrm{kHz}\right)^n.
```

Sweep содержит `alpha_0=0--3 dB/cm`, `n=1` и optimistic `n=2`, а также
perfect co-cure, `0.10 mm` и `0.25 mm` контактные слои. Для source-informed
сценария `alpha_0=2.2 dB/cm`, `n=1`, bond `0.10 mm`:

| Metric | Value |
|---|---:|
| weighted one-way power transmission | 0.2105 |
| weighted reflection | 0.0931 |
| pulse broadening | 1.000 |
| pulse correlation | 0.99969 |
| postcursor/main | 0.0656 |
| center/edge group-delay difference near carrier | about 1.27 us |

Импульс сохраняется хорошо, но сохранять уже почти нечего: доминирует
поглощение в 18 mm bulk core. Разница между `T+R` и единицей является именно
поглощением, а не скрытой нормировкой.

Число `2.2 dB/cm` не названо измерением на 242 kHz. Это линейная экстраполяция
примерно `9 dB/cm @ 1 MHz` из независимого Araldite/12-um-W состава, то есть
conservative cross-recipe bracket, а не подмена отсутствующего измерения
DP190/250-um-W. При квадратичном законе она дала бы около
`0.53 dB/cm` и прошла бы дешёвый gate. Для break-even нужен показатель степени
приблизительно `n>=1.85`; именно поэтому coupon sweep обязателен.

## G2 -- aperture без дефекта

Source-informed `alpha_0=2.2 dB/cm`, bond `0.10 mm`:

| Metric | Composite rings | Straight footprint | Bare |
|---|---:|---:|---:|
| pulse peak at target / bare | 0.598 | 0.476 | 1.000 |
| candidate / straight | 1.258 | -- | -- |
| carrier transverse FWHM | 20.0 mm | 26.0 mm | 26.0 mm |
| carrier sidelobe/peak | 0.295 | 0.393 | 0.393 |
| weighted one-way power T | 0.210 | 0.229 | 1.000 |

Кольца действительно фокусируют: они сужают профиль и выигрывают у устройства
с тем же footprint/material inventory. Но абсолютные потери делают их хуже
голого датчика. В scalar axial scan нет изолированного максимума в `60 mm`;
phase coherence там максимальна, но геометрическое убывание сильнее слабого
axial build-up. Поэтому `DOF` оставлен `NaN`, а не вычислен вокруг ложного
глобального пика.

## G3 -- defect gate

Full paired `none`/`void_r2_focus` не запущен: G1/G2 закрывает его до измерения
низкочастотных свойств. Дешёвый reciprocal point-scatterer proxy показывает
границу решения (bond `0.10 mm`):

| Loss assumption | focus peak / bare | echo peak proxy / bare |
|---|---:|---:|
| lossless bound | 1.200 | 1.500 |
| `0.30 dB/cm`, `n=1` | 1.090 | 1.228 |
| `0.53 dB/cm`, `n=2` | 1.009 | 1.030 |
| `0.65 dB/cm`, `n=1` | 0.975 | 0.973 |
| `1.0 dB/cm`, `n=1` | 0.873 | 0.772 |
| `2.2 dB/cm`, `n=1` | 0.598 | 0.352 |
| `3.0 dB/cm`, `n=1` | 0.466 | 0.211 |

Численный root для point-scatterer break-even при `n=1` и `0.10 mm` bond:
`alpha_242 = 0.618 +/- 0.001 dB/cm`.

`echo_to_ringing` proxy лучше straight примерно на `13.8 dB` в source-informed
сценарии, но это не абсолютный echo-to-ringing: point-scatterer strength в
scalar Green model произвольно нормирован. Поэтому значение хранится только в
raw surrogate sweep, а common column оставлен `NaN`.

## Неопределённости и пропущенные gates

1. Главная неопределённость -- `alpha(f)` ниже `1 MHz`; ни один источник не
   даёт нужные рецептуры непосредственно на `242 kHz`.
2. Свойства промежуточных колец линейно интерполированы между двумя реально
   измеренными non-porous autoclave anchors. Их надо измерить после cure.
3. `0.10 mm` axial contact смоделирован; perfect bonding не принят молча.
   Радиальные seams в scalar model не рассеивают P/S waves и остаются риском.
4. 2-D scalar Huygens surrogate не является axisymmetric/full-3D elastic FEM.
5. G4 manufacturing perturbations не запускались, потому что G3 не прошёл.
   После material gate нужны `+/-0.10 mm`, offset `+/-1 mm`, Al speed `+/-1%`
   и случайные ring-property errors.
6. Общий `void_r2_focus` признан, но не рассчитан; его строка содержит `NaN`.

## Воспроизведение

Из корня worktree на commit `fa32a3a5...`:

```bash
julia --startup-file=no --project=. src/code/metamaterial/topology_composite_rings_tests.jl
julia --startup-file=no --project=. src/code/metamaterial/topology_composite_rings_benchmark.jl
```

Raw outputs намеренно находятся только в
`tmp/topology_benchmark/composite_rings/`:

- `ring_design.csv`;
- `local_transfer.csv`;
- `loss_bond_sweep.csv`;
- `common_metrics.csv`;
- `carrier_transverse_profile.csv`, `carrier_axial_profile.csv`;
- `loss_gate.png`, `carrier_profiles.png`.

Полный suite также пройден:

```bash
julia --startup-file=no --project=. test/runtests.jl
```

Исторические тесты `target-frequency library` и `minimum-gap pilot` требуют
уже существующие generated fixtures в `tmp/gap_pilot_220khz` и
`tmp/period_length_pilot_220khz`; они не являются выходами этой топологии и не
добавлены в Git.

## Решение для следующей итерации

1. Отлить три coupons: `9.92`, `~15`, `19.5/20.22 vol% W` тем же процессом.
2. Измерить `cP`, density, `S11/S21` и attenuation на `162--322 kHz`, минимум
   на трёх толщинах для отделения interface loss от bulk loss.
3. Продолжать full paired-defect FEM только если верхняя доверительная граница
   `alpha_242 <= 0.60 dB/cm`, velocity recipes удерживаются в `+/-1%`, а
   `0.10 mm` contact воспроизводим без voids.
4. Иначе закрыть DP190/W bulk lens как `stop` и запустить coupon research
   единственного PEEK/PEI/CFR-PEEK backup.
