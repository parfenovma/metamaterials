# Контракт O3 F1: broadband transient transfer и paired ensemble

> Зафиксировано 15 августа 2026 года до первого transient solve. Этот документ
> продолжает прошедший carrier gate F0. Он не открывает O4 после одного
> deterministic pulse и не разрешает менять frozen operators.

## 1. Решение F1

F1 отвечает на два вложенных вопроса:

```text
F1a  Переносится ли общий 17 x 17 физический оператор в прямой pulse-echo
     transient с конечной полостью, backwall и P/S mode conversion?

F1b  Сохраняют ли frozen static/motion классы detection/localization advantage
     на raw no-defect/defect trials с contact, material, ringing и noise?
```

F1 считается пройденным только после обоих уровней. Положительный F1a
разрешает stochastic F1b, но сам не является POD или выбором архитектуры.

## 2. Frozen physics F1a

```text
dimension                    2D plane strain, full isotropic elasticity
mesh pair                    common conforming free/void pair from F0
h_bulk / h_defect            1.0 / 0.25 mm
specimen                     125.0 x 71.6 mm, nominal aluminium 6061
array                        17 equal normal-traction/normal-receive elements
source                       five-cycle Hann-windowed 242 kHz traction burst
time integrator              average-acceleration Newmark, beta=1/4, gamma=1/2
samples per carrier cycle    30
record length                70 us
input and lateral faces      first-order P/S impedance boundary
backwall                     free reflecting boundary
defect                       free circular void r=2 mm at x=60 mm, y=0
```

Semi-discrete matrices `M`, `C`, `K` собираются один раз на mesh. Эффективная
Newmark-матрица факторизуется один раз; все 17 elemental transmit experiments
решаются одновременно как multi-RHS. Receive matrix определяется
work-conjugate boundary integral, поэтому reciprocity проверяется на полном
`R(t)[receive, transmit]`, а не на выбранной архитектуре.

## 3. Frozen operators и энергия

F1a использует те же шесть проекций, что F0:

```text
uniform A-scan
static K4,Q2
static K3,Q4
motion P12 compact
motion P17 full
FMC oracle.
```

Static projection равна `W_q^T R(t) W_q / sqrt(KQ)`, motion — выбранной
диагонали `R(t)/sqrt(P)`, FMC — `vec(R(t))/sqrt(17)`. Codebooks, motion
indices и incident-energy convention не переоптимизируются по transient.

## 4. Временные окна

Окна выводятся из геометрии, а не подбираются по сигналу:

```text
earliest void arrival     2 (x_void - r_void) / cP
void window end           2 (x_void + r_void) / cP + burst_duration
backwall arrival          2 L / cP
record end                70 us.
```

По два time steps добавляются к границам окон только как численный guard.
Paired field `Delta R(t)=R_void(t)-R_free(t)` используется в F1a для проверки
переноса. Detector в F1b получает raw trials; бесплатный defect-free reference
ему не передаётся.

## 5. Предварительно зафиксированные F1a gates

```text
all time samples finite                         yes
relative reciprocity residual over R(t)         <= 1e-8
centred paired mirror residual                  <= 3%
paired energy before earliest void arrival      <= 1%
paired energy in geometric void window          > 0
free response energy in backwall window         > 0
paired tangential response                      finite and reported
same pulse, mesh exterior and receive chain     yes
```

Architecture observables F1a: raw free energy, paired echo energy in the void
window, paired/free contrast, peak time, backwall-window energy and norm
relative to FMC. Они не являются `Pd`, `Pfa` или `alpha90`.

## 6. Минимальный F1b после F1a

F1b должен использовать независимые train/validation/test splits и добавить:

- конечный defect ensemble по depth/lateral/size/orientation;
- contact amplitude/delay variation без full recalibration на каждый контакт;
- material speed/attenuation variation;
- source-face ringing и фактический backwall response;
- correlated electronics noise и независимые raw H0 trials;
- те же `Pfa`, `Pd`, `Pdl`, `alpha90`, localization и resource-vector gates,
  которые были заморожены для O3.

Ни carrier F0, ни deterministic F1a не могут открыть O4. Если F1a не проходит
physics gate, сначала исправляется transient discretization/boundary model; ни
geometry synthesis, ни повторная оптимизация codebook не разрешены.

## 7. Frozen F1b pilot protocol

До анализа finite-defect словаря фиксируется минимальный F1b pilot:

```text
defect classes
    void r=1 mm at (60, 0) mm
    void r=2 mm at (35, 0), (60, 0), (85, 0), (60, 10) mm
    elliptical slit a=4 mm, b=0.2 mm, 30 deg at (60, 0) mm

feature record              8--70 us, every third F1a sample
background calibration      train-only mean + at most 24 nuisance PCs
PC energy target            99%
H0 validation / test        10 000 / 10 000 trials
signal validation / test    30 / 60 contexts per class, 4 noise draws each
Pfa target                  1%
amplitude grid              0.25, 0.4, 0.6, 0.8, 1.0, 1.3,
                            1.7, 2.2, 3.0, 4.0, 6.0, 8.0
base seed                   0x051b2026
```

На context применяются общие для raw free/defect signal вариации:

```text
global contact amplitude sigma       5%
per-observable gain sigma             2%
global contact delay sigma            0.20 us
material time-stretch sigma           0.5%
late-time attenuation-slope sigma     10%
additional ringing scale sigma        2% of the warped free response
```

Electronics noise имеет одинаковую variance на observation sample для всех
архитектур. Его абсолютный `sigma_e` определяется только на training split как
`median ||P_nuisance Delta_FMC|| / 5`; отдельно под кандидата noise не
масштабируется. Matched templates и nuisance projector frozen после training.

Reference amplitude выбирается как первая точка сетки, где FMC достигает
validation `Pd >= 0.95`. Threshold каждой архитектуры выбирается только по её
validation H0 при `Pfa=1%`, затем один раз применяется к test.

Pilot gates сохраняют исходную направленность O3:

```text
test Pfa upper 95% bound          <= 1.25%
aggregate Pd(reference alpha)     >= 0.90
Pdl within 0.5 lambda_P0          >= 0.80
minimum per-class Pd              >= 0.75
aggregate A-scan -> FMC closure   >= 0.50
```

Шесть классов недостаточны для промышленного POD. Даже положительный pilot
требует расширения continuous ROI перед O4; отрицательный pilot уже достаточен,
чтобы не продвигать соответствующий operator class.

## 8. Результат F1

Deterministic `F1a` прошёл после явной нулевой инициализации всех Newmark
state arrays:

```text
reciprocity residual              3.73e-14
paired mirror residual            0.549%
precursor energy fraction         2.63e-11
paired tangential/normal norm     0.662
```

Frozen `F1b` pilot затем дал отрицательный архитектурный результат. Ни один
static/motion кандидат не прошёл все gates. На FMC validation reference
`alpha=2.2` static `K4,Q2` и `K3,Q4` имеют aggregate `Pd=0.858/0.849`, а
motion `P12/P17` — `0.581/0.603`. Худший класс для всех — `void_r1_focus`;
его `Pd` равен соответственно `0.183`, `0.238`, `0.058`, `0.071`.

Одновременно не прошёл positive-control task condition: uniform A-scan имеет
`alpha90=0.635`, FMC — `1.523`. Поэтому `A-scan -> FMC` gap отрицателен и
closure не определяется. Причина — шесть изолированных дефектов хорошо
разделяются по времени, а когерентный uniform transmit сохраняет высокий SNR;
в pilot отсутствуют full-wave structural clutter и несколько одновременных
reflectors, ради которых FMC должен давать information advantage.

Следовательно, F1 не открывает O4. Текущие frozen кандидаты получают no-go на
этом pilot, а сам task prior должен быть расширен плотной lateral ROI,
multi-scatterer/background ensemble и PAUT positive-control gate до нового
архитектурного сравнения. Этот вывод не является утверждением, что A-scan
доминирует PAUT в общей дефектоскопии.
