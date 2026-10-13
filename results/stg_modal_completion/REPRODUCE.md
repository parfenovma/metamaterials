# Воспроизведение и соглашения

Запускать из корня проекта в существующей среде Julia 1.12.3. Ниже `julia`
означает исполняемый файл этой среды. В выполненном расчёте использован
`/Users/mixx3/.julia/juliaup/julia-1.12.3+0.aarch64.apple.darwin14/bin/julia`.
Зависимости берутся из существующих Project.toml и Manifest.toml.

```sh
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=. src/code/metamaterial/run_stg_completion.jl --mesh
julia --startup-file=no --project=. src/code/metamaterial/run_stg_completion.jl --convert
julia --startup-file=no --project=. src/code/metamaterial/run_stg_completion.jl --modes
julia --startup-file=no --project=. src/code/metamaterial/run_stg_completion.jl --pilot
julia --startup-file=no --project=. src/code/metamaterial/run_stg_completion.jl --sweep 0
julia --startup-file=no --project=. src/code/metamaterial/run_stg_completion.jl --sweep 1
julia --startup-file=no --project=. src/code/metamaterial/refine_stg_completion.jl --modes
julia --startup-file=no --project=. src/code/metamaterial/refine_stg_completion.jl --solve 0
julia --startup-file=no --project=. src/code/metamaterial/refine_stg_completion.jl --solve 1
julia --startup-file=no --project=. src/code/metamaterial/refine_stg_completion.jl --level2 0
julia --startup-file=no --project=. src/code/metamaterial/refine_stg_completion.jl --level2 1
julia --startup-file=no --project=. src/code/metamaterial/stg_bloch.jl
julia --startup-file=no --project=. src/code/metamaterial/stg_bloch.jl --complex
julia --startup-file=no --project=. src/code/metamaterial/validate_stg_completion.jl
julia --startup-file=no --project=. src/code/metamaterial/finish_stg_controls.jl
julia --startup-file=no --project=. src/code/metamaterial/probe_stg_modal_inputs.jl
julia --startup-file=no --project=. src/code/metamaterial/analyze_stg_completion.jl
julia --startup-file=no --project=. src/code/metamaterial/plot_stg_completion_fields.jl
julia --startup-file=no --project=. test/staggered_waveguide_geometry_tests.jl
julia --startup-file=no --project=. test/stg_packet_tests.jl
julia --startup-file=no --project=. test/multimode_prepared_tests.jl
```

Группы 0/1 одной стадии независимы. Расчёт геометрии Gmsh и конвертация
в Gridap JSON выполняются отдельными процессами; FEM читает JSON.
Предварительно собранные объёмные матрицы переиспользуются только для той же
геометрии, материала по ρ/cP/cS, порядка элемента, квадратуры и области потерь.

## Частота, фаза, время

Начальная сетка 20–360 кГц с шагом 2.5 кГц. После проверки IFFT она уточнена
до 1.25 кГц для номинальных потерь и до 0.625 кГц для STG без потерь.
Второе уточнение потребовалось из-за рождения дополнительных портовых ветвей
около 151 кГц. N=1/4 без потерь рассчитаны с шагом 5 кГц в диапазоне
80–280 кГц и с дополнительными контрольными возбуждениями.

Временное соглашение: `exp(+i omega t)`, бегущая вправо волна:
`exp(-i k x)`. Направление определяется знаком потока энергии, а не знаком k:
одна из отслеживаемых ветвей имеет участок обратного фазового движения.

Портовые формы нормированы на 1 Вт/м. Их произвольная комплексная фаза
устранена последовательным согласованием перекрытий; противоположная мода
имеет сопряжённую форму перемещений. Новые частоты при уточнении согласуются
с уже сохранёнными формами, не получают независимые случайные фазы.

Коэффициенты приведены к плоскостям концов центральной части:

```math
t_{m0}^{core}=t_{m0}^{ports}\exp[i(k_0+k_m)L_{lead}],\qquad
r_{m0}^{core}=r_{m0}^{ports}\exp[i(k_0+k_m)L_{lead}].
```

Для вещественного четырёхпериодного Hann-импульса вычисляется FFT. На
положительных частотах его спектр умножается на комплексные коэффициенты
каждой ветви; отрицательная половина задаётся сопряжением. Отсутствующие
распространяющиеся ветви имеют нулевой распространяющийся вклад. Плотность
точек FFT и длительность окна согласованы с частотным шагом; `dt=0.1 мкс`.
Импульс начинается через 50 мкс после начала окна, чтобы проверить края.
При сравнении частотных сеток обе реконструкции вычисляются на одном FFT-окне.

Модальные временные сигналы имеют мощностную спектральную нормировку.
Интегралы их квадратов дают относительные энергии через Parseval. Их квадрат
в отдельный момент времени не объявляется точным локальным физическим
потоком многомодального поля. CWT — Морле, общий масштаб относительно
падающего импульса, с исключением cone of influence.

## Файлы

- `run_manifest.json` — первоначальные параметры и критерии, хеши источников
  перед основным расчётом; `final_manifest.json` — итоговая фиксация.
- `models/`, `meshes/` — геометрии и физические метки Domain/Core/Source/Microphone.
- `port_modes_120.jld2`, `midpoint_modes.jld2`, `second_modes_group*.jld2` —
  собственные формы с согласованной фазой и идентификаторами ветвей.
- `cases/` — гармонические коэффициенты, поля в отдельных контрольных точках,
  параметры случая и SHA256 JSON-модели. Имена содержат N, длину продолжения,
  частоту, масштаб потерь, число затухающих мод и разрешение портов.
- `sweep_group*.csv`, `refined_group*.csv`, `second_refined_group*.csv` —
  исходная и уточнённые частотные серии.
- `packet_summary.csv`, `packet_modes.csv`, `packet_checks.csv` — энергии
  импульсов, состав мод и проверка сходимости. `packet_*.jld2` содержит сигналы.
- `bloch.jld2`, `bloch_checks.csv`, `complex_bloch.*` — вещественная дисперсия,
  сеточный контроль и комплексные множители. Ячейка не содержит концевые переходы.
- `validation_cases.csv`, `length_controls.csv` — уточнение FEM/портов,
  продолжения 28 мм, малое прохождение N=4 и прямые референсы равной длины.
- `modal_input_scan.csv`, `modal_input_matrix_*.jld2`, `modal_matrix_checks.csv` —
  весь блок передачи и отражения для возбуждения слева, нормированный через
  измеренную матрицу входящих амплитуд. Проверены `T' T + R' R ≈ I` и взаимность
  отражения. Это не измеренная матрица реального преобразователя.

Гармонические случаи кешируются. При изменении физики, фазового соглашения,
решателя или геометрии нужен новый каталог либо явная очистка затронутых
случаев. Не следует смешивать кэши разных постановок. Бинарные результаты,
CSV и изображения игнорируются Git по существующим правилам проекта;
Markdown-резюме и исходные программы остаются доступными для версионирования.
