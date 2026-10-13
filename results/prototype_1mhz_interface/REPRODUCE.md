# Воспроизведение проверки интерфейса

Команды выполняются из корня репозитория с окружением из `Project.toml` и
`Manifest.toml`.

```sh
julia --startup-file=no --project=. \
  src/code/metamaterial/prototype_1mhz_interface_mesher.jl

julia --startup-file=no --project=. \
  src/code/metamaterial/step1b_convert_models.jl \
  tmp/prototype_1mhz_interface/meshes \
  tmp/prototype_1mhz_interface/models

OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=. \
  src/code/metamaterial/prototype_1mhz_interface_solver.jl solve
```

Результат проверки записывается в
`tmp/prototype_1mhz_interface/VERDICT.md`, а двадцать исходных строк
частотного прогона — в `tmp/prototype_1mhz_interface/sweep.csv` и
`tmp/prototype_1mhz_interface/sweep.jld2`.

Повторный анализ уже рассчитанной таблицы:

```sh
julia --startup-file=no --project=. \
  src/code/metamaterial/prototype_1mhz_interface_solver.jl analyze
```
