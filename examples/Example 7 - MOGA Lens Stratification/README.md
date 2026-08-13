<!-- SPDX-License-Identifier: MIT -->

# Example 7 — MOGA Lens Stratification

Двухкритериальная стратификация линзы Люнеберга: минимизация RMSE(ε)
аппроксимации непрерывного закона ε(r) = 2 − (r/a)² кусочно-постоянным
профилем из L слоёв при одновременной минимизации максимального скачка
Δε_max между слоями.

Two-objective stratification of the Luneburg lens: NSGA-II-family genetic
optimization of shell radii and permittivities, with the RMSE endpoint of
the Pareto front verified against the independent Lloyd optimum of
piecewise-constant L2 approximation.

Результаты для L = 8 (канон): RMSE(ε) = 0.0320 — на 13 % лучше схемы
Fuchs et al. (IEEE TAP 2007, 0.0368) и на 23 % лучше равношаговой (0.0416).

## Запуск / Run

```bash
python3 moga_stratification.py
```

Выход: `moga_stratification.csv` (RMSE трёх схем при L = 1…16),
`moga_layers.csv` (радиусы/проницаемости слоёв), рисунки
`moga-rmse-comparison.*`, `moga-pareto-l8.*`.

Зависимости: numpy, matplotlib.
