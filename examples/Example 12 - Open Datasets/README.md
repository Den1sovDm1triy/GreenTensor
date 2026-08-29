<!-- SPDX-License-Identifier: MIT -->

# Example 12 — Открытые датасеты GreenTensor: что в них и как пользоваться

Библиотекой посчитаны два открытых датасета. Они устроены симметрично — одна
и та же физика, взятая в двух направлениях:

| Датасет | Направление | Что на входе → что на выходе | Репозиторий |
|---|---|---|---|
| **LL-Surrogate-6000** | прямое (суррогат) | `k0a, ε1..ε4` → `Q_sca`, `log10 Q_back` | [luneburg-lens-nn-dataset](https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset) |
| **RPM-Inverse-6000** | обратное (идентификация) | спектр `lg Q_back` в 32 точках → `t/R, ε′, ε″, μ′, μ″` | [rpm-coating-inverse-dataset](https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset) |

Первый заменяет полноволновой расчёт внутри цикла оптимизации, второй
восстанавливает параметры покрытия по спектру. Оба — MIT, оба воспроизводимы с
фиксированным зерном. Общий обзор с описанием схем — `DATASETS.md` в корне
репозитория; связанные примеры — 10 (обратная задача) и 11 (прямая).

<p align="center">
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/luneburg-lens-nn-dataset/main/figures/gt-surrogate-parity.png" width="49%" alt="LL-Surrogate-6000 parity"/>
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/rpm-coating-inverse-dataset/main/figures/inverse-parity.png" width="49%" alt="RPM-Inverse-6000 parity"/>
</p>

<p align="center"><em>Что стоит за этими датасетами: слева — суррогат линзы против расчёта
GreenTensor (R² = 0,992 по Q_sca и 0,804 по резкой цели log10 Q_back), справа — параметры покрытия,
восстановленные сетью по спектру (R² от 0,799 для ε″ до 0,985 для толщины).</em></p>

## Запуск / Run

```bash
git clone https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset
git clone https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset
python3 open_datasets_tour.py
```

Скрипт ищет датасеты в таком порядке: переменные окружения
`GT_LL_DATASET` / `GT_RPM_DATASET` → текущая папка → рядом с репозиторием →
рядом с примером. Если копии нет — печатает команду клонирования и идёт дальше;
ничего не скачивает сам.

```bash
# если датасеты лежат в другом месте
export GT_LL_DATASET=/path/to/luneburg-lens-nn-dataset
export GT_RPM_DATASET=/path/to/rpm-coating-inverse-dataset
python3 open_datasets_tour.py
```

## Что печатает

1. **Состав** каждого датасета: формы массивов, диапазоны входов и целей.
2. **Проверку меток — главное в этом примере.** Скрипт берёт несколько строк и
   пересчитывает их заново: для линзы Люнеберга — прямо этой библиотекой
   (`green_tensor.sphere_core.MieSphere`), для покрытия — решателем, вложенным
   в сам датасет. Совпадение до машинной точности (`~1e-14`) означает, что
   данные и код связаны, а не просто лежат рядом.
3. **Что делать дальше**: какие скрипты обучения запускать и с какими
   опорными метриками сравниваться.

Ожидаемый вывод по проверке:

```
  largest deviation on these rows: 1.15e-14  (labels reproduce exactly)
  largest deviation: 2.55e-15  (spectra reproduce to machine precision)
```

## Как этим пользоваться

1. **Взять как бенчмарк.** Опорная сеть по линзе даёт `R² ≈ 0.99` по `Q_sca` и
   `R² ≈ 0.80` по резкой цели `log10 Q_back`; вторая цель намеренно оставлена
   трудной — это ориентир для сравнения архитектур, а не решённая задача.
   Обратная сеть по покрытию проверяется дополнительно при шуме спектра 0,5 дБ —
   именно эта метрика честна для диагностики.
2. **Пересобрать данные под свою задачу.** В каждом репозитории есть
   `generate_dataset.py`: у линзы он зовёт GreenTensor, у покрытия — вложенный
   самодостаточный решатель. Меняются диапазоны параметров, число слоёв, сетка
   частот — и получается ваш датасет с той же физикой.
3. **Считать свои конфигурации.** Функция проверки из этого скрипта — готовый
   образец вызова: `MieSphere(k0=..., a=[0.25, 0.5, 0.75, 1.0], eps=[...],
   problem="diffraction").cross_sections()`.

## Требования

NumPy — обязательно; SciPy и Matplotlib нужны для расчётной части библиотеки
(проверка меток линзы). Без них скрипт всё равно отработает: он честно скажет,
что проверку пропускает, и покажет состав датасетов.

---

## English (short)

A guided tour of the two open datasets generated with GreenTensor: the forward
Luneburg-lens surrogate set (`LL-Surrogate-6000`) and the inverse coating
identification set (`RPM-Inverse-6000`). The script locates the datasets
(env vars `GT_LL_DATASET` / `GT_RPM_DATASET`, cwd, or next to the repository),
prints their schema and ranges, and — the point of the example — **verifies the
stored labels by recomputing them**: the lens set against this library, the
coating set against the solver bundled with it. Agreement is to machine
precision. Nothing is downloaded; if a dataset is absent the script prints the
`git clone` command. Overview: `DATASETS.md`; related examples: 10 and 11.
