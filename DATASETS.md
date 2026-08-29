<!-- SPDX-License-Identifier: MIT -->

# Open datasets built with GreenTensor

Two open benchmark datasets are generated with this library. Both ship data,
reference models and reproduction scripts; neither reproduces the text of the
accompanying papers. Both are MIT-licensed, both are fixed-seed reproducible,
and both can be regenerated from scratch with the code in this repository.

<p align="center">
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/luneburg-lens-nn-dataset/main/figures/gt-surrogate-parity.png" width="49%" alt="LL-Surrogate-6000: surrogate vs GreenTensor parity"/>
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/rpm-coating-inverse-dataset/main/figures/inverse-parity.png" width="49%" alt="RPM-Inverse-6000: recovered vs true coating parameters"/>
</p>

<p align="center"><em>Left — the forward set: a plain MLP against the GreenTensor labels
(R² = 0.992 for Q_sca, 0.804 for the sharp log10 Q_back). Right — the inverse set: coating
parameters recovered from a spectrum (R² from 0.799 for ε″ to 0.985 for thickness).</em></p>

| Dataset | Direction | Rows | Input → output | Repository |
|---|---|---|---|---|
| **LL-Surrogate-6000** | forward (surrogate) | 6000 | `k0a, ε1..ε4` → `Q_sca`, `log10 Q_back` | [luneburg-lens-nn-dataset](https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset) |
| **RPM-Inverse-6000** | inverse (identification) | 6000 | 32-point spectrum `lg Q_back(k0R)` → `t/R, ε′, ε″, μ′, μ″` | [rpm-coating-inverse-dataset](https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset) |

The pair is deliberately symmetric: the same physics, taken in both directions.
The first replaces a full-wave solve inside an optimization loop; the second
recovers material parameters from a measured-like spectrum.

## LL-Surrogate-6000 — forward map of a four-layer Luneburg lens

A four-layer spherical Luneburg-type lens with fixed normalized layer radii
`r_l / a = 0.25, 0.50, 0.75, 1.00`. Inputs are the electrical size
`k0a ∈ [π, 8π]` and four permittivities `ε1 ≥ ε2 ≥ ε3 ≥ ε4`, each in `[1, 2]`.
Targets are the total scattering cross-section `Q_sca` and `log10 Q_back`.

Labels are computed by the multilayer Mie / tensor Green's-function solver of
this library (`green_tensor.sphere_core.MieSphere`), so the whole set is built
in minutes on a laptop instead of hours in a mesh package.

Reference model (NumPy-only MLP, shipped in the dataset repo): held-out
`R² ≈ 0.99` for `Q_sca`, `R² ≈ 0.80` for `log10 Q_back`; 5-fold CV
`0.989 ± 0.002` and `0.798 ± 0.016`; the learning curve saturates near
`N ≈ 3600`. The sharp target is deliberately left hard — it is a benchmark, not
a solved problem.

<p align="center">
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/luneburg-lens-nn-dataset/main/figures/ll-profile-approx.png" width="49%" alt="Continuous Luneburg law vs equal-step layering"/>
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/luneburg-lens-nn-dataset/main/figures/dataset-coverage.png" width="49%" alt="Feature-space coverage of the 6000 samples"/>
</p>

<p align="center"><em>Left — what the four layers approximate: the continuous law ε(r) = 2 − (r/a)²
against equal-step stratification (L = 4 used by the dataset, L = 8 for comparison).
Right — how the 6000 samples cover the input space.</em></p>

<p align="center">
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/luneburg-lens-nn-dataset/main/figures/learning-curve.png" width="60%" alt="R2 vs training-set size"/>
</p>

<p align="center"><em>Learning curve: R² against training-set size — the smooth target saturates
early, the sharp one keeps gaining, which is why 6000 rows is the useful size.</em></p>

Related example in this repository: `examples/Example 11 - NN Luneburg Lens Surrogate`.

## RPM-Inverse-6000 — coating identification from a back-scattering spectrum

A PEC sphere of radius `R` carrying one homogeneous magnetodielectric
(radar-absorbing) coating of thickness `t`. The observable is back-scattering
efficiency on a 32-point grid in `k0R ∈ [1.5, 8.5]`; the array `S` stores
`lg Q_back` at those points. The network learns the inverse map: spectrum →
`t/R`, `ε′`, `ε″`, `μ′`, `μ″`.

The dataset repository bundles its own self-contained NumPy/SciPy solver
(`forward_solver.py`, travelling-wave transfer through the layer), so the data
can be regenerated without installing this library; the physics is the same
tensor Green's-function formulation. Each sample passes an energy-balance check
before it enters the set.

<p align="center">
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/rpm-coating-inverse-dataset/main/figures/spectra-examples.png" width="70%" alt="Back-scattering spectra of a bare PEC sphere and three coatings"/>
</p>

<p align="center"><em>The input signal of the benchmark: back-scattering of a bare PEC sphere
against three coating types. A thick lossy coating suppresses the return by tens of decibels and
reshapes the interference pattern — that reshaping is what the network reads.</em></p>

<p align="center">
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/rpm-coating-inverse-dataset/main/figures/dataset-coverage.png" width="49%" alt="Parameter-space coverage"/>
  <img src="https://raw.githubusercontent.com/Den1sovDm1triy/rpm-coating-inverse-dataset/main/figures/respectral-check.png" width="49%" alt="Respectral check of reconstructed configurations"/>
</p>

<p align="center"><em>Left — parameter-space coverage of the 6000 samples. Right — the respectral
check: spectra recomputed from the recovered parameters against the originals, the honest test for
an inverse problem where different parameter sets can produce similar spectra.</em></p>

Related example in this repository: `examples/Example 10 - NN Coating Parameter Identification`.

## Reproducing and verifying

Every label is reproducible from the parameters that sit next to it. For the
Luneburg set, the check runs against this library directly:

```python
import numpy as np
from green_tensor.sphere_core import MieSphere

d = np.load("ll_surrogate_6000.npz")
X, Y = d["X"], d["Y"]
cs = MieSphere(k0=X[0, 0], a=[0.25, 0.5, 0.75, 1.0],
               eps=list(X[0, 1:]), problem="diffraction").cross_sections()
print(Y[0], cs["q_sca"], np.log10(max(cs["q_back"], 1e-8)))   # identical
```

A ready-made tour of both datasets — load, inspect, verify against the solver —
is `examples/Example 12 - Open Datasets`.

## Getting the data

```bash
git clone https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset
git clone https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset
```

Each repository contains `.npz` and `.csv` versions of the table, the reference
model, robustness scripts (5-fold CV, seed spread, learning curve), figures and
a `CITATION.cff`.

## Licence and citation

Both datasets are MIT, data and code alike. If you use them, cite the dataset
repository (see its `CITATION.cff`) together with the accompanying paper, and —
when the labels come from this solver — GreenTensor itself (see `CITATION.cff`
here).

---

## RU — кратко

Библиотекой посчитаны два открытых датасета; они устроены симметрично, как
прямая и обратная задача одной физики.

- **LL-Surrogate-6000** — прямая (суррогатная) задача: по электрическому размеру
  `k0a` и четырём проницаемостям четырёхслойной линзы Люнеберга предсказываются
  сечения `Q_sca` и `log10 Q_back`. 6000 строк, метки посчитаны
  `green_tensor.sphere_core.MieSphere`. Опорная сеть даёт `R² ≈ 0.99` по
  `Q_sca` и `R² ≈ 0.80` по резкой цели `log10 Q_back` — вторая намеренно
  оставлена открытым ориентиром для сравнения архитектур.
  Репозиторий: <https://github.com/Den1sovDm1triy/luneburg-lens-nn-dataset>
- **RPM-Inverse-6000** — обратная задача: по спектру обратного рассеяния из 32
  точек (`k0R` от 1,5 до 8,5) восстанавливаются толщина покрытия и комплексные
  ε, μ. 6000 пар, решатель вложен в сам датасет, у каждой метки проверен
  энергобаланс.
  Репозиторий: <https://github.com/Den1sovDm1triy/rpm-coating-inverse-dataset>

Оба датасета — MIT, воспроизводимы с фиксированным зерном, пересчитываются
скриптами из своих репозиториев. Экскурсия по обоим с проверкой меток —
`examples/Example 12 - Open Datasets`; связанные примеры — 10 и 11.
