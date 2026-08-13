<!-- SPDX-License-Identifier: MIT -->

# green_tensor/legacy — archival research scripts / архивные скрипты

**EN.** Original standalone research scripts from the early GreenTensor work. They are
kept verbatim for reproducibility and provenance. They are **not** part of the public
API, are **not** imported anywhere, and are **not** maintained. For new work use the
unified solvers in [`green_tensor.solvers`](../solvers.py).

**RU.** Оригинальные самостоятельные исследовательские скрипты ранней версии GreenTensor.
Сохранены как есть — для воспроизводимости и истории. Они **не** входят в публичный API,
**нигде** не импортируются и **не** поддерживаются. Для новой работы используйте единые
решатели в [`green_tensor.solvers`](../solvers.py).

| File | Назначение / Purpose |
|------|----------------------|
| `02_cilindre.py` | Цилиндр (численный VIE/GMRES) / cylinder (numerical VIE/GMRES) |
| `03_decart_dipol.py` | Декартов диполь / Cartesian dipole |
| `03_decart_plane.py` | Декартова плоская волна / Cartesian plane wave |
| `Bistatic_RCS_lin_polar.py` | Бистатическая ЭПР, линейная поляризация / bistatic RCS, linear pol. |
| `Bistatic_RCS_lin+circle_polar.py` | Бистатическая ЭПР, линейная + круговая / bistatic RCS, linear + circular |
| `lin_polar.py` | Линейная поляризация / linear polarization |
| `examples.json` | Входные данные для скриптов выше / input data for the scripts above |

> The canonical, maintained sphere solver lives in `green_tensor/01_sphere.py`
> (imported through the `green_tensor.sphere_core` facade) — **not** here.
>
> Каноническое поддерживаемое решение для сферы — в `green_tensor/01_sphere.py`
> (импорт — через фасад `green_tensor.sphere_core`), **не** здесь.
