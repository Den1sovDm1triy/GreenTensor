# SPDX-License-Identifier: MIT
"""library_bridge — упрощённый мост к пакету green_tensor.

В версии Studio для сервера этот модуль искал пакет по GREENTENSOR_PATH / соседним
каталогам. В браузерной сборке (Pyodide) пакет green_tensor записан прямо в корень
виртуальной ФС и импортируется напрямую — поиск не нужен.
"""
from __future__ import annotations

import functools


@functools.lru_cache(maxsize=1)
def load():
    import green_tensor as gt  # noqa: WPS433 — импорт в браузерной ФС
    return gt
