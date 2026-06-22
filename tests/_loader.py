"""Безопасная загрузка расчётных модулей GreenTensor для тестов.

Имена файлов с цифровым префиксом (``01_sphere.py``) нельзя импортировать
обычным ``import``, поэтому грузим их по пути через importlib. matplotlib
переводится в неинтерактивный backend ``Agg`` ДО импорта, чтобы расчётные
модули не открывали окна.
"""
from __future__ import annotations

import importlib.util
import os

import matplotlib

matplotlib.use("Agg")

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_module(rel_path: str, name: str):
    """Загрузить модуль из green_tensor/ по относительному пути."""
    path = os.path.join(_ROOT, rel_path)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # __main__-блоки защищены guard'ом — побочных эффектов нет
    return mod


def load_sphere():
    return load_module(os.path.join("green_tensor", "01_sphere.py"), "gt_sphere")
