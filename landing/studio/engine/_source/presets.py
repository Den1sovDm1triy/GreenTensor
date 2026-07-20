# SPDX-License-Identifier: MIT

"""presets — задачи GreenTensor Studio.

Две группы:
  * «Сравнение с Ansys HFSS (МКЭ)» — кейсы слоистых сфер с угловой ДН, рассчитанной в
    Ansys HFSS (метод конечных элементов). Источники данных:
      - ``studio/data/examples.json`` (массивы DN_NORM*) — срез E-плоскости (φ=0);
      - ``studio/data/Ephi=0/90_from_HFSS.csv`` (Example 2) — ОБЕ плоскости E (φ=0) и H (φ=90).
    При выборе кейса сцена считается, а кривые HFSS накладываются на соответствующие
    плоскости диаграммы — независимая проверка аналитического ядра ``01_sphere.py`` против
    конечно-элементного решения.
  * «Демонстрации» — возможности библиотеки без данных HFSS (линза Максвелла, бесконечный
    слоистый цилиндр, кластер сфер GMM).

Слои сферы: ``a_rel`` — нормированная ВНЕШНЯЯ граница слоя (внешняя = 1.0), изнутри
наружу; внешняя среда (вакуум) учитывается ядром. Металл — большим ``eps`` (PEC-предел).
"""
from __future__ import annotations

import csv
import json
import re
import math
import os

_PI = math.pi
_DATA = os.path.join(os.path.dirname(__file__), "data")
_EXAMPLES = os.path.join(_DATA, "examples.json")
# PEC-прокси: eps=1e40 из examples.json переполняет ряд ВСВФ; eps=1e7 даёт тот же
# PEC-предел и численно устойчив (как в рабочих расчётах сферы).
_PEC_CAP = 1.0e7


def _mk_layers(layers):
    """layers: список (a_rel, eps_re, eps_im) или (..., mu_re, mu_im) -> список dict-слоёв."""
    out = []
    for L in layers:
        a_rel, eps_re, eps_im = L[0], L[1], L[2]
        mu_re = L[3] if len(L) > 3 else 1.0
        mu_im = L[4] if len(L) > 4 else 0.0
        out.append({"a_rel": a_rel, "eps_re": eps_re, "eps_im": eps_im,
                    "mu_re": mu_re, "mu_im": mu_im})
    return out


def _sphere(layers, *, position=(0.0, 0.0, 0.0), radius=1.0):
    """Радиально-слоистая сфера (внешний радиус ``radius``, слои по a_rel изнутри наружу)."""
    return {"type": "sphere", "position": list(position), "radius": radius,
            "layers": _mk_layers(layers)}


def _cyl_inf(radius, layers, theta_deg=90.0, mode="TM"):
    """Бесконечный радиально-слоистый цилиндр (ось ∥ z). theta_deg — угол падения от оси."""
    return {"type": "cylinder_inf", "position": [0.0, 0.0, 0.0],
            "radius": radius, "theta_deg": theta_deg, "mode": mode,
            "layers": _mk_layers(layers)}


def _scene(body, k, *, polarization="linear", problem="diffraction", phi_deg=0.0, toch=None,
           khat=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0)):
    bodies = body if isinstance(body, list) else [body]
    return {
        "bodies": bodies,
        "radiation": {
            "k": k, "polarization": polarization, "problem": problem,
            "phi_deg": phi_deg, "toch": toch,
            "khat": list(khat), "pol": list(pol),
        },
        "sweep": {"enable": True, "k0_min": 0.25, "k0_max": max(2.0 * k, 6.0), "points": 160},
    }


# --------------------------------------------------------------------------- #
# Кривые Ansys HFSS (МКЭ) — наложение на диаграмму по плоскостям.
# Каждый кейс несёт reference = список {plane, name, theta_deg, dB}; кривая каждой
# плоскости (E/H) ложится на соответствующую расчётную плоскость.
# --------------------------------------------------------------------------- #
def _sorted_curve(theta, db, plane, name):
    order = sorted(range(len(theta)), key=lambda j: theta[j])
    return {"plane": plane, "name": name,
            "theta_deg": [float(theta[j]) for j in order],
            "dB": [float(db[j]) for j in order]}


def _ref_dn(dn):
    """DN_NORM (examples.json): нормированная ДН в дБ на нативной сетке 01_sphere
    (alpha = 0.01 + k°) -> кривая E-плоскости на стандартной θ Studio (вперёд = 0):
    θ = (180 − alpha) mod 360."""
    dn = [float(x) for x in dn]
    alpha = [0.01 + k for k in range(len(dn))]
    theta = [(180.0 - a) % 360.0 for a in alpha]
    return [_sorted_curve(theta, dn, "E", "Ansys HFSS, МКЭ (E-плоскость)")]


def _ref_hfss_csv(fname, plane, name):
    """HFSS-CSV Example 2 (колонки Theta[deg], dB20normalize(rETotal)). Угол HFSS уже в
    конвенции «вперёд = 0» -> стандартная θ Studio напрямую: θ = Theta mod 360."""
    path = os.path.join(_DATA, fname)
    theta, db = [], []
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            theta.append(float(row["Theta [deg]"]) % 360.0)
            db.append(float(row["dB20normalize(rETotal) []"]))
    return _sorted_curve(theta, db, plane, name)


# --------------------------------------------------------------------------- #
# Кейсы examples.json (срез E-плоскости): 4 слоистых сферы
# --------------------------------------------------------------------------- #
def _eval_k0(expr):
    """k0 в examples.json задан строкой: число или произведение вида "10 * (math.pi)"."""
    text = str(expr).strip()
    try:
        return float(text)
    except ValueError:
        pass
    m = re.fullmatch(r"([0-9.]+)\s*\*\s*\(?\s*math\.pi\s*\)?", text)
    if m:
        return float(m.group(1)) * math.pi
    raise ValueError(f"неразбираемое выражение k0: {expr!r}")


def _case_layers(case):
    a, eps = case["a"], case["eps"]
    miy = case.get("miy", [1.0] * len(a))
    out = []
    for i in range(len(a)):
        e = float(eps[i])
        if e > _PEC_CAP:
            e = _PEC_CAP
        out.append((float(a[i]), e, 0.0, float(miy[i]), 0.0))
    return out


def _parse_examples(path=_EXAMPLES):
    raw = open(path, encoding="utf-8").read()
    dec = json.JSONDecoder()
    objs, i = [], 0
    while i < len(raw):
        while i < len(raw) and raw[i] in " \t\r\n":
            i += 1
        if i >= len(raw):
            break
        obj, i = dec.raw_decode(raw, i)
        objs.append(obj)
    cases = [o for o in objs if "n" in o and "DN_NORM_Ansys" in o]
    arrays = {k: v for o in objs for k, v in o.items()
              if k.startswith("DN_NORM") and isinstance(v, list)}
    return cases, arrays


_BENCH_META = [
    ("luneburg5_hfss",      "Линза Люнеберга, 5 слоёв"),
    ("metal_sphere5_hfss",  "Металлическая сфера, 5 слоёв (PEC)"),
    ("metal_sphere2_hfss",  "Металлическая сфера, 2 слоя (PEC)"),
    ("metal_coated_hfss",   "Металл + диэлектрическое покрытие"),
]


def _examples_presets():
    """4 кейса из examples.json: сцена + кривая Ansys HFSS только E-плоскости."""
    try:
        cases, arrays = _parse_examples()
    except Exception:
        return []
    out = []
    for idx, case in enumerate(cases[:len(_BENCH_META)]):
        pid, base_title = _BENCH_META[idx]
        dn_key = case["DN_NORM_Ansys"]
        if dn_key not in arrays:
            continue
        is_metal = any(float(e) > _PEC_CAP for e in case["eps"])
        desc = (f"Кейс из examples.json ({case.get('//', '').strip()}). Угловая ДН рассчитана "
                f"в Ansys HFSS (метод конечных элементов) — массив {dn_key}, СРЕЗ E-плоскости "
                f"(φ=0; H-плоскости в этом наборе нет). При выборе кривая HFSS накладывается на "
                f"E-плоскость диаграммы (дБ) рядом с аналитикой 01_sphere. k0 = {case['k0']}.")
        if is_metal:
            desc += (" PEC моделируется eps=1e7 (eps=1e40 из источника переполняет ряд ВСВФ; "
                     "предел тот же).")
        out.append({
            "id": pid,
            "title": base_title + " (Ansys HFSS, E)",
            "group": "Сравнение с Ansys HFSS (МКЭ)",
            "description": desc,
            "scene": _scene(_sphere(_case_layers(case)), k=_eval_k0(case["k0"]), toch=None),
            "reference": _ref_dn(arrays[dn_key]),
        })
    return out


def _example2_preset():
    """Линза Люнеберга 4 слоя (Example 2): HFSS по ОБЕИМ плоскостям E (φ=0) и H (φ=90)."""
    try:
        refs = [
            _ref_hfss_csv("Ephi=0_from_HFSS.csv", "E", "Ansys HFSS, МКЭ (E-плоскость)"),
            _ref_hfss_csv("Ephi=90_from_HFSS.csv", "H", "Ansys HFSS, МКЭ (H-плоскость)"),
        ]
    except Exception:
        return []
    return [{
        "id": "luneburg4_hfss_eh",
        "title": "Линза Люнеберга, 4 слоя — обе плоскости (Ansys HFSS, E+H)",
        "group": "Сравнение с Ansys HFSS (МКЭ)",
        "description": "4-слойная линза Люнеберга, k0=6π. Полная угловая ДН аналитики "
                       "01_sphere сравнивается с Ansys HFSS (метод конечных элементов) в ОБЕИХ "
                       "главных плоскостях: E (φ=0) и H (φ=90). Кривые HFSS "
                       "(dB20normalize(rETotal), Example 2) накладываются на соответствующие "
                       "плоскости диаграммы. Совпадение: медиана ~0.2 дБ (E), ~0.6 дБ (H).",
        "scene": _scene(_sphere([
            (0.25, 1.94, 0.0), (0.50, 1.75, 0.0), (0.75, 1.44, 0.0), (1.00, 1.00, 0.0),
        ]), k=6.0 * _PI, toch=40),
        "reference": refs,
    }]


# --------------------------------------------------------------------------- #
# Полусфера на проводящем экране (метод зеркальных изображений)
# --------------------------------------------------------------------------- #
def _hemisphere_preset():
    """Полусферическая ЛЛ на бесконечном PEC-экране со смещённым облучателем."""
    body = _sphere([
        (0.25, 1.94, 0.0), (0.50, 1.75, 0.0), (0.75, 1.44, 0.0), (1.00, 1.00, 0.0),
    ])
    body["hemisphere"] = {"enabled": True, "feed_offset_deg": 30.0}
    sc = _scene(body, k=6.0 * _PI, toch=40)
    sc["sweep"]["enable"] = False           # спектр к полупространственной задаче неприменим
    return [{
        "id": "hemisphere_luneburg",
        "title": "Полусфера: линза Люнеберга на PEC-экране (облучатель θ′ = 30°)",
        "group": "Демонстрации",
        "description": "Полусферическая 4-слойная линза Люнеберга (k0=6π) на бесконечном "
                       "проводящем экране. Решение — методом зеркальных изображений: "
                       "E_θ(θ)=f_θ(θ−θ′)+f_θ(π−θ−θ′), E_φ(θ)=f_φ(θ−θ′)−f_φ(π−θ−θ′); "
                       "эквивалентно отбору сферических гармоник по чётности "
                       "(τ_n и π_n имеют противоположную чётность при θ→π−θ). "
                       "Точечный облучатель (красный маркер на 3D-сцене) смещён на θ′=30° "
                       "от нормали — луч выходит зеркально, в −30°, практически без потерь "
                       "сканирования (свойство линзы Люнеберга). Ниже экрана поле "
                       "тождественно 0; граничное условие E_tan=0 выполняется автоматически. "
                       "θ′ меняется в «Свойствах тела»; ДН накладывается на 3D-модель "
                       "(чекбокс «ДН» на 3D-сцене). Мощность, уходившая у полной сферы в "
                       "нижнее полупространство, целиком переносится зеркальным членом в "
                       "верхнее (энергобаланс проверен на тождестве P_i_up = P_d_down).",
        "scene": sc,
    }]


# csc²-пресет: механизм диссертации Денисова (гл. «Модифицированные линзы») и статьи
# ВАК-2021 «Формирование ДН специального типа с помощью полусферической ЛЛ» —
# косекансная ДН задаётся ОДНИМ облучателем, наклонённым к экрану: θ′ = 83° от
# нормали (7° от плоскости экрана). Линза — семислойная ЛЛ из той же работы
# (радиус 12 см, f = 10 ГГц → k0a = 8π). Лучший экран по статье — d = 8λ = 240 мм.
_CSC_FEED_DEG = 83.0
_CSC_LENS = [
    (0.39, 1.93, 0.0), (0.56, 1.77, 0.0), (0.68, 1.61, 0.0), (0.78, 1.46, 0.0),
    (0.88, 1.31, 0.0), (0.96, 1.16, 0.0), (1.00, 1.00, 0.0),
]
# HFSS-кривая лучшего экрана из данных статьи (угол в CSV — относительно максимума луча)
_CSC_HFSS_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..",
    "_dissSource", "publicationDenisov",
    "VAK_2021_Формирование диаграммы направленности специального типа "
    "с помощью полусферической ЛЛ", "Gain Plot 1.csv")
_CSC_HFSS_D_KEY = "0.2398339664"        # столбец R_ekran: экран d = 240 мм = 8λ


def _csc_ideal_reference():
    """Идеальный косекансный закон E ∝ csc(ε), 0 дБ на элевации луча ε0 = 90°−θ′.

    θ отсчитан от нормали к экрану; луч выходит зеркально облучателю (θ = −θ′),
    сегмент — от луча до элевации 50° (сектор аппроксимации из статьи ВАК-2021)."""
    e0 = 90.0 - _CSC_FEED_DEG
    theta, db = [], []
    elev = e0
    while elev <= 50.0 + 1e-9:
        theta.append(round(360.0 - (90.0 - elev), 2))     # abs: −(90−elev) → 270..360
        db.append(round(20.0 * math.log10(math.sin(math.radians(e0)) /
                                          math.sin(math.radians(elev))), 3))
        elev += 0.5
    order = sorted(range(len(theta)), key=lambda j: theta[j])
    return {"plane": "E", "name": "Идеальный csc-закон (E ∝ csc ε)",
            "theta_deg": [theta[j] for j in order], "dB": [db[j] for j in order]}


def _csc_hfss_reference():
    """Кривая Ansys HFSS для конечного экрана d = 8λ (данные статьи ВАК-2021).

    В CSV угол отсчитан от направления максимума луча; в абсолютную сетку Studio
    кривая ставится пиком в θ = −θ′ = −83° (луч зеркален облучателю). None при
    отсутствии файла данных (regen вне репо диссертации)."""
    try:
        with open(_CSC_HFSS_CSV, encoding="utf-8") as f:
            rows = [r for r in csv.DictReader(f) if r["R_ekran []"] == _CSC_HFSS_D_KEY]
    except OSError:
        return None
    if not rows:
        return None
    pts = []
    for r in rows:
        rel = float(r["Theta [deg]"])                     # угол от максимума луча
        signed = -_CSC_FEED_DEG + rel                     # abs со знаком: пик в −83°
        if -90.0 <= signed <= 45.0:                       # верхнее полупространство, сектор статьи
            pts.append((signed, float(r["dB10normalize(GainTotal) []"])))
    # порядок — по знаковому углу (непрерывный путь через 0°); mod 360 только при выводе,
    # иначе сортировка по θ∈[0,360) порвала бы кривую хордой через центр диаграммы
    pts.sort()
    return {"plane": "E", "name": "Ansys HFSS: конечный экран d = 8λ (ВАК-2021)",
            "color": "#ff7f0e", "dash": "dot",
            "theta_deg": [round((t + 360.0) % 360.0, 2) for t, _ in pts],
            "dB": [round(v, 3) for _, v in pts]}


def _hemisphere_csc_preset():
    """Косеканс-квадратная ДН: наклон одиночного облучателя к экрану (θ′ = 83°)."""
    body = _sphere(_CSC_LENS)
    body["hemisphere"] = {"enabled": True,
                          "feeds": [{"offset_deg": _CSC_FEED_DEG, "amp": 1.0,
                                     "phase_deg": 0.0}]}
    sc = _scene(body, k=8.0 * _PI, toch=40)
    sc["sweep"]["enable"] = False
    refs = [_csc_ideal_reference()]
    hfss = _csc_hfss_reference()
    if hfss:
        refs.append(hfss)
    return [{
        "id": "hemisphere_csc",
        "title": "Полусфера: косеканс-квадратная (csc²) ДН — наклон облучателя к экрану",
        "group": "Демонстрации",
        "description": "Косеканс-квадратная ДН — классика обзорных РЛС: E(ε) ∝ csc(ε) даёт "
                       "постоянную мощность эха от цели на фиксированной высоте независимо "
                       "от дальности. На полусферической ЛЛ она задаётся УГЛОМ ОТКЛОНЕНИЯ "
                       "ОДИНОЧНОГО ОБЛУЧАТЕЛЯ ОТ ЭКРАНА (механизм диссертации Денисова и "
                       "статьи ВАК-2021): облучатель наклонён к экрану на θ′ = 83° от нормали "
                       "(7° от плоскости экрана), луч выходит зеркально (θ = −83°), прижат к "
                       "горизонту, а его склон к зениту вместе с интерференцией зеркального "
                       "образа формирует косекансный спад. Линза — семислойная ЛЛ из той же "
                       "работы (a_l = 0.39…1.0; ε_l = 1.93…1.0; k0a = 8π). Красный пунктир — "
                       "идеальный csc-закон; оранжевая точечная — Ansys HFSS для конечного "
                       "экрана d = 8λ = 240 мм (наилучшее приближение по статье). Модель "
                       "Studio — идеализация с бесконечным экраном и плосковолновым "
                       "облучателем, поэтому сравнение с HFSS качественное: положение луча и "
                       "характер косекансного склона воспроизводятся, тонкая структура "
                       "лепестков зависит от конечного экрана и реального облучателя. "
                       "θ′ редактируется в «Свойствах тела»; ДН накладывается на 3D-модель "
                       "(чекбокс «ДН»).",
        "scene": sc,
        "reference": refs,
    }]


# csc² на ПОЛНОЙ сфере группой облучателей: «второе решение» диссертации —
# MIMO-ЛЛ (гл. «Многолучевая антенна»): 12 облучателей вдоль дуги с шагом 6°,
# косекансный режим по таблице «Режимы работы первичных облучателей»:
# мощности источников 1–8 = 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1 Вт
# (амплитуды = √P), источники 9–12 выключены. Линза — шестислойная
# оптимизированная (гл. «Аппроксимация»); f = 24.45 ГГц, радиус 75 мм →
# k0a = 2π·75/12.265 ≈ 38.42. Луч каждого облучателя выходит напротив него
# (θ = θ′ − 180°); источник 1 (самый мощный) — у горизонта (депрессия 21°),
# источник 8 — депрессия 63°.
_MIMO_POWERS = [1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
_MIMO_LENS = [
    (0.43, 1.91, 0.0), (0.60, 1.73, 0.0), (0.74, 1.55, 0.0),
    (0.85, 1.36, 0.0), (0.95, 1.18, 0.0), (1.00, 1.00, 0.0),
]
_MIMO_K = 38.42


def _mimo_feeds():
    return [{"offset_deg": round(69.0 - 6.0 * i, 1),
             "amp": round(math.sqrt(p), 4), "phase_deg": 0.0}
            for i, p in enumerate(_MIMO_POWERS)]


def _mimo_csc_reference():
    """Идеальный csc-закон в секторе лучей MIMO-режима, 0 дБ на депрессии 21°.

    Депрессия δ = 90° − θ′; направление луча θ = 270° − δ (mod 360)."""
    theta, db = [], []
    d = 21.0
    while d <= 65.0 + 1e-9:
        theta.append(round(270.0 - d, 2))
        db.append(round(20.0 * math.log10(math.sin(math.radians(21.0)) /
                                          math.sin(math.radians(d))), 3))
        d += 0.5
    order = sorted(range(len(theta)), key=lambda j: theta[j])
    return [{"plane": "E", "name": "Идеальный csc-закон (E ∝ csc δ)",
             "theta_deg": [theta[j] for j in order], "dB": [db[j] for j in order]}]


def _sphere_csc_mimo_preset():
    """csc²-ДН на полной сферической ЛЛ группой облучателей (MIMO-режим)."""
    body = _sphere(_MIMO_LENS)
    body["hemisphere"] = {"enabled": True, "screen": False, "feeds": _mimo_feeds()}
    sc = _scene(body, k=_MIMO_K, toch=None)
    sc["sweep"]["enable"] = False
    return [{
        "id": "sphere_csc_mimo",
        "title": "Сфера: csc²-ДН группой облучателей (MIMO-ЛЛ, реж. 1.0…0.1 Вт)",
        "group": "Демонстрации",
        "description": "Второй способ формирования косеканс-квадратной ДН из диссертации — "
                       "на ПОЛНОЙ сферической линзе Люнеберга группой облучателей "
                       "(MIMO-ЛЛ для 5G, гл. «Многолучевая антенна»): 12 первичных "
                       "облучателей вдоль дуги с шагом 6°, косекансный режим по таблице "
                       "«Режимы работы первичных облучателей» — мощности источников 1–8: "
                       "1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1 Вт (амплитуды √P), "
                       "источники 9–12 выключены. Линза — шестислойная оптимизированная "
                       "аппроксимация (a_l = 0.43…1.0; ε_l = 1.91…1.0), f = 24.45 ГГц, "
                       "радиус 75 мм (k0a ≈ 38.4; ширина луча ≈ 6° — под шаг решётки). "
                       "Каждый облучатель даёт «карандашный» луч напротив себя "
                       "(θ = θ′ − 180°): источник 1 — у горизонта (депрессия 21°), "
                       "источник 8 — депрессия 63°; спад мощностей выстраивает огибающую "
                       "косекансного типа в секторе 21°–63° (красный пунктир — идеальный "
                       "csc-закон; края сектора совпадают с ним точно, середина на "
                       "2–3 дБ выше — режим воспроизводит таблицу диссертации без "
                       "дополнительной подгонки). В отличие от полусферного пресета здесь "
                       "нет экрана — форма ДН перестраивается электронно, без механики: "
                       "веса редактируются в «Свойствах тела» (режим «сфера с линейкой "
                       "облучателей»), ДН накладывается на 3D-модель (чекбокс «ДН»).",
        "scene": sc,
        "reference": _mimo_csc_reference(),
    }]


# --------------------------------------------------------------------------- #
# Демонстрации возможностей (без данных HFSS)
# --------------------------------------------------------------------------- #
_DEMO_PRESETS = [
    {
        "id": "maxwell7",
        "title": "Линза Максвелла («рыбий глаз») — 7 слоёв",
        "group": "Демонстрации",
        "description": "7-слойная аппроксимация линзы Максвелла «рыбий глаз», k0=π "
                       "(точное ядро слоистой сферы 01_sphere).",
        "scene": _scene(_sphere([
            (0.15, 3.90, 0.0), (0.30, 3.52, 0.0), (0.45, 2.98, 0.0), (0.63, 2.48, 0.0),
            (0.82, 1.76, 0.0), (1.00, 1.50, 0.0),
        ]), k=_PI, toch=20),
    },
    {
        "id": "cylinder_inf_coated",
        "title": "Покрытый бесконечный цилиндр (ядро/оболочка)",
        "group": "Демонстрации",
        "description": "Бесконечный радиально-слоистый круговой цилиндр (R=0.5): ядро ε=4 "
                       "(до a=0.6) под стеклоподобной оболочкой ε=2.25, k0=4, нормальное "
                       "падение. Строгий метод ТФГ/эквивалентных линий передачи "
                       "(Дайлис–Шабунин 2017): поляризации TM и TE расцеплены, точное "
                       "сохранение энергии (Q_abs≈0). Выводятся эффективности Q_sca/Q_ext "
                       "(на диаметр) и азимутальная ДН |T(φ)| для TM и TE.",
        "scene": _scene(_cyl_inf(0.5, [
            (0.60, 4.00, 0.0), (1.00, 2.25, 0.0),
        ]), k=4.0),
    },
    {
        "id": "two_spheres_gmm",
        "title": "Кластер двух сфер (GMM)",
        "group": "Демонстрации",
        "description": "Две одинаковые диэлектрические сферы (R=0.5, ε=2.25 без потерь) на оси z "
                       "(центры z=±1.5; описывающие сферы НЕ пересекаются, зазор 2.0), k0=3, "
                       "поперечное падение k̂∥x. Многократное рассеяние решается обобщённым "
                       "методом Ми (GMM) через ТОЧНУЮ теорему сложения Крузана–Стейна: каждая "
                       "сфера — точная T-матрица канонического ядра 01_sphere, трансляции сшивают "
                       "связь. Это СТРОГОЕ аналитическое решение (не приближение). Проверено "
                       "прогоном: энергобаланс |C_abs|/C_sca≈2e-3 (среда без потерь; остаток — "
                       "усечение ряда GMM).",
        "scene": _scene([
            _sphere([(1.00, 2.25, 0.0)], position=(0.0, 0.0, 1.5), radius=0.5),
            _sphere([(1.00, 2.25, 0.0)], position=(0.0, 0.0, -1.5), radius=0.5),
        ], k=3.0, khat=(1.0, 0.0, 0.0), pol=(0.0, 0.0, 1.0)),
    },
]


# Example 2 (обе плоскости) первым — это самый полный кейс сравнения с HFSS.
# Антенные пресеты (полусфера/сфера с облучателями) — первыми среди демонстраций.
PRESETS = (_example2_preset() + _examples_presets()
           + _hemisphere_preset() + _hemisphere_csc_preset()
           + _sphere_csc_mimo_preset() + _DEMO_PRESETS)


def list_presets():
    """Список пресетов для UI (id, title, group, description, scene, [reference])."""
    return PRESETS


def preset_by_id(pid: str):
    for p in PRESETS:
        if p["id"] == pid:
            return p
    return None
