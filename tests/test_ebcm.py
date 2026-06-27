# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""Калибровка EBCM-ядра (green_tensor/ebcm.py) по точному пределу сферы.

Метод нулевого поля для тела вращения обязан в пределе r(θ)=a давать диагональную
T-матрицу Ми T=diag(−a_n,−b_n) (в конвенции проекта — сопряжённой). Проверяем:
  [1] диэлектрическая сфера: диагональ EBCM == t_vector сферы (~1e-12); TE/TM расцеплены;
  [2] магнитная сфера (μ≠1): то же (импедансный множитель s_E=√(ε/μ));
  [3] сечения EBCM-сферы через GMM == аналитика Ми;
  [4] drop-in: EBCM-сфера ≡ scatterer.LayeredSphere в связанном кластере (бит-в-бит).

Запуск: python3 tests/test_ebcm.py | pytest tests/test_ebcm.py
"""
from __future__ import annotations

import cmath
import math
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _ROOT)

from green_tensor import ebcm, gmm, vswf  # noqa: E402
from green_tensor.scatterer import LayeredSphere  # noqa: E402
from analytic_mie import q_sca, q_ext  # noqa: E402


def _proj_diag(sphere, k, nmax):
    """Диагональ (по n) полной T сферы из t_vector проекта: (tM_n, tN_n)."""
    tv = np.asarray(sphere.t_vector(k, nmax))
    modes = vswf.mode_list(nmax)
    K = len(modes)
    seen, tM, tN = set(), {}, {}
    for i, (n, m) in enumerate(modes):
        if n not in seen:
            tM[n] = tv[i]
            tN[n] = tv[K + i]
            seen.add(n)
    ns = sorted(tM)
    return ns, np.array([tM[n] for n in ns]), np.array([tN[n] for n in ns])


def _cross_block_max(T, nmax):
    """Максимум модуля кросс-блоков M↔N (должны быть ~0 для сферы)."""
    K = len(vswf.mode_list(nmax))
    return max(np.max(np.abs(T[:K, K:])), np.max(np.abs(T[K:, :K])))


def _offdiag_n_max(T, nmax):
    """Максимум недиагональных по n элементов в M-блоке (должны быть ~0 для сферы)."""
    K = len(vswf.mode_list(nmax))
    MM = T[:K, :K]
    return np.max(np.abs(MM - np.diag(np.diag(MM))))


def test_ebcm_sphere_matches_mie():
    print("\n[1] EBCM диэлектрической сферы → Ми (диагональ + расцепление TE/TM):")
    for a, k, eps, nmax in [(0.5, 2.0, 2.25, 6), (1.0, 3.0, 2.5, 8)]:
        T = ebcm.tmatrix_axisym(ebcm.sphere_curve(a), k, eps, 1.0, nmax)
        ns, tM, tN = ebcm.diagonal_from_tmatrix(T, nmax)
        _, pM, pN = _proj_diag(LayeredSphere([0, 0, 0], a, [eps]), k, nmax)
        dM = np.max(np.abs(tM - pM)) / np.max(np.abs(pM))
        dN = np.max(np.abs(tN - pN)) / np.max(np.abs(pN))
        cross = _cross_block_max(T, nmax)
        offn = _offdiag_n_max(T, nmax)
        print(f"    a={a} k={k} eps={eps}: dtM={dM:.1e} dtN={dN:.1e} "
              f"кросс M↔N={cross:.1e} внедиаг.n={offn:.1e}")
        assert dM < 1e-11 and dN < 1e-11, "EBCM-сфера разошлась с Ми"
        assert cross < 1e-10, "сфера: TE/TM должны быть расцеплены (кросс-блоки ~0)"
        assert offn < 1e-10, "сфера: T должна быть диагональна по n"


def test_ebcm_magnetic_sphere_matches_mie():
    print("\n[2] EBCM магнитной сферы (μ≠1) → Ми (импеданс s_E=√(ε/μ)):")
    for a, k, eps, mu, nmax in [(0.5, 2.0, 2.25, 1.7, 6), (0.6, 2.5, 3.0, 2.0, 7)]:
        T = ebcm.tmatrix_axisym(ebcm.sphere_curve(a), k, eps, mu, nmax)
        ns, tM, tN = ebcm.diagonal_from_tmatrix(T, nmax)
        _, pM, pN = _proj_diag(LayeredSphere([0, 0, 0], a, [eps], miy=[mu]), k, nmax)
        dM = np.max(np.abs(tM - pM)) / np.max(np.abs(pM))
        dN = np.max(np.abs(tN - pN)) / np.max(np.abs(pN))
        print(f"    a={a} k={k} eps={eps} mu={mu}: dtM={dM:.1e} dtN={dN:.1e}")
        assert dM < 1e-11 and dN < 1e-11, "магнитная EBCM-сфера разошлась с Ми"


def test_ebcm_sphere_cross_sections_vs_mie():
    print("\n[3] сечения EBCM-сферы через GMM == аналитика Ми:")

    class _TW:
        def __init__(self, T):
            self._T = T
            self.position = np.zeros(3)

        def t_matrix(self, k, nmax):
            return self._T

    for a, k, eps, nmax in [(0.5, 2.0, 2.25, 6), (0.5, 2.0, 2.25 + 0.3j, 6)]:
        x = k * a
        T = ebcm.tmatrix_axisym(ebcm.sphere_curve(a), k, eps, 1.0, nmax)
        s = _TW(T)
        c = gmm.isolated_scattered(s, k, (1, 0, 0), (0, 0, 1), nmax)
        cs = gmm.cross_sections([s], [c], k, (1, 0, 0), (0, 0, 1))
        area = math.pi * a * a
        ref_sca = q_sca(cmath.sqrt(eps), x, nmax)
        ref_ext = q_ext(cmath.sqrt(eps), x, nmax)
        d_sca = abs(cs["c_sca"] / area - ref_sca) / ref_sca
        d_ext = abs(cs["c_ext"] / area - ref_ext) / abs(ref_ext)
        print(f"    eps={eps!s:>12}: dQsca={d_sca:.1e} dQext={d_ext:.1e}")
        assert d_sca < 2e-3 and d_ext < 2e-3, "сечения EBCM-сферы не совпали с Ми"


def test_ebcm_drop_in_cluster():
    print("\n[4] drop-in: EBCM-сфера ≡ LayeredSphere в связанном кластере:")

    class _TW:
        def __init__(self, T, pos):
            self._T = T
            self.position = np.asarray(pos, float)

        def t_matrix(self, k, nmax):
            return self._T

    a, k, eps, nmax = 0.5, 2.0, 2.25, 6
    pos = [[0, 0, 1.0], [0, 0, -1.0]]
    khat, pol = (1, 0, 0), (0, 0, 1)
    T = ebcm.tmatrix_axisym(ebcm.sphere_curve(a), k, eps, 1.0, nmax)
    ref = gmm.solve_cluster([LayeredSphere(p, a, [eps]) for p in pos], k, khat, pol, nmax)
    mix = gmm.solve_cluster([_TW(T, pos[0]), LayeredSphere(pos[1], a, [eps])], k, khat, pol, nmax)
    rel = np.max(np.abs(mix[1] - ref[1])) / np.max(np.abs(ref[1]))
    print(f"    max rel |Δc| = {rel:.1e}")
    assert rel < 1e-11, "EBCM-сфера должна быть взаимозаменяема со сферой проекта"


def _energy_balance(T, k, nmax, pol):
    """(C_sca, |C_abs|/C_sca) одиночного тела с T-матрицей при падении k̂=x, поляр. pol."""

    class _TW:
        def __init__(self, T):
            self._T = T
            self.position = np.zeros(3)

        def t_matrix(self, k, nmax):
            return self._T

    s = _TW(T)
    c = gmm.isolated_scattered(s, k, (1, 0, 0), pol, nmax)
    cs = gmm.cross_sections([s], [c], k, (1, 0, 0), pol)
    return cs["c_sca"], abs(cs["c_abs"]) / cs["c_sca"]


def test_ebcm_spheroid_energy_and_rayleigh():
    """Гладкий сфероид (prolate/oblate, диэлектрик и магнит): энергобаланс на уровне
    квадратуры (lossless ⇒ C_abs≈0) + рэлеевский предел == аналит. поляризуемость.
    Это тест НЕсферической машины (внедиаг. по n + кросс-блоки TE↔TM)."""
    import math
    from green_tensor import ellipsoid
    print("\n[5] EBCM сфероид — энергобаланс + Рэлей:")
    for aeq, cax, k, mu, nm, tag in [(0.4, 0.8, 2.0, 1.0, 7, "prolate"),
                                     (0.8, 0.4, 2.0, 1.0, 7, "oblate"),
                                     (0.4, 0.8, 2.0, 1.7, 7, "prolate-mag")]:
        T = ebcm.tmatrix_axisym(ebcm.spheroid_curve(aeq, cax), k, 2.25, mu, nm)
        bal = max(_energy_balance(T, k, nm, (0, 0, 1))[1], _energy_balance(T, k, nm, (0, 1, 0))[1])
        print(f"    {tag} a_eq={aeq} c_ax={cax} mu={mu}: |C_abs|/C_sca={bal:.1e}")
        assert bal < 1e-3, f"сфероид: энергобаланс нарушен ({bal:.1e}) — кросс-блоки/нормировка"
    # Рэлеевский предел: малый сфероид → аналит. поляризуемость (осевая/экваториальная)
    aeq, cax, k, nm, eps = 0.05, 0.10, 1.0, 4, 3.0
    al = ellipsoid.polarizability_homogeneous(aeq, aeq, cax, eps, 1.0)  # сфероид = эллипсоид b=a
    T = ebcm.tmatrix_axisym(ebcm.spheroid_curve(aeq, cax), k, eps, 1.0, nm)
    for pol, ai, lbl in [((0, 0, 1), al[2], "axial"), ((0, 1, 0), al[0], "equat")]:
        csca = _energy_balance(T, k, nm, pol)[0]
        cray = k ** 4 / (6.0 * math.pi) * abs(ai) ** 2
        d = abs(csca - cray) / cray
        print(f"    Рэлей {lbl}: C_sca={csca:.3e} аналит={cray:.3e} d={d:.1e}")
        assert d < 5e-3, f"Рэлеевский предел сфероида разошёлся ({d:.1e})"


def test_ebcm_cone_energy():
    """Конечный конус (кусочная образующая: боковина + дно): энергобаланс lossless≈0.
    Острая вершина/ребро — но для умеренных конусов EBCM сходится."""
    print("\n[6] EBCM конус — энергобаланс:")
    for R, L, k, nm in [(0.4, 0.8, 2.0, 9), (0.5, 0.7, 2.0, 9), (0.3, 0.9, 2.0, 9)]:
        T = ebcm.tmatrix_axisym_segments(ebcm.cone_segments(R, L), k, 2.25, 1.0, nm)
        bal = max(_energy_balance(T, k, nm, (0, 0, 1))[1], _energy_balance(T, k, nm, (0, 1, 0))[1])
        print(f"    R={R} L={L} k={k}: |C_abs|/C_sca={bal:.1e}")
        assert bal < 1e-2, f"конус: энергобаланс нарушен ({bal:.1e})"


def test_ebcm_cylinder_energy():
    """Конечный цилиндр (кусочно: крышки + боковина). Проверяем СХОДЯЩИЙСЯ режим
    (вытянутый, умеренный размер): энергобаланс lossless ≈0.

    ВАЖНО (ограничение метода, не бага): у конечного цилиндра ДВА острых ребра 90°
    (условие Мейкснера). Классический EBCM в double precision для цилиндра сходится
    лишь до ~1e-3…1e-2 и чувствителен к параметрам (для плоского/кубического аспекта
    R≳H вообще расходится — энергобаланс растёт с nmax). Высокая точность для цилиндра
    требует расширенной точности (mpmath) при обращении Q, либо метода инвариантного
    вложения (IITM), либо дискретных источников (NFM-DS). См. docstring ebcm и память
    greentensor-ebcm-spec. Здесь — устойчивый сходящийся случай."""
    print("\n[7] EBCM конечный цилиндр (сходящийся режим) — энергобаланс:")
    R, H, k = 0.3, 0.8, 2.5                               # устойчиво ~5e-3 по nmax 8..10
    for nm in (8, 9, 10):
        T = ebcm.tmatrix_axisym_segments(ebcm.cylinder_segments(R, H), k, 2.25, 1.0, nm)
        bal = max(_energy_balance(T, k, nm, (0, 0, 1))[1], _energy_balance(T, k, nm, (0, 1, 0))[1])
        print(f"    R={R} H={H} k={k} nmax={nm}: |C_abs|/C_sca={bal:.1e}")
        assert bal < 1.2e-2, f"цилиндр (сходящийся режим): энергобаланс нарушен ({bal:.1e})"


def test_ebcm_layered_coated_sphere():
    """Слоистая рекурсия (Петерсон–Стрём) на ПОКРЫТОЙ сфере → покрытые коэф. Ми
    (mie_core, Аден–Керкер). Калибрует рекурсию на точном пределе."""
    from green_tensor.mie_core import MieSphere
    print("\n[8] EBCM слоистая рекурсия: покрытая сфера → Ми:")
    for R, a1, ec, es, k in [(0.5, 0.6, 2.5, 4.0, 2.0), (0.6, 0.5, 3.0, 2.0, 2.5)]:
        nmax = 6
        T = ebcm.tmatrix_axisym_layered(ebcm.sphere_curve(R), [ec, es], [1.0, 1.0], [a1, 1.0], k, nmax)
        ref = MieSphere(k0=k * R, a=[a1, 1.0], eps=[ec, es], miy=[1.0, 1.0], toch=nmax)
        Mn, Nn = ref.coefficients()
        ns, tM, tN = ebcm.diagonal_from_tmatrix(T, nmax)
        dN = np.max(np.abs(tN - (-Mn))) / np.max(np.abs(Mn))
        dM = np.max(np.abs(tM - (-Nn))) / np.max(np.abs(Nn))
        print(f"    R={R} a1={a1} eps=[{ec},{es}] k={k}: dN={dN:.1e} dM={dM:.1e}")
        assert dN < 1e-11 and dM < 1e-11, f"покрытая сфера разошлась с Ми ({dN:.1e})"


def test_ebcm_layered_reduction():
    """Редукция: слоистое тело с ОДИНАКОВЫМИ слоями == однородное (тот же T)."""
    print("\n[9] EBCM слоистая редукция (одинаковые слои == однородное):")
    k0, nmax = 2.0, 7
    for curve, tag in [(ebcm.sphere_curve(0.5), "sphere"), (ebcm.spheroid_curve(0.4, 0.8), "spheroid")]:
        Thom = ebcm.tmatrix_axisym(curve, k0, 2.25, 1.0, nmax)
        Tlay = ebcm.tmatrix_axisym_layered(curve, [2.25, 2.25], [1.0, 1.0], [0.5, 1.0], k0, nmax)
        d = np.max(np.abs(Tlay - Thom))
        print(f"    {tag}: max|T_layered_identical − T_homog| = {d:.1e}")
        assert d < 1e-4, f"редукция слоистого к однородному нарушена ({d:.1e})"


def test_ebcm_layered_spheroid_energy():
    """Слоистый сфероид (ядро/оболочка), lossless энергобаланс. ВНИМАНИЕ: T-матричная
    рекурсия использует исходящие функции на внутр. границах ⇒ для многих слоёв / малого
    ядра обусловленность падает (тут 2 слоя, умеренное ядро); высокая точность —
    импедансная форма/расш. точность (см. tmatrix_layered docstring)."""
    print("\n[10] EBCM слоистый сфероид — энергобаланс (2 слоя):")
    k0, nmax = 2.0, 7
    curve = ebcm.spheroid_curve(0.4, 0.8)
    for eps, mu, tag in [([4.0, 2.25], [1.0, 1.0], "diel"), ([4.0, 2.25], [1.5, 1.0], "mag-core")]:
        T = ebcm.tmatrix_axisym_layered(curve, eps, mu, [0.5, 1.0], k0, nmax)
        bal = max(_energy_balance(T, k0, nmax, (0, 0, 1))[1], _energy_balance(T, k0, nmax, (0, 1, 0))[1])
        print(f"    {tag} eps={eps} mu={mu}: |C_abs|/C_sca={bal:.1e}")
        assert bal < 2.5e-2, f"слоистый сфероид: энергобаланс {bal:.1e}"


def _recip_amp(bodies, k, nmax, pairs):
    """Физическая взаимность амплитуды рассеяния (конвенция-НЕзависимая):
    S(k̂_i,p̂_i→k̂_s,p̂_s) = S(−k̂_s,p̂_s→−k̂_i,p̂_i). В отличие от симметрии элементов
    сырой T в неунитарном базисе проекта (∮|M|²=n(n+1), завышает невязку ~×4), эта
    проверка измеряет НАПРАВЛЕННУЮ корректность, которой энергобаланс НЕ покрывает.
    Возвращает (медиана, макс) относительной невязки по набору пар направлений."""
    from green_tensor import gmm as _g
    Rf = 2000.0 / k
    errs = []
    for ki, pi, ks, ps in pairs:
        _, c, _ = _g.solve_cluster(bodies, k, tuple(ki), tuple(pi), nmax)
        A = (np.asarray(ps) @ _g.scattered_field(bodies, c, k, np.array([np.asarray(ks) * Rf]))[0]) * Rf
        _, c2, _ = _g.solve_cluster(bodies, k, tuple(-np.asarray(ks)), tuple(ps), nmax)
        B = (np.asarray(pi) @ _g.scattered_field(bodies, c2, k, np.array([-np.asarray(ki) * Rf]))[0]) * Rf
        errs.append(abs(A - B) / max(abs(A), abs(B), 1e-30))
    return float(np.median(errs)), float(max(errs))


def test_ebcm_reciprocity():
    """ВЗАИМНОСТЬ EBCM-несфер — направленная проверка корректности, которой энергобаланс
    (оптическая теорема) структурно НЕ ловит (он самосогласован по построению c=T·d).
    Гладкие тела (сфероид prolate/oblate, конус) и СКРУГЛЁННЫЙ цилиндр взаимны ~1e-3;
    острые рёбра 90° цилиндра ломают взаимность (скругление рёбер — реализованная мера,
    scatterer.FiniteCylinder edge_p)."""
    from green_tensor.scatterer import Spheroid, Cone, FiniteCylinder
    print("\n[11] EBCM взаимность (физическая, амплитуда рассеяния):")
    rng = np.random.default_rng(7)

    def perp(kh):
        a = np.array([1.0, 0, 0]) if abs(kh[2]) > 0.9 else np.array([0, 0, 1.0])
        p = np.cross(kh, a)
        return p / np.linalg.norm(p)

    def rdir():
        v = rng.normal(size=3)
        return v / np.linalg.norm(v)

    pairs = []
    for _ in range(6):
        ki, ks = rdir(), rdir()
        pairs.append((ki, perp(ki), ks, perp(ks)))

    k, nmax = 2.0, 8
    cases = [
        ("spheroid prolate", [Spheroid([0, 0, 0], 0.4, 0.8, 2.25)], 1.0e-3),
        ("spheroid oblate", [Spheroid([0, 0, 0], 0.8, 0.4, 2.25)], 1.0e-3),
        ("cone", [Cone([0, 0, 0], 0.4, 0.8, 2.25)], 1.0e-3),
        ("rounded cylinder", [FiniteCylinder([0, 0, 0], 0.3, 0.8, 2.25)], 3.0e-3),
    ]
    for tag, body, tol in cases:
        med, mx = _recip_amp(body, k, nmax, pairs)
        print(f"    {tag:18s}: median={med:.2e} max={mx:.2e} (tol median<{tol:.0e})")
        assert med < tol, f"{tag}: взаимность нарушена (median={med:.1e})"

    # регрессия: скругление рёбер ДОЛЖНО улучшать взаимность острого цилиндра
    medS, _ = _recip_amp([FiniteCylinder([0, 0, 0], 0.3, 0.8, 2.25, edge_p=None)], k, nmax, pairs)
    medR, _ = _recip_amp([FiniteCylinder([0, 0, 0], 0.3, 0.8, 2.25)], k, nmax, pairs)
    print(f"    цилиндр острый median={medS:.2e}  скруглённый median={medR:.2e} (скругление лучше)")
    assert medR < medS, "скругление рёбер должно улучшать взаимность цилиндра"


def test_ebcm_layered_cylinder_guarded():
    """РЕГРЕССИЯ/ЧЕСТНОСТЬ: слоистый конечный цилиндр через EBCM численно НЕУСТОЙЧИВ
    (рекурсия Петерсона–Стрёма на острых сегментах: исходящие h_n на внутр. границах +
    рёбра ⇒ унитарность нарушена в разы). Решатель ДОЛЖЕН предупреждать (UserWarning), а не
    тихо возвращать нефизичную T. Для слоистого цилиндра используйте бесконечный ТФГ-решатель
    LayeredCylinderSolver или однослойный (скруглённый) FiniteCylinder."""
    import warnings
    from green_tensor.scatterer import FiniteCylinder
    print("\n[12] EBCM слоистый конечный цилиндр — честный guard:")
    fc = FiniteCylinder([0, 0, 0], 0.3, 0.8, [3.0, 2.0], a_norm=[0.6, 1.0])
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        T = fc.t_matrix(2.0, 8)
    modes = vswf.mode_list(8)
    nn = np.array([n * (n + 1) for (n, m) in modes], float)
    D = np.concatenate([nn, nn])
    S = (np.sqrt(D)[:, None]) * T * (1.0 / np.sqrt(D)[None, :])
    uni = np.linalg.norm(S.conj().T @ S + 0.5 * (S + S.conj().T)) / max(np.linalg.norm(S), 1e-30)
    warned = any(issubclass(x.category, UserWarning) for x in w)
    print(f"    унитарность |S†S+ReS|/|S| = {uni:.2f} (нефизично, ≫0); предупреждение={'да' if warned else 'НЕТ'}")
    assert uni > 0.5, "ожидалась явная неустойчивость слоистого острого цилиндра"
    assert warned, "слоистый конечный цилиндр (EBCM) должен предупреждать о неустойчивости"


_TESTS = [test_ebcm_sphere_matches_mie, test_ebcm_magnetic_sphere_matches_mie,
          test_ebcm_sphere_cross_sections_vs_mie, test_ebcm_drop_in_cluster,
          test_ebcm_spheroid_energy_and_rayleigh, test_ebcm_cone_energy,
          test_ebcm_cylinder_energy, test_ebcm_layered_coated_sphere,
          test_ebcm_layered_reduction, test_ebcm_layered_spheroid_energy,
          test_ebcm_reciprocity, test_ebcm_layered_cylinder_guarded]

if __name__ == "__main__":
    ok = True
    for fn in _TESTS:
        try:
            fn()
        except AssertionError as exc:
            print(f"  ❌ {fn.__name__}: {exc}")
            ok = False
    print("\n✅ EBCM проверки пройдены." if ok else "\n❌ Есть провалы.")
    sys.exit(0 if ok else 1)
