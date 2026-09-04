# SPDX-License-Identifier: MIT

"""nearfield — АБСОЛЮТНОЕ поле слоистой сферы в произвольных точках (В/м).

Внешнее разложение по исходящим ВСВФ: коэффициенты рассеяния берутся из
канонического ядра (``01_sphere.py``) через диагональную T-матрицу, возбуждение —
проекция плоской волны амплитуды :math:`E_0`. Оба режима задачи различаются
только коэффициентами (``problem="diffraction"`` / ``"antenna"``), как и в канале
дальней зоны, поэтому карты обоих режимов строятся единым кодом.

Все величины АБСОЛЮТНЫЕ: поле возвращается в В/м при заданной амплитуде
источника :math:`E_0` (ф-ла 3.18 Панченко — см. :meth:`sphere_core.MieSphere.dipole_e0`).
Нормировка на максимум не применяется.

ОГРАНИЧЕНИЕ. Разложение справедливо во ВНЕШНЕЙ области (r > R): внутри тела ряд
исходящих ВСВФ расходится (корректная физика, а не ошибка), поэтому точки r < R
по умолчанию маскируются NaN.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.special as _sp

from ._angular import pi_tau as _pi_tau
from . import gmm as _gmm
from . import tmatrix as _tmatrix
from . import vswf as _vswf

__all__ = ["near_field_sphere", "external_field_sphere",
           "magnetic_field_sphere", "antenna_pattern_sphere"]


def near_field_sphere(radius, eps, k: float, pts, *, e0: complex = 1.0,
                      problem: str = "diffraction", a_norm=None, miy=None,
                      toch: int | None = None, nmax: int | None = None,
                      khat=(0.0, 0.0, 1.0), pol=(1.0, 0.0, 0.0),
                      total: bool = False, mask_inside: bool = True,
                      mode_tol: float = 1e-12) -> np.ndarray:
    r"""Абсолютное поле слоистой сферы в точках ``pts`` [В/м].

    ``radius`` [м], ``eps`` — проницаемости слоёв изнутри наружу, ``k`` —
    ФИЗИЧЕСКОЕ волновое число [рад/м], ``pts`` — массив (N,3) декартовых точек [м].
    ``e0`` — амплитуда источника :math:`E_0` [В/м]. ``problem`` — режим возбуждения.
    ``khat``/``pol`` — направление и поляризация проецируемой волны возбуждения.

    ``total=False`` (по умолчанию) — только рассеянное/излучённое поле;
    ``total=True`` добавляет падающую волну (полное поле, как карта ``|E|`` в HFSS).
    ``mask_inside`` — ставить NaN при r < ``radius`` (внешнее разложение там неприменимо).

    Возвращает комплексный массив (N,3) — декартовы компоненты поля [В/м].
    """
    from .solvers import SphereSolver          # локальный импорт: избегаем цикла

    sph = SphereSolver(radius, eps, a_norm=a_norm, miy=miy)
    mie = sph.mie(k, toch=toch, problem=problem)
    Mn, Nn = mie.coefficients()
    T = _tmatrix.from_mie_coeffs(Mn, Nn, mie.k0)
    if nmax is None:
        nmax = len(Mn)
    nmax = int(min(nmax, len(Mn)))

    modes = _vswf.mode_list(nmax)
    K = len(modes)
    a_inc = _gmm.plane_wave_coeffs(k, khat, np.asarray(pol, complex) * e0,
                                   (0.0, 0.0, 0.0), nmax)
    aM, aN = a_inc[:K], a_inc[K:]

    # рассеянные коэффициенты: c = T·a (T диагональна по n, одинакова по m)
    cM = np.array([T.t_M[n - 1] * aM[i] for i, (n, _m) in enumerate(modes)], complex)
    cN = np.array([T.t_N[n - 1] * aN[i] for i, (n, _m) in enumerate(modes)], complex)

    # отбор значимых мод: плоская волна вдоль оси возбуждает лишь |m|=1
    peak = max(np.max(np.abs(cM)), np.max(np.abs(cN)), 1e-300)
    sel = np.nonzero((np.abs(cM) > mode_tol * peak) | (np.abs(cN) > mode_tol * peak))[0]

    pts = np.asarray(pts, dtype=float)
    E = np.zeros((pts.shape[0], 3), dtype=complex)
    for i in sel:
        n, m = modes[i]
        Mg, Ng = _vswf.mn_grid(n, m, k, pts, "out")
        E += cM[i] * Mg + cN[i] * Ng

    if total:
        kvec = k * np.asarray(khat, dtype=float)
        E = E + np.asarray(pol, complex)[None, :] * e0 * np.exp(1j * (pts @ kvec))[:, None]

    if mask_inside:
        rr = np.linalg.norm(pts, axis=1)
        E[rr < float(radius), :] = np.nan

    return E


# --------------------------------------------------------------------------- #
# Прямое радиальное продолжение ряда дальней зоны (точные h_n, без ВСВФ/GMM)
# --------------------------------------------------------------------------- #
def _hankel_out(nmax: int, rho):
    r"""h_n^{(1)}(rho) и D_n = [rho·h_n]'/rho для n=1..nmax; форма (nmax, len(rho))."""
    n = np.arange(1, nmax + 1)[:, None]
    x = np.asarray(rho, dtype=float)[None, :]
    jn = _sp.spherical_jn(n, x); yn = _sp.spherical_yn(n, x)
    jnp = _sp.spherical_jn(n, x, derivative=True)
    ynp = _sp.spherical_yn(n, x, derivative=True)
    h = jn + 1j * yn
    hp = jnp + 1j * ynp
    return h, h / x + hp


def external_field_sphere(radius, eps, k: float, pts, *, e0: complex = 1.0,
                          problem: str = "diffraction", a_norm=None, miy=None,
                          toch: int | None = None, khat=(0.0, 0.0, 1.0),
                          current: float | None = None, length: float | None = None,
                          z0: float = 376.730313668,
                          total: bool = False, mask_inside: bool = True) -> np.ndarray:
    r"""АБСОЛЮТНОЕ внешнее поле слоистой сферы [В/м] — точный радиальный ряд.

    Продолжение того же ряда, чья асимптотика даёт диаграмму: вместо асимптотики
    подставлены точные радиальные функции :math:`h_n^{(1)}(k r)` и
    :math:`[\rho h_n]'/\rho`, поэтому нормировка наследуется из дальней зоны
    аналитически:

    .. math::
        E_\theta = A\cos\varphi\sum_n i^{n}\tfrac{2n+1}{n(n+1)}
        \Bigl[i a_n\tau_n\tfrac{[\rho h_n]'}{\rho}-b_n\pi_n h_n\Bigr],\quad
        E_\varphi = -A\sin\varphi\sum_n i^{n}\tfrac{2n+1}{n(n+1)}
        \Bigl[i a_n\pi_n\tfrac{[\rho h_n]'}{\rho}-b_n\tau_n h_n\Bigr],

    где :math:`a_n=\overline{M_n}` (электрический, TM), :math:`b_n=\overline{N_n}`
    (магнитный, TE) — коэффициенты канонического ядра для выбранного ``problem``.
    При :math:`kr\to\infty` даёт :math:`|E_\theta|=A|S_2|/(kr)`.

    Амплитуда источника :math:`A`:

    * ``problem="diffraction"``: ``A = e0`` — амплитуда падающей плоской волны [В/м];
      ряд сразу даёт абсолютное рассеянное поле.
    * ``problem="antenna"``: задайте ``current`` (I, А) и ``length`` (l, м) —
      тогда :math:`A = I k^2 l Z_0/(4\pi)`. Это замкнутая форма согласования
      нормировок: ряд несёт множитель :math:`1/(kr)`, а амплитуда элементарного
      вибратора (Панченко, 3.18) — :math:`1/r`, откуда лишний множитель :math:`k`.
      При :math:`kr\to\infty` результат равен :math:`E_0(r)\,|S_2|` с
      :math:`E_0(r)=I k l Z_0/(4\pi r)`.

    ``khat`` — направление распространения возбуждения, поддерживаются ``(0,0,±1)``.
    ``total=True`` добавляет падающую плоскую волну (только для дифракции).
    Точки r < ``radius`` маскируются NaN (внешнее разложение там неприменимо).

    Возвращает комплексный массив (N,3) — ДЕКАРТОВЫ компоненты поля [В/м].
    """
    return _series(radius, eps, k, pts, e0=e0, problem=problem, a_norm=a_norm,
                   miy=miy, toch=toch, khat=khat, current=current, length=length,
                   z0=z0, total=total, mask_inside=mask_inside, want="E")


def magnetic_field_sphere(radius, eps, k: float, pts, *, e0: complex = 1.0,
                          problem: str = "diffraction", a_norm=None, miy=None,
                          toch: int | None = None, khat=(0.0, 0.0, 1.0),
                          current: float | None = None, length: float | None = None,
                          z0: float = 376.730313668,
                          total: bool = False, mask_inside: bool = True) -> np.ndarray:
    r"""АБСОЛЮТНОЕ внешнее МАГНИТНОЕ поле слоистой сферы [А/м].

    Дуальная сборка того же ряда, что и :func:`external_field_sphere`: замена
    :math:`a_n\leftrightarrow b_n`, :math:`\cos\varphi\leftrightarrow\sin\varphi`
    и общий множитель :math:`-1/Z_0`,

    .. math::
        Z_0H_\theta = -A\sin\varphi\sum_n i^{n}\tfrac{2n+1}{n(n+1)}
        \Bigl[i b_n\tau_n\tfrac{[\rho h_n]'}{\rho}-a_n\pi_n h_n\Bigr],\quad
        Z_0H_\varphi = -A\cos\varphi\sum_n i^{n}\tfrac{2n+1}{n(n+1)}
        \Bigl[i b_n\pi_n\tfrac{[\rho h_n]'}{\rho}-a_n\tau_n h_n\Bigr].

    Форма проверена независимо: совпадает с :math:`\nabla\times\mathbf E/(jk_0Z_0)`,
    вычисленным численным дифференцированием аналитического поля, до машинной
    точности, и с соотношением исходящей волны
    :math:`\mathbf H=\hat r\times\mathbf E/Z_0` в дальней зоне.

    Аргументы и соглашения — как у :func:`external_field_sphere`.
    Возвращает комплексный массив (N,3) — ДЕКАРТОВЫ компоненты [А/м].
    """
    return _series(radius, eps, k, pts, e0=e0, problem=problem, a_norm=a_norm,
                   miy=miy, toch=toch, khat=khat, current=current, length=length,
                   z0=z0, total=total, mask_inside=mask_inside, want="H")


def antenna_pattern_sphere(radius, eps, k: float, theta, *,
                           current: float = 1.0, length: float = 1.0,
                           gap: float = 0.0, r: float | None = None,
                           plane: str = "E", feed: str = "huygens",
                           a_norm=None, miy=None, toch: int | None = None,
                           z0: float = 376.730313668) -> dict:
    r"""ДН линзовой антенны с облучателем на поверхности — ТЕОРЕМА ВЗАИМНОСТИ.

    Прямая запись антенного возбуждения (Панченко, ф-ла 13.3) задаёт всем
    гармоникам равную амплитуду на границе и потому описывает не физический
    облучатель, а близкое к оптимальному сверхнаправленное возбуждение. Здесь ДН
    строится иначе — через взаимность, то есть целиком на уже верифицированном
    ДИФРАКЦИОННОМ решении:

    .. math::
        rE(\hat u) = \frac{I^e k_0 l Z_0}{4\pi}\,
        \Bigl|\,\hat p_e\!\cdot\!\mathbf E^{\text{полн}}(\mathbf r_d)
        + Z_0\,\hat p_m\!\cdot\!\mathbf H^{\text{полн}}(\mathbf r_d)\Bigr| ,

    где :math:`\mathbf E^{\text{полн}}`, :math:`\mathbf H^{\text{полн}}` --- полные
    поля в точке облучателя :math:`\mathbf r_d` при ЕДИНИЧНОЙ плоской волне,
    приходящей с направления наблюдения :math:`\hat u`; :math:`\hat p_e` ---
    орт электрического тока (тангенциальный), :math:`\hat p_m=\hat n\times\hat p_e`
    --- орт магнитного тока элемента Гюйгенса. Множитель :math:`1/2` в
    :math:`E_0^{(\text{Гюйг})}/2` следует из предельного перехода к свободному
    пространству: там комбинация в модуле равна :math:`2(1+\cos\Theta)/2`, то есть
    удвоенной кардиоиде элемента Гюйгенса.

    Параметры
    ---------
    theta : углы наблюдения [рад] в той же конвенции, что и ``pattern``
        (:math:`\theta=180^\circ` --- направление, антиподальное облучателю).
    current, length : ток :math:`I^e` [А] и длина :math:`l` [м] облучателя.
    gap : отступ облучателя от поверхности [м] (0 --- вплотную).
    r : если задано, дополнительно возвращается ``E`` [В/м] на расстоянии r.
    plane : ``"E"`` --- срез :math:`\varphi=0`, ``"H"`` --- срез :math:`\varphi=\pi/2`.
    feed : ``"huygens"`` (скрещённые токи, кардиоида) либо ``"dipole"``
        (только электрический ток).
    Реализована ЛИНЕЙНАЯ поляризация (срезы E и H). Компоненты круговой
    поляризации антенного режима не реализованы: доступный полноволновой экспорт
    для них внутренне противоречив (линейно поляризованный облучатель не может
    излучать поле с нулевой одной из круговых составляющих), поэтому проверить
    их нечем.

    Возвращает словарь с ключами ``theta``, ``rE`` [В] и, если задан ``r``, ``E`` [В/м].
    """
    th = np.atleast_1d(np.asarray(theta, dtype=float))
    rd = float(radius) + float(gap)
    st, ct = np.sin(th), np.cos(th)
    zeros = np.zeros_like(th)

    if plane.upper() == "E":
        nhat = np.stack([st, zeros, ct], axis=1)          # положение облучателя
        pe = np.stack([ct, zeros, -st], axis=1)           # электрический ток
    elif plane.upper() == "H":
        nhat = np.stack([zeros, st, ct], axis=1)
        pe = np.stack([np.ones_like(th), zeros, zeros], axis=1)
    else:
        raise ValueError("plane: 'E' или 'H'")
    pm = np.cross(nhat, pe)                               # магнитный ток
    pts = rd * nhat

    common = dict(problem="diffraction", a_norm=a_norm, miy=miy, toch=toch,
                  khat=(0.0, 0.0, -1.0), e0=1.0, z0=z0,
                  total=True, mask_inside=False)

    def weight(P, Pe, Pm):
        """Отклик взаимности при x̂-поляризованной падающей волне."""
        E = _series(radius, eps, k, P, current=None, length=None, want="E", **common)
        v = np.sum(E * Pe, axis=1)
        if feed.lower() == "huygens":
            H = _series(radius, eps, k, P, current=None, length=None, want="H", **common)
            v = v + z0 * np.sum(H * Pm, axis=1)
        elif feed.lower() != "dipole":
            raise ValueError("feed: 'huygens' или 'dipole'")
        return v

    wx = weight(pts, pe, pm)
    pref = current * k * length * z0 / (4.0 * math.pi)    # = |E_0 вибратора| · r

    out = {"theta": th, "rE": pref * np.abs(wx)}
    if r is not None:
        out["E"] = out["rE"] / float(r)
    return out


def _series(radius, eps, k, pts, *, e0, problem, a_norm, miy, toch, khat,
            current, length, z0, total, mask_inside, want):
    """Общее ядро: точный радиальный ряд для E [В/м] либо H [А/м]."""
    from .solvers import SphereSolver

    kz = float(np.asarray(khat, dtype=float)[2])
    if abs(abs(kz) - 1.0) > 1e-12 or abs(np.asarray(khat, float)[0]) > 1e-12 \
       or abs(np.asarray(khat, float)[1]) > 1e-12:
        raise ValueError("external_field_sphere: khat поддерживается только (0,0,±1)")

    mie = SphereSolver(radius, eps, a_norm=a_norm, miy=miy).mie(k, toch=toch, problem=problem)
    Mn, Nn = mie.coefficients()
    an, bn = np.conj(Mn), np.conj(Nn)          # a_n — электрический, b_n — магнитный
    nmax = len(an)

    if problem == "antenna":
        if current is None or length is None:
            raise ValueError("для антенной задачи задайте current [А] и length [м]")
        A = current * k * k * length * z0 / (4.0 * math.pi)
    else:
        A = complex(e0)

    pts = np.asarray(pts, dtype=float)
    # падение вдоль -z = зеркальное отражение задачи с падением вдоль +z
    mirror = (kz < 0.0)
    q = pts.copy()
    if mirror:
        q[:, 2] = -q[:, 2]

    r = np.linalg.norm(q, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        theta = np.arccos(np.clip(np.where(r > 0, q[:, 2] / np.where(r > 0, r, 1.0), 1.0), -1.0, 1.0))
    phi = np.arctan2(q[:, 1], q[:, 0])

    rho = k * r
    safe = r > float(radius)
    rho_s = np.where(safe, rho, np.nan)

    h, D = _hankel_out(nmax, rho_s)                    # (nmax, N)
    pin, taun = _pi_tau(nmax, theta)                   # (nmax, N)
    n = np.arange(1, nmax + 1)[:, None]
    w = (1j ** n) * (2 * n + 1) / (n * (n + 1))        # (nmax,1)

    if want == "E":
        # E: при tau стоит электрический коэффициент, азимут cos / (-sin)
        c1, c2 = an[:, None], bn[:, None]
        f_th, f_ph, scale = np.cos(phi), -np.sin(phi), 1.0
    else:
        # H: дуальная сборка — a<->b, азимут sin / cos, множитель 1/Z0
        c1, c2 = bn[:, None], an[:, None]
        f_th, f_ph, scale = np.sin(phi), np.cos(phi), 1.0 / z0

    Fth = A * scale * np.sum(w * (1j * c1 * taun * D - c2 * pin * h), axis=0) * f_th
    Fph = A * scale * np.sum(w * (1j * c1 * pin * D - c2 * taun * h), axis=0) * f_ph
    Fr = A * scale * np.sum(w * (1j * c1 * n * (n + 1) * np.sin(theta)[None, :] * pin * h / rho_s),
                            axis=0) * f_th

    st, ct = np.sin(theta), np.cos(theta)
    cp, sp_ = np.cos(phi), np.sin(phi)
    Fx = Fr * st * cp + Fth * ct * cp - Fph * sp_
    Fy = Fr * st * sp_ + Fth * ct * sp_ + Fph * cp
    Fz = Fr * ct - Fth * st
    F = np.stack([Fx, Fy, Fz], axis=1)

    if total:
        if problem != "diffraction":
            raise ValueError("total=True имеет смысл только для дифракции")
        # падающая плоская волна E = e0·x̂·e^{ikz}; H = (e0/Z0)·ŷ·e^{ikz}
        pol = np.asarray([1.0, 0.0, 0.0]) if want == "E" else np.asarray([0.0, 1.0, 0.0])
        amp = complex(e0) if want == "E" else complex(e0) / z0
        F = F + pol[None, :] * amp * np.exp(1j * k * q[:, 2])[:, None]

    if mirror:
        # обратное зеркало по z; для H дополнительно общий знак (псевдовектор)
        F[:, 2] = -F[:, 2]
        if want == "H":
            F[:, :] = -F[:, :]
    if mask_inside:
        F[~safe, :] = np.nan
    return F
