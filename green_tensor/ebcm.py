# SPDX-License-Identifier: MIT
# Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

"""ebcm — T-матрица несферических осесимметричных тел методом нулевого поля (EBCM/Waterman).

Строит ПОЛНУЮ (вообще недиагональную) сферическую T-матрицу в базисе ВСВФ проекта
(совместима с gmm.solve_cluster) для тела вращения с образующей r(θ). Метод: нуль-поле
Уотермана (Waterman 1971; Mishchenko, Travis & Lacis 2002, гл.5; Doicu, Wriedt & Eremin
2006):

    T = − RgQ · Q⁻¹ ,

где Q, RgQ — поверхностные интегралы по границе тела от произведений ВСВФ (внутренних
регулярных при k1=k√(εμ) и внешних — исходящих для Q, регулярных для RgQ). Интегралы
берутся численно теми же проверенными эвалуаторами vswf.mn_grid, поэтому нормировка
автоматически согласована с остальной библиотекой; калибровка — точный предел сферы
(T == diag(−a_n,−b_n), см. tests).

Конвенция: e^{-iωt}; исходящая ~ h_n^{(1)}; внешняя среда — вакуум (k); тело — один
однородный слой (eps, mu) [слоистость добавится рекурсией Петерсона–Стрёма позже].
Блочность по m: для тела вращения Q блочно-диагональна по азимуту m (φ-интеграл даёт
δ_{mm'}), поэтому собираем независимый блок на каждый m.
"""
from __future__ import annotations

import math

import numpy as np

from . import vswf


# --------------------------------------------------------------------------- #
# Поверхностная квадратура осесимметричного тела r(θ)
# --------------------------------------------------------------------------- #
def axisym_surface(r_of_theta, ntheta: int, nphi: int):
    """Узлы поверхности тела вращения и взвешенный вектор нормали n̂·dS.

    r_of_theta(theta_array) -> (r, dr/dθ). Возвращает (pts (Q,3), nvec_dS (Q,3)), где
    nvec_dS = n̂ dS уже включает все веса квадратуры. Для интеграла ∮ F·n̂ dS достаточно
    Σ_q nvec_dS[q]·F[q].

    Гаусс–Лежандр по x=cosθ (∫_0^π · sinθ dθ = ∫_{-1}^1 ·dx), равномерно по φ. Для
    осесимметричного r(θ): dS = r√(r²+r_θ²) sinθ dθ dφ, n̂=(r r̂ − r_θ θ̂)/√(r²+r_θ²),
    поэтому n̂ dS = r(r r̂ − r_θ θ̂) sinθ dθ dφ (корень сокращается).
    """
    x, wx = np.polynomial.legendre.leggauss(ntheta)      # x = cosθ ∈ (−1,1)
    theta = np.arccos(x)
    r, drdt = r_of_theta(theta)
    r = np.asarray(r, float)
    drdt = np.asarray(drdt, float)
    st, ct = np.sin(theta), np.cos(theta)

    phi = 2.0 * np.pi * np.arange(nphi) / nphi
    wphi = 2.0 * np.pi / nphi

    TH = np.repeat(theta, nphi)
    PH = np.tile(phi, ntheta)
    RR = np.repeat(r, nphi)
    DR = np.repeat(drdt, nphi)
    WX = np.repeat(wx, nphi)
    sT, cT = np.sin(TH), np.cos(TH)
    cP, sP = np.cos(PH), np.sin(PH)

    pts = (RR[:, None] * np.stack([sT * cP, sT * sP, cT], axis=-1))
    rhat = np.stack([sT * cP, sT * sP, cT], axis=-1)
    that = np.stack([cT * cP, cT * sP, -sT], axis=-1)
    # n̂ dS = r (r r̂ − r_θ θ̂) · (wx · wphi)   [sinθ уже в wx через ∫dx]
    nvec_dS = (RR[:, None] * (RR[:, None] * rhat - DR[:, None] * that)) * (WX * wphi)[:, None]
    return pts, nvec_dS


def sphere_curve(a: float):
    """Образующая сферы r(θ)=a (для калибровки EBCM по точному Ми)."""
    def f(theta):
        theta = np.asarray(theta, float)
        return np.full_like(theta, float(a)), np.zeros_like(theta)
    return f


def spheroid_curve(a_eq: float, c_ax: float):
    """Образующая сфероида (ось симметрии — z): экваториальная полуось a_eq, осевая c_ax.

    Поверхность (ρ/a_eq)²+(z/c_ax)²=1 в полярных координатах:
        r(θ)=[sin²θ/a_eq² + cos²θ/c_ax²]^{−1/2},
        r_θ = −½ r³ sin(2θ)(1/a_eq² − 1/c_ax²).
    prolate (вытянутый): c_ax>a_eq; oblate (сплюснутый): c_ax<a_eq; сфера: a_eq=c_ax.
    Гладкая образующая (без рёбер) — годится для одной квадратуры axisym_surface;
    идеальна для калибровки кросс-блоков TE↔TM EBCM (на сфере они нулевые).
    """
    ia2, ic2 = 1.0 / (a_eq * a_eq), 1.0 / (c_ax * c_ax)

    def f(theta):
        theta = np.asarray(theta, float)
        s2, c2 = np.sin(theta) ** 2, np.cos(theta) ** 2
        r = (s2 * ia2 + c2 * ic2) ** (-0.5)
        drdt = -0.5 * r ** 3 * np.sin(2.0 * theta) * (ia2 - ic2)
        return r, drdt
    return f


# --------------------------------------------------------------------------- #
# Кусочная поверхностная квадратура (тела с рёбрами: цилиндр, конус)
# --------------------------------------------------------------------------- #
def surface_segments(segments, n_per: int, nphi: int):
    """Узлы и n̂·dS для осесимметричной поверхности, заданной КУСОЧНО.

    segments — список (r_of_theta, θ_lo, θ_hi); на каждом гладком сегменте берётся своя
    Гаусс–Лежандр квадратура по θ∈[θ_lo,θ_hi] (n_per узлов), равномерно по φ. Нормаль
    n̂=(r r̂ − r_θ θ̂)/√(...) выходит корректной автоматически из r(θ),r_θ каждого сегмента
    (включая верную сторону для крышек/боковин), поэтому сегменты НЕ нужно сшивать и
    задавать нормаль вручную. n̂ dS = r(r r̂ − r_θ θ̂) sinθ dθ dφ (θ-Гаусс, sinθ — явно).
    """
    x0, w0 = np.polynomial.legendre.leggauss(n_per)       # на [−1,1]
    phi = 2.0 * np.pi * np.arange(nphi) / nphi
    wphi = 2.0 * np.pi / nphi
    pts_all, nv_all = [], []
    for r_of_theta, tlo, thi in segments:
        theta = 0.5 * (thi - tlo) * x0 + 0.5 * (thi + tlo)
        wth = 0.5 * (thi - tlo) * w0                       # вес θ-Гаусса на [tlo,thi]
        r, drdt = r_of_theta(theta)
        r = np.asarray(r, float)
        drdt = np.asarray(drdt, float)
        st = np.sin(theta)
        TH = np.repeat(theta, nphi)
        PH = np.tile(phi, len(theta))
        RR = np.repeat(r, nphi)
        DR = np.repeat(drdt, nphi)
        Wt = np.repeat(wth * st, nphi)                     # sinθ·w_θ
        sT, cT = np.sin(TH), np.cos(TH)
        cP, sP = np.cos(PH), np.sin(PH)
        rhat = np.stack([sT * cP, sT * sP, cT], axis=-1)
        that = np.stack([cT * cP, cT * sP, -sT], axis=-1)
        pts_all.append(RR[:, None] * rhat)
        nv_all.append((RR[:, None] * (RR[:, None] * rhat - DR[:, None] * that)) * (Wt * wphi)[:, None])
    return np.concatenate(pts_all), np.concatenate(nv_all)


def cylinder_segments(radius: float, half_length: float):
    """Образующая КОНЕЧНОГО цилиндра (ось z): радиус R, полудлина H. Три сегмента
    (верхняя крышка / боковина / нижняя крышка), разрез на ребре θ_e=arctan(R/H)."""
    R, H = float(radius), float(half_length)
    te = math.atan2(R, H)

    def cap_top(theta):                  # z=+H: r=H/cosθ
        c = np.cos(theta)
        return H / c, H * np.sin(theta) / c ** 2

    def side(theta):                     # ρ=R: r=R/sinθ
        s = np.sin(theta)
        return R / s, -R * np.cos(theta) / s ** 2

    def cap_bot(theta):                  # z=−H: r=−H/cosθ (cosθ<0 ⇒ r>0)
        c = np.cos(theta)
        return -H / c, H * np.sin(theta) / c ** 2

    return [(cap_top, 0.0, te), (side, te, math.pi - te), (cap_bot, math.pi - te, math.pi)]


def rounded_cylinder_curve(radius: float, half_length: float, *, p: float = 6.0):
    """ГЛАДКАЯ образующая «скруглённого цилиндра» — суперэллипс (ρ/R)^p+(z/H)^p=1.

    Острые рёбра 90° конечного цилиндра дают сингулярность Мейкснера, для которой
    разложение по СФЕРИЧЕСКИМ ВСВФ (EBCM) сходится крайне медленно и теряет взаимность
    (а расширенная точность/квадратура НЕ помогают — это усечение, не обусловленность).
    Скругление рёбер убирает сингулярность ⇒ EBCM сходится быстро и T-матрица взаимна.
    p→∞ — острый цилиндр; p≈6 — цилиндр со слегка скруглёнными рёбрами (рекомендуется).

    Возвращает ``r_of_theta(θ) -> (r, dr/dθ)`` для :func:`tmatrix_axisym` / :func:`axisym_surface`.
    """
    R, H, p = float(radius), float(half_length), float(p)

    def r_of_theta(theta):
        theta = np.asarray(theta, dtype=float)
        s = np.abs(np.sin(theta))
        c = np.abs(np.cos(theta))
        a = (s / R) ** p
        b = (c / H) ** p
        u = a + b
        r = u ** (-1.0 / p)
        with np.errstate(divide="ignore", invalid="ignore"):
            dadt = p * s ** (p - 1.0) / R ** p * np.sign(np.sin(theta)) * np.cos(theta)
            dbdt = p * c ** (p - 1.0) / H ** p * np.sign(np.cos(theta)) * (-np.sin(theta))
        dadt = np.where(s > 0, dadt, 0.0)          # на полюсах (s или c =0) производная 0
        dbdt = np.where(c > 0, dbdt, 0.0)
        drdt = (-1.0 / p) * u ** (-1.0 / p - 1.0) * (dadt + dbdt)
        return r, drdt

    return r_of_theta


def cone_segments(radius: float, height: float, *, z_apex: float | None = None):
    """Образующая КОНЕЧНОГО конуса (ось z, вершина вверху): базовый радиус R, высота L.
    Два сегмента — боковая (наклонная) поверхность и донная крышка. По умолчанию вершина
    в z_apex=3L/4, дно в z=−L/4 (центроид в нуле ⇒ начало координат внутри тела).
    Полуугол α=arctan(R/L); боковая r(θ)=z_a·sinα/sin(θ+α); дно z=z_b: r=z_b/cosθ."""
    R, L = float(radius), float(height)
    z_a = (3.0 * L / 4.0) if z_apex is None else float(z_apex)
    z_b = z_a - L                                          # дно (<0 при дефолте)
    alpha = math.atan2(R, L)
    th_rim = math.atan2(R, z_b)                            # угол на ребро основания
    if th_rim < 0:
        th_rim += math.pi

    def lateral(theta):                  # r(θ)=z_a sinα/sin(θ+α)
        s = np.sin(theta + alpha)
        r = z_a * math.sin(alpha) / s
        drdt = -z_a * math.sin(alpha) * np.cos(theta + alpha) / s ** 2
        return r, drdt

    def base(theta):                     # z=z_b: r=z_b/cosθ
        c = np.cos(theta)
        return z_b / c, z_b * np.sin(theta) / c ** 2

    return [(lateral, 0.0, th_rim), (base, th_rim, math.pi)]


# --------------------------------------------------------------------------- #
# J-интегралы и сборка Q, RgQ (поблочно по m)
# --------------------------------------------------------------------------- #
def _J_block(pts, nvec_dS, k: float, k1: complex, m: int, n_list, test_kind: str,
             src_kind: str = "reg"):
    """Четыре J-матрицы (J11,J12,J21,J22) для фиксированного m.

    Источник — ВСВФ внутри (k1), мода (n', m), радиальная функция src_kind ('reg' —
    регулярная для однородного тела/ядра; 'out' — исходящая, для оболочек слоистого
    тела); тест — внешние ВСВФ (k), мода (n, −m), вида test_kind ('out' → Q, 'reg' → RgQ):
        J^{XY}_{nn'} = Σ_q nvec_dS[q] · ( X_{n'm}(k1) × Y_{n,−m}(k) )[q].
    Возвращает (J11,J12,J21,J22), каждая N×N (N=len(n_list)); X,Y∈{M,N}: первый индекс —
    тип теста (строка n), второй — тип источника (столбец n').
    """
    N = len(n_list)
    # источники (k1, радиальная функция src_kind): Msrc[j], Nsrc[j] для n'=n_list[j], мода +m
    Msrc, Nsrc = [], []
    for npr in n_list:
        Mg, Ng = vswf.mn_grid(npr, m, k1, pts, src_kind)
        Msrc.append(Mg)
        Nsrc.append(Ng)
    # тест (внешние, k): Mt[i], Nt[i] для n=n_list[i], мода −m
    Mt, Nt = [], []
    for nn in n_list:
        Mg, Ng = vswf.mn_grid(nn, -m, k, pts, test_kind)
        Mt.append(Mg)
        Nt.append(Ng)

    def integ(A, B):  # Σ_q nvec_dS·(A×B)
        return np.sum(nvec_dS * np.cross(A, B))

    J11 = np.empty((N, N), complex)   # test M (row n)  × src M (col n')
    J12 = np.empty((N, N), complex)   # test N × src M
    J21 = np.empty((N, N), complex)   # test M × src N
    J22 = np.empty((N, N), complex)   # test N × src N
    for i in range(N):
        for j in range(N):
            J11[i, j] = integ(Msrc[j], Mt[i])
            J12[i, j] = integ(Msrc[j], Nt[i])
            J21[i, j] = integ(Nsrc[j], Mt[i])
            J22[i, j] = integ(Nsrc[j], Nt[i])
    return J11, J12, J21, J22


def _Q_from_J(J, k: float, k1: complex, eps_rel: complex, mu_rel: complex,
              eps_o: complex = 1.0, mu_o: complex = 1.0):
    """Собрать блочную Q (2N×2N) из (J11,J12,J21,J22), порядок блоков [M(TE); N(TM)].

    Импедансный множитель s_E умножает «магнитный» (H-непрерывность) член — тот, где тип
    ИСТОЧНИКА противоположен типу блока (для M-блока это N-источник, J21). В общем случае
    интерфейса между внешней средой (eps_o,mu_o; волн. число k) и внутренней (eps_rel,mu_rel;
    k1) s_E — относительный импеданс √[(ε_i/μ_i)/(ε_o/μ_o)]; для внешней среды-вакуума
    (eps_o=mu_o=1) сводится к √(ε/μ). Выведено из предела сферы: T^{MM}=−RgQ^{MM}/Q^{MM}=−b_n
    (Bohren–Huffman 4.53), T^{NN}→−a_n. Общий множитель −i k k1 сокращается в T.

        Q^{MM} = s_E·J21 + J12   (→ t_M = −b_n)
        Q^{NN} = s_E·J12 + J21   (→ t_N = −a_n)
        Q^{MN} = s_E·J11 + J22   (кросс TE↔TM; на сфере J11=J22=0 ⇒ 0)
        Q^{NM} = s_E·J22 + J11
    """
    J11, J12, J21, J22 = J
    s_E = np.sqrt((eps_rel / mu_rel) / (eps_o / mu_o))
    pref = -1j * k * k1
    QMM = pref * (s_E * J21 + J12)
    QNN = pref * (s_E * J12 + J21)
    QMN = pref * (s_E * J11 + J22)
    QNM = pref * (s_E * J22 + J11)
    return np.block([[QMM, QMN], [QNM, QNN]])


def tmatrix_axisym(r_of_theta, k: float, eps: complex, mu: complex, nmax: int,
                   *, ntheta: int | None = None, nphi: int | None = None) -> np.ndarray:
    """ПОЛНАЯ T-матрица (2K×2K, базис [M;N], moded vswf.mode_list) однородного
    осесимметричного тела с образующей r(θ), методом нулевого поля.

    k — волновое число снаружи (вакуум); eps,mu — относительные проницаемости тела;
    nmax — усечение мультиполей. Блоки по m собираются и инвертируются независимо.

    ВАЛИДИРОВАНО (tests/test_ebcm.py): сфера → Ми ~1e-15 (диэлектрик и магнит); гладкие
    несферические тела (сфероид prolate/oblate, магнитный) — энергобаланс на уровне
    квадратуры (~1e-5) и совпадение с рэлеевской поляризуемостью; совместимо с GMM.

    НОРМИРОВКА (ключевой момент): ВСВФ проекта ортонормированы так, что ∮|M|²dΩ=n(n+1)
    (а не 1, как у Мищенко). Поэтому собранная из них T-матрица отличается от физической
    НЕунитарным подобием diag(1/(n(n+1))); поправка T[i,j]·=n_j(n_j+1)/(n_i(n_i+1)) внизу.
    На сфере T диагональна ⇒ подобие тождественно (сфера остаётся точной); вне диагонали
    (несферические тела) подобие восстанавливает энергобаланс. Неунитарность подобия и
    есть причина, по которой без поправки энергобаланс нарушался ~10%.
    """
    if ntheta is None:
        ntheta = 2 * nmax + 6
    if nphi is None:
        nphi = 2 * nmax + 2
    pts, nvec_dS = axisym_surface(r_of_theta, ntheta, nphi)
    return _tmatrix_from_surface(pts, nvec_dS, k, eps, mu, nmax)


def tmatrix_axisym_segments(segments, k: float, eps: complex, mu: complex, nmax: int,
                            *, n_per: int | None = None, nphi: int | None = None) -> np.ndarray:
    """ПОЛНАЯ T-матрица тела вращения с РЁБРАМИ (конечный цилиндр, конус) — образующая
    задана кусочно (см. :func:`cylinder_segments`, :func:`cone_segments`). Каждый гладкий
    сегмент интегрируется своей квадратурой (не сшивать через ребро). Остальное — как в
    :func:`tmatrix_axisym` (та же сборка Q, нормировка, конвенция).

    Острые рёбра/вершина (условие Мейкснера) замедляют сходимость по nmax и портят
    обусловленность Q; для тонких конусов нужны дискретные источники (NFM-DS) — здесь
    не реализованы. Контроль качества — энергобаланс (lossless ⇒ Q_abs≈0) и сходимость.
    """
    if n_per is None:
        n_per = 2 * nmax + 6
    if nphi is None:
        nphi = 2 * nmax + 2
    pts, nvec_dS = surface_segments(segments, n_per, nphi)
    return _tmatrix_from_surface(pts, nvec_dS, k, eps, mu, nmax)


def _tmatrix_from_surface(pts, nvec_dS, k: float, eps: complex, mu: complex,
                          nmax: int) -> np.ndarray:
    """Сборка полной T-матрицы по готовой поверхностной квадратуре (pts, n̂·dS).

    Общее ядро для гладких (:func:`tmatrix_axisym`) и кусочных
    (:func:`tmatrix_axisym_segments`) тел. Поблочно по m: T_m=−RgQ_m·Q_m⁻¹; затем
    нормировочная поправка ВСВФ (∮|M|²=n(n+1)) подобием diag(1/(n(n+1))) и сопряжение
    под конвенцию проекта (см. :func:`tmatrix_axisym`).
    """
    eps = complex(eps)
    mu = complex(mu)
    k1 = k * np.sqrt(eps * mu)
    modes = vswf.mode_list(nmax)
    K = len(modes)
    idx = {nm: i for i, nm in enumerate(modes)}
    T = np.zeros((2 * K, 2 * K), complex)

    for m in range(-nmax, nmax + 1):
        n_list = list(range(max(1, abs(m)), nmax + 1))
        if not n_list:
            continue
        Jout = _J_block(pts, nvec_dS, k, k1, m, n_list, "out")
        Jreg = _J_block(pts, nvec_dS, k, k1, m, n_list, "reg")
        Q = _Q_from_J(Jout, k, k1, eps, mu)
        RgQ = _Q_from_J(Jreg, k, k1, eps, mu)
        Tm = -RgQ @ np.linalg.inv(Q)                     # 2N×2N: [[MM,MN],[NM,NN]]

        N = len(n_list)
        rows = [idx[(n, m)] for n in n_list]
        for bi, ri in enumerate(rows):
            for bj, rj in enumerate(rows):
                T[ri, rj] += Tm[bi, bj]                   # MM
                T[ri, rj + K] += Tm[bi, bj + N]           # MN
                T[ri + K, rj] += Tm[bi + N, bj]           # NM
                T[ri + K, rj + K] += Tm[bi + N, bj + N]   # NN
    # Нормировочная поправка: ВСВФ проекта имеют ∮|M|²=n(n+1), собранная T связана с
    # физической неунитарным подобием diag(1/(n(n+1))) — применяем (на диагонали/сфере
    # тождественно, вне диагонали восстанавливает энергобаланс).
    nn = np.array([n * (n + 1) for (n, m) in modes], dtype=float)
    dvec = np.concatenate([1.0 / nn, 1.0 / nn])
    T = dvec[:, None] * T * (1.0 / dvec[None, :])
    # Конвенция проекта (mie_core хранит сопряжённые коэффициенты): сопрягаем для
    # совместимости с gmm и сферой scatterer.LayeredSphere.
    return np.conj(T)


# --------------------------------------------------------------------------- #
# Слоистые тела: рекурсия Петерсона–Стрёма (T-матрица оболочка-за-оболочкой)
# --------------------------------------------------------------------------- #
def _interface_P(pts, nvec_dS, k_o, k_i, eps_o, mu_o, eps_i, mu_i, nmax):
    """Четыре полных (2K×2K) интерфейсных матрицы P[src,test], src/test ∈ {reg,out}.

    src — радиальная функция внутреннего поля (k_i): reg (j_n) или out (h_n^{(1)});
    test — внешняя пробная (k_o): reg или out. Возвращает (P_rr, P_or, P_ro, P_oo)
    (первый индекс — src, второй — test). Внешняя среда (eps_o,mu_o) может быть не вакуум.
    """
    modes = vswf.mode_list(nmax)
    K = len(modes)
    idx = {nm: i for i, nm in enumerate(modes)}
    combos = [("reg", "reg"), ("out", "reg"), ("reg", "out"), ("out", "out")]
    P = {c: np.zeros((2 * K, 2 * K), complex) for c in combos}
    for m in range(-nmax, nmax + 1):
        n_list = list(range(max(1, abs(m)), nmax + 1))
        if not n_list:
            continue
        N = len(n_list)
        rows = [idx[(n, m)] for n in n_list]
        for sk, tk in combos:
            Jb = _J_block(pts, nvec_dS, k_o, k_i, m, n_list, tk, src_kind=sk)
            Qb = _Q_from_J(Jb, k_o, k_i, eps_i, mu_i, eps_o, mu_o)
            Pm = P[(sk, tk)]
            for bi, ri in enumerate(rows):
                for bj, rj in enumerate(rows):
                    Pm[ri, rj] += Qb[bi, bj]
                    Pm[ri, rj + K] += Qb[bi, bj + N]
                    Pm[ri + K, rj] += Qb[bi + N, bj]
                    Pm[ri + K, rj + K] += Qb[bi + N, bj + N]
    return P[("reg", "reg")], P[("out", "reg")], P[("reg", "out")], P[("out", "out")]


def tmatrix_layered(surface_builder, eps, mu, a_norm, k0: float, nmax: int):
    """ПОЛНАЯ T-матрица радиально-слоистого осесимметричного тела (рекурсия Петерсона–Стрёма).

    surface_builder(a) -> (pts, n̂·dS) строит квадратуру интерфейса, гомотетично
    масштабированного на фактор a (a_norm[p]); eps,mu — проницаемости слоёв ИЗНУТРИ
    НАРУЖУ (eps[0] — ядро); a_norm — нормированные радиусы границ (внешний = 1, по
    возрастанию); k0 — волновое число снаружи (вакуум).

    Рекурсия от ядра наружу: на каждом интерфейсе T_внеш = −(P_rr+P_or·T_внутр)·
    (P_ro+P_oo·T_внутр)⁻¹, где P** — интерфейсные матрицы между смежными средами, а
    T_внутр — T-матрица всего, что внутри (нагрузка). Это матричный аналог импедансной
    рекурсии слоистой сферы (ТФГ); для сферы воспроизводит покрытые коэф. Адена–Керкера.

    ВАЛИДАЦИЯ/ПРЕДЕЛ: покрытая сфера → Ми ~1e-15; редукция (одинаковые слои → однородное
    тело) ~1e-6. ЧИСЛЕННО: рекурсия использует ИСХОДЯЩИЕ (h_n^{(1)}) функции на внутренних
    границах, что для многих слоёв / малого радиуса ядра ухудшает обусловленность
    (энергобаланс ~1e-2 у 2-слойного сфероида, хуже для 3+ слоёв) — известная особенность
    T-матричной рекурсии (устойчивый вариант — импедансная форма / расш. точность).
    """
    eps = [complex(e) for e in eps]
    mu = [complex(m) for m in mu]
    nL = len(eps)
    if not (len(mu) == len(a_norm) == nL):
        raise ValueError("eps, mu, a_norm должны быть одной длины")
    k = [k0 * np.sqrt(eps[i] * mu[i]) for i in range(nL)]
    modes = vswf.mode_list(nmax)
    K = len(modes)
    T_in = np.zeros((2 * K, 2 * K), complex)
    for p in range(nL):                                   # изнутри наружу
        eps_i, mu_i, k_i = eps[p], mu[p], k[p]
        if p < nL - 1:
            eps_o, mu_o, k_o = eps[p + 1], mu[p + 1], k[p + 1]
        else:
            eps_o, mu_o, k_o = 1.0 + 0j, 1.0 + 0j, k0     # внешняя среда — вакуум
        pts, nvec_dS = surface_builder(a_norm[p])
        P_rr, P_or, P_ro, P_oo = _interface_P(pts, nvec_dS, k_o, k_i,
                                              eps_o, mu_o, eps_i, mu_i, nmax)
        T_in = -(P_rr + P_or @ T_in) @ np.linalg.inv(P_ro + P_oo @ T_in)
    nn = np.array([n * (n + 1) for (n, m) in modes], dtype=float)
    dvec = np.concatenate([1.0 / nn, 1.0 / nn])
    T = dvec[:, None] * T_in * (1.0 / dvec[None, :])
    return np.conj(T)


def tmatrix_axisym_layered(r_of_theta, eps, mu, a_norm, k0: float, nmax: int,
                           *, ntheta: int | None = None, nphi: int | None = None) -> np.ndarray:
    """Слоистая T-матрица гладкого тела вращения (сфера через const, сфероид) — обёртка
    над :func:`tmatrix_layered` с гомотетичным масштабированием образующей r(θ)."""
    if ntheta is None:
        ntheta = 2 * nmax + 6
    if nphi is None:
        nphi = 2 * nmax + 2

    def builder(a):
        def scaled(theta):
            r, dr = r_of_theta(theta)
            return a * np.asarray(r, float), a * np.asarray(dr, float)
        return axisym_surface(scaled, ntheta, nphi)
    return tmatrix_layered(builder, eps, mu, a_norm, k0, nmax)


def diagonal_from_tmatrix(T: np.ndarray, nmax: int):
    """Извлечь диагонали (t_M_n, t_N_n) из полной T (для сверки со сферой)."""
    modes = vswf.mode_list(nmax)
    K = len(modes)
    tM = {}
    tN = {}
    for i, (n, m) in enumerate(modes):
        tM.setdefault(n, []).append(T[i, i])
        tN.setdefault(n, []).append(T[i + K, i + K])
    ns = sorted(tM)
    tM_n = np.array([np.mean(tM[n]) for n in ns])
    tN_n = np.array([np.mean(tN[n]) for n in ns])
    return np.array(ns), tM_n, tN_n
