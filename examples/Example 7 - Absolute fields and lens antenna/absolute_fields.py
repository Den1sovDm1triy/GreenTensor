# SPDX-License-Identifier: MIT
"""Ненормированные характеристики поля методом ТФГ — расчётная часть примера 7.

Тела: проводящая сфера; проводящая сфера с диэлектрическим укрытием ПТФЭ; линза
Люнеберга (пять слоёв). Электрические размеры k0a = 0.5 / 2 / 15 при f0 = 10 ГГц;
линза Люнеберга рассчитывается при k0a = 15.

Для каждой конфигурации снимаются: линейная поляризация (вертикальная и
горизонтальная), круговая (основная и кросс-поляризация), карта ближней зоны в
плоскости xz и карта дальней зоны на сфере r = 10 м.

Все величины в единицах СИ (В/м); нормировка на максимум не применяется.
Амплитуда источника задаётся явно: падающая плоская волна E_0 = 1 В/м.

Результаты пишутся в data/ (диаграммы far_*.csv и сводка resources.json).
Диаграмма направленности линзовой антенны считается отдельно — figures_antenna.py.

Запуск: python3 absolute_fields.py
"""
import os, sys, time, math, tracemalloc, json
LIB = os.environ.get("GREENTENSOR_PATH", os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
sys.path.insert(0, LIB)

import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

import green_tensor as gt
from green_tensor.sphere_core import MieSphere
from green_tensor.nearfield import external_field_sphere as extfield

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
IMG = os.path.join(HERE, "figures", "aux")
os.makedirs(DATA, exist_ok=True); os.makedirs(IMG, exist_ok=True)

# ---------------- постановка ----------------
F0   = 10e9
CLIGHT = 299792458.0
K0   = 2*math.pi*F0/CLIGHT              # физическое волновое число, рад/м
LAM  = CLIGHT/F0
PEC  = 1e7                              # модель идеального проводника

# источник для АНТЕННОЙ задачи: элементарный вибратор (3.18)
I_SRC, L_SRC, R_SRC = 1.0, LAM/50.0, 1.0      # 1 А, λ/50, отсчёт E_0 на 1 м
E0_ANT = abs(MieSphere.dipole_e0(I_SRC, L_SRC, K0, R_SRC))
E0_DIF = 1.0                                   # плоская волна 1 В/м

BODIES = {
    "pec":      dict(title="Проводящая сфера",       eps=[PEC],              a_norm=[1.0]),
    "radome":   dict(title="Сфера с укрытием ПТФЭ",  eps=[PEC, 2.1-0.0004j], a_norm=[0.9, 1.0]),
    "luneburg": dict(title="Линза Люнеберга",
                     eps=[1.96, 1.84, 1.64, 1.36, 1.00],
                     a_norm=[0.2, 0.4, 0.6, 0.8, 1.0]),
}
SIZES = {"0p5": 0.5, "2": 2.0, "15": 15.0}
# конфигурации документа: однородные тела при трёх размерах и линза Люнеберга при k0a=15
CASES = [(b, s) for b in ("pec", "radome") for s in SIZES] + [("luneburg", "15")]
PROBLEM = "diffraction"

THETA = np.arange(0.01, 360.01, 1.0)*math.pi/180.0     # сетка обходит 90/270
R_FAR = 10.0                                            # м, kr≈2096 — дальняя зона всем

def e0_for(problem):
    return E0_DIF if problem == "diffraction" else E0_ANT


def abs_prefactor(problem, r_obs):
    """Множитель перевода безразмерной |S| в АБСОЛЮТНОЕ поле [В/м] на расстоянии r_obs."""
    if problem == "diffraction":
        return E0_DIF/(K0*r_obs)                      # E_0·|S|/(k r)
    return abs(MieSphere.dipole_e0(I_SRC, L_SRC, K0, r_obs))   # E_0(r) по (3.18)

# ---------------- расчёт одной конфигурации ----------------
def run_case(bkey, skey, problem):
    b, k0a = BODIES[bkey], SIZES[skey]
    R = k0a/K0
    eps, a_norm = b["eps"], b["a_norm"]
    e0 = e0_for(problem)
    tag = f"{bkey}_k{skey}_{problem}"

    tracemalloc.start(); t0 = time.perf_counter()

    # --- дальние диаграммы, АБСОЛЮТ (В/м) на расстоянии R_FAR ---
    # pattern() даёт БЕЗРАЗМЕРНУЮ амплитуду |S|; абсолютный префактор зависит от режима:
    #   дифракция: E = E_0·|S|/(k r)   (падающая плоская волна амплитуды E_0, без 1/r)
    #   антенна  : E = E_0(r)·|S|      (у E_0 по ф-ле 3.18 множитель e^{-ik0r}/r уже внутри)
    sp = gt.SphereSolver(R, eps, a_norm=a_norm)
    p_v = sp.pattern(K0, theta=THETA, problem=problem, phi=0.0)
    p_h = sp.pattern(K0, theta=THETA, problem=problem, phi=math.pi/2)
    p_c = sp.pattern(K0, theta=THETA, problem=problem, polarization="circular")
    pref = abs_prefactor(problem, R_FAR)
    # В базисе основной/кросс-поляризации (Панченко 3.15-3.17) слот "E_theta" всегда
    # несёт ОСНОВНУЮ (co-pol) составляющую среза, а "E_phi" — КРОСС, тождественно
    # нулевую в главных плоскостях. Поэтому горизонтальная = co-pol при phi=90 град.
    lin_v = {"E_theta": pref*p_v["E_theta"]}                 # E-плоскость (phi=0)
    lin_h = {"E_theta": pref*p_h["E_theta"]}                 # H-плоскость (phi=90)
    xpol  = {"E_phi":   pref*p_v["E_phi"]}                   # кросс-поляризация (контроль ~0)
    # круговая: множитель 1/2 ВЫВЕДЕН из ф-л (4.15)-(4.16) Панченко —
    # |E_R| = E_0/(2 k0 r)·|S∥−S⊥|, |E_L| = E_0/(2 k0 r)·|S∥+S⊥|;
    # величины pattern() не нормированы на единичную круговую падающую волну.
    # Рукость (сверено с HFSS): E_op ↔ RHCP, E_kp ↔ LHCP.
    circ  = {"E_op": pref/2*p_c["E_op"], "E_kp": pref/2*p_c["E_kp"]}
    for d in (lin_v, lin_h, xpol, circ):
        for kk in list(d):
            d[kk+"_max"] = float(np.nanmax(d[kk]))

    # --- карты: ближняя (плоскость xz) + дальняя (сфера r=R_FAR), один вызов ---
    ext = R + 2.0*LAM
    gx = np.linspace(-ext, ext, 141); gz = np.linspace(-ext, ext, 141)
    GX, GZ = np.meshgrid(gx, gz)
    near_pts = np.stack([GX.ravel(), np.zeros(GX.size), GZ.ravel()], axis=1)

    th_m = np.linspace(0.5, 179.5, 60)*math.pi/180.0
    ph_m = np.linspace(0.0, 360.0, 73)*math.pi/180.0
    TH, PH = np.meshgrid(th_m, ph_m)
    far_pts = np.stack([R_FAR*np.sin(TH.ravel())*np.cos(PH.ravel()),
                        R_FAR*np.sin(TH.ravel())*np.sin(PH.ravel()),
                        R_FAR*np.cos(TH.ravel())], axis=1)

    # Поле считается точным радиальным рядом (функции h_n) — продолжением того же
    # ряда, чья асимптотика даёт диаграмму; нормировка наследуется из дальней зоны.
    #   дифракция: амплитуда = E_0 падающей плоской волны;
    #   антенна  : A = I·k²·l·Z_0/(4π) — замкнутая форма согласования с (3.18).
    # Ближняя зона: полное поле при дифракции (наблюдаемая величина), излучённое — в антенне.
    # Дальняя зона: только рассеянное/излучённое (полное тонет в падающей волне).
    KHAT = (0.0, 0.0, -1.0)          # волна вдоль -z: θ=0 обратное, θ=180° вперёд
    src = {} if problem == "diffraction" else dict(current=I_SRC, length=L_SRC)
    amp = dict(e0=e0) if problem == "diffraction" else {}
    Enear = extfield(R, eps, K0, near_pts, problem=problem, a_norm=a_norm, khat=KHAT,
                total=(problem == "diffraction"), **amp, **src)
    Efar  = extfield(R, eps, K0, far_pts,  problem=problem, a_norm=a_norm, khat=KHAT,
                total=False, **amp, **src)

    # сверка: сечение φ=0 карты дальней зоны против аналитической диаграммы (абсолют)
    thc = np.deg2rad(np.arange(10.5, 171.0, 5.0))     # сетка обходит сингулярную точку 90°
    pc_ = np.stack([R_FAR*np.sin(thc), np.zeros_like(thc), R_FAR*np.cos(thc)], axis=1)
    Ec = extfield(R, eps, K0, pc_, problem=problem, a_norm=a_norm, khat=KHAT,
             total=False, **amp, **src)
    at = np.stack([np.cos(thc), np.zeros_like(thc), -np.sin(thc)], axis=1)
    lvl_map = np.abs(np.sum(Ec*at, axis=1))
    lvl_ref = pref*np.interp(np.degrees(thc), np.degrees(THETA), p_v["E_theta"])
    m_ok = lvl_ref > 0.01*np.nanmax(lvl_ref)          # значимые лепестки
    xcheck = float(np.nanmedian(np.abs(lvl_map[m_ok]-lvl_ref[m_ok])/lvl_ref[m_ok]))
    near_map = np.linalg.norm(np.abs(Enear), axis=1).reshape(GX.shape)     # |E| В/м
    far_map  = np.linalg.norm(np.abs(Efar),  axis=1).reshape(TH.shape)     # |E| В/м

    dt = time.perf_counter()-t0
    _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()

    # --- сохранение данных (абсолют) ---
    np.savetxt(os.path.join(DATA, f"far_{tag}.csv"),
        np.column_stack([np.degrees(THETA), lin_v["E_theta"], lin_h["E_theta"], xpol["E_phi"],
                         circ["E_op"], circ["E_kp"]]),
        delimiter=",", header="theta_deg,E_co_Eplane_Vm,E_co_Hplane_Vm,E_cross_Vm,E_op_Vm,E_kp_Vm", comments="")

    res = dict(tag=tag, body=bkey, k0a=k0a, radius_m=R, problem=problem,
               e0_Vm=e0, r_obs_m=R_FAR, xcheck_rel=round(xcheck,5), nmax=int(len(gt.SphereSolver(R,eps,a_norm=a_norm).mie(K0,problem=problem).coefficients()[0])),
               time_s=round(dt,3), peak_mem_MB=round(peak/1e6,3),
               E_theta_max_Vm=lin_v["E_theta_max"], E_co_Hplane_max_Vm=lin_h["E_theta_max"], E_cross_max_Vm=xpol["E_phi_max"],
               E_op_max_Vm=circ["E_op_max"], E_kp_max_Vm=circ["E_kp_max"],
               near_max_Vm=float(np.nanmax(near_map)), far_max_Vm=float(np.nanmax(far_map)))
    return res, (lin_v, lin_h, circ), (gx, gz, near_map), (th_m, ph_m, far_map), R

# ---------------- графика (абсолютные оси) ----------------
def plot_case(res, diag, near, far, R):
    tag, ttl = res["tag"], f"{BODIES[res['body']]['title']}, $k_0a$={res['k0a']:g}, {res['problem']}"
    lin_v, lin_h, circ = diag
    th = np.degrees(THETA)

    # 1) дальние диаграммы: линейная (верт/гориз) + круговая — АБСОЛЮТ В/м
    fig, ax = plt.subplots(1, 2, figsize=(12.4, 4.6))
    ax[0].plot(th, lin_v["E_theta"], lw=1.6, color="#c0392b", label=r"лин. верт. (E-плоск., $\varphi$=0)")
    ax[0].plot(th, lin_h["E_theta"], lw=1.6, color="#1f77b4", label=r"лин. гориз. (H-плоск., $\varphi$=90°)")
    ax[1].plot(th, circ["E_op"], lw=1.6, color="#2c7a2c", label=r"круг. основная $|E_{op}|$")
    ax[1].plot(th, circ["E_kp"], lw=1.6, color="#d18b00", label=r"круг. кросс $|E_{kp}|$")
    for a in ax:
        a.set_xlabel(r"$\theta$, град"); a.set_ylabel(r"$|E|$, В/м  (абсолют)")
        a.set_yscale("log"); a.grid(alpha=0.4, which="both"); a.legend(fontsize=9)
        a.set_xlim(0, 360)
    fig.suptitle(ttl, fontsize=11)
    fig.savefig(os.path.join(IMG, f"diag_{tag}.png"), dpi=105, bbox_inches="tight"); plt.close(fig)

    # 2) ближняя зона |E| В/м
    gx, gz, nm = near
    fig, a = plt.subplots(figsize=(6.0, 5.0))
    pc = a.pcolormesh(gx*1e3, gz*1e3, nm, shading="auto", cmap="inferno")
    a.add_patch(plt.Circle((0,0), R*1e3, fill=False, ec="cyan", lw=1.0, ls="--"))
    a.set_aspect("equal"); a.set_xlabel("x, мм"); a.set_ylabel("z, мм")
    a.set_title(f"Ближняя зона |E|, В/м\n{ttl}", fontsize=9)
    fig.colorbar(pc, ax=a, label="|E|, В/м (абсолют)")
    fig.savefig(os.path.join(IMG, f"near_{tag}.png"), dpi=105, bbox_inches="tight"); plt.close(fig)

    # 3) дальняя зона |E| В/м на сфере r=R_FAR
    th_m, ph_m, fm = far
    fig, a = plt.subplots(figsize=(6.8, 4.2))
    pc = a.pcolormesh(np.degrees(th_m), np.degrees(ph_m), fm, shading="auto", cmap="viridis")
    a.set_xlabel(r"$\theta$, град"); a.set_ylabel(r"$\varphi$, град")
    a.set_title(f"Дальняя зона |E|, В/м (r={R_FAR:g} м)\n{ttl}", fontsize=9)
    fig.colorbar(pc, ax=a, label="|E|, В/м (абсолют)")
    fig.savefig(os.path.join(IMG, f"farmap_{tag}.png"), dpi=105, bbox_inches="tight"); plt.close(fig)

# ---------------- прогон всей матрицы ----------------
if __name__ == "__main__":
    print(f"f={F0/1e9:g} ГГц, λ={LAM*1e3:.3f} мм, k0={K0:.3f} рад/м")
    print(f"E_0 = {E0_DIF:g} В/м (падающая плоская волна)\n")
    rows = []
    for bkey, skey in CASES:
        res, diag, near, far, R = run_case(bkey, skey, PROBLEM)
        plot_case(res, diag, near, far, R)
        rows.append(res)
        print(f"  {res['tag']:28s} N={res['nmax']:3d} "
              f"t={res['time_s']:7.2f} с  mem={res['peak_mem_MB']:7.2f} МБ  "
              f"|E|max(ближ)={res['near_max_Vm']:.4e} В/м")
    with open(os.path.join(DATA, "resources.json"), "w", encoding="utf-8") as fh:
        json.dump(rows, fh, ensure_ascii=False, indent=2)
    # таблица ресурсов
    import csv
    with open(os.path.join(DATA, "resources.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"\nГотово: {len(rows)} конфигураций -> data/, img/")
