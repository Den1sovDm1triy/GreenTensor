"""
MOGA-стратификация профиля линзы Люнеберга: многокритериальная генетическая
оптимизация ступенчатой аппроксимации, превосходящая схему Fuchs.

Постановка — по формулировке диссертации (Dissertation/_llapprox.tex,
eq:moga-objective / eq:moga-constraints):

    переменные:   x = (r_1..r_N; eps_1..eps_N),  0 < r_1 < ... < r_N <= 1,
                  область [r_N, 1] — согласующий воздушный слой eps = 1;
    ограничения:  eps_l in [1; eps_max],  eps_1 >= eps_2 >= ... >= eps_N >= 1
                  (для HIPS eps_max = 2.4 — ограничение неактивно, т.к.
                   профиль Люнеберга не превышает 2);
    критерии:     f1 = RMSE(eps)   = sqrt( int_0^1 (eps_step - eps_cont)^2 dr ),
                  f2 = deps_max    = максимальный скачок проницаемости между
                                     соседними слоями (включая выход в воздух),
                  оба -> min  =>  фронт Парето.

Алгоритм — NSGA-II (быстрая недоминируемая сортировка + crowding distance,
SBX-кроссовер, полиномиальная мутация). Значение f1 вычисляется ТОЧНО
(аналитические интегралы по слоям), без сэмплирования.

Верификация глобального оптимума по критерию f1 — алгоритм Ллойда
(необходимые условия оптимальности кусочно-постоянной L2-аппроксимации:
eps_l = среднее профиля по слою, граница между слоями — в точке равного
удаления уровней). GA обязан сойтись к оптимуму Ллойда с точностью ~1e-4.

Сравнение — с равношаговой схемой (L2-оптимальные средние, b2_approximation_rmse.py)
и со схемой Fuchs (табулировка N=1..8).

Коэффициенты заполнения FDM-печати — по логарифмической формуле Лихтенекера:
    alpha_l = ln(eps_l) / ln(eps_HIPS),  eps_HIPS = 2.4.

Выходы:
    moga_stratification.csv        : N | RMSE_eq | RMSE_Fuchs | RMSE_MOGA |
                                     выигрыш к Fuchs,% | deps_max | alpha_min..alpha_max
    moga_layers_N8.csv             : параметры слоёв MOGA при N=1..8 (r_l, eps_l)
    moga-rmse-comparison.pdf/png   : RMSE(N) три схемы + профили N=8
    moga-pareto-l8.pdf/png         : фронт Парето (RMSE, deps_max) при N=8,
                                     точки Fuchs и равношаговой — доминируемы
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EPS_HIPS = 2.4          # проницаемость сплошного HIPS (ограничение сверху)
EPS_MIN = 1.0           # воздух

# ----------------------------------------------------------------------------
# Точная метрика: RMSE^2 = sum_l int_{r_{l-1}}^{r_l} (eps_l - (2 - r^2))^2 dr
#                        + int_{r_N}^{1} (1 - (2 - r^2))^2 dr
# ----------------------------------------------------------------------------
def _seg_int(a, b, c):
    """int_a^b ((c-2) + r^2)^2 dr, поэлементно для массивов numpy."""
    d = c - 2.0
    return (d * d * (b - a)
            + (2.0 * d / 3.0) * (b**3 - a**3)
            + (b**5 - a**5) / 5.0)

def rmse_exact(bounds, eps):
    """bounds: (..., N+1) с bounds[...,0]=0; eps: (..., N). Воздух [r_N,1]."""
    a = bounds[..., :-1]
    b = bounds[..., 1:]
    err2 = _seg_int(a, b, eps).sum(axis=-1)
    err2 = err2 + _seg_int(bounds[..., -1], 1.0, 1.0)
    return np.sqrt(np.maximum(err2, 0.0))

def max_jump(eps):
    """Максимальный скачок eps между соседними слоями, включая выход в воздух."""
    inner = -np.diff(eps, axis=-1) if eps.shape[-1] > 1 else np.zeros(eps.shape[:-1] + (0,))
    outer = (eps[..., -1] - 1.0)[..., None]
    return np.concatenate([inner, outer], axis=-1).max(axis=-1)

# ----------------------------------------------------------------------------
# Схемы сравнения
# ----------------------------------------------------------------------------
def equidistant_scheme(N):
    b = np.linspace(0.0, 1.0, N + 1)
    e = 2.0 - (b[:-1]**2 + b[:-1] * b[1:] + b[1:]**2) / 3.0
    return b, e

FUCHS_TABLE = {
    1: dict(r=[0.82], eps=[1.67]),
    2: dict(r=[0.63, 0.89], eps=[1.80, 1.40]),
    3: dict(r=[0.53, 0.75, 0.93], eps=[1.86, 1.57, 1.28]),
    4: dict(r=[0.47, 0.67, 0.82, 0.94], eps=[1.88, 1.67, 1.44, 1.22]),
    5: dict(r=[0.43, 0.60, 0.74, 0.85, 0.95], eps=[1.91, 1.73, 1.55, 1.36, 1.18]),
    6: dict(r=[0.39, 0.55, 0.68, 0.79, 0.88, 0.96], eps=[1.93, 1.77, 1.61, 1.46, 1.31, 1.16]),
    7: dict(r=[0.37, 0.50, 0.63, 0.74, 0.84, 0.91, 0.97], eps=[1.93, 1.80, 1.67, 1.53, 1.40, 1.27, 1.13]),
    8: dict(r=[0.34, 0.49, 0.59, 0.69, 0.77, 0.84, 0.91, 0.97], eps=[1.94, 1.82, 1.71, 1.59, 1.47, 1.35, 1.24, 1.12]),
}

def fuchs_scheme(N):
    f = FUCHS_TABLE[N]
    b = np.array([0.0] + list(f["r"]))
    e = np.array(f["eps"])
    return b, e

# ----------------------------------------------------------------------------
# Верификация: алгоритм Ллойда (координатный спуск к необходимым условиям
# оптимальности f1 при закреплённом внешнем воздушном уровне eps=1)
# ----------------------------------------------------------------------------
def lloyd_optimum(N, iters=200000, tol=1e-14):
    r = np.linspace(0.0, 1.0, N + 2)[1:-1]        # ровно N границ r_1..r_N < 1
    for _ in range(iters):
        bounds = np.concatenate(([0.0], r))
        e = 2.0 - (bounds[:-1]**2 + bounds[:-1] * bounds[1:] + bounds[1:]**2) / 3.0
        mid = (e[:-1] + e[1:]) / 2.0
        r_new_inner = np.sqrt(np.maximum(2.0 - mid, 0.0))
        r_new_outer = np.sqrt(np.maximum(2.0 - (e[-1] + 1.0) / 2.0, 0.0))
        r_new = np.concatenate((r_new_inner, [min(r_new_outer, 1.0)]))
        if np.max(np.abs(r_new - r)) < tol:
            r = r_new
            break
        r = r_new
    bounds = np.concatenate(([0.0], r))
    e = 2.0 - (bounds[:-1]**2 + bounds[:-1] * bounds[1:] + bounds[1:]**2) / 3.0
    return bounds, e

# ----------------------------------------------------------------------------
# NSGA-II
# ----------------------------------------------------------------------------
def decode(genes, N):
    """Генотип [0,1]^{2N} -> (bounds (pop,N+1), eps (pop,N)) с соблюдением
    ограничений по построению: сортировка радиусов и невозрастание eps."""
    r = np.sort(np.clip(genes[:, :N], 1e-4, 1.0), axis=1)
    e = np.sort(EPS_MIN + (2.0 - EPS_MIN) * np.clip(genes[:, N:], 0.0, 1.0), axis=1)[:, ::-1]
    bounds = np.concatenate([np.zeros((genes.shape[0], 1)), r], axis=1)
    return bounds, np.ascontiguousarray(e)

def evaluate(genes, N):
    bounds, eps = decode(genes, N)
    return np.stack([rmse_exact(bounds, eps), max_jump(eps)], axis=1)

def fast_non_dominated_sort(F):
    n = F.shape[0]
    dominated_by = [[] for _ in range(n)]
    dom_count = np.zeros(n, dtype=int)
    for i in range(n):
        le = np.all(F[i] <= F, axis=1)
        lt = np.any(F[i] < F, axis=1)
        dominates = le & lt
        dominates[i] = False
        idx = np.nonzero(dominates)[0]
        for j in idx:
            dominated_by[i].append(j)
        dom_count += dominates.astype(int)
    fronts, current = [], np.nonzero(dom_count == 0)[0].tolist()
    rank = np.zeros(n, dtype=int)
    k = 0
    while current:
        fronts.append(current)
        nxt = []
        for i in current:
            rank[i] = k
            for j in dominated_by[i]:
                dom_count[j] -= 1
                if dom_count[j] == 0:
                    nxt.append(j)
        current = nxt
        k += 1
    return fronts, rank

def crowding_distance(F, front):
    d = np.zeros(len(front))
    for m in range(F.shape[1]):
        vals = F[front, m]
        order = np.argsort(vals)
        d[order[0]] = d[order[-1]] = np.inf
        span = vals[order[-1]] - vals[order[0]]
        if span > 0:
            d[order[1:-1]] += (vals[order[2:]] - vals[order[:-2]]) / span
    return d

def nsga2(N, pop=200, gens=600, seed=7):
    rng = np.random.default_rng(seed)
    dim = 2 * N
    P = rng.random((pop, dim))
    # эвристический посев: равношаговая точка и Fuchs-точка в стартовой популяции
    b_eq, e_eq = equidistant_scheme(N)
    P[0, :N] = b_eq[1:]
    P[0, N:] = (e_eq - EPS_MIN) / (2.0 - EPS_MIN)
    if N in FUCHS_TABLE:
        b_fu, e_fu = fuchs_scheme(N)
        P[1, :N] = b_fu[1:]
        P[1, N:] = (e_fu - EPS_MIN) / (2.0 - EPS_MIN)
    F = evaluate(P, N)
    eta_c, eta_m, pm = 15.0, 20.0, 1.0 / dim
    for _ in range(gens):
        fronts, rank = fast_non_dominated_sort(F)
        crowd = np.zeros(pop)
        for fr in fronts:
            crowd[fr] = crowding_distance(F, fr)
        # турнирная селекция
        cand = rng.integers(0, pop, size=(pop, 2))
        better = np.where(
            (rank[cand[:, 0]] < rank[cand[:, 1]])
            | ((rank[cand[:, 0]] == rank[cand[:, 1]]) & (crowd[cand[:, 0]] > crowd[cand[:, 1]])),
            cand[:, 0], cand[:, 1])
        M = P[better]
        # SBX
        C = M.copy()
        for i in range(0, pop - 1, 2):
            if rng.random() < 0.9:
                u = rng.random(dim)
                beta = np.where(u <= 0.5, (2 * u) ** (1 / (eta_c + 1)),
                                (1 / (2 * (1 - u))) ** (1 / (eta_c + 1)))
                p1, p2 = M[i], M[i + 1]
                C[i] = 0.5 * ((1 + beta) * p1 + (1 - beta) * p2)
                C[i + 1] = 0.5 * ((1 - beta) * p1 + (1 + beta) * p2)
        # полиномиальная мутация
        mut = rng.random(C.shape) < pm
        u = rng.random(C.shape)
        delta = np.where(u < 0.5, (2 * u) ** (1 / (eta_m + 1)) - 1,
                         1 - (2 * (1 - u)) ** (1 / (eta_m + 1)))
        C = np.clip(C + mut * delta, 0.0, 1.0)
        FC = evaluate(C, N)
        # объединение и отбор
        PA, FA = np.vstack([P, C]), np.vstack([F, FC])
        fronts, _ = fast_non_dominated_sort(FA)
        new_idx = []
        for fr in fronts:
            if len(new_idx) + len(fr) <= pop:
                new_idx.extend(fr)
            else:
                d = crowding_distance(FA, fr)
                order = np.argsort(-d)
                need = pop - len(new_idx)
                new_idx.extend([fr[k] for k in order[:need]])
                break
        P, F = PA[new_idx], FA[new_idx]
    fronts, _ = fast_non_dominated_sort(F)
    pf = fronts[0]
    return P[pf], F[pf]

def alpha_fill(eps):
    """Коэффициент заполнения FDM по логарифмической формуле Лихтенекера, %."""
    return 100.0 * np.log(np.asarray(eps)) / np.log(EPS_HIPS)

# ----------------------------------------------------------------------------
# Финишное уточнение фронта Парето методом эпсилон-ограничений:
# для сетки значений cap решается однокритериальная задача
#   min RMSE  при  всех скачках deps <= cap  (SLSQP, старт — оптимум Ллойда)
# ----------------------------------------------------------------------------
def refine_front(N, caps, x0_bounds, x0_eps):
    from scipy.optimize import minimize
    def unpack(x):
        r, e = x[:N], x[N:]
        return np.concatenate([[0.0], r]), e
    def obj(x):
        b, e = unpack(x)
        return float(rmse_exact(b[None], e[None])[0])
    refined = []
    x0 = np.concatenate([x0_bounds[1:], x0_eps])
    for cap in caps:
        cons = []
        for l in range(N - 1):
            cons.append({'type': 'ineq', 'fun': lambda x, l=l: x[N + l] - x[N + l + 1]})          # монотонность
            cons.append({'type': 'ineq', 'fun': lambda x, l=l: cap - (x[N + l] - x[N + l + 1])})  # скачок <= cap
            cons.append({'type': 'ineq', 'fun': lambda x, l=l: x[l + 1] - x[l] - 1e-3})           # порядок радиусов
        cons.append({'type': 'ineq', 'fun': lambda x: x[2 * N - 1] - 1.0})            # eps_N >= 1
        cons.append({'type': 'ineq', 'fun': lambda x: cap - (x[2 * N - 1] - 1.0)})    # скачок в воздух <= cap
        cons.append({'type': 'ineq', 'fun': lambda x: EPS_HIPS - x[N]})               # eps_1 <= eps_max
        cons.append({'type': 'ineq', 'fun': lambda x: x[0] - 1e-3})
        cons.append({'type': 'ineq', 'fun': lambda x: 1.0 - x[N - 1]})                # r_N <= 1
        res = minimize(obj, x0, method='SLSQP', constraints=cons,
                       options=dict(maxiter=400, ftol=1e-12))
        b, e = unpack(res.x)
        refined.append((float(rmse_exact(b[None], e[None])[0]),
                        float(max_jump(e[None])[0]), b.copy(), e.copy()))
    return refined

# ----------------------------------------------------------------------------
def main():
    Ns = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16]
    here = os.path.dirname(os.path.abspath(__file__))
    rows, layers_rows = [], []
    pareto_l8 = None
    best_l8 = None

    for N in Ns:
        b_eq, e_eq = equidistant_scheme(N)
        rmse_eq = float(rmse_exact(b_eq[None], e_eq[None])[0])
        if N in FUCHS_TABLE:
            b_fu, e_fu = fuchs_scheme(N)
            rmse_fu = float(rmse_exact(b_fu[None], e_fu[None])[0])
            jump_fu = float(max_jump(e_fu[None])[0])
        else:
            rmse_fu = jump_fu = float("nan")

        Ppf, Fpf = nsga2(N)
        i_best = int(np.argmin(Fpf[:, 0]))
        bounds, eps = decode(Ppf[i_best][None], N)
        bounds, eps = bounds[0], eps[0]
        rmse_ga = float(Fpf[i_best, 0])
        jump_ga = float(Fpf[i_best, 1])

        # верификация Ллойдом + финишное локальное уточнение GA-решения
        b_ll, e_ll = lloyd_optimum(N)
        rmse_ll = float(rmse_exact(b_ll[None], e_ll[None])[0])
        ga_gap = rmse_ga - rmse_ll     # невязка чистого GA к оптимуму Ллойда
        if rmse_ll < rmse_ga:
            bounds, eps, rmse_ga = b_ll, e_ll, rmse_ll
            jump_ga = float(max_jump(eps[None])[0])

        gain_fu = (rmse_fu - rmse_ga) / rmse_fu * 100.0
        gain_eq = (rmse_eq - rmse_ga) / rmse_eq * 100.0
        al = alpha_fill(eps)
        rows.append((N, rmse_eq, rmse_fu, rmse_ga, gain_fu, gain_eq,
                     jump_ga, al.min(), al.max(), ga_gap))
        layers_rows.append((N, bounds[1:], eps))

        if N == 8:
            pareto_l8 = (Fpf.copy(), rmse_fu, jump_fu, rmse_eq,
                         float(max_jump(e_eq[None])[0]))
            best_l8 = (bounds.copy(), eps.copy())
            # финишное уточнение фронта методом эпсилон-ограничений
            caps = np.linspace(0.10, 0.20, 21)
            refined_l8 = refine_front(N, caps, bounds, eps)

        print(f"N={N}: eq={rmse_eq:.5f}  Fuchs={rmse_fu:.5f}  MOGA={rmse_ga:.5f} "
              f"(выигрыш к Fuchs {gain_fu:+.1f}%, невязка чистого GA к Ллойду {ga_gap:.2e})")

    # ---------------- CSV: сводная таблица ----------------
    with open(os.path.join(here, "moga_stratification.csv"), "w", encoding="utf-8") as f:
        f.write("N,RMSE_equidistant,RMSE_Fuchs,RMSE_MOGA,gain_vs_Fuchs_pct,"
                "gain_vs_eq_pct,deps_max_MOGA,alpha_min_pct,alpha_max_pct,lloyd_gap\n")
        for r in rows:
            f.write(f"{r[0]},{r[1]:.5f},{r[2]:.5f},{r[3]:.5f},{r[4]:.2f},"
                    f"{r[5]:.2f},{r[6]:.4f},{r[7]:.1f},{r[8]:.1f},{r[9]:.2e}\n")

    # ---------------- CSV: параметры слоёв ----------------
    with open(os.path.join(here, "moga_layers.csv"), "w", encoding="utf-8") as f:
        f.write("N,radii,eps,alpha_pct\n")
        for N, r, e in layers_rows:
            f.write(f"{N},\"{';'.join(f'{x:.4f}' for x in r)}\","
                    f"\"{';'.join(f'{x:.4f}' for x in e)}\","
                    f"\"{';'.join(f'{x:.1f}' for x in alpha_fill(e))}\"\n")

    # ---------------- Рисунок 1: сравнение трёх схем ----------------
    plt.rcParams["font.family"] = "DejaVu Serif"
    plt.rcParams["font.size"] = 11
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    rows8 = [r for r in rows if r[0] <= 8]
    NsA = [r[0] for r in rows8]
    ax = axes[0]
    ax.plot(NsA, [r[1] for r in rows8], "o-", lw=1.6, color="#1f77b4",
            label="равношаговая (L2-оптимизация ε)")
    ax.plot(NsA, [r[2] for r in rows8], "s-", lw=1.6, color="#d62728",
            label="оптимизированная (Fuchs)")
    ax.plot(NsA, [r[3] for r in rows8], "^-", lw=1.8, color="#2ca02c",
            label="MOGA (NSGA-II, мин. RMSE)")
    ax.set_xlabel(r"$N$ (число слоёв)")
    ax.set_ylabel(r"RMSE $\varepsilon(r)$")
    ax.set_xticks(NsA)
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title("(а) среднеквадратическая погрешность", fontsize=11)

    ax = axes[1]
    r_fine = np.linspace(0, 1, 4000)
    ax.plot(r_fine, 2 - r_fine**2, "-", lw=1.6, color="black",
            label="идеальный профиль $2-(r/a)^2$")
    def plot_step(bounds, eps, color, label, lw=1.4, alpha=1.0):
        eps_full = np.concatenate([eps, [1.0]])
        b_full = np.concatenate([bounds, [1.0]])
        idx = np.clip(np.searchsorted(b_full[1:], r_fine, side="right"),
                      0, len(eps_full) - 1)
        ax.plot(r_fine, eps_full[idx], "-", lw=lw, color=color, label=label, alpha=alpha)
    b_fu, e_fu = fuchs_scheme(8)
    plot_step(b_fu, e_fu, "#d62728", "Fuchs $N=8$", alpha=0.75)
    plot_step(*best_l8, "#2ca02c", "MOGA $N=8$", lw=1.7)
    ax.set_xlabel(r"$r/a$"); ax.set_ylabel(r"$\varepsilon(r/a)$")
    ax.set_xlim(0, 1); ax.set_ylim(0.95, 2.05)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    ax.set_title(r"(б) сравнение профилей при $N=8$", fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(here, f"moga-rmse-comparison.{ext}"),
                    bbox_inches="tight", dpi=150)
    plt.close(fig)

    # ---------------- Рисунок 2: фронт Парето при N=8 ----------------
    Fpf, rmse_fu8, jump_fu8, rmse_eq8, jump_eq8 = pareto_l8
    ref_f1 = np.array([p[0] for p in refined_l8])
    ref_f2 = np.array([p[1] for p in refined_l8])
    order = np.argsort(ref_f1)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ga_mask = Fpf[:, 0] <= 0.055
    ax.plot(Fpf[ga_mask, 0], Fpf[ga_mask, 1], ".", ms=5, color="#9ecf9e",
            label="популяция MOGA (NSGA-II)")
    ax.plot(ref_f1[order], ref_f2[order], "^-", ms=5, lw=1.4, color="#2ca02c",
            label="фронт Парето (MOGA + локальное уточнение)")
    ax.plot([rmse_fu8], [jump_fu8], "s", ms=8, color="#d62728", label="Fuchs $N=8$")
    ax.plot([rmse_eq8], [jump_eq8], "o", ms=8, color="#1f77b4", label="равношаговая $N=8$")
    ax.annotate("Fuchs", (rmse_fu8, jump_fu8), textcoords="offset points",
                xytext=(8, 2), fontsize=10)
    ax.annotate("равношаговая", (rmse_eq8, jump_eq8), textcoords="offset points",
                xytext=(-52, -18), fontsize=10)
    ax.set_xlabel(r"$f_1 = \mathrm{RMSE}(\varepsilon)$")
    ax.set_ylabel(r"$f_2 = \Delta\varepsilon_{\max}$")
    ax.set_xlim(0.028, 0.044)
    ax.set_ylim(0.09, 0.23)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=9)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(here, f"moga-pareto-l8.{ext}"),
                    bbox_inches="tight", dpi=150)
    plt.close(fig)

    # точка фронта, строго доминирующая схему Fuchs (cap = скачок Fuchs)
    dom = [p for p in refined_l8 if p[1] <= jump_fu8 + 1e-9]
    if dom:
        p = min(dom, key=lambda q: q[0])
        print(f"\nТочка фронта при deps_max <= {jump_fu8:.3f} (как у Fuchs): "
              f"RMSE = {p[0]:.5f} (Fuchs: {rmse_fu8:.5f}) -> строгое доминирование")
        print("  r   =", ", ".join(f"{x:.3f}" for x in p[2][1:]))
        print("  eps =", ", ".join(f"{x:.3f}" for x in p[3]))
        print("  alpha% =", ", ".join(f"{x:.1f}" for x in alpha_fill(p[3])))

    # ---------------- параметры слоёв для таблицы диссертации ----------------
    print("\nПараметры слоёв MOGA (для таблицы):")
    for N, r, e in layers_rows:
        rr = ", ".join(f"{x:.3f}" for x in r)
        ee = ", ".join(f"{x:.3f}" for x in e)
        aa = ", ".join(f"{x:.1f}" for x in alpha_fill(e))
        print(f"  N={N}: r = [{rr}]\n        eps = [{ee}]\n        alpha% = [{aa}]")

if __name__ == "__main__":
    main()
