// SPDX-License-Identifier: MIT
// GreenTensor Studio — модуль результатов (Plotly): ДН, тепловые карты, сечения, спектр.
(function () {
  "use strict";
  const GT = (window.GT = window.GT || {});

  const DASHES = ["solid", "dash", "dot", "dashdot", "longdash", "longdashdot"];
  let animTimer = null;

  // ---- Палитра под тему (сшито с лендингом: тёмный неон-canvas / светлый научный) ----
  function themeDark() { return (document.documentElement.getAttribute("data-theme") || "dark") !== "light"; }
  function PAL() {
    return themeDark()
      ? { paper: "#0a1411", font: "#d8efe5", grid: "rgba(0,255,156,0.10)", axis: "rgba(0,255,156,0.22)",
          polar: "#0c1815", legend: "rgba(7,9,10,0.55)", legendBorder: "#1d2d27", muted: "#8aa49b",
          accent: "#00ff9c", bwLine: "#d8efe5",
          cat: ["#00ff9c", "#ff6b6b", "#4dd2ff", "#ffb547", "#c58cff", "#5ad1a8", "#ff8ac5"] }
      : { paper: "#ffffff", font: "#1a2128", grid: "#e6eaee", axis: "#c2cbd3",
          polar: "#fbfdfe", legend: "rgba(255,255,255,0.7)", legendBorder: "#d8dee4", muted: "#7a8893",
          accent: "#0d6b8a", bwLine: "#1a2128",
          cat: ["#0d6b8a", "#b23a3a", "#2f7d5b", "#7a4fb0", "#c07a16", "#1f6fb2", "#9b1d6b"] };
  }

  function isBW() { return GT.state.colorMode === "bw"; }
  function lineStyle(i) {
    const P = PAL();
    return isBW()
      ? { color: P.bwLine, dash: DASHES[i % DASHES.length], width: 2 }
      : { color: P.cat[i % P.cat.length], dash: "solid", width: 2 };
  }
  function baseLayout(extra) {
    const P = PAL();
    return Object.assign({
      paper_bgcolor: P.paper, plot_bgcolor: P.paper,
      font: { family: "Inter, -apple-system, Segoe UI, Roboto, sans-serif", size: 13, color: P.font },
      margin: { l: 56, r: 24, t: 28, b: 48 }, showlegend: true,
      legend: { bgcolor: P.legend, bordercolor: P.legendBorder, borderwidth: 1 },
    }, extra || {});
  }
  function config(name) {
    return {
      responsive: true, displaylogo: false,
      modeBarButtonsToRemove: ["lasso2d", "select2d", "toggleSpikelines"],
      toImageButtonOptions: { format: "png", filename: name || "greentensor", scale: 2 },
    };
  }

  function stopAnim() { if (animTimer) { clearInterval(animTimer); animTimer = null; } }

  // ===================== ДН (полярная) ===================== //
  function polarAxis(rot, useDB) {
    const P = PAL();
    return {
      bgcolor: P.polar,
      angularaxis: { rotation: rot, direction: "counterclockwise", gridcolor: P.grid,
        ticksuffix: "°", linecolor: P.axis },
      radialaxis: { range: useDB ? [-60, 0] : undefined, gridcolor: P.grid,
        linecolor: P.axis, title: { text: useDB ? "дБ" : "|E|" }, angle: 90 },
    };
  }

  function impTrace(imp, idx, subplot) {
    const theta = normalizeAngle(imp.x, imp.y.length);
    const t = { type: "scatterpolar", mode: "lines", name: imp.name + " (имп.)",
      theta: theta, r: imp.y,
      line: Object.assign(lineStyle(idx), { width: 1.6, dash: isBW() ? "dot" : "dash" }) };
    if (subplot) { t.subplot = subplot; t.legendgroup = subplot; }
    return t;
  }

  // Эталонная кривая (Ansys HFSS, аналитический закон и т.п.), в дБ, по плоскости:
  // E — красная пунктирная, H — оранжевая точечная; ref.color/ref.dash — переопределение.
  function refTrace(ref, subplot) {
    const col = isBW() ? PAL().bwLine
      : (ref.color || (ref.plane === "H" ? "#ff7f0e" : "#d62728"));
    const t = { type: "scatterpolar", mode: "lines", name: ref.name,
      theta: ref.theta_deg, r: ref.dB,
      line: { color: col, width: 2, dash: ref.dash || (ref.plane === "H" ? "dot" : "dash") } };
    if (subplot) { t.subplot = subplot; t.legendgroup = subplot; }
    return t;
  }

  function renderPattern() {
    const res = GT.state.results;
    const host = GT.$("#plot");
    if (!res || !res.pattern) return showEmpty(true);
    showEmpty(false);
    const useDB = GT.$("#db-toggle").checked;
    const rot = parseFloat(GT.$("#polar-rot").value) || 0;
    const pat = res.pattern;
    const series = pat.series;
    const split = GT.state.patternLayout === "split" && series.length > 1;
    const traces = [];
    const layout = baseLayout({ margin: { l: 24, r: 24, t: 34, b: 20 }, annotations: [] });

    if (!split) {
      // одна полярная панель: все вычисленные серии + все импорты
      layout.polar = polarAxis(rot, useDB);
      series.forEach((s, i) => traces.push({ type: "scatterpolar", mode: "lines",
        name: s.name, theta: pat.theta_deg, r: useDB ? s.dB : s.E, line: lineStyle(i) }));
      GT.state.imports.forEach((imp, j) => traces.push(impTrace(imp, series.length + j)));
      // кривые Ansys HFSS (по плоскостям E/H) — на общую панель
      if (useDB && GT.state.reference) GT.state.reference.forEach((rf) => traces.push(refTrace(rf)));
    } else {
      // панель на каждую плоскость; импорты раскладываются по своей панели (panel)
      const n = series.length, cols = Math.min(n, 2), rows = Math.ceil(n / cols);
      series.forEach((s, i) => {
        const col = i % cols, row = Math.floor(i / cols);
        const pol = i === 0 ? "polar" : "polar" + (i + 1);
        const px = 0.08, py = 0.14;
        const x0 = col / cols + px / 2, x1 = (col + 1) / cols - px / 2;
        const y1 = 1 - row / rows - 0.03, y0 = 1 - (row + 1) / rows + py / 2;
        const ax = polarAxis(rot, useDB);
        ax.domain = { x: [x0, x1], y: [Math.max(0, y0), Math.min(1, y1)] };
        layout[pol] = ax;
        traces.push({ type: "scatterpolar", subplot: pol, mode: "lines", name: s.name,
          theta: pat.theta_deg, r: useDB ? s.dB : s.E, line: lineStyle(i), legendgroup: pol });
        GT.state.imports.forEach((imp, j) => {
          const panel = imp.panel == null ? "all" : imp.panel;
          if (panel !== "all" && Number(panel) !== i) return;
          traces.push(impTrace(imp, n + j, pol));
        });
        layout.annotations.push({ text: s.name, x: (x0 + x1) / 2, y: Math.min(1, y1) + 0.025,
          xref: "paper", yref: "paper", showarrow: false, font: { color: PAL().muted, size: 12 } });
      });
      // кривые Ansys HFSS по плоскостям: E -> панель E (series[0]), H -> панель H (series[1])
      if (useDB && GT.state.reference) {
        const planeIdx = { E: 0, H: 1 };
        GT.state.reference.forEach((rf) => {
          const idx = planeIdx[rf.plane];
          if (idx == null || idx >= series.length) return;
          traces.push(refTrace(rf, idx === 0 ? "polar" : "polar" + (idx + 1)));
        });
      }
    }
    Plotly.react(host, traces, layout, config("greentensor_pattern"));
  }

  function normalizeAngle(x, n) {
    // авто-нормировка размерности: любой импортированный x -> равномерно в [0,360)
    if (x && x.length === n) {
      const finite = x.filter((v) => isFinite(v));
      const mn = Math.min.apply(null, finite), mx = Math.max.apply(null, finite);
      if (mx > mn && mx <= 360.5 && mn >= -0.5) return x;        // уже в градусах
      if (mx > mn) return x.map((v) => ((v - mn) / (mx - mn)) * 360); // линейно -> [0,360]
    }
    return Array.from({ length: n }, (_, i) => (i / n) * 360);   // только y -> индекс на круг
  }

  // ===================== Тепловая карта / ближнее поле ===================== //
  const DB_HM_FLOOR = -30;            // нижний предел лог-шкалы карты, дБ отн. максимума
  const ANIM_FRAMES = 30;             // кадров на период ωt = 0…2π (как в эталоне Example 13)
  const ANIM_DT = 70;                 // мс/кадр live-анимации
  const GIF_DT = 80;                  // мс/кадр GIF (эталон: 80 мс, 30 кадров)
  const GIF_W = 800, GIF_H = 600;
  const GIF_LEVELS = 48;              // градаций палитры поля в GIF (сжатие LZW ≤ 3 МБ)
  let animPhase = 0;                  // текущая фаза ωt, рад (сохраняется между кадрами)
  let wtAnnIdx = -1;                  // индекс annotation с подписью ωt (обновляется в анимации)

  // спокойная дивергентная палитра «синий → красный» (типа coolwarm), без ярких пиков;
  // ч-б — оттенки серого (низк.→свет., выс.→тёмн.)
  function heatScale() {
    return isBW()
      ? [[0, "#eef1f4"], [1, "#16202a"]]
      : [[0.0, "#3b5bbf"], [0.18, "#6a86d8"], [0.36, "#9fb2e6"],
         [0.50, "#c4c8d8"], [0.64, "#dcb39f"], [0.82, "#cc7155"],
         [1.0, "#a52a2a"]];   // синий → приглушённый красный
  }
  // Дивергентная RdBu (как matplotlib RdBu_r на слайде диссертации): −A синий, 0 белый, +A красный.
  // Ч-б: тёмный (−) → средне-серый (0) → светлый (+).
  function divScale() {
    return isBW()
      ? [[0, "#0b0f14"], [0.5, "#8a8f94"], [1, "#f7f8f9"]]
      : [[0.0, "#053061"], [0.1, "#2166ac"], [0.2, "#4393c3"], [0.3, "#92c5de"],
         [0.4, "#d1e5f0"], [0.5, "#f7f7f7"], [0.6, "#fddbc7"], [0.7, "#f4a582"],
         [0.8, "#d6604d"], [0.9, "#b2182b"], [1.0, "#67001f"]];
  }
  function fieldMax(absE) {
    let mx = 0;
    absE.forEach((row) => row.forEach((v) => { if (v != null && v > mx) mx = v; }));
    return mx > 0 ? mx : 1;
  }
  // амплитуда комплексной компоненты: max|re+i·im| — граница всех кадров (симметричная шкала)
  function compAmp(arr) {
    let a = 0;
    arr.re.forEach((row, r) => row.forEach((v, k) => {
      if (v == null) return;
      const m = Math.hypot(v, arr.im[r][k] || 0); if (m > a) a = m;
    }));
    return a > 0 ? a : 1;
  }
  // мгновенное поле E(t) = re·cos(ωt) + im·sin(ωt): физически одно и то же в обеих
  // конвенциях — движок e^{−iωt} даёт Re[E·e^{−iωt}], диссертация e^{+jωt} даёт Re[E*·e^{+jωt}].
  // С ростом ωt волна бежит в +z (фазор движка ∝ e^{+ikz}).
  function instFrame(arr, phase, absolute) {
    const c = Math.cos(phase), s = Math.sin(phase);
    return arr.re.map((row, r) => row.map((v, k) => {
      if (v == null) return null;
      const val = v * c + (arr.im[r][k] || 0) * s;
      return absolute ? Math.abs(val) : val;
    }));
  }

  // --- выбор компоненты/режима ---
  function effComp(hm) {
    const comp = GT.state.heatComp || "abs";
    if (GT.state.heatMode === "re" && comp === "abs") return hm.dominant || "x";   // знак у нормы не определён
    return comp;
  }
  function compArrays(hm, comp) {
    if (hm.comp && hm.comp[comp]) return hm.comp[comp];
    if (comp === (hm.dominant || "x") && hm.reField) return { re: hm.reField, im: hm.imField };
    return null;
  }
  function compHtml(comp) { return comp === "abs" ? "|E|" : "E<sub>" + comp + "</sub>"; }
  function fmtRu(v, d) { return Number(v).toFixed(d == null ? 2 : d).replace(/\.?0+$/, "").replace(".", ","); }
  function fmtDeg(phase) {
    const deg = Math.round(((phase % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI) * 180 / Math.PI);
    return "ωt = " + String(deg).padStart(3, " ") + "°";
  }

  // --- оверлей геометрии: контур сферы, пунктир границ слоёв, стрелка падения, подписи ---
  function overlayColor() { return "#ffffff"; }
  function overlayShapes(hm) {
    const g = hm.geometry;
    if (!g || !g.bodies) return [];
    const shapes = [];
    const white = overlayColor();
    g.bodies.forEach((b) => {
      const circ = (r, line) => ({ type: "circle", xref: "x", yref: "y",
        x0: b.cu - r, x1: b.cu + r, y0: b.cv - r, y1: b.cv + r, line, layer: "above" });
      // границы слоёв (кроме внешней) — пунктир
      (b.layers_r || []).forEach((r, i) => {
        if (r == null || i === b.layers_r.length - 1) return;
        shapes.push(circ(r, { color: "rgba(255,255,255,0.75)", width: 1, dash: "dash" }));
      });
      // внешний контур: тёмный «ореол» + белая линия (видна на любом фоне палитры)
      shapes.push(circ(b.r, { color: "rgba(0,0,0,0.35)", width: 4 }));
      shapes.push(circ(b.r, { color: white, width: 2 }));
    });
    return shapes;
  }
  function overlayAnnotations(hm, opts) {
    const g = hm.geometry || {};
    const P = PAL();
    const ann = [];
    const hs = hm.x, vs = hm.y;
    const hmin = hs[0], hmax = hs[hs.length - 1], vmin = vs[0], vmax = vs[vs.length - 1];
    const W = hmax - hmin, H = vmax - vmin;
    const boxStyle = { bgcolor: themeDark() ? "rgba(7,9,10,0.62)" : "rgba(255,255,255,0.82)",
      bordercolor: themeDark() ? "rgba(255,255,255,0.35)" : "rgba(0,0,0,0.35)", borderwidth: 1,
      borderpad: 3, font: { color: P.font, size: 12 } };
    // заголовок
    ann.push({ x: 0.5, y: 1.06, xref: "paper", yref: "paper", showarrow: false,
      text: opts.title, font: { color: P.muted, size: 12 } });
    // стрелка направления падения (проекция k̂ на плоскость) с подписью E_пад, k
    const kd = g.k_dir || [1, 0];
    const kl = Math.hypot(kd[0], kd[1]);
    // стрелка центрирована на якоре в левом верхнем углу кадра — не выходит за рамку
    // при любом направлении k̂ в плоскости (в т. ч. вертикальном)
    const cx = hmin + 0.16 * W, cy = vmax - 0.13 * H;
    const label = g.pol_inplane === false
      ? "<b>k</b> (<b>E</b><sub>пад</sub> ⊥ плоскости)" : "<b>E</b><sub>пад</sub>, <b>k</b>";
    if (g.k_inplane !== false && kl > 1e-6) {
      const len = 0.18 * W, dx = kd[0] / kl, dy = kd[1] / kl;
      const ax0 = cx - dx * len / 2, ay0 = cy - dy * len / 2;
      const tx = cx + dx * len / 2, ty = cy + dy * len / 2;
      ann.push({ x: tx, y: ty, ax: ax0, ay: ay0, xref: "x", yref: "y", axref: "x", ayref: "y",
        showarrow: true, arrowhead: 2, arrowsize: 1.2, arrowwidth: 2.5, arrowcolor: "#ffffff",
        text: label, xanchor: "left", yshift: 14, font: { color: "#ffffff", size: 13 },
        bgcolor: "rgba(0,0,0,0.30)", borderpad: 2 });
    } else {
      ann.push({ x: hmin + 0.06 * W, y: cy, xref: "x", yref: "y", showarrow: false, xanchor: "left",
        text: "⊙ <b>k</b> ⊥ плоскости", font: { color: "#ffffff", size: 13 },
        bgcolor: "rgba(0,0,0,0.30)", borderpad: 2 });
    }
    // плашка слоёв: «4 слоя, ε = 1,94 / 1,75 / 1,44 / 1,0»
    const bodies = g.bodies || [];
    let layersText = "";
    if (bodies.length === 1) {
      const b = bodies[0];
      const eps = (b.eps_re || []).map((e, i) => fmtRu(e, 2) + ((b.eps_im || [])[i] ? "−j" + fmtRu(Math.abs(b.eps_im[i]), 2) : ""));
      const n = eps.length;
      const sl = n === 1 ? "однородная" : n + (n >= 2 && n <= 4 ? " слоя" : " слоёв");
      layersText = (hm.quantity === "scattered" ? "рассеянное поле · " : "") + sl + ", ε = " + eps.join(" / ");
    } else if (bodies.length > 1) {
      layersText = "кластер: " + bodies.length + " тел · рассеянное поле";
    }
    if (layersText) ann.push(Object.assign({ x: 0.5, y: 0.02, xref: "paper", yref: "paper",
      showarrow: false, xanchor: "center", yanchor: "bottom", text: layersText }, boxStyle));
    // счётчик ωt (только для мгновенного поля)
    wtAnnIdx = -1;
    if (opts.showPhase) {
      wtAnnIdx = ann.length;
      ann.push(Object.assign({ x: 0.99, y: 0.02, xref: "paper", yref: "paper", showarrow: false,
        xanchor: "right", yanchor: "bottom", text: fmtDeg(opts.phase || 0) },
        boxStyle, { font: { color: P.font, size: 12, family: "ui-monospace, SFMono-Regular, Menlo, monospace" } }));
    }
    return ann;
  }

  function heatLayout(hm, opts) {
    const P = PAL();
    return baseLayout({
      showlegend: false,
      margin: { l: 56, r: 24, t: 34, b: 48 },
      xaxis: { title: hm.xlabel || "x", scaleanchor: "y", constrain: "domain", zeroline: false, gridcolor: P.grid },
      yaxis: { title: hm.ylabel || "z", zeroline: false, gridcolor: P.grid },
      shapes: overlayShapes(hm),
      annotations: overlayAnnotations(hm, opts),
    });
  }

  function logCheckOnce(hm) {
    if (!hm.check || hm._checkLogged) return;
    hm._checkLogged = true;
    const c = hm.check;
    if (c.error) { console.warn("[GreenTensor] проверка ближнего поля не выполнена:", c.error); return; }
    const okQ = c.q_sca_rel_err != null && c.q_sca_rel_err < 1e-6;
    const okC = c.continuity_max_rel_jump != null && c.continuity_max_rel_jump < 1e-3;
    console.info("[GreenTensor] validate(): Q_sca = " + c.q_sca + " (ядро 01_sphere: " + c.q_sca_core +
      ", отн. ошибка " + (c.q_sca_rel_err == null ? "—" : c.q_sca_rel_err.toExponential(2)) + ") " + (okQ ? "✓" : "✗") +
      "; непрерывность E_" + c.continuity_component + " на границах слоёв вдоль k̂: max отн. скачок " +
      (c.continuity_max_rel_jump == null ? "—" : c.continuity_max_rel_jump.toExponential(2)) + " " + (okC ? "✓ (≤1e-3)" : "✗"),
      c.continuity);
  }

  function renderHeatmap() {
    const res = GT.state.results;
    const host = GT.$("#plot");
    if (!res) return showEmpty(true);
    // Карта считается ПО ЗАПРОСУ (тяжёлая часть). Ключ отсутствует → ещё не считали: кнопка.
    if (!("heatmap" in res)) {
      showEmpty(true);
      GT.$("#plot-empty").innerHTML =
        "Тепловая карта ближнего поля считается отдельно (самая тяжёлая часть расчёта — " +
        "сетка |E| по узлам). " +
        '<button class="btn" id="hm-go" type="button">Рассчитать тепловую карту</button>';
      const b = GT.$("#hm-go");
      if (b) b.onclick = () => GT.compute({ pattern: false, cross: false, sweep: false,
        heatmap: true, merge: true });
      return;
    }
    if (!res.heatmap) {                              // null → поля нет для этой задачи
      showEmpty(true);
      GT.$("#plot-empty").innerHTML =
        "Тепловая карта недоступна для этой задачи: решатель не даёт ближнее поле " +
        "(2D-цилиндр, полусфера на экране). Используйте «Диаграмма (ДН)» и другие вкладки.";
      return;
    }
    showEmpty(false);
    const hm = res.heatmap;
    logCheckOnce(hm);
    syncHeatControls(hm);
    const mode = GT.state.heatMode || "abs";
    const comp = effComp(hm);
    const arr = comp === "abs" ? null : compArrays(hm, comp);
    const planeTxt = "плоскость " + (hm.plane_title || hm.plane || "xz");
    const resTxt = hm.res ? " · " + hm.res + "×" + hm.res : "";
    const maskTxt = hm.metal_masked ? " · металл-области скрыты" : "";
    const qty = hm.quantity === "scattered" ? "рассеянное поле" : "ближнее поле";
    let trace, title;
    if (mode === "re" && arr) {
      // мгновенное поле СО ЗНАКОМ; шкала симметричная и фиксированная на всю анимацию
      const A = compAmp(arr);
      trace = { type: "heatmap", x: hm.x, y: hm.y, z: instFrame(arr, animPhase, false),
        colorscale: divScale(), zmin: -A, zmax: A, zsmooth: "best",
        colorbar: { title: { text: "Re[E<sub>" + comp + "</sub>·e<sup>jωt</sup>], отн. ед." }, outlinewidth: 0 } };
      title = "Re[E<sub>" + comp + "</sub>·e<sup>jωt</sup>] (конв. дисс.) · " + qty + " · " + planeTxt + resTxt + maskTxt;
    } else {
      let absZ;
      if (comp === "abs" || !arr) absZ = hm.absE;
      else absZ = arr.re.map((row, r) => row.map((v, k) => (v == null ? null : Math.hypot(v, arr.im[r][k] || 0))));
      const mx = fieldMax(absZ);                    // нормировка на максимум поля
      const log = GT.state.heatScale !== "lin";
      let z, cbTitle, zmin, zmax;
      if (log) {                                    // лог (дБ отн. макс): и пики, и переходы слоёв
        z = absZ.map((r) => r.map((v) => (v == null ? null : Math.max(DB_HM_FLOOR, 20 * Math.log10(Math.max(v, 1e-300) / mx)))));
        cbTitle = "дБ отн. макс"; zmin = DB_HM_FLOOR; zmax = 0;
      } else {                                      // линейная, нормированная
        z = absZ.map((r) => r.map((v) => (v == null ? null : v / mx)));
        cbTitle = compHtml(comp) + " / макс"; zmin = 0; zmax = 1;
      }
      trace = { type: "heatmap", x: hm.x, y: hm.y, z, colorscale: heatScale(),
        zmin, zmax, colorbar: { title: { text: cbTitle }, outlinewidth: 0 }, zsmooth: "best" };
      title = (comp === "abs" ? "|E|" : "|E<sub>" + comp + "</sub>|") + " · " + qty + " · " +
        (log ? "лог (дБ отн. макс)" : "линейная, норм.") + " · " + planeTxt + resTxt + maskTxt;
    }
    Plotly.react(host, [trace], heatLayout(hm, { title, showPhase: mode === "re" && !!arr, phase: animPhase }),
      config(mode === "re" ? "greentensor_field" : "greentensor_heatmap"));
  }

  // --- анимация бегущей волны: кадры Re[E·e^{jωt}] по ωt = 0…2π, шкала фиксирована ---
  function setAnimBtn(on) {
    const btn = GT.$("#anim-btn");
    if (btn) btn.textContent = on ? "⏸ Стоп" : "▶ Бегущая волна";
  }
  function setHeatMode(mode) {
    GT.state.heatMode = mode;
    GT.$$("#hm-mode button").forEach((x) => x.classList.toggle("active", x.dataset.hm === mode));
  }
  function startAnim() {
    const res = GT.state.results;
    if (!res || !res.heatmap) return false;
    const hm = res.heatmap;
    if (GT.state.heatMode !== "re") setHeatMode("re");     // бегущая волна = поле со знаком
    const arr = compArrays(hm, effComp(hm));
    if (!arr) return false;
    stopAnim();
    renderHeatmap();
    setAnimBtn(true);
    const host = GT.$("#plot");
    const step = 2 * Math.PI / ANIM_FRAMES;
    animTimer = setInterval(() => {
      if (!host._fullLayout) { stopAnim(); setAnimBtn(false); return; }   // график очищен
      animPhase = (animPhase + step) % (2 * Math.PI);
      const upd = wtAnnIdx >= 0 ? { ["annotations[" + wtAnnIdx + "].text"]: fmtDeg(animPhase) } : {};
      Plotly.update(host, { z: [instFrame(arr, animPhase, false)] }, upd, [0]).catch(() => {});
    }, ANIM_DT);
    return true;
  }
  function toggleAnim() {
    if (animTimer) { stopAnim(); setAnimBtn(false); return; }
    startAnim();
  }

  // ===================== GIF: кадры на offscreen-canvas → GT.GIFEncoder ===================== //
  function hex2rgb(h) { const v = parseInt(h.slice(1), 16); return [(v >> 16) & 255, (v >> 8) & 255, v & 255]; }
  function scaleLUT(scale, levels) {                 // Plotly-colorscale → массив [r,g,b] длиной levels
    const stops = scale.map((s) => [s[0], hex2rgb(s[1])]);
    const lut = [];
    for (let i = 0; i < levels; i++) {
      const t = levels === 1 ? 0 : i / (levels - 1);
      let j = 0; while (j < stops.length - 2 && t > stops[j + 1][0]) j++;
      const [t0, c0] = stops[j], [t1, c1] = stops[j + 1];
      const u = t1 > t0 ? Math.min(1, Math.max(0, (t - t0) / (t1 - t0))) : 0;
      lut.push([0, 1, 2].map((k) => Math.round(c0[k] + (c1[k] - c0[k]) * u)));
    }
    return lut;
  }
  function niceTicks(lo, hi, n) {
    const raw = (hi - lo) / n, p = Math.pow(10, Math.floor(Math.log10(raw)));
    let step = 10 * p;
    for (const c of [1, 2, 2.5, 5, 10]) if (raw <= c * p) { step = c * p; break; }
    const out = [];
    for (let v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }
  // текст с под-/надстрочными фрагментами: parts = [{t:"Re[E"},{t:"x",sub:1},{t:"·e"},{t:"jωt",sup:1},{t:"]"}]
  function drawRich(ctx, parts, x, y, size, color, bold) {
    ctx.save(); ctx.fillStyle = color; ctx.textBaseline = "alphabetic"; ctx.textAlign = "left";
    let cx = x;
    parts.forEach((p) => {
      const s = p.sub || p.sup ? Math.round(size * 0.7) : size;
      ctx.font = (p.bold || bold ? "bold " : "") + (p.italic ? "italic " : "") + s + "px Inter, -apple-system, Segoe UI, Roboto, sans-serif";
      const dy = p.sub ? size * 0.25 : p.sup ? -size * 0.45 : 0;
      ctx.fillText(p.t, cx, y + dy);
      cx += ctx.measureText(p.t).width;
    });
    ctx.restore();
    return cx - x;
  }
  function richWidth(ctx, parts, size, bold) {
    let w = 0;
    parts.forEach((p) => {
      const s = p.sub || p.sup ? Math.round(size * 0.7) : size;
      ctx.font = (p.bold || bold ? "bold " : "") + (p.italic ? "italic " : "") + s + "px Inter, -apple-system, Segoe UI, Roboto, sans-serif";
      w += ctx.measureText(p.t).width;
    });
    return w;
  }

  // Рисует один кадр анимации 800×600 (стиль слайда: заголовок, оси, контур, стрелка, шкала, ωt).
  // lut — палитра поля; возвращает индексы палитры для области карты не нужны: квантуется энкодером.
  function drawGifFrame(ctx, hm, arr, comp, phase, A, lut, colors) {
    const Wc = GIF_W, Hc = GIF_H;
    ctx.fillStyle = colors.paper; ctx.fillRect(0, 0, Wc, Hc);
    // область карты: квадрат
    const px0 = 118, py0 = 58, S = 490;
    const hs = hm.x, vs = hm.y, n = hs.length, m = vs.length;
    const hmin = hs[0], hmax = hs[n - 1], vmin = vs[0], vmax = vs[m - 1];
    const c = Math.cos(phase), s = Math.sin(phase);
    const img = ctx.createImageData(S, S);
    const d = img.data;
    const re = arr.re, im = arr.im;
    const bg = hex2rgb(colors.paper);
    const L = lut.length;
    for (let py = 0; py < S; py++) {
      const gv = (1 - (py + 0.5) / S) * (m - 1);            // строка 0 — верх = vmax
      const j0 = Math.min(m - 2, Math.floor(gv)), tv = gv - j0;
      for (let px = 0; px < S; px++) {
        const gh = (px + 0.5) / S * (n - 1);
        const i0 = Math.min(n - 2, Math.floor(gh)), th = gh - i0;
        const r00 = re[j0][i0], r01 = re[j0][i0 + 1], r10 = re[j0 + 1][i0], r11 = re[j0 + 1][i0 + 1];
        const o = (py * S + px) * 4;
        if (r00 == null || r01 == null || r10 == null || r11 == null) {
          d[o] = bg[0]; d[o + 1] = bg[1]; d[o + 2] = bg[2]; d[o + 3] = 255; continue;
        }
        const v00 = r00 * c + im[j0][i0] * s, v01 = r01 * c + im[j0][i0 + 1] * s;
        const v10 = r10 * c + im[j0 + 1][i0] * s, v11 = r11 * c + im[j0 + 1][i0 + 1] * s;
        const val = (v00 * (1 - th) + v01 * th) * (1 - tv) + (v10 * (1 - th) + v11 * th) * tv;
        let idx = Math.round((val + A) / (2 * A) * (L - 1));
        if (idx < 0) idx = 0; else if (idx >= L) idx = L - 1;
        const col = lut[idx];
        d[o] = col[0]; d[o + 1] = col[1]; d[o + 2] = col[2]; d[o + 3] = 255;
      }
    }
    ctx.putImageData(img, px0, py0);
    const X = (h) => px0 + (h - hmin) / (hmax - hmin) * S;
    const Y = (v) => py0 + (vmax - v) / (vmax - vmin) * S;
    const g = hm.geometry || {};
    // геометрия: пунктир слоёв + контур
    (g.bodies || []).forEach((b) => {
      ctx.save();
      (b.layers_r || []).forEach((r, i) => {
        if (r == null || i === b.layers_r.length - 1) return;
        ctx.beginPath(); ctx.setLineDash([6, 5]); ctx.lineWidth = 1.2; ctx.strokeStyle = "rgba(255,255,255,0.8)";
        ctx.arc(X(b.cu), Y(b.cv), r / (hmax - hmin) * S, 0, 2 * Math.PI); ctx.stroke();
      });
      ctx.setLineDash([]);
      ctx.beginPath(); ctx.lineWidth = 5; ctx.strokeStyle = "rgba(0,0,0,0.35)";
      ctx.arc(X(b.cu), Y(b.cv), b.r / (hmax - hmin) * S, 0, 2 * Math.PI); ctx.stroke();
      ctx.beginPath(); ctx.lineWidth = 2.5; ctx.strokeStyle = "#ffffff";
      ctx.arc(X(b.cu), Y(b.cv), b.r / (hmax - hmin) * S, 0, 2 * Math.PI); ctx.stroke();
      ctx.restore();
    });
    // рамка области и оси с делениями
    ctx.save();
    ctx.strokeStyle = colors.axis; ctx.lineWidth = 1; ctx.strokeRect(px0 + 0.5, py0 + 0.5, S, S);
    ctx.fillStyle = colors.font; ctx.font = "13px Inter, -apple-system, Segoe UI, Roboto, sans-serif";
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    niceTicks(hmin, hmax, 6).forEach((t) => {
      const x = X(t); ctx.beginPath(); ctx.moveTo(x, py0 + S); ctx.lineTo(x, py0 + S + 5); ctx.stroke();
      ctx.fillText(fmtRu(t, 2).replace("-", "−"), x, py0 + S + 8);
    });
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    niceTicks(vmin, vmax, 6).forEach((t) => {
      const y = Y(t); ctx.beginPath(); ctx.moveTo(px0 - 5, y); ctx.lineTo(px0, y); ctx.stroke();
      ctx.fillText(fmtRu(t, 2).replace("-", "−"), px0 - 9, y);
    });
    ctx.font = "italic 16px Inter, -apple-system, Segoe UI, Roboto, sans-serif";
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.fillText(hm.xlabel || "x", px0 + S / 2, py0 + S + 26);
    ctx.save(); ctx.translate(px0 - 58, py0 + S / 2); ctx.rotate(-Math.PI / 2);
    ctx.textBaseline = "middle"; ctx.fillText(hm.ylabel || "z", 0, 0); ctx.restore();
    ctx.restore();
    // заголовок
    const kA = hm.k != null && hm.radius != null ? fmtRu(hm.k * hm.radius, 2) : null;
    const titleParts = [{ t: "Ближнее поле: " + (hm.quantity === "scattered" ? "рассеянное" : "полное") + ", Re[E" },
      { t: comp, sub: 1, italic: 1 }, { t: "·e" }, { t: "jωt", sup: 1, italic: 1 }, { t: "]" },
      { t: kA ? "  ·  k₀a = " + kA : "" }, { t: "  (GreenTensor)" }];
    const tw = richWidth(ctx, titleParts, 19);
    drawRich(ctx, titleParts, px0 + S / 2 - tw / 2, 36, 19, colors.font);
    // стрелка E_пад, k
    const kd = g.k_dir || [1, 0], kl = Math.hypot(kd[0], kd[1]);
    const cxp = X(hmin + 0.16 * (hmax - hmin)), cyp = Y(vmax - 0.13 * (vmax - vmin));
    ctx.save(); ctx.strokeStyle = "#ffffff"; ctx.fillStyle = "#ffffff"; ctx.lineWidth = 3;
    if (g.k_inplane !== false && kl > 1e-6) {
      const len = 0.18 * S, dx = kd[0] / kl, dy = -kd[1] / kl;   // экран: y вниз
      const ax0 = cxp - dx * len / 2, ay0 = cyp - dy * len / 2;
      const tx = cxp + dx * len / 2, ty = cyp + dy * len / 2;
      ctx.beginPath(); ctx.moveTo(ax0, ay0); ctx.lineTo(tx, ty); ctx.stroke();
      const ang = Math.atan2(dy, dx);
      ctx.beginPath(); ctx.moveTo(tx, ty);
      ctx.lineTo(tx - 12 * Math.cos(ang - 0.45), ty - 12 * Math.sin(ang - 0.45));
      ctx.lineTo(tx - 12 * Math.cos(ang + 0.45), ty - 12 * Math.sin(ang + 0.45)); ctx.closePath(); ctx.fill();
      const lbl = g.pol_inplane === false
        ? [{ t: "k", bold: 1 }, { t: " (E", bold: 1 }, { t: "пад", sub: 1 }, { t: " ⊥ плоскости)" }]
        : [{ t: "E", bold: 1 }, { t: "пад", sub: 1 }, { t: ", " }, { t: "k", bold: 1 }];
      drawRich(ctx, lbl, Math.min(ax0, tx), Math.min(ay0, ty) - 10, 15, "#ffffff");
    } else {
      drawRich(ctx, [{ t: "⊙ " }, { t: "k", bold: 1 }, { t: " ⊥ плоскости" }],
        X(hmin + 0.06 * (hmax - hmin)), cyp, 15, "#ffffff");
    }
    ctx.restore();
    // плашка слоёв + ωt
    const bodies = g.bodies || [];
    let layersText = "";
    if (bodies.length === 1) {
      const b = bodies[0];
      const eps = (b.eps_re || []).map((e) => fmtRu(e, 2));
      const nl = eps.length;
      layersText = (nl === 1 ? "однородная" : nl + (nl >= 2 && nl <= 4 ? " слоя" : " слоёв")) + ", ε = " + eps.join(" / ");
    } else if (bodies.length > 1) layersText = "кластер: " + bodies.length + " тел · рассеянное поле";
    const box = (text, xRight, yBottom, mono) => {
      ctx.save();
      ctx.font = (mono ? "13px ui-monospace, SFMono-Regular, Menlo, monospace" : "13px Inter, -apple-system, Segoe UI, Roboto, sans-serif");
      const w = ctx.measureText(text).width + 14, h = 24;
      const x = xRight - w, y = yBottom - h;
      ctx.fillStyle = colors.box; ctx.fillRect(x, y, w, h);
      ctx.strokeStyle = colors.boxBorder; ctx.lineWidth = 1; ctx.strokeRect(x + 0.5, y + 0.5, w - 1, h - 1);
      ctx.fillStyle = colors.font; ctx.textAlign = "left"; ctx.textBaseline = "middle";
      ctx.fillText(text, x + 7, y + h / 2 + 1);
      ctx.restore();
      return w;
    };
    const wtW = box(fmtDeg(phase), px0 + S - 8, py0 + S - 8, true);
    if (layersText) box(layersText, px0 + S - 8 - wtW - 8, py0 + S - 8, false);
    // цветовая шкала
    const cbx = px0 + S + 28, cby = py0, cbw = 22, cbh = S;
    for (let i = 0; i < cbh; i++) {
      const t = 1 - i / (cbh - 1);
      const col = lut[Math.round(t * (L - 1))];
      ctx.fillStyle = "rgb(" + col[0] + "," + col[1] + "," + col[2] + ")";
      ctx.fillRect(cbx, cby + i, cbw, 1);
    }
    ctx.save();
    ctx.strokeStyle = colors.axis; ctx.strokeRect(cbx + 0.5, cby + 0.5, cbw, cbh);
    ctx.fillStyle = colors.font; ctx.font = "13px Inter, -apple-system, Segoe UI, Roboto, sans-serif";
    ctx.textAlign = "left"; ctx.textBaseline = "middle";
    niceTicks(-A, A, 5).forEach((t) => {
      const y = cby + (1 - (t + A) / (2 * A)) * cbh;
      ctx.beginPath(); ctx.moveTo(cbx + cbw, y); ctx.lineTo(cbx + cbw + 4, y); ctx.stroke();
      ctx.fillText(fmtRu(t, 2).replace("-", "−"), cbx + cbw + 7, y);
    });
    ctx.translate(cbx + cbw + 62, cby + cbh / 2); ctx.rotate(-Math.PI / 2);
    const cbParts = [{ t: "Re[E" }, { t: comp, sub: 1, italic: 1 }, { t: "·e" }, { t: "jωt", sup: 1, italic: 1 }, { t: "], отн. ед." }];
    const cw = richWidth(ctx, cbParts, 15);
    drawRich(ctx, cbParts, -cw / 2, 0, 15, colors.font);
    ctx.restore();
  }

  async function exportGIF() {
    const res = GT.state.results;
    if (!res || !res.heatmap) { GT.toast("Сначала рассчитайте тепловую карту", "error"); return; }
    if (!GT.GIFEncoder) { GT.toast("GIF-энкодер не загружен", "error"); return; }
    const hm = res.heatmap;
    if (GT.state.heatMode !== "re") setHeatMode("re");
    const comp = effComp(hm);
    const arr = compArrays(hm, comp);
    if (!arr) { GT.toast("Для этой карты нет комплексного поля (re/im)", "error"); return; }
    const wasAnim = !!animTimer;
    stopAnim(); setAnimBtn(false);
    const P = PAL();
    const colors = { paper: P.paper, font: P.font, axis: P.axis,
      box: themeDark() ? "rgba(7,9,10,0.62)" : "rgba(255,255,255,0.85)",
      boxBorder: themeDark() ? "rgba(255,255,255,0.35)" : "rgba(0,0,0,0.35)" };
    const lut = scaleLUT(divScale(), GIF_LEVELS);
    // палитра GIF: градации поля + цвета темы + переходы (сглаживание текста/линий)
    const pal = lut.slice();
    const bg = hex2rgb(P.paper), fg = hex2rgb(P.font);
    const blend = (a, b, t) => [0, 1, 2].map((k) => Math.round(a[k] + (b[k] - a[k]) * t));
    pal.push(bg, fg, [255, 255, 255], [0, 0, 0], hex2rgb(P.muted));
    for (let i = 1; i < 24; i++) pal.push(blend(bg, fg, i / 24));                  // текст на фоне
    for (let i = 1; i < 12; i++) pal.push(blend(fg, [255, 255, 255], i / 12));
    // белый контур/стрелка поверх поля: смеси белого с центральными цветами палитры
    const mids = [lut[0], lut[Math.floor(L4(lut, 0.25))], lut[Math.floor(L4(lut, 0.5))], lut[Math.floor(L4(lut, 0.75))], lut[lut.length - 1]];
    mids.forEach((c) => { for (let i = 1; i < 6; i++) pal.push(blend(c, [255, 255, 255], i / 6)); });
    while (pal.length > 256) pal.pop();
    const A = compAmp(arr);
    const canvas = document.createElement("canvas");
    canvas.width = GIF_W; canvas.height = GIF_H;
    const ctx = canvas.getContext("2d", { willReadFrequently: true });
    const enc = new GT.GIFEncoder(GIF_W, GIF_H, pal, { loop: 0, delay: GIF_DT });
    GT.busy(true, "Сборка GIF… 0/" + ANIM_FRAMES);
    const timing = { draw: 0, encode: 0, total: 0 };
    const t0 = performance.now();
    try {
      for (let f = 0; f < ANIM_FRAMES; f++) {
        const phase = 2 * Math.PI * f / ANIM_FRAMES;
        const ta = performance.now();
        drawGifFrame(ctx, hm, arr, comp, phase, A, lut, colors);
        const tb = performance.now();
        enc.addFrame(ctx.getImageData(0, 0, GIF_W, GIF_H), GIF_DT);
        timing.draw += tb - ta; timing.encode += performance.now() - tb;
        GT.busy(true, "Сборка GIF… " + (f + 1) + "/" + ANIM_FRAMES);
        await yieldToUI();
      }
      const blob = enc.finish();
      timing.total = performance.now() - t0;
      const name = "greentensor_nearfield_wave_" + (hm.plane || "xz") + "_E" + comp + ".gif";
      GT.download(name, blob, "image/gif");
      GT.toast("GIF: " + ANIM_FRAMES + " кадров, " + (blob.size / 1048576).toFixed(2) + " МБ", "ok");
      GT.lastGif = { name, size: blob.size, frames: ANIM_FRAMES, blob, timing };
      console.info("[GreenTensor] GIF " + name + ": " + ANIM_FRAMES + " кадров, " +
        (blob.size / 1024).toFixed(0) + " кБ, " + timing.total.toFixed(0) + " мс (рисование " +
        timing.draw.toFixed(0) + ", кодирование " + timing.encode.toFixed(0) + ")");
    } catch (e) {
      console.error(e);
      GT.toast("Не удалось собрать GIF: " + e.message, "error");
    } finally {
      GT.busy(false);
      if (wasAnim) startAnim();
    }
  }
  function L4(lut, t) { return t * (lut.length - 1); }
  // Уступить поток UI между кадрами. НЕ setTimeout: в фоновой/скрытой вкладке Chrome
  // троттлит таймеры до 1 раза в секунду (30 кадров → полминуты); MessageChannel не троттлится.
  function yieldToUI() {
    return new Promise((resolve) => {
      const ch = new MessageChannel();
      ch.port1.onmessage = () => { ch.port1.close(); resolve(); };
      ch.port2.postMessage(0);
    });
  }

  // --- синхронизация контролов карты с состоянием ---
  function syncHeatControls(hm) {
    const mode = GT.state.heatMode || "abs";
    GT.$$("#hm-mode button").forEach((x) => x.classList.toggle("active", x.dataset.hm === mode));
    const comp = GT.state.heatComp || "abs";
    GT.$$("#hm-comp button").forEach((x) => {
      x.classList.toggle("active", x.dataset.hc === comp);
      x.disabled = mode === "re" && x.dataset.hc === "abs";         // знак у |E| не определён
    });
    const plane = GT.state.heatPlane || "xz";
    GT.$$("#hm-plane button").forEach((x) => x.classList.toggle("active", x.dataset.hp === plane));
    const sc = GT.$("#hm-scale"); if (sc) sc.style.display = mode === "re" ? "none" : "";
    const hasComplex = !!(hm && (hm.comp || hm.reField));
    const gif = GT.$("#btn-gif"); if (gif) gif.disabled = !hasComplex;
    const anim = GT.$("#anim-btn"); if (anim) anim.disabled = !hasComplex;
  }

  // ===================== Сечения ===================== //
  function renderCross() {
    const res = GT.state.results;
    const host = GT.$("#plot");
    if (!res || !res.cross_sections) return showEmpty(true);
    showEmpty(false);
    const cs = res.cross_sections;
    const keys = Object.keys(cs);
    const labels = { q_sca: "Q_sca", q_ext: "Q_ext", q_abs: "Q_abs", q_back: "Q_back",
      c_sca: "C_sca", c_ext: "C_ext", c_abs: "C_abs",
      q_sca_TM: "Q_sca (TM)", q_ext_TM: "Q_ext (TM)", q_abs_TM: "Q_abs (TM)",
      q_sca_TE: "Q_sca (TE)", q_ext_TE: "Q_ext (TE)", q_abs_TE: "Q_abs (TE)" };
    const x = keys.map((k) => labels[k] || k);
    const y = keys.map((k) => cs[k]);
    const trace = {
      type: "bar", x, y,
      marker: { color: isBW() ? PAL().bwLine : PAL().accent, line: { color: PAL().axis, width: 1 } },
      text: y.map((v) => (v == null ? "—" : v.toFixed(4))), textposition: "outside",
    };
    Plotly.react(host, [trace], baseLayout({
      showlegend: false, yaxis: { title: "эффективность / сечение", gridcolor: PAL().grid },
      xaxis: { gridcolor: PAL().grid },
    }), config("greentensor_cross"));
  }

  // ===================== ЭПР (полярная: σ/λ², дБ) ===================== //
  function noteEmpty(html) {
    showEmpty(true);
    GT.$("#plot-empty").innerHTML = html;
  }

  // радиальный дБ-диапазон полярной ЭПР: до 10 дБ вверх, динамика ≤ 50 дБ
  function rcsDbRange(arrays) {
    const vals = [];
    arrays.forEach((a) => a.forEach((v) => { if (v != null && isFinite(v)) vals.push(v); }));
    if (!vals.length) return [-40, 0];
    const mx = Math.ceil(Math.max.apply(null, vals) / 10) * 10;
    let mn = Math.floor(Math.max(Math.min.apply(null, vals), mx - 50) / 10) * 10;
    if (mn >= mx) mn = mx - 10;
    return [mn, mx];
  }

  function rcsPolarLayout(rot, range, annText) {
    return baseLayout({
      margin: { l: 24, r: 24, t: 34, b: 20 },
      annotations: annText ? [{ x: 0.5, y: 1.06, xref: "paper", yref: "paper", showarrow: false,
        text: annText, font: { color: PAL().muted, size: 12 } }] : [],
      polar: {
        bgcolor: PAL().polar,
        angularaxis: { rotation: rot, direction: "counterclockwise", ticksuffix: "°",
          gridcolor: PAL().grid, linecolor: PAL().axis },
        radialaxis: { range: range, title: { text: "σ/λ², дБ" }, angle: 90,
          gridcolor: PAL().grid, linecolor: PAL().axis },
      },
    });
  }

  function renderBistatic() {
    const res = GT.state.results, host = GT.$("#plot");
    if (!res) return showEmpty(true);
    const bi = res.bistatic;
    if (!bi || bi.available === false) {
      return noteEmpty((bi && bi.note) || "Бистатическая ЭПР недоступна для этого тела.");
    }
    showEmpty(false);
    const rot = parseFloat(GT.$("#polar-rot").value) || 0;
    const range = rcsDbRange(bi.series.map((s) => s.dB));
    const traces = bi.series.map((s, i) => ({ type: "scatterpolar", mode: "lines",
      name: s.name, theta: bi.angle_deg, r: s.dB, line: lineStyle(i) }));
    Plotly.react(host, traces, rcsPolarLayout(rot, range,
      "Бистатическая ЭПР σ/λ² (дБ) при фиксированном падении"), config("greentensor_bistatic"));
  }

  function renderMonostatic() {
    const res = GT.state.results, host = GT.$("#plot");
    if (!res) return showEmpty(true);
    const mo = res.monostatic;
    if (mo && mo.available === false) {
      return noteEmpty(mo.note || "Моностатическая ЭПР недоступна для этого тела.");
    }
    if (!mo) {                                       // кластер: требует развёртки ракурса (по запросу)
      noteEmpty('Моностатическая ЭПР (обратное рассеяние) требует развёртки ракурса падения — ' +
        'пересчёта рассеяния под каждым углом. ' +
        '<button class="btn" id="mono-go" type="button">Рассчитать</button>');
      const b = GT.$("#mono-go");
      if (b) b.onclick = () => GT.compute({ pattern: false, heatmap: false, sweep: false,
        cross: false, monostatic: true, merge: true });
      return;
    }
    showEmpty(false);
    const rot = parseFloat(GT.$("#polar-rot").value) || 0;
    let theta = mo.aspect_deg.slice();
    let dB = mo.sigma.map((v) => (v == null ? null : 10 * Math.log10(Math.max(v, 1e-300))));
    // развёртка ракурса 0..180 → полный круг по осевой симметрии тела: σ(360−θ)=σ(θ)
    if (theta[theta.length - 1] <= 270) {
      const mt = [], mr = [];
      for (let i = theta.length - 2; i >= 1; i--) { mt.push(360 - theta[i]); mr.push(dB[i]); }
      theta = theta.concat(mt); dB = dB.concat(mr);
    }
    const range = rcsDbRange([dB]);
    const layout = rcsPolarLayout(rot, range, mo.note || "Моностатическая ЭПР (обратное рассеяние)");
    layout.showlegend = false;
    Plotly.react(host, [{ type: "scatterpolar", mode: "lines", name: "σ моност.",
      theta: theta, r: dB, line: lineStyle(0) }], layout, config("greentensor_monostatic"));
  }

  // ===================== Спектр Q(k) ===================== //
  function renderSweep() {
    const res = GT.state.results;
    const host = GT.$("#plot");
    if (!res || !res.sweep) {
      showEmpty(true);
      GT.$("#plot-empty").innerHTML =
        'Спектр не рассчитан. <button class="btn" id="sweep-go" type="button">Рассчитать спектр Q(k)</button>';
      const b = GT.$("#sweep-go");
      if (b) b.onclick = () => GT.compute({ pattern: false, heatmap: false, sweep: true });
      return;
    }
    showEmpty(false);
    const sw = res.sweep;
    // Генерик-кривая (не Q(k)): произвольные оси/серии — например η₁(θ₀)
    // режима «раскрыв облучателя» (spill-over, green_tensor/feed.py).
    if (sw.generic) {
      const traces = (sw.series || []).map((s, i) => ({
        type: "scatter", mode: s.mode || "lines", name: s.name,
        x: s.x || sw.x, y: s.y, line: lineStyle(i),
      }));
      (sw.marks || []).forEach((m, j) => traces.push({
        type: "scatter", mode: "markers", name: m.name, x: m.x, y: m.y,
        marker: { size: 9, symbol: j ? "square" : "circle" },
      }));
      Plotly.react(host, traces, baseLayout({
        xaxis: { title: sw.xlabel || "x", gridcolor: PAL().grid },
        yaxis: { title: sw.ylabel || "y", gridcolor: PAL().grid },
      }), config("greentensor_curve"));
      return;
    }
    const keys = ["q_sca", "q_ext", "q_abs", "q_back"];
    const traces = keys.filter((k) => sw[k]).map((k, i) => ({
      type: "scatter", mode: "lines", name: k.replace("q_", "Q_"),
      x: sw.k0, y: sw[k], line: lineStyle(i),
    }));
    GT.state.imports.forEach((imp, j) => {
      const x = normalizeRange(imp.x, sw.k0);
      traces.push({ type: "scatter", mode: "lines", name: imp.name + " (имп.)",
        x: x, y: imp.y, line: Object.assign(lineStyle(keys.length + j), { width: 1.6 }) });
    });
    Plotly.react(host, traces, baseLayout({
      xaxis: { title: "k₀", gridcolor: PAL().grid },
      yaxis: { title: "эффективность Q", gridcolor: PAL().grid },
    }), config("greentensor_sweep"));
  }

  function normalizeRange(x, ref) {
    const lo = Math.min.apply(null, ref), hi = Math.max.apply(null, ref);
    const fx = (x || []).filter(isFinite);
    if (!fx.length) return ref;
    const mn = Math.min.apply(null, fx), mx = Math.max.apply(null, fx);
    if (mx <= mn) return ref;
    if (mn >= lo - 1e-9 && mx <= hi + 1e-9) return x;        // уже в диапазоне
    return x.map((v) => lo + ((v - mn) / (mx - mn)) * (hi - lo)); // нормировка размерности
  }

  // ===================== общий вывод ===================== //
  function showEmpty(on) {
    GT.$("#plot-empty").classList.toggle("hidden", !on);
    if (on) Plotly.purge(GT.$("#plot"));
  }

  function syncToolbar() {
    const v = GT.state.resultView;
    const polarV = v === "pattern" || v === "bistatic" || v === "monostatic";
    GT.$("#polar-rot-wrap").style.display = polarV ? "" : "none";
    GT.$("#db-toggle-wrap").style.display = v === "pattern" ? "" : "none";
    GT.$("#pattern-layout").style.display = v === "pattern" ? "" : "none";
    const hmV = v === "heatmap";
    GT.$("#hm-scale").style.display = hmV && (GT.state.heatMode || "abs") !== "re" ? "" : "none";
    GT.$("#hm-res-wrap").style.display = hmV ? "" : "none";
    ["#hm-mode", "#hm-comp", "#hm-plane", "#btn-gif"].forEach((id) => {
      const el = GT.$(id); if (el) el.style.display = hmV ? "" : "none";
    });
    // кнопка анимации бегущей волны только в тепловой карте
    let anim = GT.$("#anim-btn");
    if (hmV) {
      if (!anim) {
        anim = document.createElement("button");
        anim.id = "anim-btn"; anim.className = "btn"; anim.type = "button";
        anim.title = "Анимация бегущей волны: кадры Re[E·e^{jωt}] по ωt = 0…2π (шкала фиксирована)";
        anim.textContent = "▶ Бегущая волна";
        anim.addEventListener("click", toggleAnim);
        const gif = GT.$("#btn-gif");
        gif.parentNode.insertBefore(anim, gif);
      }
      anim.style.display = "";
    } else if (anim) { anim.style.display = "none"; }
    GT.$("#imports-bar").hidden = !(v === "pattern" || v === "sweep") || !GT.state.imports.length;
  }

  function render() {
    stopAnim();
    setAnimBtn(false);
    syncToolbar();
    // Вкладка спектра меняет подпись, если результат несёт генерик-кривую (η₁(θ₀) и т. п.)
    const res0 = GT.state.results;
    const swBtn = GT.$('[data-rv="sweep"]');
    if (swBtn) swBtn.textContent =
      (res0 && res0.sweep && res0.sweep.generic && res0.sweep.tab_title)
        ? res0.sweep.tab_title : "Спектр Q(k)";
    const v = GT.state.resultView;
    if (v === "pattern") renderPattern();
    else if (v === "bistatic") renderBistatic();
    else if (v === "monostatic") renderMonostatic();
    else if (v === "heatmap") renderHeatmap();
    else if (v === "cross") renderCross();
    else if (v === "sweep") renderSweep();
    renderImportsList();
  }

  // ===================== CSV импорт/экспорт ===================== //
  function parseCSV(text) {
    const rows = text.split(/\r?\n/).map((l) => l.trim()).filter((l) => l.length);
    const data = [];
    let header = null;
    rows.forEach((line, idx) => {
      const parts = line.split(/[,;\t ]+/);
      const nums = parts.map(Number);
      if (nums.some((n) => !isFinite(n))) {           // строка-заголовок
        if (idx === 0) header = parts;
        return;
      }
      data.push(nums);
    });
    if (!data.length) return null;
    const ncol = data[0].length;
    if (ncol >= 2) {
      return { x: data.map((r) => r[0]), y: data.map((r) => r[1]),
        name: (header && header[1]) || "" };
    }
    return { x: data.map((_, i) => i), y: data.map((r) => r[0]),
      name: (header && header[0]) || "" };
  }

  function importCSV(file) {
    const reader = new FileReader();
    reader.onload = () => {
      const ds = parseCSV(String(reader.result));
      if (!ds) { GT.toast("CSV: не найдено числовых данных", "error"); return; }
      const id = "imp" + Date.now();
      const name = ds.name || file.name.replace(/\.csv$/i, "");
      // авто-привязка панели по имени: …E… → панель 0 (E-плоскость), …H… → панель 1
      let panel = "all";
      if (/(^|[^a-zа-я])(e|е)([^a-zа-я]|$|-?plane|-?пл)/i.test(name)) panel = 0;
      else if (/(^|[^a-zа-я])(h|н)([^a-zа-я]|$|-?plane|-?пл)/i.test(name)) panel = 1;
      GT.state.imports.push({ id, name, x: ds.x, y: ds.y, panel });
      GT.toast("Импортировано: " + file.name + " (" + ds.y.length + " точек)", "ok");
      render();
    };
    reader.onerror = () => GT.toast("Не удалось прочитать файл", "error");
    reader.readAsText(file);
  }

  function renderImportsList() {
    const ul = GT.$("#imports-list");
    ul.innerHTML = "";
    GT.state.imports.forEach((imp, i) => {
      const sw = document.createElement("span");
      sw.className = "swatch";
      sw.style.background = isBW() ? PAL().bwLine : PAL().cat[(i + 2) % PAL().cat.length];
      const inp = document.createElement("input");
      inp.type = "text"; inp.value = imp.name; inp.title = "Подпись набора (ось/легенда)";
      inp.addEventListener("input", () => { imp.name = inp.value; });
      inp.addEventListener("change", render);
      const li = document.createElement("li");
      li.appendChild(sw); li.appendChild(inp);
      // целевая панель (для раскладки ДН «По плоскостям»): «Все» или конкретная плоскость
      const res = GT.state.results;
      if (GT.state.resultView === "pattern" && res && res.pattern && res.pattern.series.length > 1) {
        const sel = document.createElement("select");
        sel.title = "Панель для сравнения (раскладка «По плоскостям»)";
        sel.appendChild(new Option("Все панели", "all"));
        res.pattern.series.forEach((s, k) => sel.appendChild(new Option("→ " + s.name, String(k))));
        sel.value = imp.panel == null ? "all" : String(imp.panel);
        sel.addEventListener("change", () => {
          imp.panel = sel.value === "all" ? "all" : Number(sel.value); render();
        });
        li.appendChild(sel);
      }
      const del = document.createElement("button");
      del.className = "imp-del"; del.textContent = "×"; del.title = "Удалить";
      del.addEventListener("click", () => {
        GT.state.imports = GT.state.imports.filter((x) => x.id !== imp.id); render();
      });
      li.appendChild(del);
      ul.appendChild(li);
    });
    GT.$("#imports-bar").hidden = !GT.state.imports.length ||
      !(GT.state.resultView === "pattern" || GT.state.resultView === "sweep");
  }

  function exportCSV() {
    const res = GT.state.results;
    const v = GT.state.resultView;
    let csv = "";
    if (v === "pattern" && res && res.pattern) {
      const useDB = GT.$("#db-toggle").checked;
      const p = res.pattern;
      csv = "theta_deg," + p.series.map((s) => '"' + s.name + '"').join(",") + "\n";
      p.theta_deg.forEach((t, r) => {
        csv += t + "," + p.series.map((s) => (useDB ? s.dB : s.E)[r]).join(",") + "\n";
      });
    } else if (v === "sweep" && res && res.sweep) {
      const sw = res.sweep;
      const ks = ["q_sca", "q_ext", "q_abs", "q_back"].filter((k) => sw[k]);
      csv = "k0," + ks.join(",") + "\n";
      sw.k0.forEach((k, r) => { csv += k + "," + ks.map((kk) => sw[kk][r]).join(",") + "\n"; });
    } else if (v === "bistatic" && res && res.bistatic && res.bistatic.series) {
      const bi = res.bistatic;
      csv = "angle_deg," + bi.series.map((s) => '"' + s.name + ' [σ/λ² дБ]"').join(",") + "\n";
      bi.angle_deg.forEach((a, r) => {
        csv += a + "," + bi.series.map((s) => s.dB[r]).join(",") + "\n";
      });
    } else if (v === "monostatic" && res && res.monostatic && res.monostatic.sigma) {
      const mo = res.monostatic;
      csv = "aspect_deg,sigma_over_lambda2\n";
      mo.aspect_deg.forEach((a, r) => { csv += a + "," + mo.sigma[r] + "\n"; });
    } else if (v === "cross" && res && res.cross_sections) {
      csv = "quantity,value\n";
      Object.entries(res.cross_sections).forEach(([k, val]) => { csv += k + "," + val + "\n"; });
    } else if (v === "heatmap" && res && res.heatmap) {
      const hm = res.heatmap;
      const comp = effComp(hm);
      const arr = comp === "abs" ? null : compArrays(hm, comp);
      const hl = hm.xlabel || "x", vl = hm.ylabel || "z";
      if (GT.state.heatMode === "re" && arr) {
        csv = "# Re/Im E_" + comp + " (engine e^{-iwt}; E(t) = re*cos(wt) + im*sin(wt)); plane " +
          (hm.plane || "xz") + "; first row = " + hl + ", first col = " + vl + "; blocks: re, im\n";
        ["re", "im"].forEach((part) => {
          csv += "# " + part + "\n," + hm.x.join(",") + "\n";
          arr[part].forEach((row, r) => { csv += hm.y[r] + "," + row.join(",") + "\n"; });
        });
      } else {
        const absZ = (comp === "abs" || !arr) ? hm.absE
          : arr.re.map((row, r) => row.map((val, k) => (val == null ? null : Math.hypot(val, arr.im[r][k] || 0))));
        csv = "# " + (comp === "abs" ? "|E|" : "|E_" + comp + "|") + " heatmap; plane " + (hm.plane || "xz") +
          "; first row = " + hl + ", first col = " + vl + "\n,";
        csv += hm.x.join(",") + "\n";
        absZ.forEach((row, r) => { csv += hm.y[r] + "," + row.join(",") + "\n"; });
      }
    } else { GT.toast("Нет данных для экспорта", "error"); return; }
    GT.download("greentensor_" + v + ".csv", csv, "text/csv;charset=utf-8");
  }

  function savePNG() {
    Plotly.downloadImage(GT.$("#plot"),
      { format: "png", scale: 2, filename: "greentensor_" + GT.state.resultView });
  }

  function onShow() {
    const host = GT.$("#plot");
    // resize отложен внутри Plotly на ~100 мс; если за это время график очищен
    // (purge при смене представления/пресета) — промис отклоняется, это не ошибка
    if (host && host.data) Plotly.Plots.resize(host).catch(() => {});
    else render();
  }

  function init() {
    if (!GT.state.patternLayout) GT.state.patternLayout = "overlay";
    if (!GT.state.heatScale) GT.state.heatScale = "db";
    if (!GT.state.heatMode) GT.state.heatMode = "abs";       // abs | re  (Re± — бегущая волна)
    if (!GT.state.heatComp) GT.state.heatComp = "abs";       // abs | x | y | z
    if (!GT.state.heatPlane) GT.state.heatPlane = "xz";      // xz | yz | xy
    GT.$$("#hm-mode button").forEach((b) =>
      b.addEventListener("click", () => {
        stopAnim(); setAnimBtn(false);
        setHeatMode(b.dataset.hm);
        if (GT.state.heatMode === "re" && (GT.state.heatComp || "abs") === "abs") {
          const hm = GT.state.results && GT.state.results.heatmap;
          GT.state.heatComp = (hm && hm.dominant) || "x";
        }
        syncToolbar();
        if (GT.state.resultView === "heatmap") renderHeatmap();
      }));
    GT.$$("#hm-comp button").forEach((b) =>
      b.addEventListener("click", () => {
        if (b.disabled) return;
        stopAnim(); setAnimBtn(false);
        GT.state.heatComp = b.dataset.hc;
        if (GT.state.resultView === "heatmap") renderHeatmap();
      }));
    GT.$$("#hm-plane button").forEach((b) =>
      b.addEventListener("click", async () => {
        if (GT.state.heatPlane === b.dataset.hp) return;
        GT.state.heatPlane = b.dataset.hp;
        GT.$$("#hm-plane button").forEach((x) => x.classList.toggle("active", x === b));
        if (!GT.state.results) return;              // применится при следующем расчёте
        stopAnim(); setAnimBtn(false);
        await GT.compute({ pattern: false, heatmap: true, sweep: false, merge: true });
      }));
    const gifBtn = GT.$("#btn-gif");
    if (gifBtn) gifBtn.addEventListener("click", exportGIF);
    GT.$$("#hm-scale button").forEach((b) =>
      b.addEventListener("click", () => {
        GT.state.heatScale = b.dataset.hs;
        GT.$$("#hm-scale button").forEach((x) => x.classList.toggle("active", x === b));
        if (GT.state.resultView === "heatmap") renderHeatmap();
      }));
    const resSel = GT.$("#hm-res");
    if (resSel) {
      resSel.value = String(GT.state.heatRes);
      resSel.addEventListener("change", async () => {
        GT.state.heatRes = parseInt(resSel.value, 10) || 64;
        if (!GT.state.results) return;              // нет расчёта — применится при следующем
        stopAnim();                                 // старые re/imField станут несовместимы по размеру
        await GT.compute({ pattern: false, heatmap: true, sweep: false, merge: true });
      });                                           // GT.compute сам перерисует активное представление
    }
    GT.$$("#result-view button").forEach((b) =>
      b.addEventListener("click", () => {
        GT.state.resultView = b.dataset.rv;
        GT.$$("#result-view button").forEach((x) => x.classList.toggle("active", x === b));
        render();
      }));
    GT.$$("#pattern-layout button").forEach((b) =>
      b.addEventListener("click", () => {
        GT.state.patternLayout = b.dataset.pl;
        GT.$$("#pattern-layout button").forEach((x) => x.classList.toggle("active", x === b));
        if (GT.state.resultView === "pattern") renderPattern();
      }));
    GT.$("#polar-rot").addEventListener("input", () => {
      GT.$("#polar-rot-val").textContent = GT.$("#polar-rot").value;
      const v = GT.state.resultView;
      if (v === "pattern") renderPattern();
      else if (v === "bistatic") renderBistatic();
      else if (v === "monostatic") renderMonostatic();
    });
    GT.$("#db-toggle").addEventListener("change", renderPattern);
    GT.$("#btn-png").addEventListener("click", savePNG);
    GT.$("#btn-export").addEventListener("click", exportCSV);
    GT.$("#btn-import").addEventListener("click", () => GT.$("#csv-file").click());
    GT.$("#csv-file").addEventListener("change", (e) => {
      const f = e.target.files[0]; if (f) importCSV(f); e.target.value = "";
    });
  }

  // Установка режима карты извне (пресеты): {heatMode, heatComp, heatPlane, resultView}
  function setHeatOptions(o) {
    o = o || {};
    if (o.heatMode) GT.state.heatMode = o.heatMode;
    if (o.heatComp) GT.state.heatComp = o.heatComp;
    if (o.heatPlane) GT.state.heatPlane = o.heatPlane;
    if (o.resultView) {
      GT.state.resultView = o.resultView;
      GT.$$("#result-view button").forEach((x) => x.classList.toggle("active", x.dataset.rv === o.resultView));
    }
    syncHeatControls(GT.state.results && GT.state.results.heatmap);
  }

  GT.results = { init, render, onShow, setHeatOptions, startAnim, stopAnim, exportGIF };
})();
