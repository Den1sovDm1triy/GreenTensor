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

  // Кривая Ansys HFSS (МКЭ), в дБ, по плоскости: E — красная пунктирная, H — оранжевая точечная.
  function refTrace(ref, subplot) {
    const col = isBW() ? PAL().bwLine : (ref.plane === "H" ? "#ff7f0e" : "#d62728");
    const t = { type: "scatterpolar", mode: "lines", name: ref.name,
      theta: ref.theta_deg, r: ref.dB,
      line: { color: col, width: 2, dash: ref.plane === "H" ? "dot" : "dash" } };
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

  // ===================== Тепловая карта ===================== //
  const DB_HM_FLOOR = -30;            // нижний предел лог-шкалы карты, дБ отн. максимума
  // спокойная дивергентная палитра «синий → красный» (типа coolwarm), без ярких пиков;
  // ч-б — оттенки серого (низк.→свет., выс.→тёмн.)
  function heatScale() {
    return isBW()
      ? [[0, "#eef1f4"], [1, "#16202a"]]
      : [[0.0, "#3b5bbf"], [0.18, "#6a86d8"], [0.36, "#9fb2e6"],
         [0.50, "#c4c8d8"], [0.64, "#dcb39f"], [0.82, "#cc7155"],
         [1.0, "#a52a2a"]];   // синий → приглушённый красный
  }
  function fieldMax(absE) {
    let mx = 0;
    absE.forEach((row) => row.forEach((v) => { if (v != null && v > mx) mx = v; }));
    return mx > 0 ? mx : 1;
  }
  function renderHeatmap() {
    const res = GT.state.results;
    const host = GT.$("#plot");
    if (res && !res.heatmap) {                       // расчёт есть, но поля нет (напр. беск. цилиндр)
      showEmpty(true);
      GT.$("#plot-empty").innerHTML =
        "Тепловая карта недоступна для этого тела: решатель не даёт ближнее поле " +
        "(2D-цилиндр). Используйте «Диаграмма (ДН)», «Сечения» или «Спектр Q(k)».";
      return;
    }
    if (!res || !res.heatmap) return showEmpty(true);
    showEmpty(false);
    const hm = res.heatmap;
    const mx = fieldMax(hm.absE);                 // нормировка на максимум поля
    const log = GT.state.heatScale !== "lin";
    let z, title, zmin, zmax;
    if (log) {                                    // лог (дБ отн. макс): и пики, и переходы слоёв
      z = hm.absE.map((r) => r.map((v) => (v == null ? null : Math.max(DB_HM_FLOOR, 20 * Math.log10(v / mx)))));
      title = "дБ отн. макс"; zmin = DB_HM_FLOOR; zmax = 0;
    } else {                                      // линейная, нормированная
      z = hm.absE.map((r) => r.map((v) => (v == null ? null : v / mx)));
      title = "|E| / |E|ₘₐₖₛ"; zmin = 0; zmax = 1;
    }
    const trace = {
      type: "heatmap", x: hm.x, y: hm.y, z, colorscale: heatScale(),
      zmin, zmax, colorbar: { title, outlinewidth: 0 }, zsmooth: "best",
    };
    const sub = (log ? "лог (дБ отн. макс)" : "линейная, норм.") +
      (hm.res ? " · " + hm.res + "×" + hm.res : "") +
      (hm.metal_masked ? " · металл-области скрыты" : "");
    const layout = baseLayout({
      showlegend: false,
      xaxis: { title: "x", scaleanchor: "y", constrain: "domain", zeroline: false, gridcolor: PAL().grid },
      yaxis: { title: "z", zeroline: false, gridcolor: PAL().grid },
      annotations: [{ x: 0.5, y: 1.06, xref: "paper", yref: "paper", showarrow: false,
        text: (hm.quantity === "scattered" ? "Рассеянное поле |E_расс| в плоскости xz · "
          : "Ближнее поле |E| в плоскости xz · ") + sub, font: { color: PAL().muted, size: 12 } }],
    });
    Plotly.react(host, [trace], layout, config("greentensor_heatmap"));
  }

  // Анимация НАПРЯЖЁННОСТИ: мгновенное |Re[E·e^{-iωt}]| = |reField·cos(ωt)+imField·sin(ωt)|
  function toggleAnim() {
    const res = GT.state.results;
    if (!res || !res.heatmap || !res.heatmap.reField) return;
    const host = GT.$("#plot");
    const re = res.heatmap.reField, im = res.heatmap.imField;
    const btn = GT.$("#anim-btn");
    if (animTimer) { stopAnim(); if (btn) btn.textContent = "▶ Напряжённость (аним.)"; renderHeatmap(); return; }
    if (btn) btn.textContent = "⏸ Стоп";
    // амплитуда (для фиксированной шкалы цвета)
    let amax = 0;
    re.forEach((row, r) => row.forEach((v, k) => {
      if (v == null) return;
      const a = Math.hypot(v, im[r][k] || 0); if (a > amax) amax = a;
    }));
    const frame = (c, s) => re.map((row, r) => row.map((v, k) =>
      (v == null ? null : Math.abs(v * c + (im[r][k] || 0) * s))));
    Plotly.react(host, [{ type: "heatmap", x: res.heatmap.x, y: res.heatmap.y,
      z: frame(1, 0), colorscale: heatScale(), zmin: 0, zmax: amax || 1,
      colorbar: { title: "|E(t)|" }, zsmooth: "best" }],
      baseLayout({ showlegend: false,
        xaxis: { scaleanchor: "y", constrain: "domain", title: "x" }, yaxis: { title: "z" },
        annotations: [{ x: 0.5, y: 1.06, xref: "paper", yref: "paper", showarrow: false,
          text: (res.heatmap.quantity === "scattered"
            ? "Мгновенная напряжённость рассеянного поля |E_расс(r,t)|"
            : "Мгновенная напряжённость поля |E(r,t)| (колебание во времени)"),
          font: { color: PAL().muted, size: 12 } }] }),
      config("greentensor_field"));
    let t = 0;
    animTimer = setInterval(() => {
      t += 0.18;
      Plotly.restyle(host, { z: [frame(Math.cos(t), Math.sin(t))] }, [0]);
    }, 60);
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
    GT.$("#hm-scale").style.display = v === "heatmap" ? "" : "none";
    GT.$("#hm-res-wrap").style.display = v === "heatmap" ? "" : "none";
    // кнопка анимации только в тепловой карте
    let anim = GT.$("#anim-btn");
    if (v === "heatmap") {
      if (!anim) {
        anim = document.createElement("button");
        anim.id = "anim-btn"; anim.className = "btn"; anim.type = "button";
        anim.textContent = "▶ Напряжённость (аним.)";
        anim.addEventListener("click", toggleAnim);
        GT.$("#btn-png").parentNode.insertBefore(anim, GT.$("#btn-png"));
      }
      anim.style.display = "";
    } else if (anim) { anim.style.display = "none"; }
    GT.$("#imports-bar").hidden = !(v === "pattern" || v === "sweep") || !GT.state.imports.length;
  }

  function render() {
    stopAnim();
    const a = GT.$("#anim-btn"); if (a) a.textContent = "▶ Напряжённость (аним.)";
    syncToolbar();
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
      csv = "# |E| heatmap; first row = x, first col = y\n,";
      csv += hm.x.join(",") + "\n";
      hm.absE.forEach((row, r) => { csv += hm.y[r] + "," + row.join(",") + "\n"; });
    } else { GT.toast("Нет данных для экспорта", "error"); return; }
    GT.download("greentensor_" + v + ".csv", csv, "text/csv;charset=utf-8");
  }

  function savePNG() {
    Plotly.downloadImage(GT.$("#plot"),
      { format: "png", scale: 2, filename: "greentensor_" + GT.state.resultView });
  }

  function onShow() {
    const host = GT.$("#plot");
    if (host && host.data) Plotly.Plots.resize(host);
    else render();
  }

  function init() {
    if (!GT.state.patternLayout) GT.state.patternLayout = "overlay";
    if (!GT.state.heatScale) GT.state.heatScale = "db";
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

  GT.results = { init, render, onShow };
})();
