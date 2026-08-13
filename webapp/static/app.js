/* SPDX-License-Identifier: MIT
 */

/* GreenTensor web — фронтенд.
   Все диаграммы строятся через ЕДИНЫЙ шаблон (шрифты, сетка, палитра, экспорт),
   цвета читаются из CSS-токенов: индиго = θ/co, янтарь = φ/cross. */
"use strict";

const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const C = {
  ink: css("--ink"), muted: css("--muted"), line: css("--line"),
  a1: css("--accent"), a2: css("--accent-2"), teal: css("--teal"), rose: css("--rose"),
  paper: css("--paper"), surface: css("--surface"),
};
const FONT = { family: css("--font") || "system-ui", size: 12.5, color: C.ink };
const DB_RANGE = [-60, 2];

/* -------------------------------------------------- единый шаблон Plotly */
function baseLayout(title) {
  return {
    title: { text: title, font: { ...FONT, size: 13.5 }, x: 0, xanchor: "left", pad: { l: 2 } },
    font: FONT,
    paper_bgcolor: "rgba(0,0,0,0)",
    plot_bgcolor: "rgba(0,0,0,0)",
    margin: { l: 56, r: 18, t: 44, b: 46 },
    legend: { orientation: "h", y: -0.16, font: { ...FONT, size: 12 } },
    hovermode: "closest",
  };
}

function rectAxes(layout, xTitle, yTitle, yRange) {
  layout.xaxis = {
    title: { text: xTitle, font: { ...FONT, size: 12 } },
    gridcolor: C.line, zeroline: false, linecolor: C.line, ticks: "outside",
    tickcolor: C.line, dtick: 30, range: [0, 360],
  };
  layout.yaxis = {
    title: { text: yTitle, font: { ...FONT, size: 12 } },
    gridcolor: C.line, zeroline: false, linecolor: C.line, ticks: "outside",
    tickcolor: C.line, range: yRange,
  };
  return layout;
}

function polarLayout(title) {
  const L = baseLayout(title);
  L.margin = { l: 36, r: 36, t: 48, b: 40 };
  L.polar = {
    bgcolor: "rgba(0,0,0,0)",
    radialaxis: {
      range: DB_RANGE, dtick: 10, angle: 90, tickangle: 90,
      gridcolor: C.line, linecolor: C.line,
      tickfont: { ...FONT, size: 10.5, color: C.muted }, ticksuffix: " дБ",
    },
    angularaxis: {
      rotation: 90, direction: "clockwise", dtick: 30,
      gridcolor: C.line, linecolor: C.line, tickfont: { ...FONT, size: 11 },
    },
  };
  return L;
}

function plotConfig(filename) {
  return {
    displaylogo: false,
    responsive: true,
    modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d", "zoomIn2d", "zoomOut2d"],
    toImageButtonOptions: { format: "png", scale: 2, filename },
  };
}

const line = (w = 2.2) => ({ width: w });

/* -------------------------------------------------- состояние и пресеты */
const state = { layers: [], data: null, problem: "diffraction" };

const PRESETS = {
  "Люнеберг · 7 слоёв": {
    k0: 10, layers: [[0.15, 3.90], [0.30, 3.52], [0.45, 2.98], [0.63, 2.48], [0.82, 1.76], [1.0, 1.50]],
  },
  "Люнеберг · 4 слоя": {
    k0: 10, layers: [[0.25, 1.96], [0.5, 1.84], [0.75, 1.64], [1.0, 1.36]],
  },
  "Металлическая сфера": {
    k0: 3, layers: [[1.0, 1.0, 1e7]],
  },
  "Металл с укрытием": {
    k0: 5, layers: [[0.5, 1.0, 1e7], [1.0, 2.25]],
  },
  "Однородный диэлектрик": {
    k0: 5, layers: [[1.0, 2.25]],
  },
};

function applyPreset(name) {
  const p = PRESETS[name];
  document.getElementById("k0").value = p.k0;
  state.layers = p.layers.map(([r, er, ei]) => ({ r, eps_re: er, eps_im: ei || 0, mu_re: 1 }));
  renderLayers();
}

function renderLayers() {
  const box = document.getElementById("layers");
  box.innerHTML = "";
  state.layers.forEach((L, i) => {
    const row = document.createElement("div");
    row.className = "layer-row";
    row.innerHTML = `
      <input type="number" step="any" value="${L.r}"      data-i="${i}" data-k="r">
      <input type="number" step="any" value="${L.eps_re}" data-i="${i}" data-k="eps_re">
      <input type="number" step="any" value="${L.eps_im}" data-i="${i}" data-k="eps_im">
      <input type="number" step="any" value="${L.mu_re}"  data-i="${i}" data-k="mu_re">
      <button class="del" title="Удалить слой" type="button">×</button>`;
    row.querySelectorAll("input").forEach((inp) =>
      inp.addEventListener("change", (e) => {
        state.layers[+e.target.dataset.i][e.target.dataset.k] = parseFloat(e.target.value) || 0;
      }));
    row.querySelector(".del").addEventListener("click", () => {
      state.layers.splice(i, 1); renderLayers();
    });
    box.appendChild(row);
  });
}

/* -------------------------------------------------- расчёт */
async function run() {
  const btn = document.getElementById("run");
  const status = document.getElementById("status");
  btn.disabled = true;
  status.className = "status"; status.textContent = "Расчёт…";

  const problems = ["diffraction"];
  if (document.getElementById("problemAntenna").checked) problems.push("antenna");

  const tochRaw = document.getElementById("toch").value;
  const payload = {
    k0: parseFloat(document.getElementById("k0").value),
    toch: tochRaw ? parseInt(tochRaw, 10) : null,
    phi_deg: parseFloat(document.getElementById("phi").value),
    layers: state.layers,
    problems,
    sweep: {
      enable: document.getElementById("sweepOn").checked,
      k0_min: parseFloat(document.getElementById("sweepMin").value),
      k0_max: parseFloat(document.getElementById("sweepMax").value),
      points: parseInt(document.getElementById("sweepN").value, 10),
    },
  };

  try {
    const r = await fetch("/api/compute", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const data = await r.json();
    if (!r.ok) throw new Error(data.error || `HTTP ${r.status}`);
    state.data = data;
    state.problem = "diffraction";
    status.textContent = `Готово · toch = ${data.meta.toch_used}` +
      (data.meta.warnings.length ? ` · ⚠ ${data.meta.warnings.join(" ")}` : "");
    renderAll();
  } catch (e) {
    status.className = "status err";
    status.textContent = e.message;
  } finally {
    btn.disabled = false;
  }
}

/* -------------------------------------------------- отрисовка */
function renderAll() {
  const d = state.data;
  document.getElementById("results").hidden = false;

  renderGeometry(d.meta.layers_norm);
  renderCards(d.results.diffraction.cross_sections, d.meta);

  const hasAntenna = !!d.results.antenna;
  const tabs = document.getElementById("tabs");
  tabs.hidden = !hasAntenna;
  tabs.querySelectorAll(".tab").forEach((t) => {
    t.classList.toggle("active", t.dataset.problem === state.problem);
  });

  renderProblem(state.problem);
  renderSweep(d.sweep);
  renderCoeffs(d.results[state.problem].coeffs);
}

function renderProblem(problem) {
  const d = state.data;
  const th = d.theta_deg;
  const R = d.results[problem];
  const tag = problem === "antenna" ? " · антенна" : "";

  const lin = R.linear;
  const phiC = lin.custom.phi_deg;
  // в главных плоскостях (φ кратно 90°) E_φ ≡ 0 в оригинальной формуле — срез не рисуем
  const isPrincipal = Math.min(((phiC % 90) + 90) % 90, 90 - (((phiC % 90) + 90) % 90)) < 0.5;
  const tracesLin = [
    { theta: th, r: lin.E_plane.dB, name: "E-плоскость · φ=0°", line: { color: C.a1, ...line() } },
    { theta: th, r: lin.H_plane.dB, name: "H-плоскость · φ=90°", line: { color: C.a2, ...line() } },
  ];
  if (!isPrincipal) {
    tracesLin.push(
      { theta: th, r: lin.custom.dB_theta, name: `E_θ · φ=${phiC}°`,
        line: { color: C.teal, width: 1.7, dash: "dash" } },
      { theta: th, r: lin.custom.dB_phi, name: `E_φ · φ=${phiC}°`,
        line: { color: C.rose, width: 1.7, dash: "dash" } },
    );
  }
  Plotly.react("plotLinPolar",
    tracesLin.map((t) => ({ ...t, type: "scatterpolar", mode: "lines" })),
    polarLayout(`Полярная ДН${tag}`), plotConfig(`lin_polar_${problem}`));

  Plotly.react("plotLinRect",
    tracesLin.map((t) => ({ x: t.theta, y: t.r, name: t.name, mode: "lines", line: t.line })),
    rectAxes(baseLayout(`Развёртка по θ${tag}`), "θ, °", "уровень, дБ", DB_RANGE),
    plotConfig(`lin_rect_${problem}`));

  const tracesCir = [
    { theta: th, r: R.circular.dB_op, name: "Основная (co)", line: { color: C.a1, ...line() } },
    { theta: th, r: R.circular.dB_kp, name: "Кросс (cross)", line: { color: C.a2, ...line() } },
  ];
  Plotly.react("plotCirPolar",
    tracesCir.map((t) => ({ ...t, type: "scatterpolar", mode: "lines" })),
    polarLayout(`Полярная ДН${tag}`), plotConfig(`cir_polar_${problem}`));

  Plotly.react("plotCirRect",
    tracesCir.map((t) => ({ x: t.theta, y: t.r, name: t.name, mode: "lines", line: t.line })),
    rectAxes(baseLayout(`Развёртка по θ${tag}`), "θ, °", "уровень, дБ", DB_RANGE),
    plotConfig(`cir_rect_${problem}`));

  renderCoeffs(R.coeffs);
}

function renderSweep(sweep) {
  const sec = document.getElementById("sweepSection");
  if (!sweep) { sec.hidden = true; return; }
  sec.hidden = false;
  const mk = (key, name, color, dash) => ({
    x: sweep.k0, y: sweep[key], name, mode: "lines",
    line: { color, width: 2.2, dash: dash || "solid" }, connectgaps: false,
  });
  const k0 = state.data.meta.k0;
  const L = rectAxes(baseLayout("Q(k₀): полное, экстинкция, поглощение, обратное"),
    "k₀·R", "Q (норм. πR²)", null);
  L.xaxis.range = null; L.xaxis.dtick = null;
  L.shapes = [{
    type: "line", x0: k0, x1: k0, yref: "paper", y0: 0, y1: 1,
    line: { color: C.muted, width: 1, dash: "dot" },
  }];
  L.annotations = [{
    x: k0, yref: "paper", y: 1.02, text: `k₀=${k0}`, showarrow: false,
    font: { ...FONT, size: 11, color: C.muted },
  }];
  Plotly.react("plotSweep", [
    mk("q_sca", "Q_sca", C.a1),
    mk("q_ext", "Q_ext", C.a2),
    mk("q_abs", "Q_abs", C.teal),
    mk("q_back", "Q_back", C.rose, "dash"),
  ], L, plotConfig("cross_sections_sweep"));
}

function renderCoeffs(cf) {
  const L = baseLayout("Модули коэффициентов рассеяния (лог. шкала)");
  rectAxes(L, "номер гармоники n", "|M_n|, |N_n|", null);
  L.xaxis.range = [0.4, cf.n.length + 0.6]; L.xaxis.dtick = cf.n.length > 30 ? 5 : 1;
  L.yaxis.type = "log"; L.yaxis.range = null;
  Plotly.react("plotCoeffs", [
    { x: cf.n, y: cf.Mn_abs, name: "|M_n|", type: "bar", marker: { color: C.a1 }, offsetgroup: 1 },
    { x: cf.n, y: cf.Nn_abs, name: "|N_n|", type: "bar", marker: { color: C.a2 }, offsetgroup: 2 },
  ], L, plotConfig("coefficients"));
}

function renderCards(cs, meta) {
  const fmt = (v) => (v == null ? "—" : v.toPrecision(4));
  const items = [
    ["Q_sca", cs.q_sca, "sca", "полное рассеяние"],
    ["Q_ext", cs.q_ext, "ext", "экстинкция"],
    ["Q_abs", cs.q_abs, "abs", "поглощение"],
    ["Q_back", cs.q_back, "back", "обратное"],
  ];
  document.getElementById("csCards").innerHTML = items.map(([k, v, c, u]) => `
    <div class="card stat" data-c="${c}">
      <span class="k">${k}</span>
      <span class="v">${fmt(v)}</span>
      <span class="u">${u} · k₀=${meta.k0} · toch=${meta.toch_used}</span>
    </div>`).join("");
}

function renderGeometry(ln) {
  // концентрические слои; цвет — по ε′ (индиго-шкала), металл — графит
  const svg = document.getElementById("geom");
  const n = ln.a.length - 1; // без внешнего воздуха
  const epsVals = ln.eps_re.slice(0, n);
  const finite = epsVals.filter((e, i) => Math.hypot(e, ln.eps_im[i]) < 1e3);
  const lo = Math.min(1, ...finite), hi = Math.max(1.0001, ...finite);
  const colorOf = (i) => {
    if (Math.hypot(ln.eps_re[i], ln.eps_im[i]) >= 1e3) return "#3c4049";
    const t = (ln.eps_re[i] - lo) / (hi - lo);
    const l = 92 - t * 50;
    return `hsl(232 55% ${l}%)`;
  };
  let out = "";
  for (let i = n - 1; i >= 0; i--) {
    out += `<circle cx="0" cy="0" r="${ln.a[i]}" fill="${colorOf(i)}" stroke="${C.surface}" stroke-width="0.012"/>`;
  }
  out += `<circle cx="0" cy="0" r="1" fill="none" stroke="${C.ink}" stroke-width="0.014"/>`;
  // подписи ε снаружи
  for (let i = 0; i < n; i++) {
    const rMid = i === 0 ? ln.a[0] / 2 : (ln.a[i - 1] + ln.a[i]) / 2;
    const big = Math.hypot(ln.eps_re[i], ln.eps_im[i]) >= 1e3;
    const txt = big ? "металл" : `ε=${+ln.eps_re[i].toFixed(3)}${ln.eps_im[i] ? `+${+ln.eps_im[i].toFixed(3)}i` : ""}`;
    const y = -rMid;
    out += `<line x1="0" y1="${y}" x2="${1.08}" y2="${y}" stroke="${C.muted}" stroke-width="0.006" stroke-dasharray="0.02 0.02"/>`;
    out += `<text x="1.1" y="${y}" font-size="0.085" font-family="${css("--mono")}" fill="${C.ink}" dominant-baseline="middle">${txt}</text>`;
  }
  svg.innerHTML = out;
}

/* -------------------------------------------------- экспорт */
function download(name, text, type) {
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([text], { type }));
  a.download = name;
  a.click();
  URL.revokeObjectURL(a.href);
}

function exportCsv() {
  const d = state.data;
  if (!d) return;
  const cols = ["theta_deg"];
  const series = [d.theta_deg];
  const add = (name, arr) => { if (arr) { cols.push(name); series.push(arr); } };
  for (const [p, R] of Object.entries(d.results)) {
    add(`${p}.linear.E_plane.E`, R.linear.E_plane.E);
    add(`${p}.linear.E_plane.dB`, R.linear.E_plane.dB);
    add(`${p}.linear.H_plane.E`, R.linear.H_plane.E);
    add(`${p}.linear.H_plane.dB`, R.linear.H_plane.dB);
    const phiC = R.linear.custom.phi_deg;
    add(`${p}.linear.phi${phiC}.E_theta`, R.linear.custom.E_theta);
    add(`${p}.linear.phi${phiC}.E_phi`, R.linear.custom.E_phi);
    for (const k of ["E_op", "E_kp", "dB_op", "dB_kp"]) add(`${p}.circular.${k}`, R.circular[k]);
  }
  const rows = [cols.join(",")];
  for (let i = 0; i < d.theta_deg.length; i++) {
    rows.push(series.map((s) => (s[i] == null ? "" : s[i])).join(","));
  }
  download("greentensor_patterns.csv", rows.join("\n"), "text/csv");
}

/* -------------------------------------------------- инициализация */
document.getElementById("presets").append(...Object.keys(PRESETS).map((name) => {
  const b = document.createElement("button");
  b.type = "button"; b.textContent = name;
  b.addEventListener("click", () => applyPreset(name));
  return b;
}));

document.getElementById("addLayer").addEventListener("click", () => {
  const last = state.layers[state.layers.length - 1];
  state.layers.push({ r: last ? +(last.r + 0.1).toFixed(3) : 0.5, eps_re: 2, eps_im: 0, mu_re: 1 });
  renderLayers();
});

document.getElementById("run").addEventListener("click", run);
document.getElementById("dlJson").addEventListener("click", () =>
  state.data && download("greentensor_data.json", JSON.stringify(state.data), "application/json"));
document.getElementById("dlCsv").addEventListener("click", exportCsv);

document.getElementById("tabs").addEventListener("click", (e) => {
  const t = e.target.closest(".tab");
  if (!t || !state.data) return;
  state.problem = t.dataset.problem;
  document.querySelectorAll(".tab").forEach((x) => x.classList.toggle("active", x === t));
  renderProblem(state.problem);
});

applyPreset("Люнеберг · 7 слоёв");
