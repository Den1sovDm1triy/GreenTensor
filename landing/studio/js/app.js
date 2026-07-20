// SPDX-License-Identifier: MIT
// GreenTensor Studio — ядро приложения: состояние сцены, оркестрация модулей.
(function () {
  "use strict";
  const GT = (window.GT = window.GT || {});

  // ----- Состояние ----- //
  function defaultScene() {
    return {
      bodies: [{
        type: "sphere", position: [0, 0, 0], radius: 1.0,
        layers: [{ a_rel: 1.0, eps_re: 2.25, eps_im: 0, mu_re: 1, mu_im: 0 }],
        euler: [0, 0, 0],
      }],
      radiation: {
        k: 3.0, polarization: "linear", problem: "diffraction",
        phi_deg: 0, toch: null, khat: [0, 0, 1], pol: [1, 0, 0],
      },
      sweep: { enable: true, k0_min: 0.25, k0_max: 8.0, points: 160 },
    };
  }

  GT.state = {
    scene: defaultScene(),
    selected: 0,           // индекс выбранного тела
    results: null,         // последний ответ /api/compute
    imports: [],           // импортированные наборы CSV для сравнения
    reference: null,       // кривые Ansys HFSS (МКЭ) пресета: список {plane,name,theta_deg,dB}
    colorMode: "color",    // 'color' | 'bw'
    resultView: "pattern", // pattern | heatmap | cross | sweep
    heatRes: 64,           // разрешение тепловой карты (точек на сторону)
  };

  // ----- Хелперы ----- //
  const $ = (GT.$ = (sel, root) => (root || document).querySelector(sel));
  const $$ = (GT.$$ = (sel, root) => Array.from((root || document).querySelectorAll(sel)));

  let toastTimer = null;
  GT.toast = function (msg, kind) {
    const el = $("#toast");
    el.textContent = msg;
    el.className = "toast show" + (kind ? " " + kind : "");
    clearTimeout(toastTimer);
    toastTimer = setTimeout(() => (el.className = "toast"), 3200);
  };

  GT.busy = function (on, text) {
    const el = $("#busy");
    if (text) $("#busy-text").textContent = text;
    el.classList.toggle("hidden", !on);
  };

  GT.download = function (filename, text, mime) {
    const blob = new Blob([text], { type: mime || "text/plain;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = filename; document.body.appendChild(a); a.click();
    document.body.removeChild(a); URL.revokeObjectURL(url);
  };

  GT.selectedBody = () => GT.state.scene.bodies[GT.state.selected] || null;

  // ----- Перерисовка зависимых представлений ----- //
  GT.refresh = function () {
    if (GT.input) GT.input.render();
    if (GT.editor) GT.editor.syncFromState();
  };

  // ----- Расчёт ----- //
  GT.compute = async function (opts) {
    const o = opts || {};
    const want = Object.assign(
      { pattern: true, heatmap: true, sweep: false, heatmap_res: GT.state.heatRes }, o);
    const payload = Object.assign({}, GT.state.scene, { compute: want });
    const busyText = o.monostatic ? "Моностатическая ЭПР (развёртка ракурса)…"
      : o.merge ? "Пересчёт карты…" : (want.sweep ? "Свип спектра…" : "Расчёт…");
    // Движок считается в браузере (Pyodide). Пока он грузится — покажем это в статусе.
    if (GT.py && !GT.py.isReady()) GT.busy(true, "Инициализация движка (Pyodide)…");
    else GT.busy(true, busyText);
    try {
      const data = await GT.py.compute(payload);
      if (data.error) throw new Error(data.error);
      // частичный пересчёт (свип / только карта) — догружаем в результат, не затирая остальное
      if (o.merge && GT.state.results) {
        if (data.heatmap) GT.state.results.heatmap = data.heatmap;
        if (data.pattern) GT.state.results.pattern = data.pattern;
        if (data.bistatic) GT.state.results.bistatic = data.bistatic;
        if (data.monostatic) GT.state.results.monostatic = data.monostatic;
        if (data.sweep) GT.state.results.sweep = data.sweep;
        GT.state.results.meta = data.meta || GT.state.results.meta;
      } else if (want.sweep && GT.state.results) {
        GT.state.results.sweep = data.sweep;
        GT.state.results.meta = data.meta;
      } else {
        GT.state.results = data;
      }
      GT.toast("Готово · ядро: " + (data.meta && data.meta.core || "—"), "ok");
      switchTab("results");
      if (GT.results) GT.results.render();
      if (GT.editor && GT.editor.syncOverlay) GT.editor.syncOverlay();  // ДН поверх 3D-модели
      return data;
    } catch (e) {
      GT.toast(e.message, "error");
      return null;
    } finally {
      GT.busy(false);
    }
  };

  // ----- Вкладки ----- //
  function switchTab(name) {
    $$("#main-tabs button").forEach((b) => b.classList.toggle("active", b.dataset.tab === name));
    $("#tab-editor").classList.toggle("hidden", name !== "editor");
    $("#tab-results").classList.toggle("hidden", name !== "results");
    if (name === "editor" && GT.editor) GT.editor.onShow();
    if (name === "results" && GT.results) GT.results.onShow();
  }
  GT.switchTab = switchTab;

  // ----- Инициализация ----- //
  function init() {
    // модули
    if (GT.editor) GT.editor.init();
    if (GT.input) GT.input.init();
    if (GT.results) GT.results.init();
    if (GT.presets) GT.presets.load();

    // вкладки
    $$("#main-tabs button").forEach((b) =>
      b.addEventListener("click", () => switchTab(b.dataset.tab)));

    // цветовой режим (глобальный)
    $$("#colormode button").forEach((b) =>
      b.addEventListener("click", () => {
        GT.state.colorMode = b.dataset.mode;
        $$("#colormode button").forEach((x) => x.classList.toggle("active", x === b));
        if (GT.results) GT.results.render();
      }));

    // кнопка расчёта
    $("#btn-compute").addEventListener("click", () => GT.compute());

    GT.refresh();
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
