// SPDX-License-Identifier: MIT
// GreenTensor Studio — эталонные задачи: загрузка и применение всех предустановок.
(function () {
  "use strict";
  const GT = (window.GT = window.GT || {});
  let PRESETS = [];

  async function load() {
    const sel = GT.$("#preset-select");
    try {
      const resp = await fetch("studio/engine/presets.json?v=20260721b");
      const data = await resp.json();
      PRESETS = data.presets || [];
    } catch (e) {
      GT.toast("Не удалось загрузить эталонные задачи", "error");
      return;
    }
    // группировка по group
    const groups = {};
    PRESETS.forEach((p) => { (groups[p.group] = groups[p.group] || []).push(p); });
    sel.innerHTML = '<option value="">— выбрать —</option>';
    Object.keys(groups).forEach((g) => {
      const og = document.createElement("optgroup");
      og.label = g;
      groups[g].forEach((p) => {
        const o = document.createElement("option");
        o.value = p.id; o.textContent = p.title;
        o.title = p.description || "";
        og.appendChild(o);
      });
      sel.appendChild(og);
    });
    sel.addEventListener("change", () => apply(sel.value));

    // Глубокая ссылка: studio.html#preset=<id> (или ?preset=<id>) — авто-применить эталон.
    // Работает и при переходе с галереи примеров лендинга.
    applyFromLocation();
    window.addEventListener("hashchange", applyFromLocation);
  }

  function applyFromLocation() {
    const m = (location.hash + " " + location.search).match(/preset=([A-Za-z0-9_]+)/);
    if (!m) return;
    const id = m[1];
    if (!PRESETS.find((x) => x.id === id)) return;
    const sel = GT.$("#preset-select");
    if (sel) sel.value = id;
    apply(id);
  }

  function apply(id) {
    if (!id) return;
    const p = PRESETS.find((x) => x.id === id);
    if (!p) return;
    // глубокая копия сцены пресета -> состояние (все предустановки разом)
    GT.state.scene = JSON.parse(JSON.stringify(p.scene));
    // кривые Ansys HFSS (МКЭ) пресета — список по плоскостям; накладываются на диаграмму
    GT.state.reference = p.reference ? JSON.parse(JSON.stringify(p.reference)) : null;
    if (GT.state.reference && GT.state.reference.length) {   // кривые HFSS в дБ — включаем дБ
      const dbt = GT.$("#db-toggle"); if (dbt) dbt.checked = true;
    }
    // полусфера на экране: нормаль (θ=0) — вверх, экран — горизонтально
    const hemi = (p.scene.bodies || []).some((b) => b.hemisphere && b.hemisphere.enabled);
    const rot = GT.$("#polar-rot");
    if (rot) {
      rot.value = hemi ? 90 : 0;
      const out = GT.$("#polar-rot-val"); if (out) out.textContent = rot.value;
    }
    GT.state.selected = 0;
    GT.state.results = null;
    GT.refresh();
    GT.toast("Эталон: " + p.title, "ok");
    GT.compute();   // сразу считаем + (если есть) накладываем эталон Ansys
  }

  GT.presets = { load, apply, all: () => PRESETS };
})();
