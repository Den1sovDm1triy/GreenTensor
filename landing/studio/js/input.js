// SPDX-License-Identifier: MIT
// GreenTensor Studio — модуль ввода данных (геометрия, материалы, излучение).
(function () {
  "use strict";
  const GT = (window.GT = window.GT || {});

  const ICONS = { sphere: "●", cylinder_inf: "▥" };
  const NAMES = { sphere: "Сфера", cylinder_inf: "Цилиндр ∞" };
  const GEOM = {
    sphere: [["radius", "Радиус R", 1.0]],
    cylinder_inf: [["radius", "Радиус R", 0.5]],
  };

  // ---- маленькие конструкторы DOM ---- //
  function el(tag, attrs, children) {
    const e = document.createElement(tag);
    if (attrs) for (const k in attrs) {
      if (k === "class") e.className = attrs[k];
      else if (k === "html") e.innerHTML = attrs[k];
      else e.setAttribute(k, attrs[k]);
    }
    (children || []).forEach((c) => e.appendChild(typeof c === "string" ? document.createTextNode(c) : c));
    return e;
  }

  function numField(label, value, onInput, step) {
    const inp = el("input", { type: "number", step: step || "any", value: value });
    inp.addEventListener("input", () => onInput(parseFloat(inp.value)));
    return el("label", { class: "field" }, [el("span", {}, [label]), inp]);
  }

  // ---- дерево сцены ---- //
  function renderTree() {
    const ul = GT.$("#scene-tree");
    ul.innerHTML = "";
    const bodies = GT.state.scene.bodies;
    if (!bodies.length) {
      ul.appendChild(el("li", { class: "empty" }, ["Сцена пуста — добавьте примитив."]));
      return;
    }
    bodies.forEach((b, i) => {
      const li = el("li", { class: i === GT.state.selected ? "selected" : "" }, [
        el("span", { class: "body-icon" }, [ICONS[b.type] || "?"]),
        el("span", { class: "body-name" }, [NAMES[b.type] || b.type + " " + (i + 1)]),
        el("button", { class: "body-del", title: "Удалить" }, ["×"]),
      ]);
      li.addEventListener("click", (ev) => {
        if (ev.target.classList.contains("body-del")) {
          bodies.splice(i, 1);
          GT.state.selected = Math.max(0, Math.min(GT.state.selected, bodies.length - 1));
        } else {
          GT.state.selected = i;
        }
        GT.refresh();
      });
      ul.appendChild(li);
    });
  }

  // ---- свойства тела ---- //
  function renderBodyForm() {
    const host = GT.$("#body-form");
    host.innerHTML = "";
    const b = GT.selectedBody();
    if (!b) { host.className = "muted"; host.textContent = "Тело не выбрано."; return; }
    host.className = "";
    const grid = el("div", { class: "form-grid" });

    (GEOM[b.type] || []).forEach(([key, label]) => {
      grid.appendChild(numField(label, b[key] != null ? b[key] : 0,
        (v) => { b[key] = isFinite(v) ? v : 0; GT.editor.syncFromState(); }));
    });

    ["X", "Y", "Z"].forEach((ax, i) => {
      grid.appendChild(numField("Позиция " + ax, b.position[i],
        (v) => { b.position[i] = isFinite(v) ? v : 0; GT.editor.syncFromState(); }));
    });

    if (b.type === "sphere") {
      // полусфера на проводящем экране (метод зеркальных изображений)
      const h = b.hemisphere || (b.hemisphere = { enabled: false, feed_offset_deg: 0 });
      const chk = el("input", { type: "checkbox" });
      chk.checked = !!h.enabled;
      chk.addEventListener("change", () => { h.enabled = chk.checked; GT.refresh(); });
      grid.appendChild(el("label", { class: "field inline full",
        title: "Полусферическая линза на бесконечном PEC-экране: расчёт методом " +
               "зеркальных изображений (отбор гармоник по чётности)" },
        [chk, el("span", {}, ["Полусфера на PEC-экране"])]));
      if (h.enabled) {
        // линейка облучателей: одиночный луч — одна строка; линейка с весами —
        // профилированные ДН (например, csc²). Луч выходит зеркально: θ = −θ′.
        if (!h.feeds || !h.feeds.length) {
          h.feeds = [{ offset_deg: h.feed_offset_deg || 0, amp: 1, phase_deg: 0 }];
        }
        const wrap = el("div", { class: "full" });
        wrap.appendChild(el("span", { class: "muted" },
          ["Облучатели (луч выходит зеркально смещению: θ = −θ′)"]));
        const table = el("table", { class: "layers" });
        table.appendChild(el("thead", {}, [el("tr", {},
          ["θ′, °", "ампл.", "фаза, °", ""].map((t) => el("th", {}, [t])))]));
        const tbody = el("tbody");
        h.feeds.forEach((f, i) => {
          const tr = el("tr");
          [["offset_deg", 0], ["amp", 1], ["phase_deg", 0]].forEach(([key, dflt]) => {
            const inp = el("input", { type: "number", step: "any",
              value: f[key] != null ? f[key] : dflt });
            inp.addEventListener("input", () => {
              f[key] = parseFloat(inp.value);
              GT.editor.syncFromState();
            });
            tr.appendChild(el("td", {}, [inp]));
          });
          const del = el("button", { class: "row-del", title: "Удалить облучатель" }, ["×"]);
          del.addEventListener("click", () => {
            if (h.feeds.length > 1) { h.feeds.splice(i, 1); GT.refresh(); }
            else GT.toast("Нужен хотя бы один облучатель", "error");
          });
          tr.appendChild(el("td", {}, [del]));
          tbody.appendChild(tr);
        });
        table.appendChild(tbody);
        wrap.appendChild(table);
        const add = el("button", { class: "chip" }, ["＋ облучатель"]);
        add.addEventListener("click", () => {
          const last = h.feeds[h.feeds.length - 1] || { offset_deg: 0, amp: 1, phase_deg: 0 };
          h.feeds.push({ offset_deg: (last.offset_deg || 0) + 10, amp: last.amp || 1,
            phase_deg: last.phase_deg || 0 });
          GT.refresh();
        });
        wrap.appendChild(add);
        grid.appendChild(wrap);
      }
    }

    if (b.type === "cylinder_inf") {
      // бесконечный цилиндр (2D): падение задаётся углом от оси и поляризацией
      grid.appendChild(numField("Угол падения от оси, °", b.theta_deg != null ? b.theta_deg : 90,
        (v) => { b.theta_deg = isFinite(v) ? v : 90; }));
      const modeSel = el("select", {}, [
        el("option", { value: "TM" }, ["TM (E∥оси)"]),
        el("option", { value: "TE" }, ["TE (H∥оси)"]),
      ]);
      modeSel.value = b.mode || "TM";
      modeSel.addEventListener("change", () => { b.mode = modeSel.value; });
      grid.appendChild(el("label", { class: "field" },
        [el("span", {}, ["Поляризация (свип)"]), modeSel]));
    }
    host.appendChild(grid);
  }

  // ---- материалы (слои) ---- //
  function renderLayers() {
    const host = GT.$("#layers-form");
    host.innerHTML = "";
    const b = GT.selectedBody();
    if (!b) { host.className = "muted"; host.textContent = "Тело не выбрано."; return; }
    host.className = "";
    b.layers = b.layers || [{ a_rel: 1, eps_re: 1, eps_im: 0, mu_re: 1, mu_im: 0 }];

    const table = el("table", { class: "layers" });
    const head = el("tr", {}, ["a_rel", "ε′", "ε″", "μ′", "μ″", ""].map((h) => el("th", {}, [h])));
    table.appendChild(el("thead", {}, [head]));
    const tbody = el("tbody");
    b.layers.forEach((L, i) => {
      const tr = el("tr");
      [["a_rel", 0.01], ["eps_re"], ["eps_im"], ["mu_re"], ["mu_im"]].forEach(([key]) => {
        const inp = el("input", { type: "number", step: "any", value: L[key] != null ? L[key] : 0 });
        inp.addEventListener("input", () => { L[key] = parseFloat(inp.value); });
        tr.appendChild(el("td", {}, [inp]));
      });
      const del = el("button", { class: "row-del", title: "Удалить слой" }, ["×"]);
      del.addEventListener("click", () => {
        if (b.layers.length > 1) { b.layers.splice(i, 1); renderLayers(); }
        else GT.toast("Нужен хотя бы один слой", "error");
      });
      tr.appendChild(el("td", {}, [del]));
      tbody.appendChild(tr);
    });
    table.appendChild(tbody);
    host.appendChild(table);

    const add = el("button", { class: "chip" }, ["＋ слой"]);
    add.addEventListener("click", () => {
      const last = b.layers[b.layers.length - 1];
      // новый слой = новая внешняя граница (a_rel чуть больше), перенормируем к 1
      b.layers.forEach((L) => { L.a_rel = L.a_rel * 0.85; });
      b.layers.push({ a_rel: 1.0, eps_re: 1.0, eps_im: 0, mu_re: 1, mu_im: 0 });
      renderLayers();
    });
    host.appendChild(el("div", { class: "layers-actions" }, [
      el("span", { class: "muted" }, ["изнутри → наружу, внешняя a_rel = 1"]), add,
    ]));
  }

  // ---- излучение ---- //
  function renderRadiation() {
    const host = GT.$("#radiation-form");
    host.innerHTML = "";
    const r = GT.state.scene.radiation;
    const grid = el("div", { class: "form-grid" });

    grid.appendChild(numField("k (= k₀·R)", r.k, (v) => { r.k = isFinite(v) ? v : r.k; }));
    grid.appendChild(numField("toch (точность, опц.)", r.toch || "",
      (v) => { r.toch = isFinite(v) && v > 0 ? Math.round(v) : null; }));

    const polSel = el("select", {}, [
      el("option", { value: "linear" }, ["линейная"]),
      el("option", { value: "circular" }, ["круговая"]),
    ]);
    polSel.value = r.polarization;
    polSel.addEventListener("change", () => { r.polarization = polSel.value; renderRadiation(); });
    grid.appendChild(el("label", { class: "field" }, [el("span", {}, ["Поляризация"]), polSel]));

    const probSel = el("select", {}, [
      el("option", { value: "diffraction" }, ["дифракция"]),
      el("option", { value: "antenna" }, ["антенна (источник на сфере)"]),
    ]);
    probSel.value = r.problem || "diffraction";
    probSel.addEventListener("change", () => { r.problem = probSel.value; });
    grid.appendChild(el("label", { class: "field" }, [el("span", {}, ["Тип задачи"]), probSel]));

    if (r.polarization === "linear") {
      grid.appendChild(numField("φ среза, °", r.phi_deg, (v) => { r.phi_deg = isFinite(v) ? v : 0; }));
    }

    // падение: направление k̂ и поляризация (для несферических/кластеров)
    r.khat = r.khat || [0, 0, 1];
    r.pol = r.pol || [1, 0, 0];
    const k3 = el("div", { class: "full" }, [el("span", { class: "muted" }, ["направление падения k̂"])]);
    grid.appendChild(k3);
    ["x", "y", "z"].forEach((ax, i) => {
      grid.appendChild(numField("k̂·" + ax, r.khat[i], (v) => { r.khat[i] = isFinite(v) ? v : 0; }));
    });
    const p3 = el("div", { class: "full" }, [el("span", { class: "muted" }, ["вектор поляризации ê"])]);
    grid.appendChild(p3);
    ["x", "y", "z"].forEach((ax, i) => {
      grid.appendChild(numField("ê·" + ax, r.pol[i], (v) => { r.pol[i] = isFinite(v) ? v : 0; }));
    });

    host.appendChild(grid);
  }

  function render() {
    renderTree();
    renderBodyForm();
    renderLayers();
    renderRadiation();
  }

  function newBody(type) {
    const defaults = {
      sphere: { type: "sphere", position: [0, 0, 0], radius: 1.0,
        layers: [{ a_rel: 1, eps_re: 2.25, eps_im: 0, mu_re: 1, mu_im: 0 }], euler: [0, 0, 0] },
      cylinder_inf: { type: "cylinder_inf", position: [0, 0, 0], radius: 0.5,
        theta_deg: 90, mode: "TM",
        layers: [{ a_rel: 1, eps_re: 2.25, eps_im: 0, mu_re: 1, mu_im: 0 }] },
    };
    return JSON.parse(JSON.stringify(defaults[type]));
  }

  function init() {
    GT.$$("#input-module .add-bar [data-add]").forEach((btn) =>
      btn.addEventListener("click", () => {
        GT.state.scene.bodies.push(newBody(btn.dataset.add));
        GT.state.selected = GT.state.scene.bodies.length - 1;
        GT.refresh();
      }));
  }

  GT.input = { init, render };
})();
