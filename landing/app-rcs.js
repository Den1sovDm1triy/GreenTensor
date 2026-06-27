/* GreenTensor — web solver bootstrap + UI logic.
 * Loads Pyodide, imports the same green_tensor/calc.py used by the desktop
 * library, drives RCSCalculator / ScatteringCalculator from form inputs,
 * renders results with Plotly, and shows the live problem geometry in a
 * Three.js viewport. */

(function () {
  'use strict';

  /* ====================== State & presets ====================== */
  const PRESETS = {
    luneburg4: {
      name: '4-слойный Люнеберг', k0: 5, toch: 10, phi: Math.PI / 2,
      layers: [
        { a: 0.53, er: 1.86, ei: 0, mr: 1, mi: 0 },
        { a: 0.75, er: 1.57, ei: 0, mr: 1, mi: 0 },
        { a: 0.93, er: 1.28, ei: 0, mr: 1, mi: 0 },
        { a: 1.00, er: 1.00, ei: 0, mr: 1, mi: 0 },
      ],
    },
    luneburg9: {
      name: '9-слойный Люнеберг', k0: 5, toch: 10, phi: Math.PI / 2,
      layers: [
        { a: 0.34, er: 1.94, ei: 0, mr: 1, mi: 0 },
        { a: 0.49, er: 1.82, ei: 0, mr: 1, mi: 0 },
        { a: 0.59, er: 1.71, ei: 0, mr: 1, mi: 0 },
        { a: 0.69, er: 1.59, ei: 0, mr: 1, mi: 0 },
        { a: 0.77, er: 1.47, ei: 0, mr: 1, mi: 0 },
        { a: 0.84, er: 1.35, ei: 0, mr: 1, mi: 0 },
        { a: 0.91, er: 1.24, ei: 0, mr: 1, mi: 0 },
        { a: 0.97, er: 1.12, ei: 0, mr: 1, mi: 0 },
        { a: 1.00, er: 1.00, ei: 0, mr: 1, mi: 0 },
      ],
    },
    metal: {
      name: 'Металлическая сфера', k0: 0.5, toch: 10, phi: Math.PI / 2,
      layers: [
        { a: 0.1, er: 0, ei: -1.7e7, mr: 1, mi: 0 },
        { a: 0.2, er: 0, ei: -1.7e7, mr: 1, mi: 0 },
        { a: 1.0, er: 1, ei: 0,      mr: 1, mi: 0 },
      ],
    },
    head: {
      name: 'Модель головы', k0: 2.0, toch: 10, phi: Math.PI / 2,
      layers: [
        { a: 0.85, er: 45, ei: -13, mr: 1, mi: 0 },
        { a: 0.92, er: 11, ei: -2,  mr: 1, mi: 0 },
        { a: 1.00, er: 1,  ei: 0,   mr: 1, mi: 0 },
      ],
    },
  };

  const state = {
    layers: PRESETS.luneburg4.layers.map(l => ({ ...l })),
    pyodide: null,
    pyReady: false,
    activeTab: 'rcs',
    running: false,
    lastResult: null,
    lastResultKind: null,
  };

  /* ====================== DOM helpers ====================== */
  const $ = id => document.getElementById(id);

  function setStatus(cls, text, progress) {
    const el = $('solver-status'); if (!el) return;
    el.classList.remove('loading', 'ready', 'error');
    el.classList.add(cls);
    $('solver-status-text').textContent = text;
    $('solver-status-progress').textContent = progress || '';
  }
  function setRunning(on) {
    state.running = on;
    const btn = $('solver-run');
    btn.classList.toggle('running', on);
    btn.disabled = on || !state.pyReady;
    $('solver-run-text').textContent = on ? 'Расчёт…' : 'Рассчитать';
  }

  /* ====================== Color picker for ε ====================== */
  // Returns hex color reflecting layer permittivity.
  function epsColor(er, ei) {
    const mag = Math.hypot(er, ei);
    if (mag > 1e3) return '#5a6d68';        // metal: dark grey
    const re = er;
    if (re >= 2.5) return '#ff6b6b';        // high ε → red
    if (re >= 1.7) return '#ffb547';        // medium-high → orange
    if (re >= 1.2) return '#ffe066';        // medium → yellow
    if (re >= 1.05) return '#21d68a';       // low → green
    return '#00ff9c';                       // vacuum-like
  }

  /* ====================== Layer table ====================== */
  function renderLayers() {
    const tbody = $('solver-layers-body'); if (!tbody) return;
    tbody.innerHTML = '';
    state.layers.forEach((layer, idx) => {
      const tr = document.createElement('tr');
      tr.innerHTML = `
        <td>${idx + 1}</td>
        <td><input type="number" step="0.01" min="0" data-key="a"  data-idx="${idx}" value="${layer.a}"></td>
        <td><input type="number" step="0.1"           data-key="er" data-idx="${idx}" value="${layer.er}"></td>
        <td><input type="number" step="0.1"           data-key="ei" data-idx="${idx}" value="${layer.ei}"></td>
        <td><input type="number" step="0.1"           data-key="mr" data-idx="${idx}" value="${layer.mr}"></td>
        <td><input type="number" step="0.1"           data-key="mi" data-idx="${idx}" value="${layer.mi}"></td>
      `;
      tbody.appendChild(tr);
    });
    tbody.querySelectorAll('input').forEach(inp => {
      inp.addEventListener('input', e => {
        const i = +e.target.dataset.idx;
        const k = e.target.dataset.key;
        const v = parseFloat(e.target.value);
        if (!Number.isNaN(v)) {
          state.layers[i][k] = v;
          updateScene();
          updateOverlay();
        }
      });
    });
    $('solver-del-layer').disabled = state.layers.length <= 2;
    updateScene();
    updateOverlay();
    renderLegend();
  }

  function applyPreset(key) {
    const p = PRESETS[key]; if (!p) return;
    state.layers = p.layers.map(l => ({ ...l }));
    $('solver-k0').value = p.k0;
    $('solver-toch').value = p.toch;
    $('solver-phi').value = p.phi.toFixed(4);
    renderLayers();
  }

  /* ====================== Tabs ====================== */
  function setTab(name) {
    state.activeTab = name;
    document.querySelectorAll('.tb-tab, .solver-tab').forEach(t => {
      const on = t.dataset.tab === name;
      t.classList.toggle('active', on);
      t.setAttribute('aria-selected', on);
    });
    document.body.classList.toggle('tab-sweep', name === 'sweep');
  }

  /* ====================== 3D scene (Three.js) ====================== */
  const scene3 = {
    inited: false,
    scene: null, camera: null, renderer: null,
    root: null, layerMeshes: [], wave: null,
    width: 0, height: 0,
  };

  function initScene() {
    if (scene3.inited || typeof THREE === 'undefined') return;
    const canvas = $('scene-canvas'); if (!canvas) return;
    const wrap = canvas.parentElement;
    const w = wrap.clientWidth, h = wrap.clientHeight;
    scene3.width = w; scene3.height = h;

    scene3.scene = new THREE.Scene();
    scene3.camera = new THREE.PerspectiveCamera(40, w / Math.max(h, 1), 0.1, 100);
    scene3.camera.position.set(2.6, 1.5, 3.4);
    scene3.camera.lookAt(0, 0, 0);

    scene3.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    scene3.renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    scene3.renderer.setSize(w, h, false);
    scene3.renderer.setClearColor(0x000000, 0);

    const ambient = new THREE.AmbientLight(0xffffff, 0.55);
    scene3.scene.add(ambient);
    const dir = new THREE.DirectionalLight(0xffffff, 0.8);
    dir.position.set(3, 4, 5); scene3.scene.add(dir);
    const sc0 = themeColors();
    const rim = new THREE.DirectionalLight(new THREE.Color(sc0.accent), 0.45);
    rim.position.set(-3, -2, -4); scene3.scene.add(rim);
    scene3.rim = rim;

    scene3.root = new THREE.Group();
    scene3.scene.add(scene3.root);

    // Coordinate axes (subtle)
    const axesMat = new THREE.LineBasicMaterial({ color: new THREE.Color(sc0.muted), transparent: true, opacity: 0.5 });
    scene3.axesMat = axesMat;
    const axesGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(-1.4, 0, 0), new THREE.Vector3(1.4, 0, 0),
      new THREE.Vector3(0, -1.4, 0), new THREE.Vector3(0, 1.4, 0),
      new THREE.Vector3(0, 0, -1.4), new THREE.Vector3(0, 0, 1.4),
    ]);
    const axesObj = new THREE.LineSegments(axesGeo, axesMat); scene3.root.add(axesObj);

    // Incident wave: arrow and parallel rings approaching the sphere from -z
    const waveGroup = new THREE.Group();
    const waveColor = new THREE.Color(sc0.accent);
    scene3.waveMats = [];
    for (let i = 0; i < 4; i++) {
      const g = new THREE.RingGeometry(0.6, 0.61, 64);
      const m = new THREE.MeshBasicMaterial({
        color: waveColor.clone(), transparent: true, opacity: 0.35 - i * 0.06, side: THREE.DoubleSide,
      });
      scene3.waveMats.push(m);
      const ring = new THREE.Mesh(g, m);
      ring.position.z = -2.0 - i * 0.45;
      ring.userData.kind = 'wave';
      ring.userData.phase = i * 0.45;
      waveGroup.add(ring);
    }
    // arrow
    const arrowMat = new THREE.LineBasicMaterial({ color: new THREE.Color(sc0.accent) });
    scene3.arrowMat = arrowMat;
    const arrowGeo = new THREE.BufferGeometry().setFromPoints([
      new THREE.Vector3(0, 0, -3.2), new THREE.Vector3(0, 0, -1.2),
    ]);
    waveGroup.add(new THREE.Line(arrowGeo, arrowMat));
    const headGeo = new THREE.ConeGeometry(0.05, 0.18, 12);
    const headMat = new THREE.MeshBasicMaterial({ color: new THREE.Color(sc0.accent) });
    scene3.headMat = headMat;
    const head = new THREE.Mesh(headGeo, headMat);
    head.rotation.x = Math.PI / 2; head.position.set(0, 0, -1.05);
    waveGroup.add(head);
    scene3.wave = waveGroup; scene3.root.add(waveGroup);

    // Resize handling
    const ro = new ResizeObserver(() => {
      const w2 = wrap.clientWidth, h2 = wrap.clientHeight;
      if (!w2 || !h2) return;
      scene3.width = w2; scene3.height = h2;
      scene3.camera.aspect = w2 / h2;
      scene3.camera.updateProjectionMatrix();
      scene3.renderer.setSize(w2, h2, false);
    });
    ro.observe(wrap);

    scene3.inited = true;
    animate();
  }

  function disposeMesh(m) {
    if (!m) return;
    if (m.geometry) m.geometry.dispose();
    if (m.material) {
      if (Array.isArray(m.material)) m.material.forEach(x => x.dispose());
      else m.material.dispose();
    }
  }

  function updateScene() {
    if (!scene3.inited) return;
    // Remove old layer meshes
    scene3.layerMeshes.forEach(m => { scene3.root.remove(m); disposeMesh(m); });
    scene3.layerMeshes = [];

    const layers = state.layers;
    if (!layers.length) return;
    const aMax = layers[layers.length - 1].a || 1;

    // Build outermost first so transparency stacks correctly: render order from outer → inner.
    layers.slice().reverse().forEach((layer, idxFromEnd) => {
      const idx = layers.length - 1 - idxFromEnd;
      const r = (layer.a / aMax) * 1.0;
      const color = new THREE.Color(epsColor(layer.er, layer.ei));
      const isLast = (idx === layers.length - 1);
      const isMetal = Math.hypot(layer.er, layer.ei) > 1e3;

      // Wireframe
      const wireGeo = new THREE.SphereGeometry(r, 48, 32);
      const wireMat = new THREE.MeshBasicMaterial({
        color, wireframe: true, transparent: true,
        opacity: isLast ? 0.18 : (isMetal ? 0.55 : 0.32),
      });
      const wire = new THREE.Mesh(wireGeo, wireMat);
      scene3.root.add(wire); scene3.layerMeshes.push(wire);

      // Solid translucent shell (skip for outermost vacuum to keep view clean)
      if (!(isLast && layer.er <= 1.05 && Math.abs(layer.ei) < 0.01)) {
        const meshGeo = new THREE.SphereGeometry(r, 48, 32);
        const meshMat = new THREE.MeshPhongMaterial({
          color, transparent: true,
          opacity: isMetal ? 0.85 : Math.min(0.18 + 0.04 * layer.er, 0.45),
          depthWrite: false,
          shininess: isMetal ? 90 : 28,
          emissive: color, emissiveIntensity: isMetal ? 0.05 : 0.12,
          side: THREE.FrontSide,
        });
        const mesh = new THREE.Mesh(meshGeo, meshMat);
        mesh.renderOrder = idx; // inner shells drawn first
        scene3.root.add(mesh); scene3.layerMeshes.push(mesh);
      }
    });
  }

  function renderLegend() {
    const el = $('scene-legend'); if (!el) return;
    el.innerHTML = '';
    state.layers.forEach((layer, i) => {
      const div = document.createElement('div');
      div.className = 'item';
      const sw = document.createElement('div');
      sw.className = 'swatch';
      sw.style.background = epsColor(layer.er, layer.ei);
      const txt = document.createElement('span');
      const epsText = layer.ei !== 0 ? `${layer.er}${layer.ei < 0 ? '' : '+'}${layer.ei}j` : `${layer.er}`;
      txt.textContent = `L${i + 1} · a=${layer.a} · ε=${epsText}`;
      div.appendChild(sw); div.appendChild(txt);
      el.appendChild(div);
    });
  }

  function updateOverlay() {
    if ($('ov-n')) $('ov-n').textContent = state.layers.length;
    if ($('ov-k0')) $('ov-k0').textContent = ($('solver-k0') ? $('solver-k0').value : '—');
    if ($('ov-eps') && state.layers[0]) {
      const l = state.layers[0];
      $('ov-eps').textContent = l.ei !== 0 ? `${l.er}${l.ei < 0 ? '' : '+'}${l.ei}j` : `${l.er}`;
    }
  }

  function animate() {
    if (!scene3.inited) return;
    requestAnimationFrame(animate);
    const t = performance.now() * 0.001;
    if (scene3.root) scene3.root.rotation.y = t * 0.15;
    if (scene3.wave) {
      scene3.wave.children.forEach((c, i) => {
        if (c.userData && c.userData.kind === 'wave') {
          c.position.z = -2.0 - ((t + c.userData.phase) % 1.4);
        }
      });
    }
    scene3.renderer.render(scene3.scene, scene3.camera);
  }

  // Re-skin the persistent scene materials (wave/arrow/axes/rim) on theme toggle.
  function applySceneTheme() {
    if (!scene3.inited || typeof THREE === 'undefined') return;
    const c = themeColors();
    const accent = new THREE.Color(c.accent);
    if (scene3.rim) scene3.rim.color.set(accent);
    if (scene3.axesMat) scene3.axesMat.color.set(new THREE.Color(c.muted));
    (scene3.waveMats || []).forEach(m => m.color.set(accent));
    if (scene3.arrowMat) scene3.arrowMat.color.set(accent);
    if (scene3.headMat) scene3.headMat.color.set(accent);
  }

  /* ====================== Pyodide bootstrap ====================== */
  async function bootstrap() {
    try {
      if (typeof loadPyodide !== 'function') throw new Error('Pyodide CDN не загрузился');
      setStatus('loading', 'Загружается CPython runtime…', 'шаг 1/3');
      state.pyodide = await loadPyodide({ indexURL: 'https://cdn.jsdelivr.net/pyodide/v0.27.7/full/' });

      setStatus('loading', 'Устанавливаются numpy + scipy…', 'шаг 2/3');
      await state.pyodide.loadPackage(['numpy', 'scipy']);

      setStatus('loading', 'Подключается green_tensor/calc.py…', 'шаг 3/3');
      const src = await fetch('calc.py').then(r => {
        if (!r.ok) throw new Error('calc.py не найден (' + r.status + ')');
        return r.text();
      });
      state.pyodide.FS.writeFile('calc.py', src);

      // Stub matplotlib so calc.py imports cleanly without dragging MPL into wasm.
      state.pyodide.runPython(`
import sys, types
class _Stub(types.ModuleType):
    def __init__(self, name='stub'): super().__init__(name)
    def __getattr__(self, name): return _Stub(name)
    def __call__(self, *a, **kw): return self
    def __getitem__(self, k): return _Stub('item')
for m in ('matplotlib', 'matplotlib.pyplot', 'matplotlib.ticker'):
    sys.modules[m] = _Stub(m)
import calc
`);

      state.pyReady = true;
      setStatus('ready', 'Движок готов', '');
      $('solver-run').disabled = false;
      $('solver-export').disabled = false;
      $('solver-run-text').textContent = 'Рассчитать';
    } catch (e) {
      console.error(e);
      setStatus('error', 'Ошибка инициализации: ' + (e && e.message ? e.message : e), '');
    }
  }

  /* ====================== Param collection ====================== */
  function collectParams() {
    const k0 = parseFloat($('solver-k0').value);
    const toch = parseInt($('solver-toch').value, 10);
    const phi = parseFloat($('solver-phi').value);
    const layers = state.layers;
    return {
      k0, toch, phi, n: layers.length,
      a: layers.map(l => l.a),
      eps_pairs: layers.map(l => [l.er, l.ei]),
      miy_pairs: layers.map(l => [l.mr, l.mi]),
    };
  }

  function pyNum(x) { return Number.isFinite(x) ? x.toString() : '0'; }
  function pyList(arr) { return '[' + arr.map(pyNum).join(',') + ']'; }
  function pyPairs(pairs) { return '[' + pairs.map(p => '[' + pyNum(p[0]) + ',' + pyNum(p[1]) + ']').join(',') + ']'; }

  /* ====================== Plotly themes ====================== */
  // Colors are read live from the active theme's CSS variables, so light/dark
  // toggling re-skins every chart on the next render (see gt-theme-change).
  function cssVar(name, fallback) {
    try {
      const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
      return v || fallback;
    } catch (e) { return fallback; }
  }
  function themeColors() {
    return {
      accent: cssVar('--accent', '#00ff9c'),
      warn:   cssVar('--warn', '#ffb547'),
      green2: cssVar('--accent-2', '#21d68a'),
      danger: cssVar('--danger', '#ff6b6b'),
      text:   cssVar('--text', '#d8efe5'),
      muted:  cssVar('--muted', '#6f8a82'),
      paper:  cssVar('--code-bg', '#0a1411'),
      grid:   cssVar('--accent-dim', 'rgba(0,255,156,0.12)'),
    };
  }
  function polarLayout() {
    const c = themeColors();
    return {
      polar: {
        radialaxis: {
          range: [0, 60], color: c.muted, gridcolor: c.grid,
          tickfont: { family: 'JetBrains Mono, monospace', size: 10, color: c.muted },
          tickvals: [10, 20, 30, 40, 50, 60],
          ticktext: ['-50', '-40', '-30', '-20', '-10', '0'],
        },
        angularaxis: {
          color: c.muted, gridcolor: c.grid,
          tickfont: { family: 'JetBrains Mono, monospace', size: 10, color: c.muted },
          rotation: 90, direction: 'counterclockwise',
        },
        bgcolor: c.paper,
      },
      paper_bgcolor: c.paper, plot_bgcolor: c.paper,
      font: { color: c.text, family: 'Inter, sans-serif' },
      margin: { l: 20, r: 20, t: 30, b: 20 }, showlegend: false,
    };
  }
  function cartesianLayout(xlabel, ylabel) {
    const c = themeColors();
    return {
      xaxis: {
        title: { text: xlabel, font: { color: c.muted, size: 11 } },
        color: c.muted, gridcolor: c.grid,
        tickfont: { family: 'JetBrains Mono, monospace', size: 10, color: c.muted },
      },
      yaxis: {
        title: { text: ylabel, font: { color: c.muted, size: 11 } },
        color: c.muted, gridcolor: c.grid,
        tickfont: { family: 'JetBrains Mono, monospace', size: 10, color: c.muted },
      },
      paper_bgcolor: c.paper, plot_bgcolor: c.paper,
      font: { color: c.text, family: 'Inter, sans-serif' },
      margin: { l: 50, r: 14, t: 30, b: 40 },
      legend: { font: { color: c.text, size: 10 }, bgcolor: 'rgba(0,0,0,0)', orientation: 'h', y: 1.12 },
    };
  }
  const PLOTLY_CFG = { responsive: true, displaylogo: false, modeBarButtonsToRemove: ['lasso2d', 'select2d'] };

  function clearEmpty(id) {
    const host = $(id); if (!host) return;
    host.querySelectorAll('.viewport-plot-empty, .solver-plot-empty').forEach(el => el.remove());
  }

  /* ====================== RCS ====================== */
  function runRcs() {
    const p = collectParams();
    if (p.n < 2) throw new Error('Нужно минимум 2 слоя');
    const code = `
import calc, numpy as np
rcs = calc.RCSCalculator(
    k0=${p.k0}, toch=${p.toch}, n=${p.n}, phi=${p.phi},
    a=${pyList(p.a)},
    eps=[complex(re, im) for re, im in ${pyPairs(p.eps_pairs)}],
    miy=[complex(re, im) for re, im in ${pyPairs(p.miy_pairs)}],
)
rcs.run_calculation()
def _safe(arr):
    arr = np.asarray(arr, dtype=float)
    return np.where(np.isfinite(arr), arr, -60.0).tolist()
{
  'teta_deg': (rcs.teta * 180.0 / np.pi).tolist(),
  'E_teta_dB': _safe(rcs.DN_NORM_lin_dB_teta),
  'E_phi_dB': _safe(rcs.DN_NORM_lin_dB_phi),
  'circle_op_dB': _safe(rcs.DN_NORM_circle_dB_op),
  'circle_kp_dB': _safe(rcs.DN_NORM_circle_dB_kp),
  'n': ${p.n}, 'toch': ${p.toch},
}
`.trim();
    const t0 = performance.now();
    const proxy = state.pyodide.runPython(code);
    const res = proxy.toJs({ dict_converter: Object.fromEntries });
    if (proxy.destroy) proxy.destroy();
    const dt = performance.now() - t0;
    state.lastResult = res; state.lastResultKind = 'rcs'; state.lastDt = dt;
    plotRcs(res, dt);
  }
  function plotRcs(res, dt) {
    const c = themeColors();
    const teta = res.teta_deg;
    const shift = v => Math.max(v, -60) + 60;
    clearEmpty('rcs-plot-teta'); clearEmpty('rcs-plot-phi');
    Plotly.react('rcs-plot-teta', [{
      type: 'scatterpolar', r: res.E_teta_dB.map(shift), theta: teta, mode: 'lines',
      line: { color: c.accent, width: 2 }, name: 'E_θ',
    }], polarLayout(), PLOTLY_CFG);
    Plotly.react('rcs-plot-phi', [{
      type: 'scatterpolar', r: res.E_phi_dB.map(shift), theta: teta, mode: 'lines',
      line: { color: c.warn, width: 2 }, name: 'E_φ',
    }], polarLayout(), PLOTLY_CFG);
    $('rcs-m-time').textContent = (dt / 1000).toFixed(2) + 'с';
    $('rcs-m-n').textContent = res.n;
    $('rcs-m-toch').textContent = res.toch;
  }

  /* ====================== Sweep ====================== */
  function runSweep() {
    const p = collectParams();
    const k0a_start = parseFloat($('sweep-start').value);
    const k0a_stop = parseFloat($('sweep-stop').value);
    const k0a_step = parseFloat($('sweep-step').value);
    if (k0a_stop <= k0a_start) throw new Error('k₀a stop должен быть > start');
    if (k0a_step <= 0) throw new Error('k₀a step должен быть > 0');
    const npts = Math.ceil((k0a_stop - k0a_start) / k0a_step) + 1;
    if (npts > 500) throw new Error('Слишком много точек (' + npts + '), увеличьте шаг.');
    const code = `
import calc, numpy as np
sc = calc.ScatteringCalculator(toch=${p.toch}, k0a_start=${k0a_start}, k0a_stop=${k0a_stop}, k0a_step=${k0a_step})
sc.rcs_calculator = calc.RCSCalculator(
    toch=${p.toch}, n=${p.n}, phi=${p.phi},
    a=${pyList(p.a)},
    eps=[complex(re, im) for re, im in ${pyPairs(p.eps_pairs)}],
    miy=[complex(re, im) for re, im in ${pyPairs(p.miy_pairs)}],
)
k0a, sigma_s, sigma_r, sigma_p, sigma_theta, sigma_phi = sc.calculate()
def _f(a):
    a = np.asarray(a, dtype=float)
    return np.where(np.isfinite(a), a, 0.0).tolist()
{
  'k0a': _f(k0a), 'sigma_s': _f(sigma_s), 'sigma_r': _f(sigma_r),
  'sigma_p': _f(sigma_p), 'sigma_theta': _f(sigma_theta), 'sigma_phi': _f(sigma_phi),
}
`.trim();
    const t0 = performance.now();
    const proxy = state.pyodide.runPython(code);
    const res = proxy.toJs({ dict_converter: Object.fromEntries });
    if (proxy.destroy) proxy.destroy();
    const dt = performance.now() - t0;
    state.lastResult = res; state.lastResultKind = 'sweep'; state.lastDt = dt;
    plotSweep(res, dt);
  }
  function plotSweep(res, dt) {
    const c = themeColors();
    clearEmpty('sweep-plot-sigma'); clearEmpty('sweep-plot-rp');
    Plotly.react('sweep-plot-sigma', [
      { x: res.k0a, y: res.sigma_s,     mode: 'lines', name: 'σ_s',  line: { color: c.accent, width: 2 } },
      { x: res.k0a, y: res.sigma_theta, mode: 'lines', name: 'σ_θ',  line: { color: c.danger, width: 1.5, dash: 'dash' } },
      { x: res.k0a, y: res.sigma_phi,   mode: 'lines', name: 'σ_φ',  line: { color: c.green2, width: 1.5, dash: 'dash' } },
    ], Object.assign(cartesianLayout('k₀a', 'σ'), { showlegend: true }), PLOTLY_CFG);
    Plotly.react('sweep-plot-rp', [
      { x: res.k0a, y: res.sigma_r, mode: 'lines', name: 'σ_r (back)',    line: { color: c.danger, width: 2 } },
      { x: res.k0a, y: res.sigma_p, mode: 'lines', name: 'σ_p (forward)', line: { color: c.warn, width: 2 } },
    ], Object.assign(cartesianLayout('k₀a', 'σ'), { showlegend: true }), PLOTLY_CFG);
    const maxS = Math.max.apply(null, res.sigma_s.filter(Number.isFinite));
    $('sweep-m-time').textContent = (dt / 1000).toFixed(2) + 'с';
    $('sweep-m-pts').textContent = res.k0a.length;
    $('sweep-m-max').textContent = maxS.toFixed(2);
  }

  /* ====================== Run dispatch ====================== */
  async function run() {
    if (state.running || !state.pyReady) return;
    setRunning(true);
    await new Promise(r => setTimeout(r, 16));
    try {
      if (state.activeTab === 'sweep') runSweep(); else runRcs();
    } catch (e) {
      console.error(e);
      setStatus('error', 'Ошибка расчёта: ' + (e && e.message ? e.message : e), '');
    } finally {
      setRunning(false);
      if (state.pyReady) setStatus('ready', 'Движок готов', '');
    }
  }

  /* ====================== CSV export ====================== */
  function exportCsv() {
    const r = state.lastResult; if (!r) return;
    let csv, name;
    if (state.lastResultKind === 'rcs') {
      csv = 'theta_deg,E_theta_dB,E_phi_dB,circle_op_dB,circle_kp_dB\n';
      for (let i = 0; i < r.teta_deg.length; i++) {
        csv += [r.teta_deg[i], r.E_teta_dB[i], r.E_phi_dB[i], r.circle_op_dB[i], r.circle_kp_dB[i]].join(',') + '\n';
      }
      name = 'greentensor_rcs.csv';
    } else {
      csv = 'k0a,sigma_s,sigma_r,sigma_p,sigma_theta,sigma_phi\n';
      for (let i = 0; i < r.k0a.length; i++) {
        csv += [r.k0a[i], r.sigma_s[i], r.sigma_r[i], r.sigma_p[i], r.sigma_theta[i], r.sigma_phi[i]].join(',') + '\n';
      }
      name = 'greentensor_sweep.csv';
    }
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = name; a.click();
    URL.revokeObjectURL(url);
  }

  /* ====================== Hash deeplinks ====================== */
  // solver.html#preset=luneburg4  → apply preset on load
  // solver.html#preset=metal&tab=sweep → apply preset and switch tab
  function applyHashDeeplink() {
    const h = (window.location.hash || '').replace(/^#/, '');
    if (!h) return;
    const params = Object.fromEntries(h.split('&').map(p => p.split('=')));
    if (params.preset && PRESETS[params.preset]) applyPreset(params.preset);
    if (params.tab === 'sweep' || params.tab === 'rcs') setTab(params.tab);
  }

  /* ====================== Init ====================== */
  function init() {
    if (!$('solver-status')) return; // not on a solver page

    initScene();
    renderLayers();
    setTab('rcs');
    applyHashDeeplink();
    window.addEventListener('hashchange', applyHashDeeplink);

    document.querySelectorAll('.tb-preset, .solver-preset').forEach(b => b.addEventListener('click', () => applyPreset(b.dataset.preset)));
    document.querySelectorAll('.tb-tab, .solver-tab').forEach(t => t.addEventListener('click', () => setTab(t.dataset.tab)));
    $('solver-add-layer').addEventListener('click', () => {
      const last = state.layers[state.layers.length - 1] || { a: 1, er: 1, ei: 0, mr: 1, mi: 0 };
      state.layers.splice(state.layers.length - 1, 0, { a: last.a * 0.9, er: last.er, ei: last.ei, mr: last.mr, mi: last.mi });
      renderLayers();
    });
    $('solver-del-layer').addEventListener('click', () => {
      if (state.layers.length > 2) { state.layers.splice(state.layers.length - 2, 1); renderLayers(); }
    });
    $('solver-run').addEventListener('click', run);
    $('solver-export').addEventListener('click', exportCsv);
    $('solver-k0').addEventListener('input', updateOverlay);

    // Light/dark toggle (dispatched by solver.html) → re-skin scene + redraw charts.
    window.addEventListener('gt-theme-change', () => {
      applySceneTheme();
      if (!state.lastResult) return;
      if (state.lastResultKind === 'rcs') plotRcs(state.lastResult, state.lastDt || 0);
      else if (state.lastResultKind === 'sweep') plotSweep(state.lastResult, state.lastDt || 0);
    });

    bootstrap();
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
