// SPDX-License-Identifier: MIT
// GreenTensor Studio — модуль 3D-редактора (Three.js r128, собственные orbit-controls).
(function () {
  "use strict";
  const GT = (window.GT = window.GT || {});

  let renderer, scene, camera, bodiesGroup, incidenceArrow, grid, overlayGroup, raf = null;
  // сферические координаты камеры вокруг target
  const cam = { radius: 6, theta: 1.1, phi: 0.9, target: null };
  // Палитра сцены под тему (сшито с лендингом): тёмный неон-canvas / светлый научный.
  const COL = {};
  function paletteFor() {
    const dark = (document.documentElement.getAttribute("data-theme") || "dark") !== "light";
    return dark
      ? { body: 0x21d68a, sel: 0xffb547, edge: 0x2a4038, bg: 0x07090a, grid1: 0x223530, grid2: 0x16211c }
      : { body: 0x0d6b8a, sel: 0xe08a1e, edge: 0x2a3a44, bg: 0xf4f6f8, grid1: 0xc2cbd3, grid2: 0xe6eaee };
  }
  function buildGrid() {
    if (grid) { scene.remove(grid); grid.geometry.dispose(); grid.material.dispose(); }
    grid = new THREE.GridHelper(12, 24, COL.grid1, COL.grid2);
    grid.rotation.x = Math.PI / 2;          // сетка в плоскости z=0 (XY)
    scene.add(grid);
  }
  function applyTheme() {
    Object.assign(COL, paletteFor());
    if (!renderer) return;
    renderer.setClearColor(COL.bg, 1);
    buildGrid();
    syncFromState();                        // пересобрать тела с новыми цветами материалов
    syncOverlay();                          // перекрасить наложенную ДН
  }

  function init() {
    Object.assign(COL, paletteFor());        // цвета сцены по текущей теме
    const host = GT.$("#viewport");
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setPixelRatio(window.devicePixelRatio || 1);
    renderer.setClearColor(COL.bg, 1);
    host.appendChild(renderer.domElement);

    scene = new THREE.Scene();
    camera = new THREE.PerspectiveCamera(45, 1, 0.01, 1000);
    cam.target = new THREE.Vector3(0, 0, 0);

    scene.add(new THREE.AmbientLight(0xffffff, 0.65));
    const d1 = new THREE.DirectionalLight(0xffffff, 0.6); d1.position.set(4, 6, 8); scene.add(d1);
    const d2 = new THREE.DirectionalLight(0xffffff, 0.25); d2.position.set(-6, -3, -4); scene.add(d2);

    buildGrid();
    scene.add(new THREE.AxesHelper(2.2));

    bodiesGroup = new THREE.Group();
    scene.add(bodiesGroup);

    setupControls(host);
    setupViewButtons();
    const dnChk = GT.$("#dn-overlay");
    if (dnChk) dnChk.addEventListener("change", syncOverlay);
    window.addEventListener("resize", resize);
    resize();
    syncFromState();
    loop();
  }

  // ---- собственные orbit-controls ---- //
  function setupControls(host) {
    let dragging = 0, lx = 0, ly = 0;
    const el = renderer.domElement;
    el.addEventListener("contextmenu", (e) => e.preventDefault());
    el.addEventListener("mousedown", (e) => {
      dragging = e.button === 2 ? 2 : 1; lx = e.clientX; ly = e.clientY; e.preventDefault();
    });
    window.addEventListener("mouseup", () => (dragging = 0));
    window.addEventListener("mousemove", (e) => {
      if (!dragging) return;
      const dx = e.clientX - lx, dy = e.clientY - ly; lx = e.clientX; ly = e.clientY;
      if (dragging === 1) {                 // вращение
        cam.theta = Math.min(Math.PI - 0.01, Math.max(0.01, cam.theta - dy * 0.008));
        cam.phi -= dx * 0.008;
      } else {                              // панорама (ПКМ)
        const f = cam.radius * 0.0016;
        const right = new THREE.Vector3().subVectors(camera.position, cam.target).cross(camera.up).normalize();
        const up = camera.up.clone();
        cam.target.addScaledVector(right, -dx * f).addScaledVector(up, dy * f);
      }
    });
    el.addEventListener("wheel", (e) => {
      cam.radius = Math.min(200, Math.max(0.3, cam.radius * (1 + Math.sign(e.deltaY) * 0.1)));
      e.preventDefault();
    }, { passive: false });
  }

  function setupViewButtons() {
    GT.$$(".hud-axes [data-view]").forEach((b) =>
      b.addEventListener("click", () => setView(b.dataset.view)));
  }

  function setView(v) {
    if (v === "fit") return fit();
    if (v === "x") { cam.theta = Math.PI / 2; cam.phi = 0; }
    else if (v === "y") { cam.theta = Math.PI / 2; cam.phi = Math.PI / 2; }
    else if (v === "z") { cam.theta = 0.01; cam.phi = 0; }
    else { cam.theta = 0.95; cam.phi = 0.9; }   // iso
  }

  function fit() {
    let rmax = 1;
    GT.state.scene.bodies.forEach((b) => {
      const e = b.radius || 0;
      const d = Math.hypot.apply(null, b.position || [0, 0, 0]);
      rmax = Math.max(rmax, e + d);
    });
    cam.radius = rmax * 3.2; cam.target.set(0, 0, 0);
  }

  // ---- построение тел из состояния ---- //
  function rotationMatrix(euler) {
    const [a, b, g] = euler || [0, 0, 0];
    const m = new THREE.Matrix4();
    const rz = (t) => new THREE.Matrix4().makeRotationZ(t);
    const ry = (t) => new THREE.Matrix4().makeRotationY(t);
    return m.multiply(rz(a)).multiply(ry(b)).multiply(rz(g)); // z-y-z (как в библиотеке)
  }

  function isHemisphere(b) {
    return b.type === "sphere" && b.hemisphere && b.hemisphere.enabled;
  }

  function geomFor(b) {
    if (isHemisphere(b)) {
      // купол вверх (+z): полусфера три-джей строится вокруг +y → поворот
      const g = new THREE.SphereGeometry(b.radius || 1, 48, 24, 0, Math.PI * 2, 0, Math.PI / 2);
      g.rotateX(Math.PI / 2);
      return g;
    }
    if (b.type === "sphere") return new THREE.SphereGeometry(b.radius || 1, 48, 32);
    if (b.type === "cylinder_inf") {
      const r = b.radius || 0.5;
      const g = new THREE.CylinderGeometry(r, r, 10 * r, 48, 1, true); // длинный, без крышек — «бесконечный»
      g.rotateX(Math.PI / 2);                                          // ось ∥ z
      return g;
    }
    return new THREE.SphereGeometry(0.5, 16, 12);
  }

  // Атрибуты полусферы на PEC-экране: основание, экран, точечный облучатель.
  function hemisphereExtras(b, sel) {
    const r = b.radius || 1;
    const g = new THREE.Group();
    // плоское основание купола (z=0)
    const base = new THREE.Mesh(
      new THREE.CircleGeometry(r, 48),
      new THREE.MeshStandardMaterial({ color: sel ? COL.sel : COL.body, metalness: 0.1,
        roughness: 0.55, transparent: true, opacity: 0.35, side: THREE.DoubleSide }));
    base.position.z = 0.001;
    g.add(base);
    // проводящий экран (PEC) — большой полупрозрачный диск в плоскости z=0
    const scr = new THREE.Mesh(
      new THREE.CircleGeometry(3 * r, 64),
      new THREE.MeshStandardMaterial({ color: 0x9aa4ad, metalness: 0.75, roughness: 0.35,
        transparent: true, opacity: 0.28, side: THREE.DoubleSide }));
    scr.position.z = -0.002;
    g.add(scr);
    const rim = new THREE.LineLoop(
      new THREE.BufferGeometry().setFromPoints(Array.from({ length: 96 }, (_, i) => {
        const a = (i / 96) * Math.PI * 2;
        return new THREE.Vector3(3 * r * Math.cos(a), 3 * r * Math.sin(a), -0.002);
      })),
      new THREE.LineBasicMaterial({ color: 0x9aa4ad, transparent: true, opacity: 0.6 }));
    g.add(rim);
    // точечный облучатель: красный конус на куполе под углом θ′ от нормали, остриём к центру
    const tp = ((b.hemisphere.feed_offset_deg || 0) * Math.PI) / 180;
    const dir = new THREE.Vector3(Math.sin(tp), 0, Math.cos(tp));
    const feed = new THREE.Mesh(
      new THREE.ConeGeometry(0.06 * r, 0.16 * r, 20),
      new THREE.MeshStandardMaterial({ color: 0xd23c3c, emissive: 0x551212,
        metalness: 0.2, roughness: 0.4 }));
    feed.position.copy(dir.clone().multiplyScalar(1.09 * r));
    feed.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.clone().negate());
    g.add(feed);
    const feedDot = new THREE.Mesh(
      new THREE.SphereGeometry(0.035 * r, 16, 12),
      new THREE.MeshStandardMaterial({ color: 0xd23c3c, emissive: 0x551212 }));
    feedDot.position.copy(dir.clone().multiplyScalar(1.19 * r));
    g.add(feedDot);
    return g;
  }

  function syncFromState() {
    if (!bodiesGroup) return;
    for (let i = bodiesGroup.children.length - 1; i >= 0; i--) {
      const c = bodiesGroup.children[i];
      c.geometry && c.geometry.dispose();
      bodiesGroup.remove(c);
    }
    GT.state.scene.bodies.forEach((b, i) => {
      const geom = geomFor(b);
      const sel = i === GT.state.selected;
      const mat = new THREE.MeshStandardMaterial({
        color: sel ? COL.sel : COL.body, metalness: 0.1, roughness: 0.55,
        transparent: true, opacity: sel ? 0.55 : 0.4,
      });
      const mesh = new THREE.Mesh(geom, mat);
      const edges = new THREE.LineSegments(
        new THREE.EdgesGeometry(geom, 18),
        new THREE.LineBasicMaterial({ color: sel ? COL.sel : COL.edge, transparent: true, opacity: 0.6 }));
      const wrap = new THREE.Group();
      wrap.add(mesh); wrap.add(edges);
      if (isHemisphere(b)) wrap.add(hemisphereExtras(b, sel));
      wrap.position.set.apply(wrap.position, b.position || [0, 0, 0]);
      if (b.type !== "sphere") wrap.setRotationFromMatrix(rotationMatrix(b.euler));
      bodiesGroup.add(wrap);
    });
    drawIncidence();
  }

  function drawIncidence() {
    if (incidenceArrow) { scene.remove(incidenceArrow); incidenceArrow = null; }
    const bodies = GT.state.scene.bodies;
    const ci = bodies.length === 1 && bodies[0].type === "cylinder_inf" ? bodies[0] : null;
    const hemi = bodies.length === 1 && isHemisphere(bodies[0]) ? bodies[0] : null;
    // у бесконечного цилиндра падение задаётся углом от оси z (90° = нормальное, ⊥ оси)
    let k = GT.state.scene.radiation.khat || [0, 0, 1];
    if (ci) {
      const th = (ci.theta_deg != null ? ci.theta_deg : 90) * Math.PI / 180;
      k = [Math.sin(th), 0, Math.cos(th)];
    }
    if (hemi) {
      // полусфера: освещение со стороны облучателя (θ′ от нормали), волна идёт к линзе
      const tp = ((hemi.hemisphere.feed_offset_deg || 0) * Math.PI) / 180;
      k = [-Math.sin(tp), 0, -Math.cos(tp)];
    }
    const dir = new THREE.Vector3(k[0], k[1], k[2]);
    if (dir.length() < 1e-9) return;
    dir.normalize();
    const start = dir.clone().multiplyScalar(-3.4);
    incidenceArrow = new THREE.ArrowHelper(dir, start, 2.2, 0xb23a3a, 0.45, 0.26);
    scene.add(incidenceArrow);
  }

  // ---- наложение ДН поверх модели (полярная кривая в плоскости среза xz) ---- //
  function clearOverlay() {
    if (!overlayGroup) return;
    overlayGroup.traverse((o) => {
      if (o.geometry) o.geometry.dispose();
      if (o.material) o.material.dispose();
    });
    scene.remove(overlayGroup);
    overlayGroup = null;
  }

  function syncOverlay() {
    if (!scene) return;
    clearOverlay();
    const chk = GT.$("#dn-overlay");
    const res = GT.state.results;
    if (!chk || !chk.checked || !res || !res.pattern || !res.pattern.polar) return;
    const pat = res.pattern;
    if (!pat.theta_deg || !pat.series || !pat.series.length) return;
    let R = 0;
    GT.state.scene.bodies.forEach((b) => { R = Math.max(R, b.radius || 1); });
    R = R || 1;
    const hemi = !!pat.hemisphere;
    // базис плоскости среза: полусфера — 0° = нормаль (+z); иначе 0° = направление падения
    let axis0, axis90;
    if (hemi) {
      axis0 = new THREE.Vector3(0, 0, 1); axis90 = new THREE.Vector3(1, 0, 0);
    } else {
      const kh = GT.state.scene.radiation.khat || [0, 0, 1];
      axis0 = new THREE.Vector3(kh[0], kh[1], kh[2]);
      if (axis0.length() < 1e-9) axis0.set(0, 0, 1);
      axis0.normalize();
      const pol = GT.state.scene.radiation.pol || [1, 0, 0];
      axis90 = new THREE.Vector3(pol[0], pol[1], pol[2])
        .addScaledVector(axis0, -new THREE.Vector3(pol[0], pol[1], pol[2]).dot(axis0));
      if (axis90.length() < 1e-9) axis90.set(1, 0, 0);
      axis90.normalize();
    }
    const FLOOR = -60;
    const dark = (document.documentElement.getAttribute("data-theme") || "dark") !== "light";
    const colors = dark ? [0x00ff9c, 0xff6b6b, 0x4dd2ff] : [0x0d6b8a, 0xb23a3a, 0x1f6fb2];
    const toPoint = (thetaDeg, dB) => {
      const t = (thetaDeg * Math.PI) / 180;
      const rr = R * (0.22 + 2.1 * Math.max(0, dB - FLOOR) / (0 - FLOOR));
      return new THREE.Vector3()
        .addScaledVector(axis0, rr * Math.cos(t))
        .addScaledVector(axis90, rr * Math.sin(t));
    };
    overlayGroup = new THREE.Group();
    pat.series.forEach((s, i) => {
      let entries = pat.theta_deg.map((td, j) => [td, s.dB[j]])
        .filter(([, dB]) => dB != null && isFinite(dB));
      let closed = true;
      if (hemi) {
        // только верхнее полупространство, разомкнутая дуга −90°..+90° через 0°
        entries = entries
          .map(([td, dB]) => [td > 180 ? td - 360 : td, dB])
          .filter(([td]) => Math.abs(td) <= 90)
          .sort((a, b) => a[0] - b[0]);
        closed = false;
      }
      if (!entries.length) return;
      const pts = entries.map(([td, dB]) => toPoint(td, dB));
      if (closed) pts.push(pts[0].clone());
      const line = new THREE.Line(
        new THREE.BufferGeometry().setFromPoints(pts),
        new THREE.LineBasicMaterial({ color: colors[i % colors.length],
          transparent: true, opacity: 0.95 }));
      overlayGroup.add(line);
    });
    scene.add(overlayGroup);
  }

  function resize() {
    const host = GT.$("#viewport");
    if (!host) return;
    const w = host.clientWidth || 1, h = host.clientHeight || 1;
    renderer.setSize(w, h, false);
    camera.aspect = w / h; camera.updateProjectionMatrix();
  }

  function loop() {
    const st = Math.sin(cam.theta);
    camera.position.set(
      cam.target.x + cam.radius * st * Math.cos(cam.phi),
      cam.target.y + cam.radius * st * Math.sin(cam.phi),
      cam.target.z + cam.radius * Math.cos(cam.theta));
    camera.up.set(0, 0, 1);
    camera.lookAt(cam.target);
    renderer.render(scene, camera);
    raf = requestAnimationFrame(loop);
  }

  function onShow() { resize(); }

  GT.editor = { init, syncFromState, setView, onShow, applyTheme, syncOverlay };
})();
