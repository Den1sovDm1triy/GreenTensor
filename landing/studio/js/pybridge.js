// SPDX-License-Identifier: MIT
// GreenTensor Studio — браузерный мост к движку: Pyodide вместо серверного /api/compute.
// Загружает Pyodide (numpy+scipy), пишет пакет green_tensor + studio в виртуальную ФС,
// затем выполняет studio.compute.compute(scene) прямо в браузере. Сервер не нужен.
(function () {
  "use strict";
  const GT = (window.GT = window.GT || {});

  const ENGINE_BASE = "studio/engine/";
  const PYODIDE_CDN = "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/";

  let pyodide = null;
  let readyPromise = null;

  function setBoot(text, done) {
    const el = document.getElementById("boot-status");
    if (el) el.textContent = text;
    const wrap = document.getElementById("boot-overlay");
    if (wrap && done) wrap.classList.add("hidden");
  }

  // Стаб matplotlib: ядро 01_sphere.py импортирует его лениво, но в браузере он не нужен
  // (никакой отрисовки на стороне Python — графики строит Plotly). Ставится ДО первого compute.
  const MPL_STUB = `
import sys, types
for _m in ('matplotlib', 'matplotlib.pyplot', 'matplotlib.ticker', 'matplotlib.colors'):
    if _m not in sys.modules:
        sys.modules[_m] = types.ModuleType(_m)
def _noop(*a, **k):
    return None
for _name in ('plot','show','figure','subplots','savefig','close','xlabel','ylabel',
              'title','legend','grid','polar','tight_layout','xlim','ylim'):
    setattr(sys.modules['matplotlib.pyplot'], _name, _noop)
`;

  async function loadEngineIntoFS(py) {
    const manifest = await fetch(ENGINE_BASE + "manifest.json").then((r) => r.json());
    const files = manifest.files || [];
    // создаём каталоги пакетов и пишем каждый .py в виртуальную ФС Pyodide
    const dirs = new Set();
    files.forEach((f) => {
      const parts = f.split("/");
      let acc = "";
      for (let i = 0; i < parts.length - 1; i++) {
        acc += (i ? "/" : "") + parts[i];
        dirs.add(acc);
      }
    });
    dirs.forEach((d) => {
      try { py.FS.mkdir(d); } catch (e) { /* уже существует */ }
    });
    let n = 0;
    for (const f of files) {
      const src = await fetch(ENGINE_BASE + f).then((r) => {
        if (!r.ok) throw new Error("не найден " + f);
        return r.text();
      });
      py.FS.writeFile(f, src);
      setBoot("Загрузка движка… " + (++n) + "/" + files.length);
    }
    // корень ФС (где лежат green_tensor/ и studio/) — в путь импорта
    py.runPython("import sys\nif '' not in sys.path: sys.path.insert(0, '')");
  }

  async function boot() {
    setBoot("Загрузка среды Python (Pyodide)…");
    // pyodide.js подключается тегом <script> в studio.html → глобальная loadPyodide
    pyodide = await window.loadPyodide({ indexURL: PYODIDE_CDN });
    setBoot("Загрузка numpy и scipy…");
    await pyodide.loadPackage(["numpy", "scipy"]);
    setBoot("Стаб matplotlib…");
    pyodide.runPython(MPL_STUB);
    await loadEngineIntoFS(pyodide);
    setBoot("Инициализация движка…");
    // прогрев: импорт compute (тянет green_tensor); ошибки импорта видны сразу
    pyodide.runPython("from studio.compute import compute, ComputeError");
    setBoot("Движок готов.", true);
  }

  GT.py = {
    // единый промис готовности (идемпотентно)
    ready() {
      if (!readyPromise) readyPromise = boot();
      return readyPromise;
    },
    isReady() {
      return !!pyodide;
    },
    // Расчёт сцены. scene — обычный JS-объект; возвращает результат-объект (как /api/compute).
    async compute(scene) {
      await this.ready();
      pyodide.globals.set("_scene_json", JSON.stringify(scene));
      const out = pyodide.runPython(`
import json
from studio.compute import compute, ComputeError
try:
    _res = compute(json.loads(_scene_json))
    _out = json.dumps(_res, ensure_ascii=False)
except ComputeError as _e:
    _out = json.dumps({"error": str(_e)}, ensure_ascii=False)
except Exception as _e:
    import traceback; traceback.print_exc()
    _out = json.dumps({"error": "Внутренняя ошибка расчёта: %s" % _e}, ensure_ascii=False)
_out
`);
      return JSON.parse(out);
    },
  };

  // Старт загрузки сразу при готовности DOM (чтобы Pyodide грузился, пока юзер настраивает сцену)
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", () => GT.py.ready());
  } else {
    GT.py.ready();
  }
})();
