# Регенерация presets.json

`../presets.json` — статический дамп эталонных задач Studio (сцены + эталонные
кривые Ansys HFSS), потребляемый браузерным расчётчиком (`studio/js/presets.js`).

Источник: `presets.py` (+ `data/` — examples.json и HFSS-CSV). Пересобрать:

```bash
cd greentenor_landing/studio/engine
PYTHONPATH="$PWD:_source" python3 - <<'PY'
import sys, json, pathlib
sys.path.insert(0, "_source")
import presets                      # _source/presets.py
pathlib.Path("presets.json").write_text(
    json.dumps({"presets": presets.list_presets()}, ensure_ascii=False), encoding="utf-8")
print("presets.json regenerated:", len(presets.list_presets()), "presets")
PY
```

`presets.py` использует публичный API `green_tensor` только для сцен — при импорте
достаточно, чтобы пакет был доступен (он лежит рядом в `../green_tensor`).
Провенанс кривых HFSS — см. комментарии в `presets.py` и `data/examples.json`.
