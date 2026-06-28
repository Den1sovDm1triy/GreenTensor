# SPDX-License-Identifier: MIT

"""Локальный сервер веб-приложения GreenTensor (только stdlib, без зависимостей).

Запуск:
    python3 webapp/server.py [--port 8765]

Открыть в браузере: http://127.0.0.1:8765
"""
from __future__ import annotations

import argparse
import json
import os
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

from compute import ComputeError, compute

STATIC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")
MAX_BODY = 1_000_000  # 1 МБ — с запасом для формы параметров


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=STATIC_DIR, **kwargs)

    def log_message(self, fmt, *args):  # компактный лог
        print(f"[web] {self.address_string()} {fmt % args}")

    def _json(self, code: int, obj: dict):
        body = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        if self.path != "/api/compute":
            self._json(404, {"error": "Неизвестный маршрут."})
            return
        try:
            length = int(self.headers.get("Content-Length", 0))
        except ValueError:
            length = 0
        if not (0 < length <= MAX_BODY):
            self._json(400, {"error": "Пустое или слишком большое тело запроса."})
            return
        try:
            payload = json.loads(self.rfile.read(length).decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError):
            self._json(400, {"error": "Некорректный JSON."})
            return
        try:
            self._json(200, compute(payload))
        except ComputeError as exc:
            self._json(400, {"error": str(exc)})
        except Exception as exc:  # noqa: BLE001 — сообщение в UI, трейс в консоль
            import traceback
            traceback.print_exc()
            self._json(500, {"error": f"Внутренняя ошибка расчёта: {exc}"})


def main():
    parser = argparse.ArgumentParser(description="GreenTensor web app")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    srv = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"GreenTensor web: http://127.0.0.1:{args.port}  (Ctrl+C — остановить)")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
