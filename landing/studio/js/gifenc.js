// SPDX-License-Identifier: MIT
// GreenTensor Studio — минимальный GIF89a-энкодер (LZW) без зависимостей.
// Кадры принимаются как ImageData (RGBA); палитра ГЛОБАЛЬНАЯ, 256 цветов, задаётся
// снаружи (в Studio — LUT палитры поля + цвета темы/оверлея), пиксели квантуются
// по ближайшему цвету палитры с кэшем. Цикл — NETSCAPE2.0 (loop=0 — бесконечно).
(function () {
  "use strict";
  const GT = (window.GT = window.GT || {});

  // ---- LZW (GIF-вариант, без early-change), результат — массив байт кодового потока ----
  function lzwEncode(pixels, minCodeSize) {
    const CLEAR = 1 << minCodeSize, EOI = CLEAR + 1;
    let codeSize = minCodeSize + 1, next = EOI + 1;
    const out = [];
    let bitBuf = 0, bitCnt = 0;
    const emit = (code) => {
      bitBuf |= code << bitCnt; bitCnt += codeSize;
      while (bitCnt >= 8) { out.push(bitBuf & 255); bitBuf >>>= 8; bitCnt -= 8; }
    };
    const table = new Int32Array(4096 * 256);        // ключ = (префикс<<8)|байт
    table.fill(-1);
    emit(CLEAR);
    let prefix = pixels[0];
    for (let i = 1; i < pixels.length; i++) {
      const c = pixels[i];
      const key = (prefix << 8) | c;
      const e = table[key];
      if (e >= 0) { prefix = e; continue; }
      emit(prefix);
      if (next < 4096) {
        // разрядность растёт ПЕРЕД добавлением записи 2^codeSize (декодер отстаёт на одну
        // запись: он заводит строку при чтении СЛЕДУЮЩЕГО кода) — как в классическом compress
        if (next >= (1 << codeSize) && codeSize < 12) codeSize++;
        table[key] = next++;
      } else {                                         // словарь полон → сброс
        emit(CLEAR); table.fill(-1); next = EOI + 1; codeSize = minCodeSize + 1;
      }
      prefix = c;
    }
    emit(prefix); emit(EOI);
    if (bitCnt > 0) out.push(bitBuf & 255);
    return out;
  }

  function GIFEncoder(width, height, palette, opts) {
    // palette — массив [r,g,b] (≤256); opts: {loop:0, delay:80 (мс)}
    this.w = width; this.h = height;
    this.opts = Object.assign({ loop: 0, delay: 80 }, opts || {});
    this.pal = palette.slice(0, 256);
    while (this.pal.length < 256) this.pal.push([0, 0, 0]);
    this.cache = new Map();                            // rgb(24 бит) → индекс палитры
    this.chunks = [];
    this._header();
  }

  GIFEncoder.prototype._u16 = function (v) { return [v & 255, (v >> 8) & 255]; };

  GIFEncoder.prototype._header = function () {
    const b = [0x47, 0x49, 0x46, 0x38, 0x39, 0x61];   // "GIF89a"
    b.push(...this._u16(this.w), ...this._u16(this.h));
    b.push(0xF7, 0, 0);                                // GCT 256 цветов, 8 бит/канал
    this.pal.forEach((c) => b.push(c[0] & 255, c[1] & 255, c[2] & 255));
    // NETSCAPE2.0 — число повторов (0 = бесконечно)
    b.push(0x21, 0xFF, 0x0B, 0x4E, 0x45, 0x54, 0x53, 0x43, 0x41, 0x50, 0x45, 0x32, 0x2E, 0x30,
      0x03, 0x01, ...this._u16(this.opts.loop | 0), 0x00);
    this.chunks.push(Uint8Array.from(b));
  };

  // индекс ближайшего цвета палитры (евклид в RGB) с кэшем
  GIFEncoder.prototype.nearest = function (r, g, b) {
    const key = (r << 16) | (g << 8) | b;
    let idx = this.cache.get(key);
    if (idx !== undefined) return idx;
    let best = 0, bd = 1e9;
    const P = this.pal;
    for (let i = 0; i < P.length; i++) {
      const dr = P[i][0] - r, dg = P[i][1] - g, db = P[i][2] - b;
      const d = dr * dr + dg * dg + db * db;
      if (d < bd) { bd = d; best = i; if (d === 0) break; }
    }
    this.cache.set(key, best);
    return best;
  };

  // Кадр: ImageData (RGBA, размер w×h) либо готовый Uint8Array индексов.
  GIFEncoder.prototype.addFrame = function (frame, delayMs) {
    const n = this.w * this.h;
    let idx;
    if (frame instanceof Uint8Array && frame.length === n) {
      idx = frame;
    } else {
      const d = frame.data;
      idx = new Uint8Array(n);
      for (let i = 0, p = 0; i < n; i++, p += 4) idx[i] = this.nearest(d[p], d[p + 1], d[p + 2]);
    }
    const delay = Math.max(1, Math.round((delayMs == null ? this.opts.delay : delayMs) / 10));
    const b = [0x21, 0xF9, 0x04, 0x00, ...this._u16(delay), 0x00, 0x00];   // GCE, без прозрачности
    b.push(0x2C, 0, 0, 0, 0, ...this._u16(this.w), ...this._u16(this.h), 0x00); // дескриптор
    b.push(0x08);                                                              // min code size
    const stream = lzwEncode(idx, 8);
    for (let i = 0; i < stream.length; i += 255) {
      const len = Math.min(255, stream.length - i);
      b.push(len);
      for (let j = 0; j < len; j++) b.push(stream[i + j]);
    }
    b.push(0x00);
    this.chunks.push(Uint8Array.from(b));
  };

  GIFEncoder.prototype.finish = function () {
    this.chunks.push(Uint8Array.from([0x3B]));
    return new Blob(this.chunks, { type: "image/gif" });
  };

  GT.GIFEncoder = GIFEncoder;
  GT.lzwEncode = lzwEncode;
})();
