<!-- SPDX-License-Identifier: MIT -->

# Example 8 — RPM Coated PEC Sphere

Проводящая (PEC) сфера радиуса 100 мм с однослойным радиопоглощающим
покрытием ZIK 51-2190 (диэлектрическое, 4 мм) или ZIK 51-2191
(магнетодиэлектрическое, 2 мм); измеренные ε(f), μ(f) на частотах
3 / 5,5 / 7,5 / 10 ГГц. Вычисляется эффективность подавления обратного
рассеяния относительно голой сферы.

PEC sphere with a single-layer radar-absorbing coating (lossy
magnetodielectric). For strongly absorbing layers the computation
automatically switches to the arbitrary-precision branch of
`green_tensor.multilayer_mie` (mpmath, 60 digits) once
|Im(m·k0·r)| > 8 — double precision loses ~1 dB in Q_back by
|Im z| ≈ 18 and collapses (NaN) beyond ≈ 30.

Reference values: ZIK 51-2190 gives +4.4 dB suppression at 3 GHz fading
to −2.8 dB at 10 GHz; ZIK 51-2191 produces an RCS *enhancement* of
13.2/13.4/8.4/6.9 dB at 3/5.5/7.5/10 GHz.

## Запуск / Run

```bash
python3 zik_rpm_coatings.py
```

Зависимости: numpy, scipy; mpmath — для сильно поглощающих покрытий.
