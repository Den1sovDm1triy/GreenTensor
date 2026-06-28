.. SPDX-License-Identifier: MIT
.. Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

GreenTensor
===========

**GreenTensor** — библиотека аналитического анализа рассеяния электромагнитных
волн на гетерогенных структурах. Область охвата — только **точные аналитические
решения** метода тензорных функций Грина (ТФГ): радиально-слоистая сфера (Ми/ТФГ),
бесконечный слоистый/косой цилиндр (точная 2D-аналитика) и их строгая суперпозиция —
кластер невзаимопересекающихся сфер (Generalized Multiparticle Mie, GMM, теорема
сложения Крузана–Стейна). Приближённые и чисто численные методы (EBCM для тел с
рёбрами, квазистатика Рэлея, разложение тел в кластер сфер) в библиотеку намеренно
НЕ входят — только верифицируемая аналитика.

*EN.* **GreenTensor** is a library for analytic electromagnetic-scattering analysis
on heterogeneous structures. Its scope is restricted to **exact analytic** tensor
Green's function (TGF) solutions: the radially-layered sphere (Mie/TGF), the infinite
layered/oblique cylinder (exact 2D analytics) and their rigorous superposition — a
cluster of non-overlapping spheres (Generalized Multiparticle Mie, GMM, Cruzan–Stein
addition theorem). Approximate and purely numerical methods (EBCM for edged bodies,
Rayleigh quasi-statics, decomposition of bodies into a sphere cluster) are
intentionally NOT included — only verifiable analytics.

The library underpins two IEEE publications:

* `GreenTensor: Efficient Electromagnetic Scattering Analysis on Heterogeneous
  Spherical Structures <https://ieeexplore.ieee.org/document/10615009>`_
* `GreenTensor: A Python Library for Electromagnetic Wave Propagation Analysis
  <https://ieeexplore.ieee.org/document/10584060>`_

The full mathematical derivation is in ``GreenTensor_Theory.tex`` at the repository root.

.. note::

   В библиотеке остаются только точные аналитические ТФГ-решения (сфера,
   бесконечный цилиндр) и их строгая суперпозиция (кластер невзаимопересекающихся
   сфер, GMM) — все проверены на машинную точность против независимых эталонов и
   законов сохранения/взаимности. Приближённые/численные методы исключены из
   области охвата.

Contents
--------

.. toctree::
   :maxdepth: 2

   intro
   usage
   solvers
   api
