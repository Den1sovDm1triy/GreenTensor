.. SPDX-License-Identifier: MIT
.. Scientific scope: scientific research and engineering modeling in classical electrodynamics, antenna theory, microwave devices, and electromagnetic scattering.

GreenTensor
===========

**GreenTensor** — библиотека аналитического анализа рассеяния электромагнитных
волн на гетерогенных структурах. Математическим ядром служит точное решение для
радиально-слоистой сферы методом тензорных функций Грина (ТФГ / Ми). Поверх ядра
построено семейство аналитических решателей (сфероид, эллипсоид, цилиндр, конус) с
единым T-матричным интерфейсом в сферическом базисе векторных волновых функций и
единый движок сборки сложной геометрии (Generalized Multiparticle Mie, GMM).

*EN.* **GreenTensor** is a library for analytic electromagnetic-scattering analysis
on heterogeneous structures. Its mathematical core is the exact radially-layered-sphere
solution by the tensor Green's function (TGF / Mie) method. On top of the core, a
family of analytic solvers (spheroid, ellipsoid, cylinder, cone) shares one T-matrix
interface in the spherical vector-wave-function basis, composed by a single
complex-geometry assembly engine (Generalized Multiparticle Mie, GMM).

The library underpins two IEEE publications:

* `GreenTensor: Efficient Electromagnetic Scattering Analysis on Heterogeneous
  Spherical Structures <https://ieeexplore.ieee.org/document/10615009>`_
* `GreenTensor: A Python Library for Electromagnetic Wave Propagation Analysis
  <https://ieeexplore.ieee.org/document/10584060>`_

The full mathematical derivation is in ``GreenTensor_Theory.tex`` at the repository root.

.. note::

   Точные ветви (сфера, бесконечный цилиндр) и квазистатические решатели
   (эллипсоид, сфероид) проверены на машинную точность против независимых
   эталонов. Полноволновые ветви сфероида/конечного цилиндра/конуса честно
   поднимают :class:`NotImplementedError` — заглушек нет.

Contents
--------

.. toctree::
   :maxdepth: 2

   intro
   usage
   solvers
   api
