"""
GreenTensor — open-source Python library for electromagnetic scattering
analysis on multilayer spherical structures based on Tensor Green's
Functions (TGF).

This module is the public entry point. The two workhorse classes are
re-exported here so that user code can simply do::

    from greentensor import RCSCalculator, ScatteringCalculator

Both classes live in :mod:`green_tensor.calc`. Single-frequency RCS
analysis is performed via :class:`RCSCalculator`; parametric sweep over
electrical size :math:`k_0 a` is performed via :class:`ScatteringCalculator`.

Authors:
    * Dmitriy V. Denisov   (UrFU, dv.denisov@urfu.ru)
    * Vladislav Ya. Noskov (UrFU, v.y.noskov@urfu.ru)
    * Ilia O. Skumatenko   (UrFU, ilya.skumatenko@urfu.ru)

References:
    [1] Panchenko B.A. Scattering and Absorption of EM Waves by
        Inhomogeneous Spherical Bodies. Moscow: Radiotekhnika, 2013.
    [2] Denisov D.V., Skumatenko I.O., Fadeev V.O., Shesterov M.A.
        GreenTensor Electromagnetic Library for Calculate Scattering on
        Multilayer Spherical Structures. IEEE EDM 2024.
        DOI 10.1109/EDM61683.2024.10615009
    [3] Denisov D.V., Skumatenko I.O., Shesterov M.A. GreenTensor
        Electromagnetic Library for Calculating Scattering Patterns on
        a Human Head Model. IEEE USBEREIT 2024.
        DOI 10.1109/USBEREIT61901.2024.10584060
"""
from green_tensor.calc import RCSCalculator, ScatteringCalculator

__version__ = "0.2.0"
__all__ = ["RCSCalculator", "ScatteringCalculator"]
