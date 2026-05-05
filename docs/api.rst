API Reference
=============

GreenTensor exposes a class-based API in :mod:`green_tensor.calc`. The
public entry point :mod:`greentensor` re-exports the two main classes for
convenience.

Public entry point
------------------

.. automodule:: greentensor
   :members:
   :undoc-members:
   :show-inheritance:

Class-based API (single-frequency RCS, parametric sweep)
--------------------------------------------------------

.. automodule:: green_tensor.calc
   :members: RCSCalculator, ScatteringCalculator
   :undoc-members:
   :show-inheritance:


Legacy monolithic scripts
-------------------------

The following scripts in :mod:`green_tensor` are preserved as historical
references — they were the source from which the class-based API in
``calc.py`` was refactored. New user code should prefer
:class:`green_tensor.calc.RCSCalculator` and
:class:`green_tensor.calc.ScatteringCalculator`.

* ``green_tensor.lin_polar`` — diffraction at linear polarization,
  superseded by :class:`RCSCalculator`.
* ``green_tensor.Bistatic_RCS_lin_polar`` — bistatic RCS, linear polarization.
* ``green_tensor.Bistatic_RCS_lin_circle_polar`` (file is named
  ``Bistatic_RCS_lin+circle_polar.py``) — bistatic RCS, both polarizations.


Reference datasets
------------------

``green_tensor/examples.json`` is a multi-document JSON file with parameter
sets and reference scattering patterns from Ansys HFSS. See
:doc:`usage` for description of the available reference cases.
