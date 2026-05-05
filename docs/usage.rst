Usage
=====

.. _installation:

Installation
------------

Clone the repository and install the dependencies:

.. code-block:: bash

   git clone https://github.com/Den1sovDm1triy/GreenTensor.git
   cd GreenTensor
   pip install -r requirements.txt

GreenTensor requires **Python 3.9 or newer**. The runtime dependencies
are NumPy, SciPy and Matplotlib.

To verify the installation, run a one-line import test:

.. code-block:: bash

   python -c "from greentensor import RCSCalculator, ScatteringCalculator; print('OK')"


Quick start
-----------

GreenTensor exposes two main classes from the public entry point
``greentensor``:

* :class:`green_tensor.calc.RCSCalculator` — single-frequency analysis of
  bistatic radar cross-section (RCS) for a multilayer sphere.
* :class:`green_tensor.calc.ScatteringCalculator` — parametric sweep of
  integrated cross-sections over the electrical size :math:`k_0 a`.

Single-frequency RCS calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compute the bistatic RCS pattern of a four-layer Luneburg lens at
:math:`k_0 a = 2\pi` (i.e. :math:`a = \lambda`):

.. code-block:: python

   from greentensor import RCSCalculator
   import math

   calc = RCSCalculator(
       k0=2 * math.pi,
       toch=10,                     # series truncation
       n=4,                         # number of layers (last = air)
       a=[0.25, 0.5, 0.75, 1.0],    # normalized layer radii
       eps=[1.94, 1.75, 1.44, 1.0], # eps_r per layer
       miy=[1, 1, 1, 1],            # mu_r per layer
   )
   calc.run_calculation()
   calc.plot_results()              # 4-panel matplotlib figure

After ``run_calculation()`` the following arrays are populated on the
calculator instance:

* ``Mn``, ``Nn`` — modal coefficients :math:`M_n,\,N_n`
* ``E_op``, ``E_kp`` — main and cross polarization (circular case)
* ``E_teta``, ``E_phi`` — field components (linear case)
* ``DN_NORM_lin_dB_*``, ``DN_NORM_circle_dB_*`` — normalized patterns in dB

Parametric sweep over electrical size
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compute the integrated scattering cross-section :math:`\sigma_s(k_0 a)`
and the radar cross-section :math:`\sigma_r(k_0 a)` for a default sphere:

.. code-block:: python

   from greentensor import ScatteringCalculator
   import matplotlib.pyplot as plt

   sweep = ScatteringCalculator(
       toch=15,
       k0a_start=0.25,
       k0a_stop=10.0,
       k0a_step=0.05,
   )
   k0a, sigma_s, sigma_r, sigma_p, sigma_theta, sigma_phi = sweep.calculate()
   sweep.plot_results()


Reference datasets
------------------

The file ``green_tensor/examples.json`` ships canonical parameter sets
together with reference scattering patterns produced in Ansys HFSS:

* ``DN_NORM1`` — Luneburg lens in vacuum, 5 equispaced layers, :math:`k_0 a = 10\pi`.
* ``DN_NORM4`` — metallic spheres of 2 and 5 homogeneous shells, :math:`k_0 a = 4\pi`.
* ``DN_NORM5`` — metallic sphere with a dielectric coating, :math:`k_0 a = 4\pi`.

Worked examples are in the ``examples/`` folder:

* ``examples/Example 1 — Luneburg Lens Bistatic RCS`` reproduces Fig. 11
  of Greenwood & Jin (1999).
* ``examples/Example 2 — Luneburg Lens Scattering Diagram`` provides a
  ready-to-open Ansys Electronics Desktop project
  (``LuneburgLens_Layer4_ravnoshagApprox.aedtz``) and CSV exports of the
  HFSS pattern (``Ephi=0_from_HFSS.csv``, ``Ephi=90_from_HFSS.csv``) for
  side-by-side validation against GreenTensor.
* ``examples/Example 3 — Luneburg Lens for UAV Scattering Diagram``
  computes the diagram for a UAV-borne Luneburg reflector.


References
----------

#. Panchenko B.A. *Scattering and Absorption of Electromagnetic Waves by
   Inhomogeneous Spherical Bodies.* Moscow: Radiotekhnika, 2013.
#. Denisov D.V. *Antenna and Diffraction Characteristics of Luneburg
   Lenses under Circular Polarization.* PhD thesis, NNGTU, 2015.
#. Denisov D.V., Skumatenko I.O., Fadeev V.O., Shesterov M.A.
   *GreenTensor Electromagnetic Library for Calculate Scattering on
   Multilayer Spherical Structures.* IEEE EDM 2024,
   `DOI 10.1109/EDM61683.2024.10615009
   <https://doi.org/10.1109/EDM61683.2024.10615009>`_.
#. Denisov D.V., Skumatenko I.O., Shesterov M.A. *GreenTensor
   Electromagnetic Library for Calculating Scattering Patterns on a
   Human Head Model.* IEEE USBEREIT 2024,
   `DOI 10.1109/USBEREIT61901.2024.10584060
   <https://doi.org/10.1109/USBEREIT61901.2024.10584060>`_.
