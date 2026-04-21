Theory of Parametric Blade Generation
======================================

BladeGen implements the parametric turbomachinery blade design method
described by Koller *et al.* [1]_ and further developed at DLR [2]_.
The central idea is to perform all aerodynamic design work in a *conformal
design plane* — a flat Euclidean coordinate system in which angles between
curves exactly equal the physical flow angles on the three-dimensional blade
surface.

The pipeline from aerodynamic intent to OCCT solid consists of three
conceptually independent layers:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Layer
     - Responsibility
   * - :doc:`coordinate_system`
     - Defines the conformal :math:`(m', \theta)` design plane and its
       relationship to the meridional streamline
   * - :doc:`profile`
     - Parameterises the 2-D blade section as a B-spline camber line with
       a superimposed thickness distribution
   * - :doc:`blade_assembly`
     - Stacks 2-D sections along the span and maps them to 3-D Cartesian
       coordinates via the inverse conformal map

.. toctree::
   :maxdepth: 2
   :hidden:

   coordinate_system
   profile
   blade_assembly

----

.. rubric:: References

.. [1] Köller, U., Mönig, R., Küsters, B., and Schreiber, H. A. (1999).
       *Development of Advanced Compressor Airfoils for Heavy-Duty Gas
       Turbines — Part I: Design and Optimization.*
       ASME Journal of Turbomachinery, 121(3), 397–405.

.. [2] Aulich, M., and Siller, U. (2011).
       *High-Dimensional Constrained Multiobjective Optimization of a
       Fan Stage.* ASME Paper GT2011-45415.
