The Conformal :math:`(m', \theta)` Design Plane
================================================

Motivation
----------

A turbomachinery blade surface is a curved surface of revolution.  Designing
a blade profile directly in three-dimensional Cartesian space is cumbersome
because the physical angle between two curves on the surface depends on the
local radius :math:`r`.  The conformal coordinate system eliminates this
radius dependence: in the :math:`(m', \theta)` plane, the angle between any
two curves is identical to the physical angle those curves subtend on the 3-D
surface.

Meridional Geometry
-------------------

Consider a body of revolution with axis of symmetry :math:`x`.  Any point on
the blade surface can be specified by:

- :math:`x` — axial coordinate
- :math:`r` — radial coordinate (distance from the axis)
- :math:`\theta` — circumferential angle [rad]

The **meridional arc length** element along a streamline in the :math:`(x,r)`
half-plane is

.. math::
   \mathrm{d}m = \sqrt{\mathrm{d}x^2 + \mathrm{d}r^2}.

The metric on the blade surface expressed in :math:`(m, \theta)` coordinates
is *not* flat:

.. math::
   \mathrm{d}s^2 = \mathrm{d}m^2 + r^2 \, \mathrm{d}\theta^2.

The factor :math:`r^2` means that equal increments of :math:`\theta` cover
different arc lengths at different radii.

The Conformal Coordinate :math:`m'`
-------------------------------------

The conformal meridional coordinate :math:`m'` is defined by the line
integral

.. math::
   m'(s) = \int_0^s \frac{\mathrm{d}m}{r(m)},

where the integration runs along the construction streamline from the leading
edge (:math:`s = 0`) to the trailing edge.  Substituting this change of
variables transforms the surface metric to

.. math::
   \mathrm{d}s^2 = r^2 \left( \mathrm{d}m'^{\,2} + \mathrm{d}\theta^2 \right).

The factor :math:`r^2` is now a *conformal factor* that multiplies an
isotropic Euclidean metric.  Because conformal maps preserve angles, the
angle between any two curves in the flat :math:`(m', \theta)` plane equals
the physical angle between those curves on the 3-D surface.

This is the key property that makes the design plane useful:

.. admonition:: Design-plane angle correspondence

   The inlet metal angle :math:`\beta_1`, the exit metal angle
   :math:`\beta_2`, and the blade-to-blade flow angles are *all* measured
   in the :math:`(m', \theta)` plane with no correction for radius.

Normalisation
-------------

In practice :math:`m'` is normalised by its total value at the trailing edge:

.. math::
   \bar{m}' = \frac{m'}{m'_{\max}}, \qquad \bar{m}' \in [0, 1].

The normalised coordinate is what BladeGen uses internally;
:math:`\bar{m}' = 0` at the leading edge and :math:`\bar{m}' = 1` at the
trailing edge.

Numerical Integration via RK4
------------------------------

The construction streamline is represented as a B-spline curve in the
:math:`(x, r)` half-plane parameterised by :math:`t \in [0,1]`.  The
ordinary differential equation for :math:`m'(t)` is

.. math::
   \frac{\mathrm{d}m'}{\mathrm{d}t} = \frac{\|\dot{\mathbf{x}}(t)\|}{r(t)},
   \qquad m'(0) = 0,

where :math:`\|\dot{\mathbf{x}}\| = \sqrt{\dot{x}^2 + \dot{r}^2}` is the
speed of the B-spline parameterisation.  Because the right-hand side is
independent of :math:`m'` (it depends only on :math:`t` through the
streamline geometry), the equation is pure quadrature.  The classical
fixed-step RK4 scheme achieves 4th-order convergence; 200 steps yields a
relative error below :math:`10^{-8}` for smooth turbine meridional channels.

The simultaneous integration of the physical arc length :math:`m(t)` uses

.. math::
   \frac{\mathrm{d}m}{\mathrm{d}t} = \|\dot{\mathbf{x}}(t)\|,

sharing the same RK4 loop at no additional cost.

Inverse Map: :math:`m' \rightarrow t \rightarrow (x, r)`
----------------------------------------------------------

Given a target design-plane point :math:`(\bar{m}'^*, \theta^*)`, the 3-D
Cartesian position is recovered in two steps:

1. **Invert** :math:`m'(t)` to find :math:`t^*` such that
   :math:`m'(t^*) = \bar{m}'^* \cdot m'_{\max}`.  This is done by binary
   search followed by linear interpolation in the pre-built tabulation of
   :math:`m'(t_i)`.

2. **Evaluate** the streamline: :math:`(x, r) = \mathbf{x}(t^*)`, then
   convert to Cartesian:

.. math::
   (x,\, y,\, z) = \bigl(x,\; r\sin\theta,\; r\cos\theta\bigr).

The binary search uses primal (``double``) values for index selection so that
no branch enters the CppAD automatic-differentiation tape, while the linear
interpolation is carried out with AD-typed :math:`T` values so that gradients
flow through the result.

Implementation
--------------

The relevant source files are:

- :file:`src/blade/StreamlineMPrime.hpp` — ``BSplineCurve2D<T>``,
  ``integrate_mprime()``, ``invert_mprime()``, ``conformal_to_3d()``
- :file:`src/blade/BSplineBasis.hpp` — Cox–de Boor evaluation and derivative
  used by the streamline B-spline
