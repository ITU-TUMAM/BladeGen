2-D Blade Profile Parameterisation
===================================

All profile design work takes place in the flat :math:`(\bar{m}', \theta)`
design plane (see :doc:`coordinate_system`).  The profile is split into two
independent B-spline curves:

1. A **camber line** :math:`\theta_c(\bar{m}')` that encodes the aerodynamic
   turning, and
2. A **thickness distribution** :math:`\Delta\theta(\bar{m}')` that is added
   symmetrically on either side of the camber line to form the suction and
   pressure surfaces.

Primary Design Variables
------------------------

The table below lists all user-facing parameters of
``TurbineProfileParams<T>`` and their physical meaning.

.. list-table::
   :widths: 20 15 65
   :header-rows: 1

   * - Parameter
     - Units
     - Description
   * - :math:`\beta_1`
     - rad
     - **Inlet metal angle** — tangent of the camber line at the leading
       edge with respect to the :math:`\bar{m}'` axis.  Positive values
       lean toward the direction of rotation.  Typical axial-turbine range:
       40°–70°.
   * - :math:`\beta_2`
     - rad
     - **Exit metal angle** — tangent of the camber line at the trailing
       edge.  Typically negative (counter-rotation direction): −60° to −75°.
   * - :math:`\theta_\mathrm{LE}`
     - rad
     - Circumferential position of the leading edge (:math:`\bar{m}' = 0`).
   * - :math:`\theta_\mathrm{TE}`
     - rad
     - Circumferential position of the trailing edge (:math:`\bar{m}' = 1`).
       The total camber angle is :math:`\theta_\mathrm{TE} - \theta_\mathrm{LE}`.
   * - :math:`t_\mathrm{max}`
     - :math:`\theta`-units
     - Maximum thickness (in the design plane).
   * - :math:`t_\mathrm{max,loc}`
     - —
     - Chordwise position of maximum thickness, normalised to :math:`[0,1]`.
       Typical value: 0.30–0.40.
   * - :math:`t_\mathrm{TE}`
     - :math:`\theta`-units
     - Trailing-edge thickness (0 = sharp TE).
   * - :math:`\psi_\mathrm{TE}`
     - rad
     - TE wedge half-angle; sets
       :math:`\mathrm{d}\Delta\theta/\mathrm{d}\bar{m}'|_1 = -2\tan(\psi_\mathrm{TE}/2)`.
   * - :math:`r_\mathrm{LE}`
     - :math:`\theta`-units
     - Leading-edge radius (reserved for future blending; currently stored
       but not yet used to modify the thickness control points).

Camber Line
-----------

**Curve family.**  A clamped cubic B-spline on :math:`\bar{m}' \in [0, 1]`
with 4 control points and knot vector
:math:`\{0, 0, 0, 0, 1, 1, 1, 1\}`.

**Endpoint tangent conditions.**  For a clamped B-spline of degree :math:`p`
with :math:`n` control points :math:`Q_0, \ldots, Q_{n-1}`, the standard
endpoint tangent formulas are

.. math::
   \left. \frac{\mathrm{d}\theta_c}{\mathrm{d}\bar{m}'} \right|_0
   = p \cdot \frac{Q_1 - Q_0}{\xi_{p+1} - \xi_1},
   \qquad
   \left. \frac{\mathrm{d}\theta_c}{\mathrm{d}\bar{m}'} \right|_1
   = p \cdot \frac{Q_{n-1} - Q_{n-2}}{\xi_n - \xi_{n-p-1}},

where :math:`\xi_i` denotes the :math:`i`-th knot.  With the clamped cubic
knot vector all interior knot spans equal 1, so the formulas simplify to

.. math::
   \left. \frac{\mathrm{d}\theta_c}{\mathrm{d}\bar{m}'} \right|_0
   = 3 (Q_1 - Q_0),
   \qquad
   \left. \frac{\mathrm{d}\theta_c}{\mathrm{d}\bar{m}'} \right|_1
   = 3 (Q_3 - Q_2).

**Enforcement of metal angles.**  Setting

.. math::
   Q_0 &= \theta_\mathrm{LE}, \\
   Q_1 &= \theta_\mathrm{LE} + \tfrac{1}{3}\tan\beta_1, \\
   Q_2 &= \theta_\mathrm{TE} - \tfrac{1}{3}\tan\beta_2, \\
   Q_3 &= \theta_\mathrm{TE}

enforces the metal angles *exactly* — no penalty term, no iteration.  The
two free interior control points (:math:`Q_1`, :math:`Q_2`) are entirely
determined by the aerodynamic angles.

Thickness Distribution
----------------------

**Curve family.**  A clamped degree-4 B-spline with 6 control points
:math:`T_0, \ldots, T_5` and interior knot at :math:`\bar{m}' = 0.5`,
giving the knot vector

.. math::
   \{0, 0, 0, 0, 0,\; 0.5,\; 1, 1, 1, 1, 1\}.

The single interior knot provides one additional inflection degree of
freedom while keeping the overall curve smooth (:math:`C^3` everywhere).

**Boundary conditions.**  Five conditions are imposed:

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Condition
     - Physical meaning
   * - :math:`\Delta\theta(0) = 0`
     - Zero thickness at the leading edge
   * - :math:`\Delta\theta'(0) = 0`
     - Smooth LE blend (zero slope — blunt LE approach)
   * - :math:`\Delta\theta \approx t_\mathrm{max}` near :math:`t_\mathrm{max,loc}`
     - Aerodynamic peak thickness target
   * - :math:`\Delta\theta(1) = t_\mathrm{TE}`
     - Prescribed trailing-edge thickness
   * - :math:`\Delta\theta'(1) = -2\tan(\psi_\mathrm{TE}/2)`
     - TE wedge angle

**Control-point mapping.**  The endpoint derivative of a clamped degree-4
spline with 6 control points satisfies

.. math::
   \Delta\theta'(0) = \frac{4(T_1 - T_0)}{\xi_5},
   \qquad
   \Delta\theta'(1) = \frac{4(T_5 - T_4)}{1 - \xi_5},

where :math:`\xi_5 = 0.5`.  Setting :math:`T_0 = T_1 = 0` simultaneously
enforces zero thickness and zero slope at the LE.  The remaining control
points are assigned as:

.. math::
   T_2 &= 0.70\; t_\mathrm{max}, \\
   T_3 &= t_\mathrm{max}, \\
   T_4 &= t_\mathrm{TE} + 0.2 \cdot 2\tan(\psi_\mathrm{TE}/2), \\
   T_5 &= t_\mathrm{TE}.

The factor 0.70 on :math:`T_2` shapes the rise toward the peak; the
factor 0.2 on :math:`T_4` positions the wedge slope over the last
~20 % of chord.

Surface Construction
--------------------

**Camber-line normal.**  In the Euclidean :math:`(\bar{m}', \theta)` plane
the unit tangent to the camber line is

.. math::
   \hat{\mathbf{t}} = \frac{1}{\sqrt{1 + \theta_c'^{\,2}}}
   \begin{pmatrix} 1 \\ \theta_c' \end{pmatrix},

and the inward unit normal (pointing toward the suction side) is

.. math::
   \hat{\mathbf{n}} = \frac{1}{\sqrt{1 + \theta_c'^{\,2}}}
   \begin{pmatrix} -\theta_c' \\ 1 \end{pmatrix}.

**Suction and pressure sides.**  At each chordwise station :math:`\bar{m}'_c`:

.. math::
   \text{Suction side:} \quad
   (\bar{m}', \theta)_\mathrm{SS} &= (\bar{m}'_c, \theta_c) + \tfrac{1}{2}\Delta\theta \cdot \hat{\mathbf{n}}, \\
   \text{Pressure side:} \quad
   (\bar{m}', \theta)_\mathrm{PS} &= (\bar{m}'_c, \theta_c) - \tfrac{1}{2}\Delta\theta \cdot \hat{\mathbf{n}}.

The sign convention is such that the suction side is on the high-:math:`\theta`
side for a turbine rotor blade with positive camber.

Point Clustering
----------------

Points are distributed along the chord using **cosine clustering**:

.. math::
   u_i = \frac{1}{2}\left(1 - \cos\frac{\pi i}{N-1}\right), \quad
   i = 0, 1, \ldots, N-1,

where :math:`N` is the number of points per side.  Cosine spacing
concentrates points near :math:`\bar{m}' = 0` (LE) and :math:`\bar{m}' = 1`
(TE), where curvature is highest, improving the quality of the downstream
OCCT loft.

The closed profile loop is assembled as:

1. Suction side LE → TE (:math:`N` points, cosine-spaced).
2. Pressure side TE → LE (:math:`N` points, same cosine distribution in
   reverse), closing the loop at the LE.

Implementation
--------------

- :file:`src/blade/ProfileSection.hpp` — ``TurbineProfileParams<T>``,
  ``CamberLine<T>``, ``ThicknessDistribution<T>``, ``ProfileSection<T>``
- :file:`src/blade/BSplineBasis.hpp` — ``basis_nonzero()``,
  ``basis_deriv_nonzero()``, ``clamped_knots()``
