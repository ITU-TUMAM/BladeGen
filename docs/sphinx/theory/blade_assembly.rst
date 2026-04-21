3-D Blade Assembly
==================

Once the conformal design plane and the 2-D profile parameterisation are in
place (see :doc:`coordinate_system` and :doc:`profile`), the three-dimensional
blade solid is built in three steps:

1. **Section generation** — each spanwise station produces a set of 3-D
   surface points via the inverse conformal map.
2. **Span interpolation** — intermediate sections are inserted between
   user-defined stations by linearly blending all profile parameters.
3. **OCCT skinning** — the ordered point loops are converted to OCCT wire
   loops and lofted into a closed shell solid.

Spanwise Stations and Interpolation
------------------------------------

The user provides at least two *design stations*: the hub (:math:`\eta = 0`)
and the tip (:math:`\eta = 1`), where :math:`\eta` is the normalised span
coordinate.  Additional stations at intermediate span fractions may be
supplied to capture spanwise twist, thickness taper, or custom streamline
shapes.

``BladeParams::n_interp_sections`` controls how many uniformly-spaced
intermediate stations are inserted between adjacent user-defined stations.
All ten profile parameters (:math:`\beta_1, \beta_2, \theta_\mathrm{LE},
\theta_\mathrm{TE}, t_\mathrm{max}, \ldots`) are interpolated linearly in
:math:`\eta` between the bounding design stations.

Section-to-3D Mapping
----------------------

For each spanwise station the mapping from the design plane to 3-D Cartesian
space proceeds as follows.

**Step 1 — Build the** :math:`m'` **table.**
Integrate the meridional conformal coordinate along the construction
streamline (RK4, 200 steps by default):

.. math::
   m'(t_i) = \int_0^{t_i} \frac{\|\dot{\mathbf{x}}(\tau)\|}{r(\tau)} \, \mathrm{d}\tau,
   \qquad i = 0, \ldots, N_\mathrm{steps}.

**Step 2 — Sample the design-plane profile.**
Generate the closed :math:`2N` point loop
:math:`\{(\bar{m}'_j, \theta_j)\}_{j=0}^{2N-1}` using cosine clustering
(see :doc:`profile`).

**Step 3 — Map each point to 3-D.**
For each design-plane point :math:`(\bar{m}', \theta)`:

a. Scale to absolute: :math:`m'_\mathrm{abs} = \bar{m}' \cdot m'_\mathrm{max}`.
b. Invert the table: find :math:`t^*` such that
   :math:`m'(t^*) = m'_\mathrm{abs}` (binary search + linear interpolation).
c. Evaluate the streamline: :math:`(x, r) = \mathbf{x}(t^*)`.
d. Apply circumferential stacking offset: :math:`\theta_\mathrm{total} = \theta + \Delta\theta`.
e. Convert to Cartesian:

.. math::
   \mathbf{p} = \begin{pmatrix} x + \Delta x \\ r \sin\theta_\mathrm{total} \\ r \cos\theta_\mathrm{total} \end{pmatrix}.

Stacking
---------

*Stacking* refers to the alignment of corresponding spanwise points across
different sections.  Pure *radial stacking* aligns the sections such that a
reference point (commonly the centroid or the leading edge) lies on a radial
line.  BladeGen supports two additional stacking degrees of freedom:

.. list-table::
   :widths: 20 80
   :header-rows: 1

   * - Offset
     - Physical effect
   * - :math:`\Delta x`
     - **Sweep** — axial shift of the section relative to the
       radially-stacked position.  Forward sweep reduces shock losses on
       transonic blades; aft sweep shifts the axial load distribution.
   * - :math:`\Delta\theta`
     - **Lean** — circumferential shift of the section.  Lean modifies the
       spanwise pressure gradient and can be used to control secondary flows
       near endwalls.

Both offsets are applied in step (d)/(e) above after the conformal-to-3D
mapping, so they act as rigid body translations in the physical space and do
not distort the profile shape.

OCCT Solid Construction
------------------------

The ordered 3-D point loops for all spanwise sections are handed to
``OcctBlade``, which:

1. Converts each point loop to an OCCT ``BRepBuilderAPI_MakePolygon`` wire.
2. Fills each wire to a planar face (``BRepBuilderAPI_MakeFace``) to create
   section caps.
3. Lofts the ordered wire sequence using ``BRepOffsetAPI_ThruSections`` with
   ruled (linear) interpolation between adjacent sections to produce a
   closed shell.
4. Runs ``BRepBuilderAPI_Sewing`` to merge coincident edges between the
   lofted shell and the hub/tip cap faces.
5. Passes the result through ``BRepBuilderAPI_MakeSolid`` to obtain a
   fully closed, orientable solid suitable for export.

The solid can be exported to IGES, STEP (AP203 or AP214), or binary/ASCII
STL via the :file:`src/io/` exporters.

AD-Compatible Formulation
--------------------------

Every function in the geometry stack — B-spline evaluation, RK4 integration,
conformal inversion, profile sampling — is templated on the scalar type
:math:`T`.  Substituting ``T = CppAD::AD<double>`` records a CppAD tape of
the complete map

.. math::
   \boldsymbol{\xi} \mapsto \{\mathbf{p}_j\}_{j=0}^{2N n_\mathrm{sec}-1}

from the design-variable vector :math:`\boldsymbol{\xi}` (all
:math:`\beta_1, \beta_2, t_\mathrm{max}, \ldots` across all stations) to the
3-D surface-point coordinates.  The reverse-mode Jacobian
:math:`\partial \mathbf{p} / \partial \boldsymbol{\xi}` can then be computed
at the cost of a single additional pass, enabling gradient-based optimisation
of any objective that depends on the surface geometry (e.g., aerodynamic loss,
mechanical stress, or passage area distribution).

Implementation
--------------

- :file:`src/blade/BladeSection.hpp` — ``BladeSpanSection<T>``,
  ``BladeGeometry<T>``, ``StackingOffset<T>``
- :file:`src/blade/OcctBlade.hpp` / :file:`.cpp` — ``OcctBlade``
- :file:`src/params/BladeParams.hpp` / :file:`.cpp` — ``BladeParams``,
  ``StationParams``, interpolation and span assembly
- :file:`src/pipeline/BladeRunner.hpp` / :file:`.cpp` — ``BladeRunner``,
  ``BladeRunnerConfig``, ``BladeRunnerResult``
