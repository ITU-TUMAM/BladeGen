/// @file blade/BladeSection.hpp
/// @brief 3D blade section and multi-span blade geometry.
///
/// A BladeSpanSection<T> combines a meridional streamline, a 2D profile in the
/// conformal (m′, θ) design plane, and a stacking offset (Δx, Δθ) to produce
/// a set of 3D surface points for one spanwise station.
///
/// BladeGeometry<T> collects sections ordered hub → tip, applies stacking laws,
/// and exposes the 3D point arrays needed by the OCCT skinning layer.

#pragma once

#include "ProfileSection.hpp"
#include "StreamlineMPrime.hpp"

#include <array>
#include <vector>

namespace PCAD::Blade {

// ─── Stacking offset ─────────────────────────────────────────────────────────
/// Rigid shift applied to a section after conformal mapping to 3D.
/// Lean (Δθ) and sweep (Δx) are the classical turbine stacking degrees of
/// freedom; both are zero for a pure radially-stacked blade.
template <typename T>
struct StackingOffset
{
    T dx     = T(0.0);   ///< axial offset [same length unit as streamline]
    T dtheta = T(0.0);   ///< circumferential offset [rad]
};

// ─── Single span section ─────────────────────────────────────────────────────
/// One spanwise station of the blade.
///
/// Workflow per section:
///   1.  integrate_mprime() on the streamline → MprimeTable
///   2.  profile.generate(n_pts) → (m′_norm, θ) loop in the design plane
///   3.  conformal_to_3d() for each point → 3D loop before stacking
///   4.  Apply stacking offset
template <typename T>
struct BladeSpanSection
{
    BSplineCurve2D<T>    streamline;  ///< meridional construction streamline
    ProfileSection<T>    profile;     ///< 2D profile in the design plane
    StackingOffset<T>    stacking;    ///< lean / sweep offsets

    /// Pre-build the m′ integration table.  Call this once after constructing
    /// the section and before calling points3d().
    [[nodiscard]] MprimeTable<T> build_table(int n_steps = 200) const
    {
        return integrate_mprime<T>(streamline, n_steps);
    }

    /// Evaluate the 3D point loop for this section.
    ///
    /// @param tbl     pre-integrated m′ table (from build_table())
    /// @param n_pts   number of points per side (total loop = 2·n_pts)
    /// @returns       ordered 3D points (x, y, z) around the profile
    [[nodiscard]] std::vector<std::array<double, 3>>
    points3d(const MprimeTable<T>& tbl, int n_pts = 64) const
    {
        const auto loop_2d = profile.generate(n_pts);
        std::vector<std::array<double, 3>> pts;
        pts.reserve(loop_2d.size());

        for (const auto& pp : loop_2d) {
            // mp is already normalised to [0,1] by ProfileSection::generate
            const T theta_total = pp.theta + stacking.dtheta;
            const auto [x3, y3, z3] =
                conformal_to_3d<T>(pp.mp, theta_total, streamline, tbl);

            const double x_out = to_double(x3) + to_double(stacking.dx);
            pts.push_back({ x_out, to_double(y3), to_double(z3) });
        }
        return pts;
    }

    /// Convenience overload — builds the m′ table internally.
    [[nodiscard]] std::vector<std::array<double, 3>>
    points3d(int n_pts = 64, int n_rk4_steps = 200) const
    {
        const auto tbl = build_table(n_rk4_steps);
        return points3d(tbl, n_pts);
    }
};

// ─── Multi-span blade geometry ────────────────────────────────────────────────
/// Collection of spanwise sections ordered hub (index 0) → tip (index back).
///
/// For AD-based optimisation: make sections with T = CppAD::AD<double> and
/// call all_sections_3d() inside a recorded function to obtain
/// ∂(surface points) / ∂(design parameters) via reverse-mode Jacobian.
template <typename T>
struct BladeGeometry
{
    std::vector<BladeSpanSection<T>> sections;

    /// Generate 3D point loops for all span sections.
    ///
    /// @param n_pts        profile sampling density (points per side)
    /// @param n_rk4_steps  RK4 integration steps for m′
    /// @returns            outer vector = sections (hub→tip),
    ///                     inner vector = 3D loop points
    [[nodiscard]]
    std::vector<std::vector<std::array<double, 3>>>
    all_sections_3d(int n_pts = 64, int n_rk4_steps = 200) const
    {
        std::vector<std::vector<std::array<double, 3>>> result;
        result.reserve(sections.size());

        for (const auto& sec : sections)
            result.push_back(sec.points3d(n_pts, n_rk4_steps));

        return result;
    }
};

// ─── Convenience factory ─────────────────────────────────────────────────────
/// Build a BladeGeometry<double> for a simple annular streamline
/// (straight hub–tip passage) given a hub radius, tip radius, and axial chord.
///
/// Both hub and tip use the same aerodynamic profile (common starting point
/// for preliminary design); the stacking is radially pure.
///
/// @param r_hub       hub radius [m]
/// @param r_tip       tip radius [m]
/// @param axial_chord axial extent of the blade [m]
/// @param params_hub  profile parameters at hub
/// @param params_tip  profile parameters at tip
/// @param n_spans     number of intermediate span sections (total = n_spans + 2)
[[nodiscard]] inline BladeGeometry<double> make_annular_blade(
    double r_hub, double r_tip, double axial_chord,
    const TurbineProfileParams<double>& params_hub,
    const TurbineProfileParams<double>& params_tip,
    int n_spans = 3)
{
    BladeGeometry<double> geom;
    const int n_total = n_spans + 2;

    for (int i = 0; i < n_total; ++i) {
        const double span_frac = static_cast<double>(i) / (n_total - 1);
        const double r = r_hub + span_frac * (r_tip - r_hub);

        // Straight axial streamline at this radius:
        //   x(t) = t * axial_chord,   r(t) = r  (constant radius)
        BSplineCurve2D<double> sl;
        sl.degree = 1;                            // linear — straight line
        sl.knots  = clamped_knots(1, 2);          // {0, 0, 1, 1}
        sl.Px     = {0.0, axial_chord};
        sl.Pr     = {r,   r};

        // Linearly interpolate profile parameters hub → tip
        const double s = span_frac;
        TurbineProfileParams<double> p;
        p.beta1     = (1-s)*params_hub.beta1     + s*params_tip.beta1;
        p.beta2     = (1-s)*params_hub.beta2     + s*params_tip.beta2;
        p.theta_LE  = (1-s)*params_hub.theta_LE  + s*params_tip.theta_LE;
        p.theta_TE  = (1-s)*params_hub.theta_TE  + s*params_tip.theta_TE;
        p.t_max     = (1-s)*params_hub.t_max     + s*params_tip.t_max;
        p.t_max_loc = (1-s)*params_hub.t_max_loc + s*params_tip.t_max_loc;
        p.t_TE      = (1-s)*params_hub.t_TE      + s*params_tip.t_TE;
        p.wedge_TE  = (1-s)*params_hub.wedge_TE  + s*params_tip.wedge_TE;
        p.r_LE      = (1-s)*params_hub.r_LE      + s*params_tip.r_LE;

        BladeSpanSection<double> sec;
        sec.streamline = std::move(sl);
        sec.profile    = ProfileSection<double>::from_params(p);
        // sec.stacking left at default (radial stacking)

        geom.sections.push_back(std::move(sec));
    }
    return geom;
}

} // namespace PCAD::Blade
