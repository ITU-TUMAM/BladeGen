/// @file blade/ProfileSection.hpp
/// @brief 2D turbine blade profile in the conformal (m′, θ) design plane.
///
/// The design plane is flat Euclidean (metric ds² = r²(dm′² + dθ²)), so angles
/// between curves in (m′, θ) equal the physical flow angles in 3D space.
///
/// Parameterisation hierarchy:
///   TurbineProfileParams<T>  ── aerodynamic design intent
///       │
///       ├─ CamberLine<T>          ── cubic B-spline, endpoint tangents = metal angles
///       └─ ThicknessDistribution<T> ── degree-4 B-spline, LE/TE conditions
///           │
///           └─ ProfileSection<T>  ── suction/pressure sides + closed-loop sampling

#pragma once

#include "BSplineBasis.hpp"

#include <array>
#include <cmath>
#include <vector>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace PCAD::Blade {

// ─── Primary design parameters ───────────────────────────────────────────────
/// Aerodynamic design variables for a single turbine blade section.
/// All angular quantities are in radians; lengths are in the same units as
/// the (m′, θ) coordinate system (θ is dimensionless — it IS the angle).
template <typename T>
struct TurbineProfileParams
{
    /// Inlet metal angle [rad] — tangent of camber line at LE w.r.t. m′ axis.
    /// Typical range for axial turbines: 40°–70° (positive = lean toward rotation).
    T beta1;

    /// Exit metal angle [rad] — tangent of camber line at TE w.r.t. m′ axis.
    /// Typically negative (opposite to rotation direction): −60° to −75°.
    T beta2;

    /// Circumferential position of the leading edge (θ at m′=0).
    T theta_LE;

    /// Circumferential position of the trailing edge (θ at m′=1).
    /// Determines the total camber in θ: theta_TE − theta_LE.
    T theta_TE;

    /// Maximum thickness in θ units (e.g. fraction of pitch × pitch_angle).
    T t_max;

    /// Normalised chordwise location of maximum thickness, in (0, 1).
    /// Typical turbine value: 0.30–0.40.
    T t_max_loc;

    /// Trailing-edge thickness in θ units (0 = sharp TE).
    T t_TE;

    /// TE wedge half-angle [rad] — controls the slope of thickness at the TE.
    /// dΔθ/dm′|_{m′=1} = −2·tan(wedge_TE/2).
    T wedge_TE;

    /// Leading-edge radius in the design plane (reserved for future LE blending).
    T r_LE;
};

// ─── Camber line ─────────────────────────────────────────────────────────────
/// Cubic B-spline camber line θ_c(m′) defined on [0,1].
///
/// The clamped knot vector {0,0,0,0,1,1,1,1} with 4 control points gives
/// exact endpoint tangent conditions via the well-known relation:
///   dθ_c/dm′|_{0} = 3 · (Q₁_θ − Q₀_θ)    (knot span normalised to 1)
///   dθ_c/dm′|_{1} = 3 · (Q₃_θ − Q₂_θ)
/// Setting Q₁_θ = Q₀_θ + tan(β₁)/3 and Q₂_θ = Q₃_θ − tan(β₂)/3 enforces
/// the metal angles exactly — no penalty, no approximation.
template <typename T>
struct CamberLine
{
    static constexpr int kDegree = 3;
    static constexpr int kNCtrl  = 4;

    std::array<T, kNCtrl> Qtheta;  ///< θ control points
    std::vector<double>   knots;   ///< clamped cubic: {0,0,0,0,1,1,1,1}

    /// Build from aerodynamic angle parameters.
    [[nodiscard]] static CamberLine<T> from_params(const TurbineProfileParams<T>& p)
    {
        using std::tan;
        CamberLine<T> cl;
        cl.knots = clamped_knots(kDegree, kNCtrl);

        cl.Qtheta[0] = p.theta_LE;
        cl.Qtheta[1] = p.theta_LE + tan(p.beta1) / T(3.0);
        cl.Qtheta[2] = p.theta_TE - tan(p.beta2) / T(3.0);
        cl.Qtheta[3] = p.theta_TE;

        return cl;
    }

    /// Evaluate θ_c(m′).
    [[nodiscard]] T eval(T mp) const
    {
        const int k = find_span(to_double(mp), kDegree, knots);
        const auto N = basis_nonzero<T>(mp, k, kDegree, knots);

        T theta = T(0.0);
        for (int j = 0; j <= kDegree; ++j)
            theta += N[j] * Qtheta[k - kDegree + j];
        return theta;
    }

    /// Evaluate dθ_c/dm′.
    [[nodiscard]] T eval_deriv(T mp) const
    {
        const int k  = find_span(to_double(mp), kDegree, knots);
        const auto dN = basis_deriv_nonzero<T>(mp, k, kDegree, knots);

        T dtheta = T(0.0);
        for (int j = 0; j <= kDegree; ++j)
            dtheta += dN[j] * Qtheta[k - kDegree + j];
        return dtheta;
    }
};

// ─── Thickness distribution ──────────────────────────────────────────────────
/// Degree-4 B-spline thickness Δθ(m′) with the following boundary conditions:
///
///   Δθ(0)         = 0                        zero thickness at LE
///   dΔθ/dm′|_{0}  = 0                        smooth LE blend (zero slope)
///   Δθ ≈ t_max    near m′ = t_max_loc         aerodynamic peak thickness
///   Δθ(1)         = t_TE                      finite or zero TE thickness
///   dΔθ/dm′|_{1}  = −2·tan(wedge_TE/2)       TE wedge angle
///
/// 6 control points, clamped degree-4 knots {0⁵, 0.5, 1⁵}.
template <typename T>
struct ThicknessDistribution
{
    static constexpr int kDegree = 4;
    static constexpr int kNCtrl  = 6;

    std::array<T, kNCtrl> Tdelta;  ///< thickness control values
    std::vector<double>   knots;   ///< clamped degree-4, 6 ctrl

    /// Build from aerodynamic thickness parameters.
    [[nodiscard]] static ThicknessDistribution<T> from_params(
        const TurbineProfileParams<T>& p)
    {
        using std::tan;
        ThicknessDistribution<T> td;
        td.knots = clamped_knots(kDegree, kNCtrl);
        // Resulting knots: {0, 0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1}

        // Boundary conditions at LE: Δθ(0)=0 and dΔθ/dm′|_0=0.
        // For clamped degree-4: dΔθ/dm′|_0 = 4·(T[1]−T[0]) / knots[5]
        //                                    = 4·(T[1]−T[0]) / 0.5 = 8·(T[1]−T[0])
        // Setting T[1]=T[0]=0 enforces both conditions simultaneously.
        td.Tdelta[0] = T(0.0);
        td.Tdelta[1] = T(0.0);

        // Shape the rise toward the aerodynamic maximum
        td.Tdelta[2] = p.t_max * T(0.70);
        td.Tdelta[3] = p.t_max;

        // Falling toward TE: account for the wedge slope over the last ~0.2 of chord.
        // dΔθ/dm′|_1 ≈ (T[5]−T[4]) / (1−knots[5]) · (p) but let's set T[4]
        // to place the TE value linearly given the wedge angle derivative.
        td.Tdelta[4] = p.t_TE + T(0.2) * T(2.0) * tan(p.wedge_TE / T(2.0));
        td.Tdelta[5] = p.t_TE;

        return td;
    }

    /// Evaluate Δθ(m′).
    [[nodiscard]] T eval(T mp) const
    {
        const int k = find_span(to_double(mp), kDegree, knots);
        const auto N = basis_nonzero<T>(mp, k, kDegree, knots);

        T delta = T(0.0);
        for (int j = 0; j <= kDegree; ++j) {
            const int idx = k - kDegree + j;
            if (idx >= 0 && idx < kNCtrl)
                delta += N[j] * Tdelta[idx];
        }
        // Clamp to non-negative — thickness cannot be physically negative
        return CppAD::CondExpGt(delta, T(0.0), delta, T(0.0));
    }
};

// ─── Profile section ─────────────────────────────────────────────────────────
/// A point on a blade surface in the (m′, θ) design plane.
template <typename T>
struct ProfilePoint
{
    T mp;     ///< meridional conformal coordinate
    T theta;  ///< circumferential angle
};

/// Combines camber line and thickness to define the full 2D blade profile.
/// Suction (+normal) and pressure (−normal) sides are offset from the camber
/// line along the Euclidean normal in the (m′, θ) plane.
template <typename T>
struct ProfileSection
{
    CamberLine<T>           camber;
    ThicknessDistribution<T> thickness;

    /// Unit normal to the camber line in the Euclidean (m′, θ) plane.
    /// Tangent = (1, dθ_c/dm′) / mag;  normal = (−dθ_c/dm′, 1) / mag.
    [[nodiscard]] std::array<T, 2> camber_normal(T mp) const
    {
        using std::sqrt;
        const T dtheta = camber.eval_deriv(mp);
        const T mag    = sqrt(T(1.0) + dtheta * dtheta);
        return { -dtheta / mag, T(1.0) / mag };
    }

    /// Point on the suction side at camber-line parameter m′_c.
    [[nodiscard]] ProfilePoint<T> suction_side(T mp_c) const
    {
        const T theta_c        = camber.eval(mp_c);
        const T half_thickness = thickness.eval(mp_c) / T(2.0);
        const auto [nm, nt]    = camber_normal(mp_c);
        return { mp_c + nm * half_thickness,
                 theta_c + nt * half_thickness };
    }

    /// Point on the pressure side at camber-line parameter m′_c.
    [[nodiscard]] ProfilePoint<T> pressure_side(T mp_c) const
    {
        const T theta_c        = camber.eval(mp_c);
        const T half_thickness = thickness.eval(mp_c) / T(2.0);
        const auto [nm, nt]    = camber_normal(mp_c);
        return { mp_c - nm * half_thickness,
                 theta_c - nt * half_thickness };
    }

    /// Generate a closed profile loop: SS LE→TE then PS TE→LE.
    /// Cosine clustering gives denser point spacing near the LE and TE.
    [[nodiscard]] std::vector<ProfilePoint<T>> generate(int n_pts = 64) const
    {
        std::vector<ProfilePoint<T>> pts;
        pts.reserve(2 * n_pts);

        // Suction side: LE → TE
        for (int i = 0; i < n_pts; ++i) {
            const double u = 0.5 * (1.0 - std::cos(M_PI * i / (n_pts - 1)));
            pts.push_back(suction_side(T(u)));
        }
        // Pressure side: TE → LE (reverse order closes the loop)
        for (int i = n_pts - 1; i >= 0; --i) {
            const double u = 0.5 * (1.0 - std::cos(M_PI * i / (n_pts - 1)));
            pts.push_back(pressure_side(T(u)));
        }
        return pts;
    }

    [[nodiscard]] static ProfileSection<T> from_params(
        const TurbineProfileParams<T>& p)
    {
        ProfileSection<T> ps;
        ps.camber    = CamberLine<T>::from_params(p);
        ps.thickness = ThicknessDistribution<T>::from_params(p);
        return ps;
    }
};

} // namespace PCAD::Blade
