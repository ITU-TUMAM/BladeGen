/// @file blade/StreamlineMPrime.hpp
/// @brief Meridional streamline as a B-spline curve and RK4 integration of the
///        conformal meridional coordinate m′ = ∫ dm/r along that streamline.
///
/// The conformal coordinate m′ produces a flat design plane (m′, θ) in which
/// angles between curves equal the physical flow angles in 3D — see the (m′, θ)
/// system derivation: ds² = r²(dm′² + dθ²).
///
/// All types are templated on T so the same code compiles for double (evaluation)
/// and CppAD::AD<double> (AD taping for gradient-based optimisation).

#pragma once

#include "BSplineBasis.hpp"

#include <array>
#include <cassert>
#include <vector>

namespace PCAD::Blade {

// ─── 2D B-spline streamline in (x, r) ────────────────────────────────────────
/// Meridional construction streamline as a B-spline curve.
/// Control point coordinates are typed T so they can be AD active variables.
template <typename T>
struct BSplineCurve2D
{
    int                 degree;
    std::vector<double> knots;  ///< Fixed knot topology — always double
    std::vector<T>      Px;     ///< x control points
    std::vector<T>      Pr;     ///< r control points

    [[nodiscard]] int n_ctrl() const noexcept
    {
        return static_cast<int>(Px.size());
    }

    /// Evaluate position (x(t), r(t)).
    [[nodiscard]] std::array<T, 2> eval(T t) const
    {
        const int k = find_span(to_double(t), degree, knots);
        const auto N = basis_nonzero<T>(t, k, degree, knots);

        T x_val = T(0.0), r_val = T(0.0);
        for (int j = 0; j <= degree; ++j) {
            const int idx = k - degree + j;
            x_val += N[j] * Px[idx];
            r_val += N[j] * Pr[idx];
        }
        return {x_val, r_val};
    }

    /// Evaluate tangent (dx/dt, dr/dt).
    [[nodiscard]] std::array<T, 2> eval_deriv(T t) const
    {
        const int k  = find_span(to_double(t), degree, knots);
        const auto dN = basis_deriv_nonzero<T>(t, k, degree, knots);

        T dx_dt = T(0.0), dr_dt = T(0.0);
        for (int j = 0; j <= degree; ++j) {
            const int idx = k - degree + j;
            dx_dt += dN[j] * Px[idx];
            dr_dt += dN[j] * Pr[idx];
        }
        return {dx_dt, dr_dt};
    }
};

// ─── ODE right-hand side ─────────────────────────────────────────────────────
/// Compute dy/dt at parameter t along the streamline:
///   y[0] = m′(t),   dy[0]/dt = ‖ẋ(t)‖ / r(t)
///   y[1] = m(t),    dy[1]/dt = ‖ẋ(t)‖
///
/// The RHS is independent of y — this is pure quadrature cast as an ODE,
/// which guarantees 4th-order convergence with no stability concerns.
template <typename T>
[[nodiscard]] std::array<T, 2> mprime_rhs(T t, const BSplineCurve2D<T>& sl)
{
    const auto [dx_dt, dr_dt] = sl.eval_deriv(t);

    using std::sqrt;
    T arc_speed = sqrt(dx_dt * dx_dt + dr_dt * dr_dt);

    const auto [x_val, r_val] = sl.eval(t);

    // Guard r → 0 at axis; turbine annulus geometry keeps r >> 1e-12 m in
    // practice, but this avoids NaN in the derivative tape.
    T r_safe = CppAD::CondExpGt(r_val, T(1e-12), r_val, T(1e-12));

    return { arc_speed / r_safe, arc_speed };
}

// ─── Integration table ───────────────────────────────────────────────────────
/// Tabulated result of integrating m′(t) and m(t) along [0, 1].
/// Stored as vectors of T so gradients flow through the table values.
template <typename T>
struct MprimeTable
{
    std::vector<T> t_vals;
    std::vector<T> mprime_vals;
    std::vector<T> m_vals;
};

// ─── RK4 integrator ──────────────────────────────────────────────────────────
/// Fixed-step classical RK4 integration of m′ and m from t=0 to t=1.
///
/// Because the ODE RHS is independent of y the method reduces to 4-point
/// Gaussian-type quadrature per step, giving exact 4th-order convergence.
/// n_steps=200 gives relative errors < 1e-8 for smooth turbine streamlines.
template <typename T>
[[nodiscard]] MprimeTable<T> integrate_mprime(const BSplineCurve2D<T>& sl,
                                               int n_steps = 200)
{
    assert(n_steps > 0);

    MprimeTable<T> tbl;
    tbl.t_vals.reserve(n_steps + 1);
    tbl.mprime_vals.reserve(n_steps + 1);
    tbl.m_vals.reserve(n_steps + 1);

    const T h = T(1.0 / n_steps);
    T t  = T(0.0);
    T mp = T(0.0);
    T m  = T(0.0);

    tbl.t_vals.push_back(t);
    tbl.mprime_vals.push_back(mp);
    tbl.m_vals.push_back(m);

    for (int i = 0; i < n_steps; ++i) {
        const auto f1 = mprime_rhs<T>(t,           sl);
        const auto f2 = mprime_rhs<T>(t + h / T(2), sl);
        const auto f3 = mprime_rhs<T>(t + h / T(2), sl);  // == f2; kept for clarity
        const auto f4 = mprime_rhs<T>(t + h,         sl);

        mp += (h / T(6.0)) * (f1[0] + T(2) * f2[0] + T(2) * f3[0] + f4[0]);
        m  += (h / T(6.0)) * (f1[1] + T(2) * f2[1] + T(2) * f3[1] + f4[1]);
        t  += h;

        tbl.t_vals.push_back(t);
        tbl.mprime_vals.push_back(mp);
        tbl.m_vals.push_back(m);
    }
    return tbl;
}

// ─── m′ inversion ────────────────────────────────────────────────────────────
/// Given a target m′* find the streamline parameter t by linear interpolation
/// in the pre-built table.
///
/// The binary search uses primal (double) values for the index selection so no
/// branch enters the AD tape, while the interpolation itself uses AD-typed T
/// values so gradients flow through the result.
template <typename T>
[[nodiscard]] T invert_mprime(T mp_target, const MprimeTable<T>& tbl)
{
    const int n  = static_cast<int>(tbl.mprime_vals.size());
    const double mp_d = to_double(mp_target);

    // Binary search on primal values — branch-free in the tape
    int lo = 0, hi = n - 2;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        (mp_d < to_double(tbl.mprime_vals[mid])) ? hi = mid : lo = mid;
    }

    const T lo_mp = tbl.mprime_vals[lo];
    const T hi_mp = tbl.mprime_vals[lo + 1];
    // Safe denominator: interval should be strictly positive, but guard anyway
    const T span  = hi_mp - lo_mp;
    const T safe_span = CppAD::CondExpGt(span, T(1e-30), span, T(1.0));
    const T alpha = CppAD::CondExpGt(span, T(1e-30),
                                     (mp_target - lo_mp) / safe_span,
                                     T(0.5));

    return tbl.t_vals[lo] + alpha * (tbl.t_vals[lo + 1] - tbl.t_vals[lo]);
}

// ─── Conformal → 3D mapping ──────────────────────────────────────────────────
/// Map a design-plane point (mp_norm, θ) to a physical 3D point (x, y, z).
///
/// @param mp_norm   normalised m′ coordinate in [0, 1]  (m′_abs = mp_norm · m′_max)
/// @param theta     circumferential angle [rad]
/// @param sl        construction streamline for this span section
/// @param tbl       pre-integrated m′ table for the same streamline
/// @returns         Cartesian point (x, r·sin θ, r·cos θ)
template <typename T>
[[nodiscard]] std::array<T, 3> conformal_to_3d(T mp_norm, T theta,
                                                const BSplineCurve2D<T>& sl,
                                                const MprimeTable<T>& tbl)
{
    const T mp_max = tbl.mprime_vals.back();
    const T mp_abs = mp_norm * mp_max;
    const T t_star = invert_mprime<T>(mp_abs, tbl);

    const auto [x, r] = sl.eval(t_star);

    using std::sin;
    using std::cos;
    return { x, r * sin(theta), r * cos(theta) };
}

} // namespace PCAD::Blade
