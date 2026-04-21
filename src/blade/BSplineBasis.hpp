// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
/// @file blade/BSplineBasis.hpp
/// @brief Templated Cox–de Boor B-spline basis and derivative evaluation.
///
/// All functions are templated on the scalar type T so they compile for both
/// plain double (geometry evaluation) and CppAD::AD<double> (AD taping).
/// OCCT is not referenced here — this is pure numerics.

#pragma once

#include <cppad/cppad.hpp>
#include <cmath>
#include <type_traits>
#include <vector>

namespace BladeGen::Blade {

// ─── Primal extraction ────────────────────────────────────────────────────────
/// Extract the underlying double from T.
/// For plain double this is a no-op; for AD<double> it returns the primal value
/// at the current evaluation point (required for branch decisions inside tapes).
template <typename T>
[[nodiscard]] double to_double(const T& x)
{
    if constexpr (std::is_floating_point_v<T>)
        return static_cast<double>(x);
    else
        return CppAD::Value(CppAD::Var2Par(x));
}

// ─── Knot-span search ─────────────────────────────────────────────────────────
/// Binary search for the knot span index k such that knots[k] ≤ t < knots[k+1].
/// Operates on double (primal) values only — never enters the AD tape.
/// @param t   parameter value in [knots[p], knots[n+1]]
/// @param p   B-spline degree
/// @param knots  knot vector of length n+p+2
/// @returns   span index k in [p, n]
[[nodiscard]] inline int find_span(double t, int p,
                                   const std::vector<double>& knots)
{
    const int n = static_cast<int>(knots.size()) - p - 2;

    if (t >= knots[static_cast<std::size_t>(n + 1)]) return n;   // clamp at upper boundary
    if (t <= knots[static_cast<std::size_t>(p)])     return p;   // clamp at lower boundary

    int lo = p, hi = n + 1;
    while (hi - lo > 1) {
        int mid = (lo + hi) / 2;
        (t < knots[static_cast<std::size_t>(mid)]) ? hi = mid : lo = mid;
    }
    return lo;
}

// ─── Non-zero basis values ────────────────────────────────────────────────────
/// Triangular (de Boor) algorithm returning the p+1 non-zero basis values
/// N_{k-p,p}(t), …, N_{k,p}(t) in result[0..p].
///
/// When a denominator is zero (repeated interior knots) the fraction is
/// taken as zero by convention; the safe-denominator trick avoids NaN in
/// the AD tape while keeping the derivative correct.
template <typename T>
[[nodiscard]] std::vector<T> basis_nonzero(T t, int k, int p,
                                            const std::vector<double>& knots)
{
    std::vector<T> N(p + 1, T(0.0));
    N[0] = T(1.0);

    std::vector<T> left (p + 1, T(0.0));
    std::vector<T> right(p + 1, T(0.0));

    for (int j = 1; j <= p; ++j) {
        left [j] = t - T(knots[k + 1 - j]);
        right[j] = T(knots[k + j]) - t;

        T saved = T(0.0);
        for (int r = 0; r < j; ++r) {
            T denom      = right[r + 1] + left[j - r];
            // safe_denom: replace near-zero denominator with 1 so division is
            // finite; CondExpGt then selects the numerically correct branch.
            T safe_denom = CppAD::CondExpGt(denom, T(1e-14), denom, T(1.0));
            T temp       = N[r] / safe_denom;
            T contrib    = CppAD::CondExpGt(denom, T(1e-14), temp, T(0.0));
            N[r]         = saved + right[r + 1] * contrib;
            saved        = left[j - r]           * contrib;
        }
        N[j] = saved;
    }
    return N;
}

// ─── Derivative of basis values ───────────────────────────────────────────────
/// Returns dN[0..p] — the derivatives of the p+1 non-zero degree-p basis
/// functions at span k, using the standard recurrence:
///   dN_{i,p}/dt = p * ( N_{i,p-1}/Δ₁  −  N_{i+1,p-1}/Δ₂ )
/// where Δ₁ = knots[i+p] − knots[i],  Δ₂ = knots[i+p+1] − knots[i+1].
template <typename T>
[[nodiscard]] std::vector<T> basis_deriv_nonzero(T t, int k, int p,
                                                  const std::vector<double>& knots)
{
    // Degree-(p-1) basis at the same span
    const std::vector<T> N = basis_nonzero<T>(t, k, p - 1, knots);
    // N[0..p-1] corresponds to N_{k-p+1, p-1} … N_{k, p-1}

    std::vector<T> dN(p + 1, T(0.0));
    for (int j = 0; j < p; ++j) {
        const double d = knots[k + j + 1] - knots[k + j + 1 - p];
        if (std::abs(d) > 1e-14) {
            const T factor = T(static_cast<double>(p) / d);
            dN[j]     -= factor * N[j];
            dN[j + 1] += factor * N[j];
        }
        // zero-knot-span: contribution is zero by convention — dN unchanged
    }
    return dN;
}

// ─── Clamped uniform knot vector ─────────────────────────────────────────────
/// Build a clamped, uniformly-spaced interior knot vector for a B-spline of
/// degree p with n control points.  Total length = n + p + 1.
[[nodiscard]] inline std::vector<double> clamped_knots(int p, int n)
{
    const int n_knots   = n + p + 1;
    const int n_interior = n - p - 1;    // may be 0

    std::vector<double> knots(static_cast<std::size_t>(n_knots), 0.0);

    for (int i = 1; i <= n_interior; ++i)
        knots[static_cast<std::size_t>(p + i)] = static_cast<double>(i) / static_cast<double>(n_interior + 1);

    for (int i = 0; i <= p; ++i)
        knots[static_cast<std::size_t>(n_knots - 1 - i)] = 1.0;

    return knots;
}

} // namespace BladeGen::Blade
