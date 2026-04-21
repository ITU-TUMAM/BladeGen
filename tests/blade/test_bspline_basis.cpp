// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "blade/BSplineBasis.hpp"

#include <numeric>

using namespace BladeGen::Blade;

// ─── find_span ────────────────────────────────────────────────────────────────
TEST_CASE("find_span returns correct span index", "[bspline][basis]")
{
    // Cubic clamped knots for 4 control points: {0,0,0,0,1,1,1,1}
    const auto knots = clamped_knots(3, 4);

    // For a clamped cubic with no interior knots, every t in [0,1] maps to span 3
    CHECK(find_span(0.0,  3, knots) == 3);
    CHECK(find_span(0.5,  3, knots) == 3);
    CHECK(find_span(0.99, 3, knots) == 3);
    CHECK(find_span(1.0,  3, knots) == 3);   // clamped upper boundary
}

TEST_CASE("find_span with interior knot", "[bspline][basis]")
{
    // Degree-4 clamped, 6 ctrl pts: {0,0,0,0,0,0.5,1,1,1,1,1}
    const auto knots = clamped_knots(4, 6);

    CHECK(find_span(0.0,   4, knots) == 4);
    CHECK(find_span(0.25,  4, knots) == 4);
    CHECK(find_span(0.499, 4, knots) == 4);
    CHECK(find_span(0.5,   4, knots) == 5);
    CHECK(find_span(0.75,  4, knots) == 5);
    CHECK(find_span(1.0,   4, knots) == 5);
}

// ─── Partition-of-unity property ─────────────────────────────────────────────
TEST_CASE("basis_nonzero sums to 1 everywhere (partition of unity)", "[bspline][basis]")
{
    // Cubic B-spline with one interior knot: {0,0,0,0,0.5,1,1,1,1}
    const auto knots = clamped_knots(3, 5);
    const int  p = 3;

    for (double t : {0.0, 0.1, 0.3, 0.49, 0.5, 0.7, 0.99, 1.0}) {
        const int k = find_span(t, p, knots);
        const auto N = basis_nonzero<double>(t, k, p, knots);

        double sum = std::accumulate(N.begin(), N.end(), 0.0);
        REQUIRE_THAT(sum, Catch::Matchers::WithinAbs(1.0, 1e-12));
    }
}

TEST_CASE("basis_nonzero non-negativity", "[bspline][basis]")
{
    const auto knots = clamped_knots(3, 5);
    const int  p = 3;

    for (double t : {0.0, 0.2, 0.5, 0.8, 1.0}) {
        const int k = find_span(t, p, knots);
        const auto N = basis_nonzero<double>(t, k, p, knots);
        for (const double n : N)
            CHECK(n >= -1e-15);   // allow tiny numerical noise
    }
}

// ─── Interpolation at endpoints ───────────────────────────────────────────────
TEST_CASE("clamped cubic interpolates endpoints", "[bspline][basis]")
{
    // 4 control points, clamped cubic: N_{0,3}(0)=1, N_{3,3}(1)=1
    const auto knots = clamped_knots(3, 4);

    {
        const int k  = find_span(0.0, 3, knots);
        const auto N = basis_nonzero<double>(0.0, k, 3, knots);
        REQUIRE_THAT(N[0], Catch::Matchers::WithinAbs(1.0, 1e-12));
        REQUIRE_THAT(N[1], Catch::Matchers::WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(N[2], Catch::Matchers::WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(N[3], Catch::Matchers::WithinAbs(0.0, 1e-12));
    }
    {
        const int k  = find_span(1.0, 3, knots);
        const auto N = basis_nonzero<double>(1.0, k, 3, knots);
        REQUIRE_THAT(N[0], Catch::Matchers::WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(N[1], Catch::Matchers::WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(N[2], Catch::Matchers::WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(N[3], Catch::Matchers::WithinAbs(1.0, 1e-12));
    }
}

// ─── Derivative via finite difference ─────────────────────────────────────────
TEST_CASE("basis_deriv_nonzero matches finite difference", "[bspline][basis]")
{
    const auto knots = clamped_knots(3, 5);
    const int  p = 3;
    const double h = 1e-6;

    for (double t0 : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        const int k = find_span(t0, p, knots);
        const auto dN_ad = basis_deriv_nonzero<double>(t0, k, p, knots);

        for (int j = 0; j <= p; ++j) {
            const int idx = k - p + j;

            // Finite-difference approximation of dN_{idx,p}/dt
            const int k1 = find_span(t0 + h, p, knots);
            const int k2 = find_span(t0 - h, p, knots);
            const auto Np = basis_nonzero<double>(t0 + h, k1, p, knots);
            const auto Nm = basis_nonzero<double>(t0 - h, k2, p, knots);

            // Extract N_{idx,p} from each evaluation
            auto get_N = [&](const std::vector<double>& N, int kk) -> double {
                const int j_local = idx - (kk - p);
                if (j_local < 0 || j_local > p) return 0.0;
                return N[j_local];
            };

            const double fd = (get_N(Np, k1) - get_N(Nm, k2)) / (2.0 * h);
            REQUIRE_THAT(dN_ad[j], Catch::Matchers::WithinAbs(fd, 1e-5));
        }
    }
}

// ─── clamped_knots structure ──────────────────────────────────────────────────
TEST_CASE("clamped_knots has correct size and boundary values", "[bspline][basis]")
{
    for (int p = 1; p <= 4; ++p) {
        for (int n = p + 1; n <= p + 4; ++n) {
            const auto k = clamped_knots(p, n);
            CHECK(static_cast<int>(k.size()) == n + p + 1);
            for (int i = 0; i <= p; ++i)   CHECK(k[i] == 0.0);
            for (int i = 0; i <= p; ++i)   CHECK(k[n + p - i] == 1.0);
        }
    }
}
