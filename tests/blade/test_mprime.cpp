#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "blade/StreamlineMPrime.hpp"

#include <cmath>

using namespace PCAD::Blade;

// ─── Straight axial streamline at constant radius ────────────────────────────
// Streamline: x(t) = L·t,  r(t) = R  (constant)
// dm/dt  = L                 →  m(1)  = L
// dm'/dt = L/R               →  m'(1) = L/R
static BSplineCurve2D<double> straight_streamline(double L, double R)
{
    BSplineCurve2D<double> sl;
    sl.degree = 1;
    sl.knots  = clamped_knots(1, 2);   // {0, 0, 1, 1}
    sl.Px     = {0.0, L};
    sl.Pr     = {R,   R};
    return sl;
}

TEST_CASE("m' integration: straight axial streamline at constant r",
          "[blade][mprime]")
{
    const double L = 0.05;   // 50 mm axial chord
    const double R = 0.30;   // 300 mm radius (hub-ish)

    const auto sl  = straight_streamline(L, R);
    const auto tbl = integrate_mprime<double>(sl, 200);

    const double mp_expected = L / R;
    const double m_expected  = L;

    REQUIRE_THAT(tbl.mprime_vals.back(),
                 Catch::Matchers::WithinRel(mp_expected, 1e-6));
    REQUIRE_THAT(tbl.m_vals.back(),
                 Catch::Matchers::WithinRel(m_expected,  1e-6));
}

TEST_CASE("m' table is monotonically increasing", "[blade][mprime]")
{
    const auto sl  = straight_streamline(0.04, 0.35);
    const auto tbl = integrate_mprime<double>(sl, 100);

    for (std::size_t i = 1; i < tbl.mprime_vals.size(); ++i)
        CHECK(tbl.mprime_vals[i] > tbl.mprime_vals[i - 1]);
}

TEST_CASE("m' starts at zero", "[blade][mprime]")
{
    const auto sl  = straight_streamline(0.04, 0.35);
    const auto tbl = integrate_mprime<double>(sl, 50);
    REQUIRE_THAT(tbl.mprime_vals.front(),
                 Catch::Matchers::WithinAbs(0.0, 1e-15));
}

TEST_CASE("invert_mprime recovers parameter t", "[blade][mprime]")
{
    const double L = 0.05, R = 0.30;
    const auto sl  = straight_streamline(L, R);
    const auto tbl = integrate_mprime<double>(sl, 200);

    // For a straight streamline m'(t) = (L/R)·t, so t = mp * R/L
    const double mp_max = tbl.mprime_vals.back();
    for (double frac : {0.1, 0.3, 0.5, 0.7, 0.9}) {
        const double mp_target = frac * mp_max;
        const double t_inv     = invert_mprime<double>(mp_target, tbl);
        REQUIRE_THAT(t_inv, Catch::Matchers::WithinAbs(frac, 1e-4));
    }
}

// ─── conformal_to_3d: straight streamline at R, x along axis ─────────────────
TEST_CASE("conformal_to_3d maps to expected 3D coordinates", "[blade][mprime]")
{
    const double L = 0.05, R = 0.30;
    const auto sl  = straight_streamline(L, R);
    const auto tbl = integrate_mprime<double>(sl, 200);

    // At mp_norm=0.5 on the straight streamline: x = L/2
    const double theta = 0.0;
    const auto [x, y, z] = conformal_to_3d<double>(0.5, theta, sl, tbl);

    REQUIRE_THAT(x, Catch::Matchers::WithinRel(L / 2.0, 1e-4));
    REQUIRE_THAT(y, Catch::Matchers::WithinAbs(0.0, 1e-10));    // sin(0)=0
    REQUIRE_THAT(z, Catch::Matchers::WithinRel(R,   1e-4));     // cos(0)=1
}

// ─── RK4 convergence: halving steps reduces error by ~16× ───────────────────
TEST_CASE("RK4 integration converges at 4th order", "[blade][mprime]")
{
    const double L = 0.05, R = 0.30;
    const auto sl = straight_streamline(L, R);

    const double ref = integrate_mprime<double>(sl, 1000).mprime_vals.back();
    const double e50 = std::abs(integrate_mprime<double>(sl,  50).mprime_vals.back() - ref);
    const double e100= std::abs(integrate_mprime<double>(sl, 100).mprime_vals.back() - ref);

    // Ratio should be ≈ 2^4 = 16 for RK4
    if (e100 > 1e-15)
        CHECK((e50 / e100) > 8.0);   // conservative: allow for double-precision floor
}

// ─── CppAD: derivatives through RK4 integration ──────────────────────────────
TEST_CASE("CppAD records m' integral and delivers exact dmp/dPr", "[blade][mprime][cppad]")
{
    // Parameterise a straight streamline at constant R.
    // m'_max = L/R.  dm'_max/dR = -L/R^2.

    using AD = CppAD::AD<double>;

    const double L = 0.05;
    const double R0 = 0.30;

    std::vector<AD> params{AD(R0)};   // single design variable: R
    CppAD::Independent(params);

    BSplineCurve2D<AD> sl;
    sl.degree = 1;
    sl.knots  = clamped_knots(1, 2);
    sl.Px     = {AD(0.0), AD(L)};
    sl.Pr     = {params[0], params[0]};   // R is the AD variable

    const auto tbl = integrate_mprime<AD>(sl, 100);

    std::vector<AD> mp_out{tbl.mprime_vals.back()};
    CppAD::ADFun<double> f(params, mp_out);

    const std::vector<double> jac = f.Jacobian(std::vector<double>{R0});

    const double analytical = -L / (R0 * R0);
    REQUIRE_THAT(jac[0], Catch::Matchers::WithinRel(analytical, 1e-5));
}
