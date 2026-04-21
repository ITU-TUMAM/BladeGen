// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "blade/ProfileSection.hpp"
#include "blade/BladeSection.hpp"

#include <cmath>

using namespace BladeGen::Blade;

// ─── Camber line endpoint conditions ─────────────────────────────────────────
TEST_CASE("CamberLine interpolates LE and TE positions", "[blade][profile]")
{
    TurbineProfileParams<double> p;
    p.beta1     = 0.9;       // ~52°
    p.beta2     = -1.1;      // ~-63°
    p.theta_LE  = 0.0;
    p.theta_TE  = 0.05;
    p.t_max     = 0.015;
    p.t_max_loc = 0.35;
    p.t_TE      = 0.0;
    p.wedge_TE  = 0.05;
    p.r_LE      = 0.002;

    const auto cl = CamberLine<double>::from_params(p);

    REQUIRE_THAT(cl.eval(0.0), Catch::Matchers::WithinAbs(p.theta_LE, 1e-12));
    REQUIRE_THAT(cl.eval(1.0), Catch::Matchers::WithinAbs(p.theta_TE, 1e-12));
}

TEST_CASE("CamberLine satisfies inlet metal angle", "[blade][profile]")
{
    TurbineProfileParams<double> p;
    p.beta1     = 0.9;
    p.beta2     = -1.1;
    p.theta_LE  = 0.0;
    p.theta_TE  = 0.05;
    p.t_max     = 0.015;
    p.t_max_loc = 0.35;
    p.t_TE      = 0.0;
    p.wedge_TE  = 0.05;
    p.r_LE      = 0.002;

    const auto cl = CamberLine<double>::from_params(p);

    REQUIRE_THAT(cl.eval_deriv(0.0),
                 Catch::Matchers::WithinRel(std::tan(p.beta1), 1e-10));
}

TEST_CASE("CamberLine satisfies exit metal angle", "[blade][profile]")
{
    TurbineProfileParams<double> p;
    p.beta1     = 0.9;
    p.beta2     = -1.1;
    p.theta_LE  = 0.0;
    p.theta_TE  = 0.05;
    p.t_max     = 0.015;
    p.t_max_loc = 0.35;
    p.t_TE      = 0.0;
    p.wedge_TE  = 0.05;
    p.r_LE      = 0.002;

    const auto cl = CamberLine<double>::from_params(p);

    REQUIRE_THAT(cl.eval_deriv(1.0),
                 Catch::Matchers::WithinRel(std::tan(p.beta2), 1e-10));
}

// ─── Thickness distribution ───────────────────────────────────────────────────
TEST_CASE("ThicknessDistribution is zero at LE", "[blade][profile]")
{
    TurbineProfileParams<double> p;
    p.beta1 = 0.9; p.beta2 = -1.1;
    p.theta_LE = 0.0; p.theta_TE = 0.05;
    p.t_max = 0.015; p.t_max_loc = 0.35;
    p.t_TE = 0.0; p.wedge_TE = 0.05; p.r_LE = 0.002;

    const auto td = ThicknessDistribution<double>::from_params(p);
    REQUIRE_THAT(td.eval(0.0), Catch::Matchers::WithinAbs(0.0, 1e-12));
}

TEST_CASE("ThicknessDistribution equals t_TE at trailing edge", "[blade][profile]")
{
    TurbineProfileParams<double> p;
    p.beta1 = 0.9; p.beta2 = -1.1;
    p.theta_LE = 0.0; p.theta_TE = 0.05;
    p.t_max = 0.015; p.t_max_loc = 0.35;
    p.t_TE = 0.003; p.wedge_TE = 0.05; p.r_LE = 0.002;

    const auto td = ThicknessDistribution<double>::from_params(p);
    REQUIRE_THAT(td.eval(1.0), Catch::Matchers::WithinAbs(p.t_TE, 1e-10));
}

TEST_CASE("ThicknessDistribution is non-negative everywhere", "[blade][profile]")
{
    TurbineProfileParams<double> p;
    p.beta1 = 0.9; p.beta2 = -1.1;
    p.theta_LE = 0.0; p.theta_TE = 0.05;
    p.t_max = 0.015; p.t_max_loc = 0.35;
    p.t_TE = 0.0; p.wedge_TE = 0.05; p.r_LE = 0.002;

    const auto td = ThicknessDistribution<double>::from_params(p);

    for (int i = 0; i <= 100; ++i) {
        const double t = i / 100.0;
        CHECK(td.eval(t) >= -1e-12);
    }
}

// ─── ProfileSection: SS/PS symmetry ──────────────────────────────────────────
TEST_CASE("ProfileSection: SS and PS meet at LE and TE", "[blade][profile]")
{
    TurbineProfileParams<double> p;
    p.beta1 = 0.9; p.beta2 = -1.1;
    p.theta_LE = 0.0; p.theta_TE = 0.05;
    p.t_max = 0.015; p.t_max_loc = 0.35;
    p.t_TE = 0.0; p.wedge_TE = 0.05; p.r_LE = 0.002;

    const auto ps = ProfileSection<double>::from_params(p);

    // At LE (m'=0): zero thickness → SS == PS
    const auto ss_LE = ps.suction_side (0.0);
    const auto ps_LE = ps.pressure_side(0.0);
    REQUIRE_THAT(ss_LE.mp,    Catch::Matchers::WithinAbs(ps_LE.mp,    1e-10));
    REQUIRE_THAT(ss_LE.theta, Catch::Matchers::WithinAbs(ps_LE.theta, 1e-10));

    // At TE (m'=1): zero TE thickness → SS == PS
    const auto ss_TE = ps.suction_side (1.0);
    const auto ps_TE = ps.pressure_side(1.0);
    REQUIRE_THAT(ss_TE.mp,    Catch::Matchers::WithinAbs(ps_TE.mp,    1e-10));
    REQUIRE_THAT(ss_TE.theta, Catch::Matchers::WithinAbs(ps_TE.theta, 1e-10));
}

TEST_CASE("ProfileSection: suction side is above camber line", "[blade][profile]")
{
    TurbineProfileParams<double> p;
    p.beta1 = 0.9; p.beta2 = -1.1;
    p.theta_LE = 0.0; p.theta_TE = 0.05;
    p.t_max = 0.015; p.t_max_loc = 0.35;
    p.t_TE = 0.0; p.wedge_TE = 0.05; p.r_LE = 0.002;

    const auto ps = ProfileSection<double>::from_params(p);

    // At mid-chord the suction side theta should differ from camber theta
    const double mp_mid = 0.5;
    const double theta_c  = ps.camber.eval(mp_mid);
    const auto   ss_mid   = ps.suction_side(mp_mid);
    const auto   pres_mid = ps.pressure_side(mp_mid);

    // Suction and pressure sides are on opposite sides of the camber line
    // (their theta offsets are equal-and-opposite)
    const double delta_ss = ss_mid.theta   - theta_c;
    const double delta_ps = pres_mid.theta - theta_c;
    REQUIRE_THAT(delta_ss + delta_ps, Catch::Matchers::WithinAbs(0.0, 1e-10));
    CHECK(std::abs(delta_ss) > 1e-6);   // non-zero thickness at mid-chord
}

// ─── BladeGeometry factory ────────────────────────────────────────────────────
TEST_CASE("make_annular_blade produces correct section count", "[blade][geometry]")
{
    TurbineProfileParams<double> phub, ptip;
    phub.beta1=0.9; phub.beta2=-1.1; phub.theta_LE=0; phub.theta_TE=0.05;
    phub.t_max=0.015; phub.t_max_loc=0.35; phub.t_TE=0; phub.wedge_TE=0.05; phub.r_LE=0.002;
    ptip = phub;
    ptip.t_max = 0.010;

    const auto geom = make_annular_blade(0.25, 0.40, 0.05, phub, ptip, 3);
    CHECK(static_cast<int>(geom.sections.size()) == 5);   // 3 + 2
}

TEST_CASE("make_annular_blade: all sections produce 3D point loops", "[blade][geometry]")
{
    TurbineProfileParams<double> p;
    p.beta1=0.9; p.beta2=-1.1; p.theta_LE=0; p.theta_TE=0.05;
    p.t_max=0.015; p.t_max_loc=0.35; p.t_TE=0; p.wedge_TE=0.05; p.r_LE=0.002;

    const auto geom   = make_annular_blade(0.25, 0.40, 0.05, p, p, 2);
    const auto all_3d = geom.all_sections_3d(32);

    CHECK(all_3d.size() == geom.sections.size());
    for (const auto& sec : all_3d)
        CHECK(sec.size() == 64);   // 2 × 32 points
}

// ─── CppAD: gradient through camber + profile ─────────────────────────────────
TEST_CASE("CppAD records camber line and delivers dtheta/d_beta1", "[blade][profile][cppad]")
{
    // theta_TE is a fixed constant; beta2 controls where Q2_theta is placed.
    // For the 4-point cubic, theta_c(1) = theta_TE always (exact endpoint interpolation).
    // Test instead: dθ_c(0.5)/d_beta1 — sensitivity of mid-chord angle to inlet angle.

    using AD = CppAD::AD<double>;

    const double beta1_0 = 0.9;

    std::vector<AD> vars{AD(beta1_0)};
    CppAD::Independent(vars);

    TurbineProfileParams<AD> p;
    p.beta1     = vars[0];
    p.beta2     = AD(-1.1);
    p.theta_LE  = AD(0.0);
    p.theta_TE  = AD(0.05);
    p.t_max     = AD(0.015);
    p.t_max_loc = AD(0.35);
    p.t_TE      = AD(0.0);
    p.wedge_TE  = AD(0.05);
    p.r_LE      = AD(0.002);

    const auto cl = CamberLine<AD>::from_params(p);
    std::vector<AD> out{cl.eval(AD(0.5))};

    CppAD::ADFun<double> f(vars, out);
    const auto jac = f.Jacobian(std::vector<double>{beta1_0});

    // Finite-difference check
    const double h    = 1e-5;
    TurbineProfileParams<double> pd;
    pd.beta1=beta1_0+h; pd.beta2=-1.1; pd.theta_LE=0; pd.theta_TE=0.05;
    pd.t_max=0.015; pd.t_max_loc=0.35; pd.t_TE=0; pd.wedge_TE=0.05; pd.r_LE=0.002;
    TurbineProfileParams<double> pm;
    pm.beta1=beta1_0-h; pm.beta2=-1.1; pm.theta_LE=0; pm.theta_TE=0.05;
    pm.t_max=0.015; pm.t_max_loc=0.35; pm.t_TE=0; pm.wedge_TE=0.05; pm.r_LE=0.002;

    const auto clp = CamberLine<double>::from_params(pd);
    const auto clm = CamberLine<double>::from_params(pm);
    const double fd = (clp.eval(0.5) - clm.eval(0.5)) / (2.0 * h);

    REQUIRE_THAT(jac[0], Catch::Matchers::WithinRel(fd, 1e-6));
}
