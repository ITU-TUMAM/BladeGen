// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cppad/cppad.hpp>
#include <vector>

using CppAD::AD;

// Record f(x) = x^2 and verify f'(3) == 6 via reverse-mode Jacobian.
TEST_CASE("CppAD differentiates x^2", "[cppad][autodiff]") {
    std::vector<AD<double>> ax{3.0};
    CppAD::Independent(ax);

    std::vector<AD<double>> ay{ax[0] * ax[0]};
    CppAD::ADFun<double> f(ax, ay);

    std::vector<double> jac = f.Jacobian(std::vector<double>{3.0});

    REQUIRE_THAT(jac[0], Catch::Matchers::WithinRel(6.0, 1e-12));
}
