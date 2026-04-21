# ------------------------------------------------------------------------------
# Project: BladeGen
# Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
# ------------------------------------------------------------------------------
# tests/math/CMakeLists.txt

add_executable(bladegen_tests_math
    test_cppad.cpp
)

target_link_libraries(bladegen_tests_math
    PRIVATE
        bladegen::cppad
        Eigen3::Eigen       # Windows cppad/configure.hpp has cppad_has_eigen=1;
                            # sparse2eigen.hpp then requires <Eigen/Sparse>.
        bladegen::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(bladegen_tests_math
    # On Linux cppad_lib is SHARED; provide its directory so the test runner
    # can find the .so even before install.
    DL_PATHS "$<TARGET_FILE_DIR:cppad_lib>"
)
