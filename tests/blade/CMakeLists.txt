# ------------------------------------------------------------------------------
# Project: BladeGen
# Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
# ------------------------------------------------------------------------------
# tests/blade/CMakeLists.txt

add_executable(bladegen_tests_blade
    test_bspline_basis.cpp
    test_mprime.cpp
    test_profile.cpp
)

target_link_libraries(bladegen_tests_blade
    PRIVATE
        bladegen::lib            # brings in CppAD, OCCT, Eigen
        bladegen::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(bladegen_tests_blade
    DL_PATHS "$<TARGET_FILE_DIR:cppad_lib>"
)
