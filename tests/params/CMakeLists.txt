# ------------------------------------------------------------------------------
# Project: BladeGen
# Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
# ------------------------------------------------------------------------------
# tests/params/CMakeLists.txt

add_executable(bladegen_tests_params
    test_blade_params.cpp
)

target_link_libraries(bladegen_tests_params
    PRIVATE
        bladegen::lib
        bladegen::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(bladegen_tests_params
    DL_PATHS "$<TARGET_FILE_DIR:cppad_lib>"
)
