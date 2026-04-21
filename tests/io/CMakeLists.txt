# ------------------------------------------------------------------------------
# Project: BladeGen
# Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
# ------------------------------------------------------------------------------
# tests/io/CMakeLists.txt

add_executable(bladegen_tests_io
    test_step_export.cpp
    test_stl_export.cpp
    test_iges_export.cpp
)

target_link_libraries(bladegen_tests_io
    PRIVATE
        bladegen::lib
        bladegen::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(bladegen_tests_io
    DL_PATHS "$<TARGET_FILE_DIR:TKernel>"
)
