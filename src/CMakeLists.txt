# ------------------------------------------------------------------------------
# Project: BladeGen
# Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
# ------------------------------------------------------------------------------
# src/CMakeLists.txt

# ── Core library ─────────────────────────────────────────────────────────────
# All geometry, parametric, and IO code compiles into a static library.
# The CLI executable and the test suite both link against it, so tests
# can exercise every module without going through main().

add_library(bladegen_lib STATIC
    geometry/OcctUtils.cpp
    io/StepExporter.cpp
    io/StlExporter.cpp
    io/IgesExporter.cpp
    blade/OcctBlade.cpp
    params/BladeParams.cpp
    pipeline/BladeRunner.cpp
)

add_library(bladegen::lib ALIAS bladegen_lib)

target_include_directories(bladegen_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(bladegen_lib
    PUBLIC
        bladegen::occt
        Eigen3::Eigen
        bladegen::cppad
        bladegen::json
    PRIVATE
        bladegen::compiler_flags
)

# ── CLI executable ────────────────────────────────────────────────────────────
add_executable(bladegen_app main.cpp)
add_executable(bladegen::app ALIAS bladegen_app)

set_target_properties(bladegen_app PROPERTIES OUTPUT_NAME "bladegen")

target_link_libraries(bladegen_app
    PRIVATE
        bladegen::lib
        bladegen::compiler_flags
)