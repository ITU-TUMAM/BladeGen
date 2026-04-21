# ------------------------------------------------------------------------------
# Project: BladeGen
# Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE
# ------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.25)

project(ParametricCAD
    VERSION 0.1.0
    DESCRIPTION "Parametric CAD template built on OpenCASCADE"
    LANGUAGES CXX
)

# ── C++ standard ────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Output directories ───────────────────────────────────────────────────────
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ── Module path ──────────────────────────────────────────────────────────────
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ── Build options ────────────────────────────────────────────────────────────
option(BLADEGEN_BUILD_TESTS   "Build the Catch2 test suite"  ON)
option(BLADEGEN_BUILD_DOCS    "Build Sphinx + Doxygen docs"  OFF)
option(BLADEGEN_ENABLE_ASAN   "Enable AddressSanitizer"      OFF)

# ── Compiler flags ───────────────────────────────────────────────────────────
include(CompilerFlags)

# ── Vendor dependencies ──────────────────────────────────────────────────────
include(VendorOCCT)

# CppAD — automatic differentiation (vendored on Windows, FetchContent on Linux)
include(VendorCppAD)

# nlohmann/json — header-only JSON, FetchContent on all platforms
include(VendorJson)

# Eigen — header-only, just expose the include path
add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/vendor/eigen"
)

# Catch2 — only needed when tests are enabled
if(BLADEGEN_BUILD_TESTS)
    set(CATCH_INSTALL_DOCS    OFF CACHE BOOL "" FORCE)
    set(CATCH_INSTALL_EXTRAS  OFF CACHE BOOL "" FORCE)
    add_subdirectory(vendor/catch2)
endif()

# ── Main library + application ───────────────────────────────────────────────
add_subdirectory(src)

# ── Tests ────────────────────────────────────────────────────────────────────
if(BLADEGEN_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# ── Docs ─────────────────────────────────────────────────────────────────────
if(BLADEGEN_BUILD_DOCS)
    include(Docs)
endif()