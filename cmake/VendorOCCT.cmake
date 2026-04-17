# cmake/VendorOCCT.cmake
# Configures OpenCASCADE Technology (OCCT) as a vendored subdirectory.
#
# Strategy
# --------
# Pull in only the geometry kernel and data exchange modules.
# All visualisation, scripting, and inspector modules are disabled —
# this keeps the Windows build fast and removes FreeType/Tcl/Tk deps.
#
# Modules kept ON:
#   FoundationClasses  — basic types, collections, math
#   ModelingData       — BRep, geometry primitives
#   ModelingAlgorithms — Boolean ops, fillets, sweeps
#   DataExchange       — STEP, IGES read/write
#
# Everything else is OFF.

set(OCCT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/vendor/occt")

if(NOT EXISTS "${OCCT_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
        "OpenCASCADE source not found at ${OCCT_SOURCE_DIR}.\n"
        "Did you run bootstrap.sh (Linux) or bootstrap.ps1 (Windows)?\n"
        "  git submodule update --init --recursive"
    )
endif()

# ── Disable all modules first ─────────────────────────────────────────────────
foreach(_mod IN ITEMS
    Draw Visualization ApplicationFramework DataExchange
    Inspector Samples Tests
)
    set(BUILD_MODULE_${_mod} OFF CACHE BOOL "" FORCE)
endforeach()

# ── Enable only what we need ──────────────────────────────────────────────────
set(BUILD_MODULE_FoundationClasses  ON  CACHE BOOL "" FORCE)
set(BUILD_MODULE_ModelingData       ON  CACHE BOOL "" FORCE)
set(BUILD_MODULE_ModelingAlgorithms ON  CACHE BOOL "" FORCE)
set(BUILD_MODULE_DataExchange       ON  CACHE BOOL "" FORCE)

# ── Strip optional third-party deps ──────────────────────────────────────────
set(USE_FREETYPE  OFF CACHE BOOL "" FORCE)
set(USE_TK        OFF CACHE BOOL "" FORCE)
set(USE_TCL       OFF CACHE BOOL "" FORCE)
set(USE_OPENGL    OFF CACHE BOOL "" FORCE)
set(USE_D3D       OFF CACHE BOOL "" FORCE)
set(USE_FREEIMAGE OFF CACHE BOOL "" FORCE)
set(USE_FFMPEG    OFF CACHE BOOL "" FORCE)
set(USE_OPENVR    OFF CACHE BOOL "" FORCE)
set(USE_RAPIDJSON OFF CACHE BOOL "" FORCE)
set(USE_DRACO     OFF CACHE BOOL "" FORCE)

# ── Redirect OCCT install targets to a throwaway subdirectory ────────────────
# Empty strings cause "install FILES given no DESTINATION!" on Linux/Ninja.
# Pointing to a build-tree subdirectory keeps OCCT's install() calls valid
# while keeping them out of the project's own install prefix.
set(_occt_install "${CMAKE_BINARY_DIR}/occt_install")
set(INSTALL_DIR        "${_occt_install}"           CACHE PATH "" FORCE)
set(INSTALL_DIR_BIN    "${_occt_install}/bin"       CACHE PATH "" FORCE)
set(INSTALL_DIR_LIB    "${_occt_install}/lib"       CACHE PATH "" FORCE)
set(INSTALL_DIR_INC    "${_occt_install}/inc"       CACHE PATH "" FORCE)
set(INSTALL_DIR_SCRIPT "${_occt_install}"           CACHE PATH "" FORCE)
set(BUILD_DOC_Overview OFF CACHE BOOL "" FORCE)

# ── Silence OCCT's own warnings so /WX doesn't trip on their headers ─────────
if(MSVC)
    add_compile_options(/experimental:external /external:anglebrackets /external:W0)
endif()

# ── Expose OCCT's generated headers before compiling its own targets ─────────
# OCCT copies all public headers to ${CMAKE_BINARY_DIR}/inc/ at configure time.
# Standard_Version.hxx (configure_file) lands in ${CMAKE_BINARY_DIR}/occt/inc/.
# Both paths are needed: the first for all source headers, the second for that
# one generated header. The global include_directories ensures OCCT's own
# toolkits find them; the INTERFACE on pcad_occt_kernel propagates to consumers.
include_directories(BEFORE
    "${CMAKE_BINARY_DIR}/inc"                        # Windows: flat copy by OCCT
    "${CMAKE_BINARY_DIR}/occt/inc"                   # Windows: Standard_Version.hxx (INSTALL_DIR_INCLUDE="inc")
    "${CMAKE_BINARY_DIR}/include/opencascade"        # Linux: collected source headers
    "${CMAKE_BINARY_DIR}/occt/include/opencascade"   # Linux: Standard_Version.hxx (INSTALL_DIR_INCLUDE="include/opencascade")
)

# ── Pull in the source tree ───────────────────────────────────────────────────
add_subdirectory("${OCCT_SOURCE_DIR}" occt EXCLUDE_FROM_ALL)

# ── Convenience alias target ──────────────────────────────────────────────────
# Consumer targets link against pcad::occt and get the right include paths
# and link libraries automatically.
add_library(pcad_occt_kernel INTERFACE)
add_library(pcad::occt ALIAS pcad_occt_kernel)

target_include_directories(pcad_occt_kernel INTERFACE
    "${CMAKE_BINARY_DIR}/inc"                        # Windows flat copy
    "${CMAKE_BINARY_DIR}/occt/inc"                   # Windows Standard_Version.hxx
    "${CMAKE_BINARY_DIR}/include/opencascade"        # Linux collected source headers
    "${CMAKE_BINARY_DIR}/occt/include/opencascade"   # Linux Standard_Version.hxx
)

target_link_libraries(pcad_occt_kernel INTERFACE
    # Foundation
    TKernel
    TKMath
    # Modeling
    TKBRep
    TKG2d
    TKG3d
    TKGeomBase
    TKGeomAlgo
    TKTopAlgo
    TKPrim
    TKShHealing
    TKBool
    TKBO
    TKFillet
    TKOffset
    # Data exchange (OCCT 7.8 names)
    TKXSBase
    TKDE
    TKDECascade
    TKDESTEP
    TKDEIGES
    TKDESTL
)