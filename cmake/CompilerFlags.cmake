# cmake/CompilerFlags.cmake
# Applies strict compiler flags to every target in the build.
# Call this once from the root CMakeLists — it sets an INTERFACE library
# that src/ and tests/ both link against.

add_library(pcad_compiler_flags INTERFACE)
add_library(pcad::compiler_flags ALIAS pcad_compiler_flags)

# ── Detect compiler family ────────────────────────────────────────────────────
set(IS_MSVC  $<CXX_COMPILER_ID:MSVC>)
set(IS_GNU   $<CXX_COMPILER_ID:GNU>)
set(IS_CLANG $<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>>)
set(IS_GCC_LIKE $<OR:${IS_GNU},${IS_CLANG}>)

# ── MSVC flags ────────────────────────────────────────────────────────────────
target_compile_options(pcad_compiler_flags INTERFACE
    $<${IS_MSVC}:
        /W4           # High warning level
        /permissive-  # Strict standards conformance
        /Zc:__cplusplus  # Report correct __cplusplus value
        /utf-8        # Source and execution charset
        /MP           # Multi-processor compilation
        /EHsc         # Standard C++ exception handling
        # Disable some noisy MSVC warnings that fire inside OCCT headers
        /wd4100       # unreferenced formal parameter (OCCT)
        /wd4127       # conditional expression is constant (OCCT)
        /wd4251       # dll-interface (OCCT export macros)
        /wd4275       # non dll-interface base class (OCCT)
    >
)

# ── GCC / Clang flags ────────────────────────────────────────────────────────
target_compile_options(pcad_compiler_flags INTERFACE
    $<${IS_GCC_LIKE}:
        -Wall
        -Wextra
        -Wshadow
        -Wcast-align
        -Wunused
        -Wconversion
        -Wsign-conversion
        -Wnull-dereference
        -Wdouble-promotion
        -Wformat=2
    >
    $<${IS_GNU}:
        -Wmisleading-indentation
        -Wduplicated-cond
        -Wduplicated-branches
        -Wlogical-op
    >
)

# ── AddressSanitizer ──────────────────────────────────────────────────────────
if(PCAD_ENABLE_ASAN)
    if(MSVC)
        message(FATAL_ERROR "ASan preset is not supported on MSVC — use linux-asan.")
    endif()
    target_compile_options(pcad_compiler_flags INTERFACE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
    )
    target_link_options(pcad_compiler_flags INTERFACE
        -fsanitize=address,undefined
    )
endif()

# ── Link-time optimisation (release builds only) ──────────────────────────────
include(CheckIPOSupported)
check_ipo_supported(RESULT _ipo_supported OUTPUT _ipo_msg)
if(_ipo_supported)
    target_compile_options(pcad_compiler_flags INTERFACE
        $<$<CONFIG:Release,RelWithDebInfo>:>
    )
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELWITHDEBINFO ON)
endif()

# ── Windows: ensure Unicode API layer ─────────────────────────────────────────
target_compile_definitions(pcad_compiler_flags INTERFACE
    $<${IS_MSVC}:
        UNICODE
        _UNICODE
        NOMINMAX          # Prevent windows.h defining min/max macros
        WIN32_LEAN_AND_MEAN
    >
)