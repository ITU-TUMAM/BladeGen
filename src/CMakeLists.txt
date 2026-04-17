# src/CMakeLists.txt

# ── Core library ─────────────────────────────────────────────────────────────
# All geometry, parametric, and IO code compiles into a static library.
# The CLI executable and the test suite both link against it, so tests
# can exercise every module without going through main().

add_library(pcad_lib STATIC
    geometry/OcctUtils.cpp
    io/StepExporter.cpp
    io/StlExporter.cpp
    io/IgesExporter.cpp
    blade/OcctBlade.cpp
    params/BladeParams.cpp
)

add_library(pcad::lib ALIAS pcad_lib)

target_include_directories(pcad_lib
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(pcad_lib
    PUBLIC
        pcad::occt
        Eigen3::Eigen
        pcad::cppad
        pcad::json
    PRIVATE
        pcad::compiler_flags
)

# ── CLI executable ────────────────────────────────────────────────────────────
add_executable(pcad_app main.cpp)
add_executable(pcad::app ALIAS pcad_app)

set_target_properties(pcad_app PROPERTIES OUTPUT_NAME "pcad")

target_link_libraries(pcad_app
    PRIVATE
        pcad::lib
        pcad::compiler_flags
)