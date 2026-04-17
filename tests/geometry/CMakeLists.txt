# tests/geometry/CMakeLists.txt

add_executable(pcad_tests_geometry
    test_primitives.cpp
)

target_link_libraries(pcad_tests_geometry
    PRIVATE
        pcad::lib
        pcad::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(pcad_tests_geometry
    DL_PATHS "$<TARGET_FILE_DIR:TKernel>"
)
