# tests/io/CMakeLists.txt

add_executable(pcad_tests_io
    test_step_export.cpp
    test_stl_export.cpp
)

target_link_libraries(pcad_tests_io
    PRIVATE
        pcad::lib
        pcad::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(pcad_tests_io
    DL_PATHS "$<TARGET_FILE_DIR:TKernel>"
)
