# tests/pipeline/CMakeLists.txt

add_executable(pcad_tests_pipeline
    test_blade_runner.cpp
)

target_link_libraries(pcad_tests_pipeline
    PRIVATE
        pcad::lib
        pcad::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(pcad_tests_pipeline
    DL_PATHS "$<TARGET_FILE_DIR:TKernel>;$<TARGET_FILE_DIR:cppad_lib>"
)
