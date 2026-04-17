# tests/blade/CMakeLists.txt

add_executable(pcad_tests_blade
    test_bspline_basis.cpp
    test_mprime.cpp
    test_profile.cpp
)

target_link_libraries(pcad_tests_blade
    PRIVATE
        pcad::lib            # brings in CppAD, OCCT, Eigen
        pcad::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(pcad_tests_blade
    DL_PATHS "$<TARGET_FILE_DIR:cppad_lib>"
)
