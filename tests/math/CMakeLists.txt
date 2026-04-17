# tests/math/CMakeLists.txt

add_executable(pcad_tests_math
    test_cppad.cpp
)

target_link_libraries(pcad_tests_math
    PRIVATE
        pcad::cppad
        Eigen3::Eigen       # Windows cppad/configure.hpp has cppad_has_eigen=1;
                            # sparse2eigen.hpp then requires <Eigen/Sparse>.
        pcad::compiler_flags
        Catch2::Catch2WithMain
)

include(Catch)
catch_discover_tests(pcad_tests_math
    # On Linux cppad_lib is SHARED; provide its directory so the test runner
    # can find the .so even before install.
    DL_PATHS "$<TARGET_FILE_DIR:cppad_lib>"
)
