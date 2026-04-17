# vendor/CMakeLists.txt
#
# This directory is NOT added via add_subdirectory from the root.
# Each dependency is wired individually:
#
#   occt/    — configured by cmake/VendorOCCT.cmake
#              exposes alias target: pcad::occt
#
#   catch2/  — added via add_subdirectory in root CMakeLists.txt
#              when PCAD_BUILD_TESTS=ON
#              exposes targets: Catch2::Catch2WithMain
#
#   cppad/   — Windows only: pre-built static binaries (x64-Debug / x64-Release)
#              configured by cmake/VendorCppAD.cmake
#              On Linux the same module fetches the source via FetchContent.
#              exposes alias target: pcad::cppad
#
#   eigen/   — header-only, exposed as INTERFACE IMPORTED target
#              Eigen3::Eigen in root CMakeLists.txt
#
# To upgrade a dependency:
#   1. cd vendor/<name>
#   2. git fetch --tags
#   3. git checkout <new-tag>
#   4. cd ../..
#   5. git add vendor/<name>
#   6. Update the tag reference in .gitmodules and bootstrap.sh/.ps1
#   7. Build and test on both Linux and Windows before committing
#
# Never modify files inside vendor/ directly.
# Never commit a submodule pointing at a branch tip — always a tag.