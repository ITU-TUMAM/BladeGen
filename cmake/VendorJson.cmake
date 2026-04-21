# cmake/VendorJson.cmake
# nlohmann/json — header-only JSON library, v3.11.3
#
# Downloaded via FetchContent on all platforms (no quirks; purely header-only).
# Tests and install rules are suppressed.
#
# Exposes alias target: bladegen::json

include(FetchContent)

FetchContent_Declare(nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
    GIT_SHALLOW    TRUE
)

set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
set(JSON_Install    OFF CACHE BOOL "" FORCE)
set(JSON_CI         OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(nlohmann_json)

add_library(bladegen::json ALIAS nlohmann_json)
