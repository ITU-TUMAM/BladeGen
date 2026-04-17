# cmake/Docs.cmake
# Documentation pipeline:
#   1. Doxygen  → XML (in build/docs/doxygen/xml)
#   2. Sphinx   → HTML via Breathe (reads Doxygen XML)
#
# Requirements (installed in the Python env):
#   pip install sphinx breathe sphinx-rtd-theme
# Or for a nicer theme:
#   pip install sphinx breathe furo
#
# Usage:
#   cmake --preset docs
#   cmake --build build/docs --target docs

find_package(Doxygen REQUIRED)

# Prefer the project-local .venv so CI and contributors get identical
# Sphinx/Breathe versions from docs/requirements.txt.
set(_venv_python_linux "${CMAKE_SOURCE_DIR}/.venv/bin/python3")
set(_venv_python_win   "${CMAKE_SOURCE_DIR}/.venv/Scripts/python.exe")

if(EXISTS "${_venv_python_linux}")
    set(Python3_EXECUTABLE "${_venv_python_linux}")
elseif(EXISTS "${_venv_python_win}")
    set(Python3_EXECUTABLE "${_venv_python_win}")
else()
    message(WARNING
        "Project .venv not found. Falling back to system Python.\n"
        "For a reproducible build run:\n"
        "  Linux:   ./docs/venv-setup.sh\n"
        "  Windows: .\\docs\\venv-setup.ps1"
    )
endif()

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# ── Doxygen configuration ────────────────────────────────────────────────────
set(DOXYGEN_PROJECT_NAME    "${PROJECT_NAME}")
set(DOXYGEN_PROJECT_VERSION "${PROJECT_VERSION}")
set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/doxygen")
set(DOXYGEN_GENERATE_HTML   NO)   # Sphinx will render HTML instead
set(DOXYGEN_GENERATE_XML    YES)  # Breathe reads this
set(DOXYGEN_GENERATE_LATEX  NO)
set(DOXYGEN_EXTRACT_ALL     YES)
set(DOXYGEN_EXTRACT_PRIVATE NO)
set(DOXYGEN_QUIET           YES)
set(DOXYGEN_WARN_AS_ERROR   YES)

# Do not recurse into vendor/
set(DOXYGEN_EXCLUDE_PATTERNS
    "*/vendor/*"
    "*/build/*"
    "*/tests/*"
)

# Enable Markdown support and brief member docs
set(DOXYGEN_MARKDOWN_SUPPORT     YES)
set(DOXYGEN_AUTOLINK_SUPPORT     YES)
set(DOXYGEN_BRIEF_MEMBER_DESC    YES)
set(DOXYGEN_JAVADOC_AUTOBRIEF    YES)

doxygen_add_docs(doxygen_xml
    "${CMAKE_SOURCE_DIR}/src"
    COMMENT "Generating Doxygen XML for Breathe"
)

# ── Sphinx configuration ──────────────────────────────────────────────────────
set(SPHINX_SOURCE_DIR "${CMAKE_SOURCE_DIR}/docs/sphinx")
set(SPHINX_BUILD_DIR  "${CMAKE_BINARY_DIR}/sphinx")

# Check sphinx-build is available
execute_process(
    COMMAND ${Python3_EXECUTABLE} -m sphinx --version
    RESULT_VARIABLE _sphinx_check
    OUTPUT_QUIET ERROR_QUIET
)
if(NOT _sphinx_check EQUAL 0)
    message(FATAL_ERROR
        "sphinx-build not found.\n"
        "Install documentation dependencies:\n"
        "  pip install sphinx breathe furo"
    )
endif()

add_custom_target(docs
    COMMAND ${Python3_EXECUTABLE} -m sphinx
        -b html
        -D "breathe_projects.${PROJECT_NAME}=${DOXYGEN_OUTPUT_DIRECTORY}/xml"
        "${SPHINX_SOURCE_DIR}"
        "${SPHINX_BUILD_DIR}/html"
    DEPENDS doxygen_xml
    COMMENT "Building Sphinx HTML documentation"
    VERBATIM
)

message(STATUS "Docs target configured — run: cmake --build . --target docs")