# docs/sphinx/conf.py
# Sphinx configuration for ParametricCAD.
#
# breathe_projects is injected at build time by cmake/Docs.cmake:
#   sphinx -D "breathe_projects.ParametricCAD=<path/to/doxygen/xml>" ...
# so it is left empty here and populated by CMake.

project   = "ParametricCAD"
author    = "Onur Tuncer"
release   = "0.1.0"
copyright = "2024, Onur Tuncer"

extensions = [
    "breathe",
    "sphinx_copybutton",
]

# ── Breathe ───────────────────────────────────────────────────────────────────
breathe_projects        = {}          # filled by CMake -D flag
breathe_default_project = "ParametricCAD"

# ── HTML output ───────────────────────────────────────────────────────────────
html_theme = "furo"
html_title = "ParametricCAD"

html_theme_options = {
    "source_repository": "https://github.com/onurtuncer/ParametricCAD/",
    "source_branch": "main",
    "source_directory": "docs/sphinx/",
}
