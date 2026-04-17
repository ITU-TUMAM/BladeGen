#!/usr/bin/env bash
# docs/venv-setup.sh
# Creates a reproducible Python virtual environment for the documentation
# pipeline (Sphinx + Breathe + Furo).
#
# Usage:
#   chmod +x docs/venv-setup.sh
#   ./docs/venv-setup.sh          # creates .venv/ in the repo root
#   source .venv/bin/activate
#   cmake --preset docs && cmake --build build/docs --target docs

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="${REPO_ROOT}/.venv"
REQUIREMENTS="${REPO_ROOT}/docs/requirements.txt"

echo "==> Checking Python 3 availability..."
if ! command -v python3 &>/dev/null; then
    echo "ERROR: python3 not found. Install Python 3.10 or newer." >&2
    exit 1
fi

PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo "    Found Python ${PYTHON_VERSION}"

echo "==> Creating virtual environment at ${VENV_DIR} ..."
python3 -m venv "${VENV_DIR}"

echo "==> Upgrading pip, setuptools, wheel ..."
"${VENV_DIR}/bin/pip" install --quiet --upgrade pip setuptools wheel

echo "==> Installing pinned documentation dependencies ..."
"${VENV_DIR}/bin/pip" install --quiet -r "${REQUIREMENTS}"

echo ""
echo "Done. Activate the environment with:"
echo "  source .venv/bin/activate"
echo ""
echo "Then build docs:"
echo "  cmake --preset docs"
echo "  cmake --build build/docs --target docs"