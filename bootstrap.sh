#!/usr/bin/env bash
# bootstrap.sh
# Initialises all git submodules for ParametricCAD.
# Run this once after cloning if you did not use --recurse-submodules.
#
# Usage:
#   chmod +x bootstrap.sh
#   ./bootstrap.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_ROOT}"

echo "╔══════════════════════════════════════════════╗"
echo "║         ParametricCAD — Bootstrap            ║"
echo "╚══════════════════════════════════════════════╝"
echo ""

# ── Git check ────────────────────────────────────────────────────────────────
if ! command -v git &>/dev/null; then
    echo "ERROR: git not found on PATH." >&2
    exit 1
fi

# ── Submodules ────────────────────────────────────────────────────────────────
echo "==> Syncing submodule URLs from .gitmodules ..."
git submodule sync
echo ""
echo "==> Initialising and updating git submodules ..."
git submodule update --init --recursive --progress
echo "    Done."
echo ""

# ── Summary ──────────────────────────────────────────────────────────────────
echo "==> Submodule status:"
git submodule status
echo ""

echo "Bootstrap complete. Next steps:"
echo ""
echo "  Build (debug):"
echo "    cmake --preset linux-debug"
echo "    cmake --build build/linux-debug"
echo ""
echo "  Run tests:"
echo "    ctest --preset linux-debug"
echo ""
echo "  Build docs (optional):"
echo "    ./docs/venv-setup.sh"
echo "    source .venv/bin/activate"
echo "    cmake --preset docs"
echo "    cmake --build build/docs --target docs"