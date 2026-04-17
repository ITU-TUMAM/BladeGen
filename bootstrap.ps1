# bootstrap.ps1
# Initialises all git submodules for ParametricCAD.
# Run this once after cloning if you did not use --recurse-submodules.
#
# Usage (from repo root in PowerShell):
#   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
#   .\bootstrap.ps1

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RepoRoot = Resolve-Path "$PSScriptRoot"
Set-Location $RepoRoot

Write-Host "================================================"
Write-Host "=         ParametricCAD - Bootstrap            ="
Write-Host "================================================"
Write-Host ""

# ── Git check ────────────────────────────────────────────────────────────────
if (-not (Get-Command git -ErrorAction SilentlyContinue)) {
    Write-Error "git not found on PATH. Install Git for Windows: https://git-scm.com"
    exit 1
}
Write-Host "==> Found $(git --version)"

# ── Submodules ────────────────────────────────────────────────────────────────
Write-Host ""
Write-Host "==> Syncing submodule URLs from .gitmodules ..."
git submodule sync
Write-Host ""
Write-Host "==> Initialising and updating git submodules ..."
git submodule update --init --recursive --progress
Write-Host "    Done."
Write-Host ""

# ── Summary ──────────────────────────────────────────────────────────────────
Write-Host "==> Submodule status:"
git submodule status
Write-Host ""

Write-Host "Bootstrap complete. Next steps:"
Write-Host ""
Write-Host "  Build (MSVC Debug):"
Write-Host "    cmake --preset windows-msvc-debug"
Write-Host "    cmake --build build/windows-msvc-debug"
Write-Host ""
Write-Host "  Run tests:"
Write-Host "    ctest --preset windows-msvc-debug"
Write-Host ""
Write-Host "  Build docs (optional):"
Write-Host "    .\docs\venv-setup.ps1"
Write-Host "    .\.venv\Scripts\Activate.ps1"
Write-Host "    cmake --preset docs"
Write-Host "    cmake --build build/docs --target docs"