# docs/venv-setup.ps1
# Creates a reproducible Python virtual environment for the documentation
# pipeline (Sphinx + Breathe + Furo) on Windows.
#
# Usage (from repo root in PowerShell):
#   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
#   .\docs\venv-setup.ps1
#   .\.venv\Scripts\Activate.ps1
#   cmake --preset docs
#   cmake --build build/docs --target docs

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$RepoRoot    = Resolve-Path "$PSScriptRoot\.."
$VenvDir     = Join-Path $RepoRoot ".venv"
$Requirements = Join-Path $RepoRoot "docs\requirements.txt"

Write-Host "==> Checking Python 3 availability..."
try {
    $PyVersion = python --version 2>&1
    Write-Host "    Found $PyVersion"
} catch {
    Write-Error "Python 3 not found. Install Python 3.10+ from https://python.org and ensure it is on PATH."
    exit 1
}

Write-Host "==> Creating virtual environment at $VenvDir ..."
python -m venv $VenvDir

$Pip = Join-Path $VenvDir "Scripts\pip.exe"

Write-Host "==> Upgrading pip, setuptools, wheel ..."
& $Pip install --quiet --upgrade pip setuptools wheel

Write-Host "==> Installing pinned documentation dependencies ..."
& $Pip install --quiet -r $Requirements

Write-Host ""
Write-Host "Done. Activate the environment with:"
Write-Host "  .\.venv\Scripts\Activate.ps1"
Write-Host ""
Write-Host "Then build docs:"
Write-Host "  cmake --preset docs"
Write-Host "  cmake --build build/docs --target docs"