# BladeGen

> Parametric turbine blade geometry generation in modern C++23,
> built on [OpenCASCADE Technology (OCCT)](https://dev.opencascade.org/).

[![CI · Linux](https://github.com/onurtuncer/BladeGen/actions/workflows/ci-linux.yml/badge.svg)](https://github.com/onurtuncer/BladeGen/actions/workflows/ci-linux.yml)
[![CI · Windows](https://github.com/onurtuncer/BladeGen/actions/workflows/ci-windows.yml/badge.svg)](https://github.com/onurtuncer/BladeGen/actions/workflows/ci-windows.yml)
[![Docs](https://github.com/onurtuncer/BladeGen/actions/workflows/docs.yml/badge.svg)](https://onurtuncer.github.io/BladeGen/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

---

## What BladeGen does

BladeGen takes a compact JSON description of a turbine (or compressor) blade
and produces an OCCT solid that can be exported to IGES, STEP, or STL.

The full pipeline is:

```
JSON parameters
    │
    ▼
BladeParams          ← per-station inlet/exit angles, thickness, LE radius, …
    │
    ▼
ProfileSection       ← B-spline camber + thickness envelope in the 2-D (m′, θ) plane
    │
    ▼
StreamlineMPrime     ← RK4 integration of the meridional conformal coordinate m′
    │
    ▼
BladeSection / BladeGeometry   ← 3-D point arrays (hub → tip)
    │
    ▼
OcctBlade            ← OCCT lofted/skinned solid
    │
    ▼
IgesExporter / StepExporter / StlExporter
```

### Key modelling features

| Feature | Detail |
|---|---|
| Profile parameterisation | B-spline camber line + symmetric thickness envelope; parameters include `t_max`, `t_max_loc`, `t_TE`, TE wedge angle, and LE radius |
| Conformal coordinate | Meridional streamline is mapped to the `m′` conformal coordinate via RK4 integration, giving a proper area-rule-consistent design plane |
| Multi-span stacking | Arbitrary number of spanwise stations; intermediate sections are interpolated linearly; lean (Δθ) and sweep (Δx) stacking offsets are supported |
| Automatic differentiation | CppAD templates thread through the geometry layer, enabling gradient-based optimisation of any blade parameter |
| Export | IGES, STEP (AP203/AP214), and binary/ASCII STL via OCCT `DataExchange` |

---

## Repository layout

```
BladeGen/
├── cmake/
│   ├── CompilerFlags.cmake   # Strict cross-platform warning flags
│   ├── VendorOCCT.cmake      # OCCT subdirectory + bladegen::occt alias target
│   ├── VendorCppAD.cmake     # CppAD: pre-built (Windows) / FetchContent (Linux)
│   └── Docs.cmake            # Doxygen → Sphinx pipeline
├── vendor/
│   ├── occt/                 # git submodule — OpenCASCADE V7_8_0
│   ├── catch2/               # git submodule — Catch2 v3.6.0
│   ├── eigen/                # git submodule — Eigen 3.4.0
│   └── cppad/                # pre-built CppAD (Windows only; x64-Debug / x64-Release)
├── src/
│   ├── blade/                # ProfileSection, BladeSection, OcctBlade, B-spline basis
│   ├── geometry/             # OCCT primitive wrappers and bounding-box utilities
│   ├── io/                   # IGES / STEP / STL exporters
│   ├── params/               # BladeParams — JSON I/O and validation
│   ├── pipeline/             # BladeRunner — single-call orchestrator
│   └── main.cpp              # CLI entry point
├── tests/
│   ├── blade/                # Profile, B-spline, and m′ unit tests
│   ├── geometry/             # Primitive creation and bounding-box tests
│   ├── io/                   # IGES, STEP, and STL export tests
│   ├── math/                 # CppAD automatic differentiation tests
│   ├── params/               # BladeParams JSON round-trip tests
│   └── pipeline/             # BladeRunner end-to-end tests
├── docs/
│   ├── sphinx/               # Sphinx source (conf.py, index.rst, api.rst)
│   ├── requirements.txt      # Pinned Python doc dependencies
│   ├── venv-setup.sh         # Create .venv on Linux
│   └── venv-setup.ps1        # Create .venv on Windows
├── assets/                   # Logos and diagrams used in documentation
├── LICENSES/                 # Third-party license texts
├── bootstrap.sh              # Submodule init (Linux / macOS)
├── bootstrap.ps1             # Submodule init (Windows PowerShell)
├── CMakeLists.txt
├── CMakePresets.json
└── .gitmodules
```

---

## Prerequisites

### Linux

| Tool | Minimum |
|---|---|
| CMake | 3.25 |
| Ninja | any |
| GCC or Clang | GCC 12 / Clang 15 |
| Python 3 | 3.10 (docs only) |
| Doxygen | any (docs only) |

```bash
# Ubuntu 22.04 / 24.04
sudo apt install cmake ninja-build gcc g++ python3 python3-venv doxygen
```

### Windows

| Tool | Notes |
|---|---|
| Visual Studio 2022 | Desktop C++ workload required |
| CMake 3.25+ | Bundled with VS, or from cmake.org |
| Python 3.10+ | Docs only — from python.org |
| Doxygen | Docs only — from doxygen.nl |

> **PowerShell execution policy** — run once before the bootstrap script:
> ```powershell
> Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
> ```

---

## Quick start

### 1 — Clone

```bash
git clone https://github.com/onurtuncer/BladeGen
cd BladeGen
```

### 2 — Bootstrap (pulls all submodules)

```bash
# Linux
./bootstrap.sh

# Windows (PowerShell)
.\bootstrap.ps1
```

### 3 — Configure and build

```bash
# Linux — debug build with tests
cmake --preset linux-debug
cmake --build build/linux-debug

# Windows — MSVC debug build with tests
cmake --preset windows-msvc-debug
cmake --build build/windows-msvc-debug
```

### 4 — Run tests

```bash
# Linux
ctest --preset linux-debug --output-on-failure

# Windows
ctest --preset windows-msvc-debug --output-on-failure
```

### 5 — Generate a blade

Provide a JSON parameter file and point `BladeRunner` at it:

```cpp
#include "params/BladeParams.hpp"
#include "pipeline/BladeRunner.hpp"

BladeGen::Pipeline::BladeRunnerConfig cfg;
cfg.output_dir    = "output/";
cfg.export_iges   = true;
cfg.export_step   = true;
cfg.profile_points = 128;

auto result = BladeGen::Pipeline::BladeRunner(cfg).run("my_blade.json");
// result.written  — list of files produced
// result.n_sections, result.pts_per_section — geometry metadata
```

Or from the command line:

```bash
./build/linux-release/bladegen my_blade.json
```

---

## Available presets

| Preset | Platform | Compiler | Config | ASan |
|---|---|---|---|---|
| `linux-debug` | Linux | GCC / Clang | Debug | off |
| `linux-release` | Linux | GCC / Clang | RelWithDebInfo | off |
| `linux-asan` | Linux | GCC / Clang | Debug | on |
| `windows-msvc-debug` | Windows | MSVC x64 | Debug | off |
| `windows-msvc-release` | Windows | MSVC x64 | Release | off |
| `windows-ninja-debug` | Windows | clang-cl | Debug | off |
| `docs` | Linux | — | Release | off |

---

## Vendored dependencies

| Library | Version | License | Vendoring strategy | Purpose |
|---|---|---|---|---|
| [OpenCASCADE](https://github.com/Open-Cascade-SAS/OCCT) | `V7_8_0` | LGPL 2.1 | git submodule | Geometry kernel, STEP/IGES/STL I/O |
| [CppAD](https://github.com/coin-or/CppAD) | `20250000.3` | EPL-2.0 / GPL-2.0+ | pre-built (Windows) · FetchContent (Linux) | Automatic differentiation |
| [Catch2](https://github.com/catchorg/Catch2) | `v3.6.0` | BSL-1.0 | git submodule | Unit and integration testing |
| [Eigen](https://gitlab.com/libeigen/eigen) | `3.4.0` | MPL-2.0 | git submodule | Header-only linear algebra |

Only the OCCT geometry kernel modules are compiled:
`FoundationClasses`, `ModelingData`, `ModelingAlgorithms`, `DataExchange`.
Visualisation, scripting (Tcl/Tk), and inspector modules are all disabled.

To upgrade a submodule:

```bash
git -C vendor/<name> fetch --tags
git -C vendor/<name> checkout <new-tag>
git add vendor/<name>
git commit -m "chore(vendor): bump <name> to <new-tag>"
```

---

## Building documentation

```bash
# Linux
./docs/venv-setup.sh
source .venv/bin/activate
cmake --preset docs
cmake --build build/docs --target docs
# Output: build/docs/sphinx/html/index.html

# Windows (PowerShell)
.\docs\venv-setup.ps1
.\.venv\Scripts\Activate.ps1
cmake --preset docs
cmake --build build/docs --target docs
```

---

## License

GPL v3 — see [LICENSE](LICENSE).

Third-party licenses (OCCT, CppAD, Catch2, Eigen) are collected under
[LICENSES/](LICENSES/).

---

## Author

**Prof. Dr. Onur Tuncer**
Aerospace Engineer · Istanbul Technical University
Email: onur.tuncer@itu.edu.tr

<p align="left">
  <img src="assets/itu_logo.png" width="180" alt="Istanbul Technical University"/>
</p>
