# Contributing to ParametricCAD

Thank you for taking the time to contribute. This document covers everything
you need to know before opening a PR.

---

## Table of contents

- [Getting started](#getting-started)
- [Branch naming](#branch-naming)
- [Commit style](#commit-style)
- [Code style](#code-style)
- [Building and testing locally](#building-and-testing-locally)
- [Pull request checklist](#pull-request-checklist)
- [Reporting bugs](#reporting-bugs)
- [Requesting features](#requesting-features)

---

## Getting started

1. Fork the repository and clone your fork:
   ```bash
   git clone --recurse-submodules https://github.com/YOUR_USERNAME/ParametricCAD
   cd ParametricCAD
   ```

2. If you forgot `--recurse-submodules`:
   ```bash
   ./bootstrap.sh        # Linux
   .\bootstrap.ps1       # Windows (PowerShell)
   ```

3. Create a branch for your work (see [Branch naming](#branch-naming) below).

4. Make your changes, write tests, run the full suite locally before pushing.

---

## Branch naming

| Type | Pattern | Example |
|---|---|---|
| New feature | `feat/<short-description>` | `feat/step-importer` |
| Bug fix | `fix/<short-description>` | `fix/brep-null-deref` |
| Documentation | `docs/<short-description>` | `docs/io-module-api` |
| Refactor | `refactor/<short-description>` | `refactor/occt-wrapper` |
| CI / tooling | `ci/<short-description>` | `ci/windows-ninja-preset` |
| Chore | `chore/<short-description>` | `chore/bump-catch2-v3.7` |

Branch names must be lowercase, hyphen-separated, no spaces or underscores.

---

## Commit style

This project follows [Conventional Commits](https://www.conventionalcommits.org/).

```
<type>(<scope>): <short summary in imperative mood>

[optional body — wrap at 72 chars]

[optional footer — e.g. Closes #42]
```

**Types:** `feat`, `fix`, `docs`, `refactor`, `test`, `ci`, `chore`, `perf`

**Scope:** the module affected — `geometry`, `parametric`, `io`, `cmake`,
`vendor`, `tests`, `docs`, `ci`

**Examples:**

```
feat(io): add IGES export via STEPControl_Writer

fix(geometry): handle null TopoDS_Shape in OcctUtils::boundingBox

docs(parametric): add Doxygen examples to DesignVariable

ci: add windows-ninja-debug preset to matrix

chore(vendor): bump Catch2 to v3.7.0
```

Keep the summary line under 72 characters. Use the body to explain *why*,
not *what* — the diff shows what changed.

---

## Code style

### C++ formatting

All C++ code is formatted with `clang-format` using the config in
`.clang-format` at the repo root. The CI format check will fail if any
file is not correctly formatted.

Format before committing:

```bash
# Format all src/ and tests/ files in-place
find src tests -name '*.cpp' -o -name '*.hpp' | xargs clang-format -i
```

Or configure your editor to format on save using the project's
`.clang-format` file.

### C++ linting

`clang-tidy` is configured in `.clang-tidy`. Fix all warnings before
opening a PR — the CI tidy check runs with `-warnings-as-errors`.

```bash
# Run tidy on a single file (requires compile_commands.json)
clang-tidy src/geometry/OcctUtils.cpp -p build/linux-debug
```

### General rules

- C++20 — use modern features (`std::span`, ranges, concepts) where they
  improve clarity, not just for novelty.
- No raw owning pointers — use `std::unique_ptr` or `std::shared_ptr`.
- Prefer `const` by default. Mark every method that does not mutate state
  as `const`.
- No `using namespace` in headers — ever.
- All public API symbols must have a Doxygen comment (`///` style).
- Keep includes minimal — forward-declare where possible, especially in
  headers that pull in OCCT types.

---

## Building and testing locally

Always run the full test suite on at least the `linux-debug` preset before
pushing. For changes touching Windows-specific code, also verify
`windows-msvc-debug` builds clean.

```bash
# Configure
cmake --preset linux-debug

# Build
cmake --build build/linux-debug

# Test
ctest --preset linux-debug

# ASan run (catches memory errors)
cmake --preset linux-asan
cmake --build build/linux-asan
ctest --preset linux-asan
```

If you are adding or changing documentation:

```bash
./docs/venv-setup.sh          # first time only
source .venv/bin/activate
cmake --preset docs
cmake --build build/docs --target docs
# Review build/docs/sphinx/html/index.html
```

---

## Pull request checklist

Before marking your PR ready for review, confirm all of the following:

- [ ] `cmake --preset linux-debug && cmake --build build/linux-debug` succeeds
- [ ] `ctest --preset linux-debug` — all tests pass
- [ ] `cmake --preset linux-asan && ctest --preset linux-asan` — no sanitizer errors
- [ ] `clang-format` — no formatting diff (`git diff --exit-code` after formatting)
- [ ] `clang-tidy` — no warnings on changed files
- [ ] New public API has Doxygen `///` comments
- [ ] New behaviour is covered by at least one Catch2 test
- [ ] `CONTRIBUTING.md` branch naming and commit style followed
- [ ] PR description explains *why* the change is needed, not just *what* it does

For changes to the CMake build system, also verify:

- [ ] `cmake --preset windows-msvc-debug` configures without error
  (use a Windows machine or the CI run)

---

## Reporting bugs

Open an issue using the **Bug report** template. Please include:

- Operating system and compiler version
- CMake preset used
- The exact error message or unexpected behaviour
- A minimal reproducer if possible

---

## Requesting features

Open an issue using the **Feature request** template. Describe the use case
first — what are you trying to achieve and why is it not possible today?
Implementation ideas are welcome but optional.

---

## Questions

For questions that are not bugs or feature requests, open a
[Discussion](../../discussions) rather than an issue.