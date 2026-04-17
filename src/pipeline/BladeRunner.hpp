/// @file pipeline/BladeRunner.hpp
/// @brief Single-object orchestrator for the full blade generation pipeline.
///
/// BladeRunner owns the configuration; calling run() executes the chain:
///   BladeParams → BladeGeometry → OCCT solid → export files
///
/// It is intentionally stateless between runs — run() takes all inputs as
/// arguments and returns all outputs in BladeRunnerResult, so it is safe to
/// call from multiple threads with different params objects.

#pragma once

#include <filesystem>
#include <string>
#include <vector>

// Forward declarations — OCCT and geometry types stay out of this header.
namespace PCAD::Params { struct BladeParams; }

namespace PCAD::Pipeline {

// ─── Configuration ────────────────────────────────────────────────────────────
/// Controls how BladeRunner samples geometry and which files it produces.
struct BladeRunnerConfig
{
    /// Points sampled per side of the profile (total loop = 2 × profile_points).
    int profile_points = 64;

    /// RK4 integration steps for the m′ conformal coordinate per section.
    int rk4_steps = 200;

    /// Directory into which output files are written.  Created if absent.
    std::filesystem::path output_dir = ".";

    /// Base filename stem (no extension).  e.g. "blade" → blade.iges, blade.step
    std::string stem = "blade";

    bool export_iges = true;   ///< Write an IGES file
    bool export_step = false;  ///< Write a STEP file
    bool export_stl  = false;  ///< Write an STL  file
};

// ─── Result ───────────────────────────────────────────────────────────────────
/// What run() produced.
struct BladeRunnerResult
{
    std::size_t n_sections        = 0;  ///< Span sections in the solid
    std::size_t pts_per_section   = 0;  ///< Points in each closed loop
    std::vector<std::filesystem::path> written;  ///< Files successfully written
};

// ─── Runner ───────────────────────────────────────────────────────────────────
class BladeRunner
{
public:
    explicit BladeRunner(BladeRunnerConfig cfg = {});

    /// Execute the full pipeline from an already-loaded parameter set.
    /// @throws std::invalid_argument if params fails validation.
    /// @throws std::runtime_error    if OCCT skinning or any export fails.
    [[nodiscard]] BladeRunnerResult run(const Params::BladeParams& params) const;

    /// Convenience overload: load params from a JSON file, then run.
    /// @throws std::runtime_error on file / parse errors.
    [[nodiscard]] BladeRunnerResult run(const std::filesystem::path& params_file) const;

    [[nodiscard]] const BladeRunnerConfig& config() const noexcept { return m_cfg; }

private:
    BladeRunnerConfig m_cfg;
};

} // namespace PCAD::Pipeline
