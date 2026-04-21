// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
/// @file pipeline/BladeRunner.cpp
/// @brief BladeRunner implementation.

#include "pipeline/BladeRunner.hpp"

#include "blade/OcctBlade.hpp"
#include "io/IgesExporter.hpp"
#include "io/StepExporter.hpp"
#include "io/StlExporter.hpp"
#include "params/BladeParams.hpp"

#include <stdexcept>

namespace BladeGen::Pipeline {

BladeRunner::BladeRunner(BladeRunnerConfig cfg)
    : m_cfg(std::move(cfg))
{}

BladeRunnerResult BladeRunner::run(const Params::BladeParams& params) const
{
    params.validate();

    // ── Step 1: geometry ──────────────────────────────────────────────────────
    const auto geom     = params.to_blade_geometry();
    const auto sections = geom.all_sections_3d(m_cfg.profile_points,
                                               m_cfg.rk4_steps);

    // ── Step 2: OCCT solid ────────────────────────────────────────────────────
    const auto solid = Blade::MakeBladeSolid(sections);
    if (!solid)
        throw std::runtime_error("BladeRunner: OCCT skinning failed — " +
                                 solid.error());

    // ── Step 3: export ────────────────────────────────────────────────────────
    if (!m_cfg.export_iges && !m_cfg.export_step && !m_cfg.export_stl)
        throw std::runtime_error("BladeRunner: no export format enabled");

    std::filesystem::create_directories(m_cfg.output_dir);

    const auto stem = m_cfg.output_dir / m_cfg.stem;

    BladeRunnerResult result;
    result.n_sections      = sections.size();
    result.pts_per_section = sections.empty() ? 0 : sections.front().size();

    if (m_cfg.export_iges) {
        const auto path = std::filesystem::path(stem.string() + ".iges");
        IO::IgesExporter exp;
        exp.AddShape(*solid, m_cfg.stem);
        exp.Write(path);
        result.written.push_back(path);
    }

    if (m_cfg.export_step) {
        const auto path = std::filesystem::path(stem.string() + ".step");
        IO::StepExporter exp;
        exp.AddShape(*solid, m_cfg.stem);
        exp.Write(path);
        result.written.push_back(path);
    }

    if (m_cfg.export_stl) {
        const auto path = std::filesystem::path(stem.string() + ".stl");
        IO::StlExporter exp;
        exp.Write(*solid, path);
        result.written.push_back(path);
    }

    return result;
}

BladeRunnerResult BladeRunner::run(const std::filesystem::path& params_file) const
{
    return run(Params::BladeParams::from_file(params_file));
}

} // namespace BladeGen::Pipeline
