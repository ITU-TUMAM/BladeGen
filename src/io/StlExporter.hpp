// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
/// @file io/StlExporter.hpp
/// @brief STL file export via OpenCASCADE StlAPI_Writer.
///
/// Meshes a shape with BRepMesh_IncrementalMesh before writing — this step
/// is required and is a common source of empty STL files when skipped.

#pragma once

#include <filesystem>

class TopoDS_Shape;

namespace BladeGen::IO {

class StlExporter
{
public:
    StlExporter() = default;

    /// Linear deflection controls triangle density (smaller = finer mesh).
    /// Default: 0.1
    void SetLinearDeflection(double deflection) noexcept
    {
        m_linearDeflection = deflection;
    }

    /// Angular deflection in radians. Default: 0.5
    void SetAngularDeflection(double radians) noexcept
    {
        m_angularDeflection = radians;
    }

    /// Mesh and write @p shape to @p path.
    /// @throws std::runtime_error on OCCT write failure.
    void Write(const TopoDS_Shape& shape, const std::filesystem::path& path) const;

private:
    double m_linearDeflection  = 0.1;
    double m_angularDeflection = 0.5;
};

} // namespace BladeGen::IO
