// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
/// @file blade/OcctBlade.hpp
/// @brief OCCT-based blade solid construction from spanwise section loops.
///
/// Takes the 3D point arrays produced by BladeGeometry::all_sections_3d() and
/// skins them into a closed solid using OpenCASCADE ThruSections lofting.
///
/// OCCT internals stay in the .cpp — only TopoDS_Shape crosses the boundary.

#pragma once

#include <TopoDS_Shape.hxx>

#include <array>
#include <expected>
#include <vector>

namespace BladeGen::Blade {

/// Create a blade solid by skinning ordered spanwise profile loops.
///
/// @param sections   outer: span stations ordered hub→tip
///                   inner: closed 3D point loop for that section
///                   All loops must have the same point count.
/// @returns          a closed BRep solid, or an error string on OCCT failure.
[[nodiscard]] std::expected<TopoDS_Shape, std::string> MakeBladeSolid(
    const std::vector<std::vector<std::array<double, 3>>>& sections);

} // namespace BladeGen::Blade
