/// @file geometry/OcctUtils.hpp
/// @brief Thin wrappers around OpenCASCADE geometry primitives.
///
/// All OCCT includes live in the corresponding .cpp — this header forward-
/// declares only what the rest of the codebase needs to know about, keeping
/// compile times sane and isolating the OCCT dependency.

#pragma once

#include <optional>
#include <string>
#include <tuple>

// TopoDS_Shape must be a complete type because it is stored inside std::optional.
#include <TopoDS_Shape.hxx>

namespace PCAD::Geometry {

// ── BoundingBox ──────────────────────────────────────────────────────────────

/// Axis-aligned bounding box with named accessors.
struct BoundingBox
{
    double xMin, yMin, zMin;
    double xMax, yMax, zMax;

    /// Returns (dx, dy, dz) — positive extents in each axis.
    [[nodiscard]] std::tuple<double, double, double> Extents() const noexcept
    {
        return { xMax - xMin, yMax - yMin, zMax - zMin };
    }
};

// ── Primitives ───────────────────────────────────────────────────────────────

/// Create a sphere of the given radius centred at the origin.
/// @returns the shape, or std::nullopt if OCCT reports a build failure.
[[nodiscard]] std::optional<TopoDS_Shape> MakeSphere(double radius);

/// Create an axis-aligned box with one corner at the origin.
/// @returns the shape, or std::nullopt if OCCT reports a build failure.
[[nodiscard]] std::optional<TopoDS_Shape> MakeBox(double dx, double dy, double dz);

// ── Queries ──────────────────────────────────────────────────────────────────

/// Compute the axis-aligned bounding box of a shape.
/// @returns nullopt if the shape is null or degenerate.
[[nodiscard]] std::optional<BoundingBox> GetBoundingBox(const TopoDS_Shape& shape);

/// Human-readable name for the top-level shape type (e.g. "Solid", "Shell").
[[nodiscard]] std::string ShapeTypeString(const TopoDS_Shape& shape);

} // namespace PCAD::Geometry
