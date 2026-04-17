/// @file geometry/OcctUtils.cpp
/// @brief Implementation of OCCT geometry wrappers.
///
/// All OCCT headers are confined here — nothing outside this translation unit
/// should need to include OCCT directly for basic primitive creation/queries.

#include "geometry/OcctUtils.hpp"

// OCCT — BRep primitives
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeBox.hxx>

// OCCT — bounding box
#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>

// OCCT — shape type
#include <TopAbs_ShapeEnum.hxx>

namespace PCAD::Geometry {

// ── Primitives ───────────────────────────────────────────────────────────────

std::optional<TopoDS_Shape> MakeSphere(double radius)
{
    try {
        BRepPrimAPI_MakeSphere builder(radius);
        builder.Build();
        if (!builder.IsDone())
            return std::nullopt;
        return builder.Shape();
    }
    catch (...) {
        return std::nullopt;
    }
}

std::optional<TopoDS_Shape> MakeBox(double dx, double dy, double dz)
{
    try {
        BRepPrimAPI_MakeBox builder(dx, dy, dz);
        builder.Build();
        if (!builder.IsDone())
            return std::nullopt;
        return builder.Shape();
    }
    catch (...) {
        return std::nullopt;
    }
}

// ── Queries ──────────────────────────────────────────────────────────────────

std::optional<BoundingBox> GetBoundingBox(const TopoDS_Shape& shape)
{
    if (shape.IsNull())
        return std::nullopt;

    Bnd_Box bbox;
    BRepBndLib::Add(shape, bbox);

    if (bbox.IsVoid())
        return std::nullopt;

    double xMin, yMin, zMin, xMax, yMax, zMax;
    bbox.Get(xMin, yMin, zMin, xMax, yMax, zMax);

    return BoundingBox{ xMin, yMin, zMin, xMax, yMax, zMax };
}

std::string ShapeTypeString(const TopoDS_Shape& shape)
{
    switch (shape.ShapeType())
    {
        case TopAbs_COMPOUND:  return "Compound";
        case TopAbs_COMPSOLID: return "CompSolid";
        case TopAbs_SOLID:     return "Solid";
        case TopAbs_SHELL:     return "Shell";
        case TopAbs_FACE:      return "Face";
        case TopAbs_WIRE:      return "Wire";
        case TopAbs_EDGE:      return "Edge";
        case TopAbs_VERTEX:    return "Vertex";
        case TopAbs_SHAPE:     return "Shape";
        default:               return "Unknown";
    }
}

} // namespace PCAD::Geometry
