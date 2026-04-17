/// @file blade/OcctBlade.cpp
/// @brief OCCT skinning implementation for the blade solid.
///
/// Pipeline:
///   1. For each spanwise section: fit a BSpline curve through the 3D points
///      using GeomAPI_PointsToBSpline.
///   2. Build a TopoDS_Wire from each curve edge.
///   3. Skin all wires with BRepOffsetAPI_ThruSections (solid + ruled=false).
///   4. Return the resulting TopoDS_Shape.

#include "OcctBlade.hpp"

// Foundation / geometry
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Geom_Curve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <NCollection_Array1.hxx>
#include <gp_Pnt.hxx>

// Shape
#include <TopoDS_Wire.hxx>

#include <stdexcept>

namespace PCAD::Blade {

std::optional<TopoDS_Shape> MakeBladeSolid(
    const std::vector<std::vector<std::array<double, 3>>>& sections)
{
    if (sections.size() < 2)
        return std::nullopt;

    const int n_pts = static_cast<int>(sections.front().size());
    for (const auto& sec : sections)
        if (static_cast<int>(sec.size()) != n_pts)
            return std::nullopt;   // all loops must have equal point count

    // ── Step 1+2: build a wire for each spanwise section ──────────────────────
    BRepOffsetAPI_ThruSections skinner(
        /*isSolid=*/true,
        /*isRuled=*/false);

    for (const auto& sec : sections) {
        // Pack 3D points into OCCT array (1-based indexing)
        NCollection_Array1<gp_Pnt> pts(1, n_pts);
        for (int i = 0; i < n_pts; ++i)
            pts.SetValue(i + 1, gp_Pnt(sec[i][0], sec[i][1], sec[i][2]));

        // Fit a BSpline through the loop; degree 3, continuity C2
        GeomAPI_PointsToBSpline fitter;
        fitter.Init(pts,
                    /*DegMin=*/3,
                    /*DegMax=*/8,
                    /*Continuity=*/GeomAbs_C2,
                    /*Tol3D=*/1.0e-6);

        if (!fitter.IsDone())
            return std::nullopt;

        const Handle(Geom_Curve) curve = fitter.Curve();

        BRepBuilderAPI_MakeEdge edge_maker(curve);
        if (!edge_maker.IsDone())
            return std::nullopt;

        BRepBuilderAPI_MakeWire wire_maker(edge_maker.Edge());
        if (!wire_maker.IsDone())
            return std::nullopt;

        skinner.AddWire(wire_maker.Wire());
    }

    // ── Step 3: skin all wires ─────────────────────────────────────────────────
    skinner.CheckCompatibility(false);
    skinner.Build();

    if (!skinner.IsDone())
        return std::nullopt;

    return skinner.Shape();
}

} // namespace PCAD::Blade
