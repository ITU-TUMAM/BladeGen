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

#include <cmath>
#include <expected>
#include <stdexcept>

namespace PCAD::Blade {

// Remove consecutive duplicate points (within tol) from a closed loop.
// The BSpline fitter requires strictly positive chord-length increments.
static std::vector<std::array<double, 3>> deduplicate(
    const std::vector<std::array<double, 3>>& in, double tol = 1e-9)
{
    std::vector<std::array<double, 3>> out;
    out.reserve(in.size());
    for (const auto& p : in) {
        if (out.empty()) { out.push_back(p); continue; }
        const auto& q = out.back();
        const double d = std::sqrt((p[0]-q[0])*(p[0]-q[0]) +
                                   (p[1]-q[1])*(p[1]-q[1]) +
                                   (p[2]-q[2])*(p[2]-q[2]));
        if (d > tol) out.push_back(p);
    }
    return out;
}

std::expected<TopoDS_Shape, std::string> MakeBladeSolid(
    const std::vector<std::vector<std::array<double, 3>>>& sections)
{
    if (sections.size() < 2)
        return std::unexpected("need at least 2 sections");

    // ── Step 1+2: build a wire for each spanwise section ──────────────────────
    BRepOffsetAPI_ThruSections skinner(
        /*isSolid=*/true,
        /*isRuled=*/false);

    for (const auto& sec : sections) {
        // Remove consecutive duplicates (arise when LE/TE thickness == 0)
        const auto pts_vec = deduplicate(sec);
        const int n_pts = static_cast<int>(pts_vec.size());
        if (n_pts < 4) return std::unexpected("section has fewer than 4 unique points");

        // Pack 3D points into OCCT array (1-based indexing)
        NCollection_Array1<gp_Pnt> pts(1, n_pts);
        for (int i = 0; i < n_pts; ++i)
            pts.SetValue(i + 1, gp_Pnt(pts_vec[i][0], pts_vec[i][1], pts_vec[i][2]));

        // Fit a BSpline through the loop; degree 3, continuity C2
        GeomAPI_PointsToBSpline fitter;
        fitter.Init(pts,
                    /*DegMin=*/3,
                    /*DegMax=*/8,
                    /*Continuity=*/GeomAbs_C2,
                    /*Tol3D=*/1.0e-6);

        if (!fitter.IsDone())
            return std::unexpected("GeomAPI_PointsToBSpline failed");

        const Handle(Geom_Curve) curve = fitter.Curve();

        BRepBuilderAPI_MakeEdge edge_maker(curve);
        if (!edge_maker.IsDone())
            return std::unexpected("BRepBuilderAPI_MakeEdge failed");

        BRepBuilderAPI_MakeWire wire_maker(edge_maker.Edge());
        if (!wire_maker.IsDone())
            return std::unexpected("BRepBuilderAPI_MakeWire failed");

        skinner.AddWire(wire_maker.Wire());
    }

    // ── Step 3: skin all wires ─────────────────────────────────────────────────
    skinner.CheckCompatibility(false);
    skinner.Build();

    if (!skinner.IsDone())
        return std::unexpected("BRepOffsetAPI_ThruSections failed");

    return skinner.Shape();
}

} // namespace PCAD::Blade
