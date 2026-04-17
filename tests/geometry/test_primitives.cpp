/// @file tests/geometry/test_primitives.cpp
/// @brief Unit tests for PCAD::Geometry primitives and bounding box queries.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "geometry/OcctUtils.hpp"

#include <TopoDS_Shape.hxx>

using namespace Catch::Matchers;

// ── MakeSphere ────────────────────────────────────────────────────────────────

TEST_CASE("MakeSphere returns a valid shape", "[geometry][sphere]")
{
    const auto shape = PCAD::Geometry::MakeSphere(5.0);
    REQUIRE(shape.has_value());
    REQUIRE_FALSE(shape->IsNull());
}

TEST_CASE("MakeSphere shape type is Solid", "[geometry][sphere]")
{
    const auto shape = PCAD::Geometry::MakeSphere(5.0);
    REQUIRE(shape.has_value());
    CHECK(PCAD::Geometry::ShapeTypeString(*shape) == "Solid");
}

TEST_CASE("MakeSphere bounding box extents match diameter", "[geometry][sphere]")
{
    const double r = 3.0;
    const auto shape = PCAD::Geometry::MakeSphere(r);
    REQUIRE(shape.has_value());

    const auto bb = PCAD::Geometry::GetBoundingBox(*shape);
    REQUIRE(bb.has_value());

    const auto [dx, dy, dz] = bb->Extents();
    CHECK_THAT(dx, WithinAbs(2.0 * r, 1e-6));
    CHECK_THAT(dy, WithinAbs(2.0 * r, 1e-6));
    CHECK_THAT(dz, WithinAbs(2.0 * r, 1e-6));
}

// ── MakeBox ───────────────────────────────────────────────────────────────────

TEST_CASE("MakeBox returns a valid shape", "[geometry][box]")
{
    const auto shape = PCAD::Geometry::MakeBox(10.0, 6.0, 4.0);
    REQUIRE(shape.has_value());
    REQUIRE_FALSE(shape->IsNull());
}

TEST_CASE("MakeBox shape type is Solid", "[geometry][box]")
{
    const auto shape = PCAD::Geometry::MakeBox(10.0, 6.0, 4.0);
    REQUIRE(shape.has_value());
    CHECK(PCAD::Geometry::ShapeTypeString(*shape) == "Solid");
}

TEST_CASE("MakeBox bounding box extents match dimensions", "[geometry][box]")
{
    const auto shape = PCAD::Geometry::MakeBox(10.0, 6.0, 4.0);
    REQUIRE(shape.has_value());

    const auto bb = PCAD::Geometry::GetBoundingBox(*shape);
    REQUIRE(bb.has_value());

    const auto [dx, dy, dz] = bb->Extents();
    CHECK_THAT(dx, WithinAbs(10.0, 1e-6));
    CHECK_THAT(dy, WithinAbs(6.0,  1e-6));
    CHECK_THAT(dz, WithinAbs(4.0,  1e-6));
}

// ── GetBoundingBox edge cases ─────────────────────────────────────────────────

TEST_CASE("GetBoundingBox returns nullopt for null shape", "[geometry][bbox]")
{
    TopoDS_Shape nullShape;
    const auto bb = PCAD::Geometry::GetBoundingBox(nullShape);
    CHECK_FALSE(bb.has_value());
}
