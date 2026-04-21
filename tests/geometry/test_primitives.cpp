// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
/// @file tests/geometry/test_primitives.cpp
/// @brief Unit tests for BladeGen::Geometry primitives and bounding box queries.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "geometry/OcctUtils.hpp"

#include <TopoDS_Shape.hxx>

using namespace Catch::Matchers;

// ── MakeSphere ────────────────────────────────────────────────────────────────

TEST_CASE("MakeSphere returns a valid shape", "[geometry][sphere]")
{
    const auto shape = BladeGen::Geometry::MakeSphere(5.0);
    REQUIRE(shape.has_value());
    REQUIRE_FALSE(shape->IsNull());
}

TEST_CASE("MakeSphere shape type is Solid", "[geometry][sphere]")
{
    const auto shape = BladeGen::Geometry::MakeSphere(5.0);
    REQUIRE(shape.has_value());
    CHECK(BladeGen::Geometry::ShapeTypeString(*shape) == "Solid");
}

TEST_CASE("MakeSphere bounding box extents match diameter", "[geometry][sphere]")
{
    const double r = 3.0;
    const auto shape = BladeGen::Geometry::MakeSphere(r);
    REQUIRE(shape.has_value());

    const auto bb = BladeGen::Geometry::GetBoundingBox(*shape);
    REQUIRE(bb.has_value());

    const auto [dx, dy, dz] = bb->Extents();
    CHECK_THAT(dx, WithinAbs(2.0 * r, 1e-6));
    CHECK_THAT(dy, WithinAbs(2.0 * r, 1e-6));
    CHECK_THAT(dz, WithinAbs(2.0 * r, 1e-6));
}

// ── MakeBox ───────────────────────────────────────────────────────────────────

TEST_CASE("MakeBox returns a valid shape", "[geometry][box]")
{
    const auto shape = BladeGen::Geometry::MakeBox(10.0, 6.0, 4.0);
    REQUIRE(shape.has_value());
    REQUIRE_FALSE(shape->IsNull());
}

TEST_CASE("MakeBox shape type is Solid", "[geometry][box]")
{
    const auto shape = BladeGen::Geometry::MakeBox(10.0, 6.0, 4.0);
    REQUIRE(shape.has_value());
    CHECK(BladeGen::Geometry::ShapeTypeString(*shape) == "Solid");
}

TEST_CASE("MakeBox bounding box extents match dimensions", "[geometry][box]")
{
    const auto shape = BladeGen::Geometry::MakeBox(10.0, 6.0, 4.0);
    REQUIRE(shape.has_value());

    const auto bb = BladeGen::Geometry::GetBoundingBox(*shape);
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
    const auto bb = BladeGen::Geometry::GetBoundingBox(nullShape);
    CHECK_FALSE(bb.has_value());
}
