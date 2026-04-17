/// @file tests/io/test_stl_export.cpp
/// @brief Integration tests for StlExporter.

#include <catch2/catch_test_macros.hpp>

#include "geometry/OcctUtils.hpp"
#include "io/StlExporter.hpp"

#include <filesystem>

namespace fs = std::filesystem;

static fs::path TempStlPath()
{
    return fs::temp_directory_path() / "pcad_test_export.stl";
}

TEST_CASE("StlExporter writes a non-empty STL file", "[io][stl]")
{
    const auto shape = PCAD::Geometry::MakeSphere(5.0);
    REQUIRE(shape.has_value());

    const auto path = TempStlPath();

    PCAD::IO::StlExporter exp;
    REQUIRE_NOTHROW(exp.Write(*shape, path));

    REQUIRE(fs::exists(path));
    CHECK(fs::file_size(path) > 0);

    fs::remove(path);
}

TEST_CASE("StlExporter finer deflection produces larger file", "[io][stl]")
{
    const auto shape = PCAD::Geometry::MakeSphere(10.0);
    REQUIRE(shape.has_value());

    const auto coarsePath = fs::temp_directory_path() / "pcad_coarse.stl";
    const auto finePath   = fs::temp_directory_path() / "pcad_fine.stl";

    PCAD::IO::StlExporter coarse;
    coarse.SetLinearDeflection(1.0);
    coarse.Write(*shape, coarsePath);

    PCAD::IO::StlExporter fine;
    fine.SetLinearDeflection(0.01);
    fine.Write(*shape, finePath);

    CHECK(fs::file_size(finePath) > fs::file_size(coarsePath));

    fs::remove(coarsePath);
    fs::remove(finePath);
}

TEST_CASE("StlExporter throws for null shape", "[io][stl]")
{
    TopoDS_Shape nullShape;
    PCAD::IO::StlExporter exp;
    CHECK_THROWS_AS(exp.Write(nullShape, TempStlPath()), std::runtime_error);
}
