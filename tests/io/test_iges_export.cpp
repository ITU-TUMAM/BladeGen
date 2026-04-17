/// @file tests/io/test_iges_export.cpp
/// @brief Integration tests for IgesExporter — writes real IGES files to a
///        temporary directory and checks they are non-empty and well-formed.

#include <catch2/catch_test_macros.hpp>

#include "blade/BladeSection.hpp"
#include "blade/OcctBlade.hpp"
#include "geometry/OcctUtils.hpp"
#include "io/IgesExporter.hpp"

#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

static fs::path TempIgesPath()
{
    return fs::temp_directory_path() / "pcad_test_export.iges";
}

// ─── Basic OCCT primitive ─────────────────────────────────────────────────────

TEST_CASE("IgesExporter writes a non-empty IGES file", "[io][iges]")
{
    const auto shape = PCAD::Geometry::MakeSphere(5.0);
    REQUIRE(shape.has_value());

    const auto path = TempIgesPath();

    PCAD::IO::IgesExporter exp;
    exp.AddShape(*shape, "TestSphere");
    REQUIRE_NOTHROW(exp.Write(path));

    REQUIRE(fs::exists(path));
    CHECK(fs::file_size(path) > 0);

    fs::remove(path);
}

TEST_CASE("IgesExporter writes multiple shapes into one file", "[io][iges]")
{
    const auto sphere = PCAD::Geometry::MakeSphere(3.0);
    const auto box    = PCAD::Geometry::MakeBox(4.0, 4.0, 4.0);
    REQUIRE(sphere.has_value());
    REQUIRE(box.has_value());

    const auto path = TempIgesPath();

    PCAD::IO::IgesExporter exp;
    exp.AddShape(*sphere, "Sphere");
    exp.AddShape(*box,    "Box");
    REQUIRE_NOTHROW(exp.Write(path));

    REQUIRE(fs::exists(path));
    CHECK(fs::file_size(path) > 0);

    fs::remove(path);
}

TEST_CASE("IgesExporter throws when no shapes staged", "[io][iges]")
{
    PCAD::IO::IgesExporter exp;
    CHECK_THROWS_AS(exp.Write(TempIgesPath()), std::runtime_error);
}

TEST_CASE("IGES file starts with correct magic header", "[io][iges]")
{
    const auto shape = PCAD::Geometry::MakeBox(1.0, 1.0, 1.0);
    REQUIRE(shape.has_value());

    const auto path = TempIgesPath();

    PCAD::IO::IgesExporter exp;
    exp.AddShape(*shape);
    exp.Write(path);

    // IGES files begin with a "S" (Start) section flag in column 73
    std::string line;
    {
        std::ifstream f(path);
        REQUIRE(f.is_open());
        std::getline(f, line);
    }  // f closed before remove

    // Column 73 (0-indexed: 72) holds the section code; first line is always 'S'
    REQUIRE(line.size() >= 73);
    CHECK(line[72] == 'S');

    fs::remove(path);
}

// ─── Full blade → IGES pipeline ───────────────────────────────────────────────

TEST_CASE("IgesExporter exports a generated blade solid", "[io][iges][blade]")
{
    using namespace PCAD::Blade;

    TurbineProfileParams<double> p;
    p.beta1     =  0.9;
    p.beta2     = -1.1;
    p.theta_LE  =  0.0;
    p.theta_TE  =  0.05;
    p.t_max     =  0.015;
    p.t_max_loc =  0.35;
    p.t_TE      =  0.0;
    p.wedge_TE  =  0.05;
    p.r_LE      =  0.002;

    const auto geom  = make_annular_blade(0.25, 0.40, 0.05, p, p, 2);
    const auto sec3d = geom.all_sections_3d(48);

    const auto blade = MakeBladeSolid(sec3d);
    REQUIRE(blade.has_value());

    const auto path = TempIgesPath();

    PCAD::IO::IgesExporter exp;
    exp.AddShape(*blade, "TurbineBlade");
    REQUIRE_NOTHROW(exp.Write(path));

    REQUIRE(fs::exists(path));
    CHECK(fs::file_size(path) > 0);

    fs::remove(path);
}
