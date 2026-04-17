/// @file tests/io/test_step_export.cpp
/// @brief Integration tests for StepExporter — writes real STEP files to a
///        temporary directory and checks they are non-empty.

#include <catch2/catch_test_macros.hpp>

#include "geometry/OcctUtils.hpp"
#include "io/StepExporter.hpp"

#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

// ── Helpers ───────────────────────────────────────────────────────────────────

static fs::path TempStepPath()
{
    return fs::temp_directory_path() / "pcad_test_export.step";
}

// ── Tests ─────────────────────────────────────────────────────────────────────

TEST_CASE("StepExporter writes a non-empty STEP file", "[io][step]")
{
    const auto shape = PCAD::Geometry::MakeSphere(5.0);
    REQUIRE(shape.has_value());

    const auto path = TempStepPath();

    PCAD::IO::StepExporter exp;
    exp.AddShape(*shape, "TestSphere");
    REQUIRE_NOTHROW(exp.Write(path));

    REQUIRE(fs::exists(path));
    CHECK(fs::file_size(path) > 0);

    fs::remove(path);
}

TEST_CASE("StepExporter writes multiple shapes into one file", "[io][step]")
{
    const auto sphere = PCAD::Geometry::MakeSphere(3.0);
    const auto box    = PCAD::Geometry::MakeBox(4.0, 4.0, 4.0);
    REQUIRE(sphere.has_value());
    REQUIRE(box.has_value());

    const auto path = TempStepPath();

    PCAD::IO::StepExporter exp;
    exp.AddShape(*sphere, "Sphere");
    exp.AddShape(*box,    "Box");
    REQUIRE_NOTHROW(exp.Write(path));

    REQUIRE(fs::exists(path));
    CHECK(fs::file_size(path) > 0);

    fs::remove(path);
}

TEST_CASE("StepExporter throws when no shapes staged", "[io][step]")
{
    PCAD::IO::StepExporter exp;
    CHECK_THROWS_AS(exp.Write(TempStepPath()), std::runtime_error);
}

TEST_CASE("StepExporter AP203 protocol produces a file", "[io][step]")
{
    const auto shape = PCAD::Geometry::MakeBox(2.0, 2.0, 2.0);
    REQUIRE(shape.has_value());

    const auto path = TempStepPath();

    PCAD::IO::StepExporter exp(PCAD::IO::StepProtocol::AP203);
    exp.AddShape(*shape, "AP203Box");
    REQUIRE_NOTHROW(exp.Write(path));

    REQUIRE(fs::exists(path));
    CHECK(fs::file_size(path) > 0);

    fs::remove(path);
}
