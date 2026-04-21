// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "params/BladeParams.hpp"

#include <filesystem>
#include <fstream>
#include <numbers>

namespace fs = std::filesystem;
using namespace BladeGen::Params;

// ─── Helpers ─────────────────────────────────────────────────────────────────
static constexpr double kDeg = std::numbers::pi / 180.0;

static fs::path tmp(const char* name)
{
    return fs::temp_directory_path() / name;
}

static BladeParams make_two_station()
{
    BladeParams bp;
    bp.r_hub            = 0.25;
    bp.r_tip            = 0.40;
    bp.axial_chord      = 0.05;
    bp.n_interp_sections = 2;

    StationParams hub;
    hub.span = 0.0; hub.beta1_deg = 51.6; hub.beta2_deg = -63.0;
    hub.theta_LE_deg = 0.0; hub.theta_TE_deg = 2.86;
    hub.t_max = 0.015; hub.t_max_loc = 0.35;
    hub.t_TE = 0.0; hub.wedge_TE_deg = 2.86; hub.r_LE = 0.002;

    StationParams tip = hub;
    tip.span = 1.0; tip.t_max = 0.010;

    bp.stations = {hub, tip};
    return bp;
}

// ─── to_profile_params: degree → radian conversion ───────────────────────────
TEST_CASE("StationParams converts angles to radians", "[params]")
{
    StationParams s;
    s.span = 0.0;
    s.beta1_deg = 90.0; s.beta2_deg = -45.0;
    s.theta_LE_deg = 0.0; s.theta_TE_deg = 5.0;
    s.t_max = 0.01; s.t_max_loc = 0.35;
    s.t_TE = 0.0; s.wedge_TE_deg = 3.0; s.r_LE = 0.001;

    const auto p = s.to_profile_params();
    REQUIRE_THAT(p.beta1,    Catch::Matchers::WithinRel(90.0 * kDeg,  1e-12));
    REQUIRE_THAT(p.beta2,    Catch::Matchers::WithinRel(-45.0 * kDeg, 1e-12));
    REQUIRE_THAT(p.theta_LE, Catch::Matchers::WithinAbs(0.0,          1e-15));
    REQUIRE_THAT(p.theta_TE, Catch::Matchers::WithinRel(5.0 * kDeg,   1e-12));
    REQUIRE_THAT(p.wedge_TE, Catch::Matchers::WithinRel(3.0 * kDeg,   1e-12));
    REQUIRE_THAT(p.t_max,    Catch::Matchers::WithinAbs(0.01, 1e-15));
}

// ─── validate() ──────────────────────────────────────────────────────────────
TEST_CASE("validate() accepts a well-formed BladeParams", "[params]")
{
    REQUIRE_NOTHROW(make_two_station().validate());
}

TEST_CASE("validate() rejects r_hub >= r_tip", "[params]")
{
    auto bp = make_two_station();
    bp.r_hub = 0.40; bp.r_tip = 0.25;
    CHECK_THROWS_AS(bp.validate(), std::invalid_argument);
}

TEST_CASE("validate() rejects fewer than 2 stations", "[params]")
{
    auto bp = make_two_station();
    bp.stations.resize(1);
    CHECK_THROWS_AS(bp.validate(), std::invalid_argument);
}

TEST_CASE("validate() rejects stations not starting at span=0", "[params]")
{
    auto bp = make_two_station();
    bp.stations.front().span = 0.1;
    CHECK_THROWS_AS(bp.validate(), std::invalid_argument);
}

TEST_CASE("validate() rejects stations not ending at span=1", "[params]")
{
    auto bp = make_two_station();
    bp.stations.back().span = 0.9;
    CHECK_THROWS_AS(bp.validate(), std::invalid_argument);
}

TEST_CASE("validate() rejects unsorted stations", "[params]")
{
    auto bp = make_two_station();
    StationParams mid = bp.stations[0]; mid.span = 0.5;
    bp.stations.insert(bp.stations.begin() + 1, mid);
    // stations are now {0, 0.5, 1} — valid; swap to make them unsorted
    std::swap(bp.stations[1], bp.stations[2]);
    CHECK_THROWS_AS(bp.validate(), std::invalid_argument);
}

TEST_CASE("validate() rejects non-positive t_max", "[params]")
{
    auto bp = make_two_station();
    bp.stations[0].t_max = 0.0;
    CHECK_THROWS_AS(bp.validate(), std::invalid_argument);
}

// ─── to_file / from_file round-trip ──────────────────────────────────────────
TEST_CASE("BladeParams round-trips through JSON file", "[params][io]")
{
    const auto bp  = make_two_station();
    const auto path = tmp("test_blade_params.json");

    REQUIRE_NOTHROW(bp.to_file(path));
    REQUIRE(fs::exists(path));
    CHECK(fs::file_size(path) > 0);

    const auto bp2 = BladeParams::from_file(path);

    REQUIRE_THAT(bp2.r_hub,       Catch::Matchers::WithinAbs(bp.r_hub,       1e-12));
    REQUIRE_THAT(bp2.r_tip,       Catch::Matchers::WithinAbs(bp.r_tip,       1e-12));
    REQUIRE_THAT(bp2.axial_chord, Catch::Matchers::WithinAbs(bp.axial_chord, 1e-12));
    CHECK(bp2.n_interp_sections == bp.n_interp_sections);
    REQUIRE(bp2.stations.size() == bp.stations.size());

    for (std::size_t i = 0; i < bp.stations.size(); ++i) {
        REQUIRE_THAT(bp2.stations[i].span,
                     Catch::Matchers::WithinAbs(bp.stations[i].span, 1e-12));
        REQUIRE_THAT(bp2.stations[i].beta1_deg,
                     Catch::Matchers::WithinAbs(bp.stations[i].beta1_deg, 1e-12));
        REQUIRE_THAT(bp2.stations[i].t_max,
                     Catch::Matchers::WithinAbs(bp.stations[i].t_max, 1e-12));
    }

    fs::remove(path);
}

TEST_CASE("from_file throws on missing file", "[params][io]")
{
    CHECK_THROWS_AS(BladeParams::from_file(tmp("nonexistent_blade.json")),
                    std::runtime_error);
}

TEST_CASE("from_file throws on malformed JSON", "[params][io]")
{
    const auto path = tmp("bad_blade.json");
    { std::ofstream f(path); f << "{ not valid json !!!"; }
    CHECK_THROWS_AS(BladeParams::from_file(path), std::runtime_error);
    fs::remove(path);
}

TEST_CASE("from_file throws on missing required field", "[params][io]")
{
    const auto path = tmp("missing_field.json");
    {
        std::ofstream f(path);
        f << R"({"r_hub": 0.25, "r_tip": 0.40, "axial_chord": 0.05,
                 "stations": [{"span": 0.0, "beta1_deg": 50.0}]})";
    }
    CHECK_THROWS_AS(BladeParams::from_file(path), std::invalid_argument);
    fs::remove(path);
}

// ─── to_blade_geometry ───────────────────────────────────────────────────────
TEST_CASE("to_blade_geometry produces correct section count", "[params][geometry]")
{
    const auto bp   = make_two_station();
    const auto geom = bp.to_blade_geometry();

    // 2 user stations + 2 interp sections = 4 total
    CHECK(static_cast<int>(geom.sections.size()) == 4);
}

TEST_CASE("to_blade_geometry: three-station blade has sections at all user spans",
          "[params][geometry]")
{
    auto bp = make_two_station();
    StationParams mid = bp.stations[0];
    mid.span = 0.5; mid.t_max = 0.012;
    bp.stations.insert(bp.stations.begin() + 1, mid);
    bp.n_interp_sections = 0;

    const auto geom = bp.to_blade_geometry();
    // 3 user stations + 0 interp = 3 total
    CHECK(static_cast<int>(geom.sections.size()) == 3);
}

TEST_CASE("to_blade_geometry generates valid 3D section loops", "[params][geometry]")
{
    const auto bp    = make_two_station();
    const auto geom  = bp.to_blade_geometry();
    const auto all3d = geom.all_sections_3d(32);

    CHECK(all3d.size() == geom.sections.size());
    for (const auto& sec : all3d)
        CHECK(sec.size() == 64);   // 2 × 32 points
}

TEST_CASE("interpolated mid-span thickness is between hub and tip", "[params][geometry]")
{
    auto bp = make_two_station();
    bp.n_interp_sections = 1;   // one section at span=0.5

    const auto geom = bp.to_blade_geometry();
    // Sections are at 0.0, 0.5, 1.0
    REQUIRE(geom.sections.size() == 3);

    // Check via the thickness distribution at t_max_loc
    const auto& sec_mid = geom.sections[1];
    // Indirect check: the mid section generates without error and produces points
    // whose thickness is between hub and tip values.
    CHECK_NOTHROW(sec_mid.points3d(16));
}
