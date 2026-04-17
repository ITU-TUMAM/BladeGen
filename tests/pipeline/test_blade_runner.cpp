#include <catch2/catch_test_macros.hpp>

#include "params/BladeParams.hpp"
#include "pipeline/BladeRunner.hpp"

#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;
using namespace PCAD;

// ─── Helpers ─────────────────────────────────────────────────────────────────
static Params::BladeParams typical_params()
{
    Params::BladeParams bp;
    bp.r_hub = 0.25; bp.r_tip = 0.40; bp.axial_chord = 0.05;
    bp.n_interp_sections = 2;

    Params::StationParams hub;
    hub.span = 0.0; hub.beta1_deg = 51.6; hub.beta2_deg = -63.0;
    hub.theta_LE_deg = 0.0; hub.theta_TE_deg = 2.86;
    hub.t_max = 0.015; hub.t_max_loc = 0.35;
    hub.t_TE = 0.0; hub.wedge_TE_deg = 2.86; hub.r_LE = 0.002;

    Params::StationParams tip = hub;
    tip.span = 1.0; tip.t_max = 0.010;

    bp.stations = {hub, tip};
    return bp;
}

static fs::path tmp_dir()
{
    auto p = fs::temp_directory_path() / "pcad_runner_test";
    fs::create_directories(p);
    return p;
}

// ─── IGES (default format) ────────────────────────────────────────────────────
TEST_CASE("BladeRunner writes an IGES file by default", "[pipeline]")
{
    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir = tmp_dir();
    cfg.stem       = "runner_default";

    const auto result = Pipeline::BladeRunner(cfg).run(typical_params());

    REQUIRE(result.written.size() == 1);
    CHECK(result.written[0].extension() == ".iges");
    CHECK(fs::exists(result.written[0]));
    CHECK(fs::file_size(result.written[0]) > 0);

    for (const auto& p : result.written) fs::remove(p);
}

// ─── Result metadata ─────────────────────────────────────────────────────────
TEST_CASE("BladeRunner result reports correct section and point counts", "[pipeline]")
{
    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir    = tmp_dir();
    cfg.stem          = "runner_meta";
    cfg.profile_points = 32;

    const auto bp     = typical_params();   // 2 stations + 2 interp = 4 sections
    const auto result = Pipeline::BladeRunner(cfg).run(bp);

    CHECK(result.n_sections      == 4);
    CHECK(result.pts_per_section == 64);  // 2 × profile_points

    for (const auto& p : result.written) fs::remove(p);
}

// ─── Multiple formats ─────────────────────────────────────────────────────────
TEST_CASE("BladeRunner writes all enabled formats", "[pipeline]")
{
    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir  = tmp_dir();
    cfg.stem        = "runner_multi";
    cfg.export_iges = true;
    cfg.export_step = true;
    cfg.export_stl  = true;

    const auto result = Pipeline::BladeRunner(cfg).run(typical_params());

    REQUIRE(result.written.size() == 3);

    std::set<std::string> exts;
    for (const auto& p : result.written) {
        CHECK(fs::exists(p));
        CHECK(fs::file_size(p) > 0);
        exts.insert(p.extension().string());
    }
    CHECK(exts.count(".iges") == 1);
    CHECK(exts.count(".step") == 1);
    CHECK(exts.count(".stl")  == 1);

    for (const auto& p : result.written) fs::remove(p);
}

TEST_CASE("BladeRunner writes only step when iges disabled", "[pipeline]")
{
    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir  = tmp_dir();
    cfg.stem        = "runner_step_only";
    cfg.export_iges = false;
    cfg.export_step = true;

    const auto result = Pipeline::BladeRunner(cfg).run(typical_params());

    REQUIRE(result.written.size() == 1);
    CHECK(result.written[0].extension() == ".step");

    for (const auto& p : result.written) fs::remove(p);
}

// ─── Output directory creation ────────────────────────────────────────────────
TEST_CASE("BladeRunner creates output_dir if it does not exist", "[pipeline]")
{
    const auto new_dir = tmp_dir() / "new_subdir_runner";
    fs::remove_all(new_dir);
    REQUIRE_FALSE(fs::exists(new_dir));

    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir = new_dir;
    cfg.stem       = "blade";

    const auto result = Pipeline::BladeRunner(cfg).run(typical_params());

    CHECK(fs::exists(new_dir));
    for (const auto& p : result.written) fs::remove(p);
    fs::remove(new_dir);
}

// ─── File-based convenience overload ─────────────────────────────────────────
TEST_CASE("BladeRunner::run(path) loads params and produces output", "[pipeline]")
{
    const auto params_path = tmp_dir() / "runner_test_params.json";
    typical_params().to_file(params_path);

    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir = tmp_dir();
    cfg.stem       = "runner_from_file";

    const auto result = Pipeline::BladeRunner(cfg).run(params_path);

    REQUIRE_FALSE(result.written.empty());
    for (const auto& p : result.written) {
        CHECK(fs::exists(p));
        fs::remove(p);
    }
    fs::remove(params_path);
}

// ─── Error paths ─────────────────────────────────────────────────────────────
TEST_CASE("BladeRunner throws when no format is enabled", "[pipeline]")
{
    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir  = tmp_dir();
    cfg.export_iges = false;
    cfg.export_step = false;
    cfg.export_stl  = false;

    CHECK_THROWS_AS(Pipeline::BladeRunner(cfg).run(typical_params()),
                    std::runtime_error);
}

TEST_CASE("BladeRunner throws on invalid params", "[pipeline]")
{
    auto bp = typical_params();
    bp.r_hub = 0.5; bp.r_tip = 0.1;  // invalid: hub > tip

    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir = tmp_dir();

    CHECK_THROWS_AS(Pipeline::BladeRunner(cfg).run(bp),
                    std::invalid_argument);
}

TEST_CASE("BladeRunner::run(path) throws on missing file", "[pipeline]")
{
    Pipeline::BladeRunnerConfig cfg;
    cfg.output_dir = tmp_dir();

    CHECK_THROWS_AS(
        Pipeline::BladeRunner(cfg).run(tmp_dir() / "nonexistent.json"),
        std::runtime_error);
}
