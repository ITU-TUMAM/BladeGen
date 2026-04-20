/// @file params/BladeParams.cpp
/// @brief JSON serialisation and geometry conversion for BladeParams.

#include "params/BladeParams.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace PCAD::Params {

// ─── Degree ↔ radian helpers ─────────────────────────────────────────────────
static constexpr double kDegToRad = std::numbers::pi / 180.0;
static constexpr double kRadToDeg = 180.0 / std::numbers::pi;

// ─── StationParams ────────────────────────────────────────────────────────────
Blade::TurbineProfileParams<double> StationParams::to_profile_params() const
{
    Blade::TurbineProfileParams<double> p;
    p.beta1     = beta1_deg    * kDegToRad;
    p.beta2     = beta2_deg    * kDegToRad;
    p.theta_LE  = theta_LE_deg * kDegToRad;
    p.theta_TE  = theta_TE_deg * kDegToRad;
    p.t_max     = t_max;
    p.t_max_loc = t_max_loc;
    p.t_TE      = t_TE;
    p.wedge_TE  = wedge_TE_deg * kDegToRad;
    p.r_LE      = r_LE;
    return p;
}

// ─── nlohmann/json from/to hooks ─────────────────────────────────────────────
static void from_json(const nlohmann::json& j, StationParams& s)
{
    s.span        = j.at("span")        .get<double>();
    s.beta1_deg   = j.at("beta1_deg")   .get<double>();
    s.beta2_deg   = j.at("beta2_deg")   .get<double>();
    s.theta_LE_deg = j.at("theta_LE_deg").get<double>();
    s.theta_TE_deg = j.at("theta_TE_deg").get<double>();
    s.t_max       = j.at("t_max")       .get<double>();
    s.t_max_loc   = j.at("t_max_loc")   .get<double>();
    s.t_TE        = j.at("t_TE")        .get<double>();
    s.wedge_TE_deg = j.at("wedge_TE_deg").get<double>();
    s.r_LE        = j.at("r_LE")        .get<double>();
}

static void to_json(nlohmann::json& j, const StationParams& s)
{
    j = nlohmann::json{
        {"span",         s.span},
        {"beta1_deg",    s.beta1_deg},
        {"beta2_deg",    s.beta2_deg},
        {"theta_LE_deg", s.theta_LE_deg},
        {"theta_TE_deg", s.theta_TE_deg},
        {"t_max",        s.t_max},
        {"t_max_loc",    s.t_max_loc},
        {"t_TE",         s.t_TE},
        {"wedge_TE_deg", s.wedge_TE_deg},
        {"r_LE",         s.r_LE}
    };
}

static void from_json(const nlohmann::json& j, BladeParams& b)
{
    b.r_hub             = j.at("r_hub")            .get<double>();
    b.r_tip             = j.at("r_tip")            .get<double>();
    b.axial_chord       = j.at("axial_chord")      .get<double>();
    b.n_interp_sections = j.value("n_interp_sections", 3);
    b.stations          = j.at("stations")         .get<std::vector<StationParams>>();
}

static void to_json(nlohmann::json& j, const BladeParams& b)
{
    j = nlohmann::json{
        {"r_hub",             b.r_hub},
        {"r_tip",             b.r_tip},
        {"axial_chord",       b.axial_chord},
        {"n_interp_sections", b.n_interp_sections},
        {"stations",          b.stations}
    };
}

// ─── Validation ───────────────────────────────────────────────────────────────
void BladeParams::validate() const
{
    auto err = [](const std::string& msg) {
        throw std::invalid_argument("BladeParams: " + msg);
    };

    if (r_hub <= 0.0)             err("r_hub must be positive");
    if (r_tip <= r_hub)           err("r_tip must be greater than r_hub");
    if (axial_chord <= 0.0)       err("axial_chord must be positive");
    if (n_interp_sections < 0)    err("n_interp_sections must be >= 0");
    if (stations.size() < 2)      err("at least 2 stations required");

    // Stations must be sorted and span [0, 1]
    for (std::size_t i = 0; i < stations.size(); ++i) {
        const double s = stations[i].span;
        if (s < 0.0 || s > 1.0) {
            std::ostringstream oss;
            oss << "station[" << i << "].span=" << s << " is out of [0,1]";
            err(oss.str());
        }
        if (i > 0 && stations[i].span <= stations[i-1].span) {
            std::ostringstream oss;
            oss << "stations must be strictly ordered by span (violation at index " << i << ")";
            err(oss.str());
        }
    }
    if (std::abs(stations.front().span - 0.0) > 1e-9)
        err("first station must have span=0 (hub)");
    if (std::abs(stations.back().span  - 1.0) > 1e-9)
        err("last station must have span=1 (tip)");

    for (std::size_t i = 0; i < stations.size(); ++i) {
        const auto& st = stations[i];
        if (st.t_max <= 0.0) {
            std::ostringstream oss;
            oss << "station[" << i << "].t_max must be positive";
            err(oss.str());
        }
        if (st.t_max_loc <= 0.0 || st.t_max_loc >= 1.0) {
            std::ostringstream oss;
            oss << "station[" << i << "].t_max_loc must be in (0,1)";
            err(oss.str());
        }
    }
}

// ─── Piecewise-linear station interpolation ───────────────────────────────────
// Assumes stations is sorted by span (validated before this is called).
static StationParams interpolate(double span,
                                 const std::vector<StationParams>& sts)
{
    // Find bounding interval
    std::size_t hi = 1;
    while (hi < sts.size() - 1 && sts[hi].span < span)
        ++hi;
    const std::size_t lo = hi - 1;

    const double t = (span - sts[lo].span) / (sts[hi].span - sts[lo].span);

    auto lerp = [t](double a, double b) { return a + t * (b - a); };

    StationParams s;
    s.span         = span;
    s.beta1_deg    = lerp(sts[lo].beta1_deg,    sts[hi].beta1_deg);
    s.beta2_deg    = lerp(sts[lo].beta2_deg,    sts[hi].beta2_deg);
    s.theta_LE_deg = lerp(sts[lo].theta_LE_deg, sts[hi].theta_LE_deg);
    s.theta_TE_deg = lerp(sts[lo].theta_TE_deg, sts[hi].theta_TE_deg);
    s.t_max        = lerp(sts[lo].t_max,        sts[hi].t_max);
    s.t_max_loc    = lerp(sts[lo].t_max_loc,    sts[hi].t_max_loc);
    s.t_TE         = lerp(sts[lo].t_TE,         sts[hi].t_TE);
    s.wedge_TE_deg = lerp(sts[lo].wedge_TE_deg, sts[hi].wedge_TE_deg);
    s.r_LE         = lerp(sts[lo].r_LE,         sts[hi].r_LE);
    return s;
}

// ─── Geometry conversion ─────────────────────────────────────────────────────
Blade::BladeGeometry<double> BladeParams::to_blade_geometry() const
{
    // Collect all span fractions: user-defined + evenly-spaced interpolated
    std::vector<double> spans;
    spans.reserve(stations.size() + n_interp_sections);

    for (const auto& st : stations)
        spans.push_back(st.span);

    for (int i = 1; i <= n_interp_sections; ++i)
        spans.push_back(static_cast<double>(i) / (n_interp_sections + 1));

    // Sort and remove near-duplicates (within 1e-6 span)
    std::ranges::sort(spans);
    spans.erase(
        std::unique(spans.begin(), spans.end(),
                    [](double a, double b) { return b - a < 1e-6; }),
        spans.end());

    Blade::BladeGeometry<double> geom;
    geom.sections.reserve(spans.size());

    for (double s : spans) {
        const double r = r_hub + s * (r_tip - r_hub);

        // Straight axial streamline at constant radius
        Blade::BSplineCurve2D<double> sl;
        sl.degree = 1;
        sl.knots  = Blade::clamped_knots(1, 2);
        sl.Px     = {0.0, axial_chord};
        sl.Pr     = {r,   r};

        const StationParams st = interpolate(s, stations);

        Blade::BladeSpanSection<double> sec;
        sec.streamline = std::move(sl);
        sec.profile    = Blade::ProfileSection<double>::from_params(
                             st.to_profile_params());

        geom.sections.push_back(std::move(sec));
    }

    return geom;
}

// ─── I/O ─────────────────────────────────────────────────────────────────────
BladeParams BladeParams::from_file(const std::filesystem::path& path)
{
    std::ifstream f(path);
    if (!f.is_open())
        throw std::runtime_error("BladeParams::from_file: cannot open '" +
                                 path.string() + "'");

    nlohmann::json j;
    try {
        f >> j;
    } catch (const nlohmann::json::parse_error& e) {
        throw std::runtime_error("BladeParams::from_file: JSON parse error in '" +
                                 path.string() + "': " + e.what());
    }

    BladeParams bp;
    try {
        from_json(j, bp);
    } catch (const nlohmann::json::exception& e) {
        throw std::invalid_argument(
            std::string("BladeParams::from_file: missing or wrong-type field: ") +
            e.what());
    }

    bp.validate();
    return bp;
}

void BladeParams::to_file(const std::filesystem::path& path) const
{
    nlohmann::json j;
    to_json(j, *this);

    std::ofstream f(path);
    if (!f.is_open())
        throw std::runtime_error("BladeParams::to_file: cannot write '" +
                                 path.string() + "'");
    f << j.dump(2) << '\n';
}

} // namespace PCAD::Params
