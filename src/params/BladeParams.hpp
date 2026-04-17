/// @file params/BladeParams.hpp
/// @brief JSON-serialisable blade design parameters.
///
/// This is the boundary between the file system and the geometry kernel.
/// All angular quantities are in degrees in the JSON file and in this struct;
/// to_profile_params() converts to radians for the blade geometry layer.
///
/// JSON schema (all lengths in metres, all angles in degrees):
/// @code
/// {
///   "r_hub": 0.25,
///   "r_tip": 0.40,
///   "axial_chord": 0.05,
///   "n_interp_sections": 3,
///   "stations": [
///     { "span": 0.0, "beta1_deg": 51.6, "beta2_deg": -63.0,
///       "theta_LE_deg": 0.0, "theta_TE_deg": 2.86,
///       "t_max": 0.015, "t_max_loc": 0.35,
///       "t_TE": 0.0, "wedge_TE_deg": 2.86, "r_LE": 0.002 },
///     { "span": 1.0, ... }
///   ]
/// }
/// @endcode
///
/// At least two stations are required (hub at span=0 and tip at span=1).
/// Additional intermediate stations are interpolated linearly.

#pragma once

#include "blade/BladeSection.hpp"

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

namespace PCAD::Params {

// ─── Per-span-station parameters ─────────────────────────────────────────────
struct StationParams
{
    double span;           ///< Normalised span position [0=hub, 1=tip]

    double beta1_deg;      ///< Inlet metal angle [deg]
    double beta2_deg;      ///< Exit  metal angle [deg]
    double theta_LE_deg;   ///< Leading-edge circumferential position [deg]
    double theta_TE_deg;   ///< Trailing-edge circumferential position [deg]

    double t_max;          ///< Maximum thickness [m]
    double t_max_loc;      ///< Chordwise location of max thickness [0–1]
    double t_TE;           ///< Trailing-edge thickness [m]
    double wedge_TE_deg;   ///< TE wedge half-angle [deg]
    double r_LE;           ///< Leading-edge radius [m]

    /// Convert to the radians-based internal struct.
    [[nodiscard]] Blade::TurbineProfileParams<double> to_profile_params() const;
};

// ─── Full blade parameter set ─────────────────────────────────────────────────
struct BladeParams
{
    double r_hub;            ///< Hub radius [m]
    double r_tip;            ///< Tip radius [m]
    double axial_chord;      ///< Axial chord (x-extent) [m]

    /// Number of evenly-spaced sections interpolated between the user-defined
    /// stations.  Total section count = n_interp_sections + stations.size().
    int n_interp_sections = 3;

    /// Spanwise stations (≥ 2, must include span=0 and span=1).
    std::vector<StationParams> stations;

    // ── Validation ───────────────────────────────────────────────────────────
    /// Throw std::invalid_argument with a descriptive message if anything is
    /// inconsistent.  Called automatically by from_file().
    void validate() const;

    // ── Geometry conversion ──────────────────────────────────────────────────
    /// Build a BladeGeometry ready for OCCT skinning.
    /// Performs piecewise-linear interpolation between stations.
    [[nodiscard]] Blade::BladeGeometry<double> to_blade_geometry() const;

    // ── I/O ──────────────────────────────────────────────────────────────────
    /// Parse from a JSON file.  Calls validate() before returning.
    /// @throws std::invalid_argument on schema violations.
    /// @throws std::runtime_error on file/parse errors.
    [[nodiscard]] static BladeParams from_file(const std::filesystem::path& path);

    /// Serialise to a JSON file (pretty-printed, 2-space indent).
    /// @throws std::runtime_error on write errors.
    void to_file(const std::filesystem::path& path) const;
};

} // namespace PCAD::Params
