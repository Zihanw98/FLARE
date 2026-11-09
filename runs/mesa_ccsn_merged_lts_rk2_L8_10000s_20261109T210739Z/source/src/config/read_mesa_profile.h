/**
 * @file read_mesa_profile.h
 * @brief Read MESA shell data and construct a radial FLARE mesh.
 */
#pragma once

#include "config.h"
#include "mesh.h"
#include "table.h"
#include <cstddef>
#include <string>
#include <vector>

namespace config
{
    /// @brief Active radial faces and MESA zone rows ordered from the innermost shell outward.
    struct MesaProfile
    {
        /// Active face radii in FLARE length units. N cells have N+1 faces.
        std::vector<double> radius_faces;
        /// Row i describes the shell between radius_faces[i] and radius_faces[i+1].
        /// Inner padding copies the innermost row. Exterior padding has no rows.
        /// Column names and values retain their original MESA units.
        Table zones;
    };

    /**
     * @brief Read a MESA profile and pad or truncate its radial grid.
     * @param[in] filename Profile path. The sixth line names the zone columns, including
     *            logR and dr. Rows run from the stellar surface toward the center.
     * @param[in] radius_inner Target inner face in FLARE length units. A negative value
     *            selects the innermost MESA face. A smaller nonnegative radius adds inner
     *            padding. A larger radius excises the inner shells, up to the stellar surface.
     *            A cut inside a shell keeps that shell's parameters for the surviving part.
     * @param[in] radius_outer Target outer face in the same units. A negative value selects
     *            the stellar surface. A cut inside a shell keeps that shell's original row.
     *            The outer face must exceed the target inner face and reach the native inner face.
     * @param[in] radius_scale Padding scale, either "linear" for equal radial widths or
     *            "log" for equal radius ratios. Logarithmic padding requires positive endpoints.
     * @param[in] inner_cells Number of inner padding cells. Zero estimates the count using
     *            the innermost shell spacing on the chosen scale, rounded to the nearest integer.
     * @param[in] outer_cells Number of outer padding cells. Zero uses the outermost shell spacing.
     * @param[in] solar_radius_cm Solar radius in cm used by the model to define logR.
     *            The default matches profile013_radius_fixed.data.
     * @pre Set physcon::units.unit_length to the number of centimeters per FLARE length unit.
     * @return Active faces and zone parameters. Use MesaProfile::radius_faces for the geometry
     *         of cut and padded cells. The logR and dr columns describe their source rows.
     * @throws std::runtime_error If the file or required columns cannot be read.
     * @throws std::invalid_argument If the units, radii, scale or resulting spacings are invalid.
     * @throws std::length_error If the grid exceeds the cell count supported by MeshBlock.
     */
    MesaProfile read_mesa_profile(const std::string& filename,
                                  double radius_inner = -1.0, double radius_outer = -1.0,
                                  const std::string& radius_scale = "linear",
                                  std::size_t inner_cells = 0, std::size_t outer_cells = 0,
                                  double solar_radius_cm = 6.96e10);

    /**
     * @brief Build a one-dimensional spherical mesh from the active faces of a MESA profile.
     * @param[in] profile Active radial faces in FLARE length units. State initialization uses
     *            the logRho, logP and velocity columns, one row per populated active cell.
     * @param[in,out] cfg Supplies the radial ghost count and boundary conditions, and gamma
     *            for state initialization. On success, receives the geometry, active cell
     *            count and radial limits derived from the profile.
     * @param[in] initialize_state If true, convert MESA density, total pressure and radial
     *            velocity to code units and construct each populated cell with cfg.gamma.
     *            If false, allocate default cells and skip the parameter table entirely.
     * @pre Build with HD_NDIM=1. Use the same unit system as when the profile was read.
     * @return A mesh with mirrored ghost faces at reflective boundaries and linearly
     *         extrapolated ghost faces elsewhere, shell volumes and face areas. Ghost and
     *         exterior padding cells retain their default state for the caller to fill.
     * @note Each populated cell takes the velocity in its own MESA row (the outer-face value).
     * @throws std::invalid_argument If the dimension, geometry data, units or gas state is invalid.
     * @throws std::runtime_error If a boundary condition or required state column is unavailable.
     * @throws std::overflow_error If the total mesh size exceeds the index range.
     */
    hd_core::MeshBlock build_mesh_from_mesa_profile(const MesaProfile& profile,
                                                   SimulationConfig& cfg,
                                                   bool initialize_state);
}
