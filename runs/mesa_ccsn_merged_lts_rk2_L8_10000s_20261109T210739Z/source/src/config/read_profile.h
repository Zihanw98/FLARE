/**
 * @file read_profile.h
 * @brief Load primitive hydrodynamic states onto an existing FLARE mesh.
 */
#pragma once

#include <string>

namespace hd_core { struct MeshBlock; }

namespace config
{
    /**
     * @brief Initialize active cells from a whitespace-delimited profile in code units.
     * @param[in] filename Profile path. The first line names rho, vx, vy, vz, pre
     *            and gamma in any order. Only the velocity components present in
     *            the compiled mesh dimension enter the hydrodynamic state.
     * @param[in,out] mesh Allocated mesh. Each data row supplies one active cell,
     *                with x varying fastest, then y, then z. Ghost cells keep their states.
     * @throws std::runtime_error If the table cannot be read, its columns or row count
     *         are invalid, or a row cannot be converted to a physical cell state.
     * @throws std::logic_error If the mesh storage does not match its dimensions.
     * @note Every row is validated before the first cell is updated.
     */
    void read_profile(const std::string& filename, hd_core::MeshBlock& mesh);
}
