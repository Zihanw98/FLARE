/**
 * @file output.h
 * @brief Simulation diagnostics and output.
 */
#pragma once
#include <cstddef>

namespace hd_core { struct MeshBlock; }
namespace config { struct SimulationConfig; }
namespace hd_user { struct UserHooks; }

namespace io
{
    /**
     * @brief Print step diagnostics and write scheduled CSV snapshots.
     * @param[in,out] mesh Allocated mesh; built-in output reads it, while the optional callback receives a mutable reference.
     * @param[in] cfg Output directory, suffix, step interval, end time, and maximum step count.
     * @param[in] hooks Optional useroutput callback, called after a snapshot is written successfully.
     * @param[in] step Number of completed steps; zero denotes the initial state.
     * @param[in] time Finite nonnegative physical time.
     * @param[in] dt Finite nonnegative timestep just taken; zero for the initial state.
     * @pre Geometry and all states being output are valid; ghost states have already been filled.
     * @note Always prints active-cell density/pressure ranges and integrated mass/energy. Writes snapshots at steps 0 and 1, interval multiples, and the end time or step limit.
     * @note CSV includes ghost cells in storage order, cell-center coordinates, an is_ghost flag, and three initial metadata lines for step/time/dt. Missing axes and velocities are zero.
     * @throws std::invalid_argument If output settings, time, geometry, or a state are invalid.
     * @throws std::runtime_error If the CSV file cannot be opened or written.
     * @throws std::filesystem::filesystem_error If the output directory cannot be created.
     */
    void Output(
        hd_core::MeshBlock& mesh,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        std::size_t step,
        double time,
        double dt
    );
}
