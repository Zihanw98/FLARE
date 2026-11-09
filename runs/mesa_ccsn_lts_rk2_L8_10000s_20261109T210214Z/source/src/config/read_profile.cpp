/**
 * @file read_profile.cpp
 * @brief Read one row per active cell in FLARE storage order.
 */
#include "read_profile.h"
#include "table.h"
#include "core/mesh.h"
#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace config
{
void read_profile(const std::string& filename, hd_core::MeshBlock& mesh)
{
    // flat_index checks every axis, the total storage size and its representability.
    mesh.flat_index(hd_core::Index{}, true);
    std::size_t active_count = 1;
    for (const int n : mesh.nactive)
        active_count *= static_cast<std::size_t>(n);

    const Table table = read_table(filename, 0, true);
    if (table.data.size() != active_count)
        throw std::runtime_error(filename + ": expected " + std::to_string(active_count)
                                 + " active-cell rows, got " + std::to_string(table.data.size()));

    // Resolve column positions once. Duplicate required names would make the mapping ambiguous.
    const auto column = [&](const char* name) -> std::size_t {
        const auto end = table.headers.end();
        const auto it = std::find(table.headers.begin(), end, name);
        if (it == end)
            throw std::runtime_error(filename + ": missing column " + name);
        if (std::find(it + 1, end, name) != end)
            throw std::runtime_error(filename + ": duplicate column " + name);
        return static_cast<std::size_t>(it - table.headers.begin());
    };
    const auto rho = column("rho");
    const auto pre = column("pre");
    const auto gamma = column("gamma");
    const std::array<std::size_t, 3> velocity{column("vx"), column("vy"), column("vz")};

    std::vector<hd_core::Cell> states;
    states.reserve(active_count);
    for (const auto& row : table.data)
    {
        hd_core::Primitives state{};
        state[hd_core::RHO] = row[rho];
        state[hd_core::PRE] = row[pre];
        for (int d = 0; d < hd_core::NDIM; ++d)
            state[hd_core::VEL(d)] = row[velocity[d]];
        try
        {
            state.validate(row[gamma]); // Also checks gamma and zero-density (vacuum) rows.
            states.emplace_back(state, row[gamma]);
        }
        catch (const std::invalid_argument& e)
        {
            throw std::runtime_error(filename + ": data row " + std::to_string(states.size() + 1)
                                     + ": " + e.what());
        }
    }

    std::size_t row_index = 0;
    for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
    {
        if (mesh.is_active(flat))
            mesh.cells[flat] = states[row_index++];
    }
}
}
