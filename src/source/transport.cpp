#include "source.h"
#include "core/mesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hd_source
{
void add_transport_flux(
    const hd_core::MeshBlock& stage, std::size_t face, double time,
    const config::SimulationConfig& cfg, const SourceList& sources,
    hd_core::Conserved& flux, std::array<double, 2>& inverse_dt)
{
    inverse_dt = {};
    for (const auto& source : sources)
    {
        if (source.kind != SourceKind::Transport) continue;
        if (!source.transport_flux)
            throw std::invalid_argument("Missing transport face function: " + source.name);
        hd_core::Conserved contribution{};
        std::array<double, 2> rates{};
        source.transport_flux(stage, face, time, cfg, contribution, rates);
        for (double value : contribution.data)
            if (!std::isfinite(value))
                throw std::invalid_argument("Transport returned a non-finite face flux: " + source.name);
        for (double value : rates)
            if (!std::isfinite(value) || value < 0.0)
                throw std::invalid_argument("Transport requires finite nonnegative face rates: " + source.name);
        flux += contribution;
        for (int side = 0; side < 2; ++side) inverse_dt[side] += rates[side];
    }
    for (double value : inverse_dt)
        if (!std::isfinite(value)) throw std::invalid_argument("Non-finite accumulated face rate");
}

double add_transport_rhs(
    const hd_core::MeshBlock& stage, double time,
    const config::SimulationConfig& cfg, const SourceList& sources,
    std::vector<hd_core::Conserved>& rhs)
{
    if (std::none_of(sources.begin(), sources.end(), [](const SourceTerm& source) {
            return source.kind == SourceKind::Transport;
        })) return 0.0;
    if (rhs.size() != stage.cells.size())
        throw std::invalid_argument("Transport RHS size must match mesh storage");
    std::vector<double> rates(stage.active_count, 0.0);
    for (std::size_t face = 0; face < stage.face_info.size(); ++face)
    {
        const auto& info = stage.face_info[face];
        hd_core::Conserved flux{};
        std::array<double, 2> rate{};
        add_transport_flux(stage, face, time, cfg, sources, flux, rate);
        const auto balance = flux * stage.face_area[info.direction].at(info.storage_coordinates);
        if (info.left_active != hd_core::MeshBlock::invalid_active_index)
        {
            const auto flat = stage.active_to_storage[info.left_active];
            rhs[flat] -= balance / stage.volume_cell[flat];
            rates[info.left_active] += rate[0];
        }
        if (info.right_active != hd_core::MeshBlock::invalid_active_index)
        {
            const auto flat = stage.active_to_storage[info.right_active];
            rhs[flat] += balance / stage.volume_cell[flat];
            rates[info.right_active] += rate[1];
        }
    }
    double maximum = 0.0;
    for (double rate : rates)
    {
        if (!std::isfinite(rate))
            throw std::invalid_argument("Non-finite accumulated transport rate");
        maximum = std::max(maximum, rate);
    }
    return maximum;
}
}
