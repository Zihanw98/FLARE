/**
 * @file source.cpp
 * @brief Physical source terms.
 */
#include "source.h"
#include "user/user_api.h"

namespace hd_source
{
    void complete_internal_source(const hd_core::MeshBlock& mesh, const SourceTerm& term,
                                  std::vector<hd_core::Conserved>& rates)
    {
        if (rates.size() != mesh.cells.size()) throw std::logic_error("Source changed RHS storage size");
        if (term.supplies_internal_energy) return;
        for (std::size_t i = 0; i < rates.size(); ++i)
        {
            const auto& U = mesh.cells[i].U;
            auto& S = rates[i];
            if (U[hd_core::RHO] == 0.0) continue;
            double thermal = S[hd_core::ENERGY];
            for (int d = 0; d < hd_core::NDIM; ++d)
            {
                const double v = U[hd_core::MOMDEN(d)] / U[hd_core::RHO];
                thermal += -v * S[hd_core::MOMDEN(d)] + 0.5 * v * v * S[hd_core::RHO];
            }
            S[hd_core::INTERNAL_ENERGY] += thermal;
        }
    }

    SourceList make_sources(
        const config::SimulationConfig& cfg,
        const hd_core::MeshBlock& mesh,
        const hd_user::UserHooks& hooks
    )
    {
        SourceList sources;
        if(cfg.constant_gravity_enabled)
        {
            sources.push_back(make_constant_gravity({cfg.gravity}));
        }
        if(mesh.geometry != hd_core::Geometry::Cartesian)
        {
            sources.push_back(make_geometric_source());
        }
        if (cfg.constant_heating_enabled)
        {
            sources.push_back(make_constant_heating(cfg.heating_rate));
        }
        if (hooks.register_sources) //User sources
        {
            hooks.register_sources(cfg, mesh, sources);
        }
        if (cfg.spherical_gravity_enabled)
        {
            if (mesh.geometry != hd_core::Geometry::Spherical)
                throw std::invalid_argument(
                    "Spherical gravity requires spherical geometry");

            sources.push_back(make_spherical_self_gravity({
                cfg.spherical_central_mass,
                cfg.spherical_self_gravity_enabled
            }));
        }
        return sources;
    }
}
