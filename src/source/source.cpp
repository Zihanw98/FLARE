/**
 * @file source.cpp
 * @brief Physical source terms.
 */
#include "source.h"
#include "user/user_api.h"

namespace hd_source
{
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
