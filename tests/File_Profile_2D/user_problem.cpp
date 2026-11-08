#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <cmath>
#include <stdexcept>

namespace
{
void initialize_units(physcon::Unit& units, const config::SimulationConfig&)
{
    units.init_unitsystem_MLT(2.0, 3.0, 5.0, 7.0);
}

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    // The framework has already populated the geometry and states from the two files.
    // This hook checks that handoff and deliberately leaves those states intact.
    if (!cfg.init_mesh_from_file || !cfg.init_profile_from_file ||
        mesh.geometry != hd_core::Geometry::Cartesian ||
        mesh.nactive[0] != 32 || mesh.nactive[1] != 32 ||
        cfg.lower[0] != 0.0 || cfg.upper[0] != 1.0 ||
        cfg.lower[1] != 0.0 || cfg.upper[1] != 0.75 ||
        physcon::units.unit_length != 3.0)
        throw std::runtime_error("File-profile initialization or unit hook was skipped");
    mesh.validate_active();
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize_units = initialize_units;
    hooks.initialize = initialize;
    return hooks;
}
