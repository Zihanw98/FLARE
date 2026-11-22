#include "user_config.h"
#include "user/user_api.h"
#include "config/config.h"
#include <stdexcept>

static_assert(hd_core::NDIM == 2 && hd_core::NCOMP == 2);
namespace
{
void initialize_units(physcon::Unit& units, const config::SimulationConfig&)
{
    units.init_unitsystem_MLT(1.0, 1.0, 1.0, 1.0); // Profiles use CGS and kelvin.
}

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (!cfg.init_profile_from_file)
        throw std::invalid_argument("Composition_Advection_2D requires a generated profile");
    // read_profile already built U through RhoP; the common initialization
    // fills ghosts and refreshes T and mu. Do not replace the profile here.
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
