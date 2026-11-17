#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include <stdexcept>

static_assert(hd_core::NDIM == 2, "Hydrogen_Profile_2D requires HD_NDIM=2");

namespace
{
void initialize_units(physcon::Unit& units, const config::SimulationConfig&)
{
    // Change the base units here if the profile represents a physical model.
    // The profile reader takes code-unit values and performs no CGS conversion.
    units.init_unitsystem_MLT(1.0, 1.0, 1.0, 1.0);
}

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (!cfg.init_profile_from_file || !cfg.enable_composition)
        throw std::invalid_argument("Hydrogen_Profile_2D requires profile input and composition");
    // read_profile has already populated all active cells. Edit profile.dat or
    // generate_profile.py to change the initial state; do not overwrite it here.
    // A future Saha initialization can be inserted here after reading the profile.
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
