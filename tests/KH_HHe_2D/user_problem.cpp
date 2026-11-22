#include "user_config.h"
#include "user/user_api.h"
#include "config/config.h"
#include "microphysics/eos/eos.h"
#include "physcon/constants.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

static_assert(hd_core::NDIM == 2 && hd_core::NCOMP == 2);

namespace
{
// Edit these constants to change the initial state; mesh/time/output settings
// remain in parameters.txt. Composition is frozen: there is no Saha or diffusion.
constexpr double initial_density = 1e-10;  // g/cm^3
constexpr double initial_pressure = 100;  // dyn/cm^2
constexpr double shear_speed = 5e5;        // cm/s, in either direction
constexpr double perturbation_speed = 1e4;
constexpr double shear_width = 0.02;      // Fraction of the box height
constexpr double perturbation_width = 0.05;
constexpr int perturbation_modes = 2;

void initialize_units(physcon::Unit& units, const config::SimulationConfig&)
{
    units.init_unitsystem_MLT(1.0, 1.0, 1.0, 1.0);
}

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    using namespace hd_core;
    if (hd_microphysics::eos_type != hd_microphysics::Equation_of_State::IdealGas_Mixed)
        throw std::invalid_argument("KH_HHe_2D uses IdealGas_Mixed with neutral H and He");
    const double lx = mesh.active_length(0), ly = mesh.active_length(1);
    for (int j = 0; j < mesh.nactive[1]; ++j)
        for (int i = 0; i < mesh.nactive[0]; ++i)
        {
            const Index index{i + mesh.nghost[0], j + mesh.nghost[1]};
            const double x = (mesh.cell_center(0, index) - cfg.lower[0]) / lx;
            const double y = (mesh.cell_center(1, index) - cfg.lower[1]) / ly;
            const double dy = mesh.cell_width(1, index) / ly;
            // Exact cell-average fraction of the initially pure H band.
            const double hydrogen = std::max(0.0, std::min(y + dy/2, 0.75)
                                                - std::max(y - dy/2, 0.25)) / dy;
            // Smooth, exactly periodic velocity shear at y/L=1/4 and 3/4.
            // Its finite thickness avoids seeding an unresolved vortex sheet.
            const double q = std::cos(2*physcon::pi*y);
            const double vx = -shear_speed * std::tanh(q/(2*physcon::pi*shear_width));
            const double envelope = std::exp(-std::pow(q/(2*physcon::pi*perturbation_width), 2));
            Primitives w{};
            w[RHO] = initial_density;
            w[PRE] = initial_pressure;
            w[VEL(0)] = vx;
            w[VEL(1)] = perturbation_speed * std::sin(2*physcon::pi*perturbation_modes*x) * envelope;
            w[COMP(0)] = std::clamp(hydrogen, 0.0, 1.0); // H I
            w[COMP(1)] = 1.0 - w[COMP(0)];             // He I
            auto& cell = mesh.at(index);
            cell.gamma = cfg.gamma;
            cell.Set_Conserved_From_Primitives(w); // RhoP; T depends on composition.
        }
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize_units = initialize_units;
    hooks.initialize = initialize;
    return hooks;
}
