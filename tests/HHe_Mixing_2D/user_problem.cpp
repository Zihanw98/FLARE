#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include "microphysics/composition/species.h"
#include "physcon/constants.h"
#include <cmath>
#include <stdexcept>

static_assert(hd_core::NDIM == 2 && hd_core::NCOMP == 2 && hd_core::NUSER == 0,
              "HHe_Mixing_2D requires NDIM=2, NCOMP=2, NUSER=0");

namespace
{
    void initialize_units(physcon::Unit& units, const config::SimulationConfig&)
    {
        // A cold, neutral monatomic ideal-gas demonstration; T is derived in Python.
        units.init_unitsystem_RHOLV(1.e-24, 1.e18, 1.e5, 1.0);
    }

    void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
    {
        if (!cfg.enable_composition || std::abs(cfg.gamma - 5.0/3.0) > 1.e-14)
            throw std::invalid_argument("H/He mixing requires composition and gamma=5/3");
        for (int d = 0; d < 2; ++d)
            if (mesh.bc_lower[d] != boundary::BoundaryCondition::Periodic ||
                mesh.bc_upper[d] != boundary::BoundaryCondition::Periodic)
                throw std::invalid_argument("H/He mixing requires periodic boundaries");
        const auto& species = hd_microphysics::species_table;
        const int hydrogen = species.find_species_from_available_table_index(1, 1, 0);
        const int helium = species.find_species_from_available_table_index(4, 2, 0);
        for (int j = 0; j < mesh.nactive[1]; ++j)
            for (int i = 0; i < mesh.nactive[0]; ++i)
            {
                const hd_core::Index index{i + mesh.nghost[0], j + mesh.nghost[1]};
                const double x = (mesh.cell_center(0, index) - cfg.lower[0]) / mesh.active_length(0);
                const double y = (mesh.cell_center(1, index) - cfg.lower[1]) / mesh.active_length(1);
                // This cosine-based layer is smooth and exactly periodic, including its derivative.
                const double layer = 0.5 * (1.0 - std::tanh(std::cos(2*physcon::pi*y) / 0.15));
                const double shear = std::cos(2*physcon::pi*y) / 0.3;
                hd_core::Primitives w{};
                w[hd_core::RHO] = 1.0 + layer;
                w[hd_core::PRE] = 2.5;
                w[hd_core::VEL(0)] = layer - 0.5;
                w[hd_core::VEL(1)] = 0.02 * std::sin(4*physcon::pi*x) * std::exp(-shear*shear);
                w[hd_core::COMP(hydrogen)] = 1.0 - layer;
                w[hd_core::COMP(helium)] = layer;
                // The framework converts X_k -> rho X_k and transports every component.
                mesh.at(index) = hd_core::Cell(w, cfg.gamma);
            }
    }
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    hooks.initialize_units = initialize_units;
    return hooks;
}
