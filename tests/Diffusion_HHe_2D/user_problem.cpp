#include "user_config.h"
#include "user/user_api.h"
#include "config/config.h"
#include "source/source.h"
#include "microphysics/eos/eos.h"
#include "microphysics/composition/species.h"
#include "physcon/constants.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

// Test-only constant diffusivity. Compile with -DHHE_DIFFUSIVITY_CGS=0
// for the otherwise identical, stationary no-diffusion control.
#ifndef HHE_DIFFUSIVITY_CGS
#define HHE_DIFFUSIVITY_CGS 1e13
#endif

static_assert(hd_core::NDIM == 2 && hd_core::NCOMP == 2);
static_assert(HHE_DIFFUSIVITY_CGS >= 0, "Diffusivity must be nonnegative");

namespace
{
constexpr double density = 1e-10; // g/cm^3
constexpr double pressure = 1.0; // dyn/cm^2
constexpr double radius = 0.20;  // Fraction of box length
constexpr double width = 0.025;
constexpr double trace_fraction = 1e-4;

void initialize_units(physcon::Unit& units, const config::SimulationConfig&)
{
    units.init_unitsystem_MLT(1.0, 1.0, 1.0, 1.0); // CGS; T in kelvin.
}

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    using namespace hd_core;
    if (hd_microphysics::eos_type != hd_microphysics::Equation_of_State::IdealGas_Mixed
        || std::abs(cfg.gamma - 5.0/3.0) > 1e-12)
        throw std::invalid_argument("Diffusion_HHe_2D requires neutral monatomic IdealGas_Mixed");
    const auto& species = hd_microphysics::species_table.species_active;
    if (species.size() != 2 || species[0].Z != 1 || species[0].A != 1
        || species[1].Z != 2 || species[1].A != 4
        || species[0].charge != 0 || species[1].charge != 0)
        throw std::invalid_argument("Diffusion_HHe_2D requires H I and He I in that order");
    for (auto flat : mesh.active_to_storage)
    {
        const auto index = mesh.coordinates(flat);
        const double x = (mesh.cell_center(0, index)-cfg.lower[0])/mesh.active_length(0)-0.5;
        const double y = (mesh.cell_center(1, index)-cfg.lower[1])/mesh.active_length(1)-0.5;
        const double hydrogen = 0.5*(1.0-std::tanh((std::hypot(x,y)-radius)/width));
        Primitives w{};
        w[RHO] = density;
        w[PRE] = pressure;
        w[COMP(0)] = trace_fraction + (1.0-2.0*trace_fraction)*hydrogen;
        w[COMP(1)] = 1.0-w[COMP(0)];
        auto& cell = mesh.cells[flat];
        cell.gamma = cfg.gamma;
        cell.Set_Conserved_From_Primitives(w);
    }
}

// NEW: binary mass diffusion, using the Transport source interface. The solver
// owns divergence, time integration and the equal/opposite exchange at a face.
void diffusion_flux(const hd_core::MeshBlock& mesh, std::size_t face, double,
                    const config::SimulationConfig& cfg, hd_core::Conserved& flux,
                    std::array<double,2>& inverse_dt)
{
    using namespace hd_core;
    const auto& info = mesh.face_info[face];
    // Periodic seams have two active neighbors. Insulate any physical boundary.
    if (info.left_active == MeshBlock::invalid_active_index
        || info.right_active == MeshBlock::invalid_active_index) return;
    const int d = info.direction;
    auto right = info.storage_coordinates, left = right;
    --left[d];
    const auto& L = mesh.at(left);
    const auto& R = mesh.at(right);
    const double distance = mesh.cell_center(d,right)-mesh.cell_center(d,left);
    const double rho_face = 0.5*(L.U[RHO]+R.U[RHO]);
    const double rho_D = rho_face*HHE_DIFFUSIVITY_CGS;
    const double J = -rho_D*(R.U[COMP(0)]/R.U[RHO]-L.U[COMP(0)]/L.U[RHO])/distance;
    flux[COMP(0)] = J;
    flux[COMP(1)] = -J; // No net barycentric mass flux; do not change RHO or momentum.

    // NEW: enthalpy carried by the species, h_k = gamma/(gamma-1)*R_g*T/A_k.
    // Frozen neutral species; no Saha, thermal conduction or radiation here.
    // Each species carries its donor cell's temperature. Centering T here can
    // drain a cold cell using its hot neighbor's enthalpy and lose positivity.
    const double T_H = J >= 0.0 ? L.temperature : R.temperature;
    const double T_He = J >= 0.0 ? R.temperature : L.temperature;
    const double h_difference = cfg.gamma/(cfg.gamma-1.0)*physcon::R_g*(T_H-0.25*T_He);
    flux[ENERGY] = flux[INTERNAL_ENERGY] = h_difference*J;

    // Include A/V but not CFL or dt. Factor 8 is a conservative explicit-step
    // allowance for this H/He + enthalpy demonstration, not a universal closure.
    const double rate = 8.0*mesh.face_area[d].at(right)*rho_D/distance;
    const auto lf = mesh.active_to_storage[info.left_active];
    const auto rf = mesh.active_to_storage[info.right_active];
    inverse_dt[0] = rate/(mesh.volume_cell[lf]*L.U[RHO]);
    inverse_dt[1] = rate/(mesh.volume_cell[rf]*R.U[RHO]);
}

void register_sources(const config::SimulationConfig&, const hd_core::MeshBlock&,
                      hd_source::SourceList& sources)
{
    std::cout << "H/He diffusion coefficient = " << HHE_DIFFUSIVITY_CGS << " cm^2/s\n";
    if (HHE_DIFFUSIVITY_CGS == 0) return;
    hd_source::SourceTerm term;
    term.name = "HHe_mass_and_enthalpy_diffusion";
    term.kind = hd_source::SourceKind::Transport;
    term.transport_flux = diffusion_flux;
    sources.push_back(term);
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize_units = initialize_units;
    hooks.initialize = initialize;
    hooks.register_sources = register_sources;
    return hooks;
}
