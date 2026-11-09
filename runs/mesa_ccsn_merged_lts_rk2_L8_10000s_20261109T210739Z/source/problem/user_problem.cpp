#include "user_config.h"
#include "user/user_api.h"
#include "core/mesh.h"
#include "config/config.h"
#include "physcon/constants.h"
#include "hydro/timestep_CFL.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

static_assert(hd_core::NDIM == 1, "MESA CCSN requires spherical 1D");

namespace
{
constexpr double solar_radius = 6.957e10;
constexpr double solar_gm = 1.3271244e26;
constexpr double injection_energy_erg = 1.6e51;
constexpr double atmosphere_density_cgs = 1e-15;
constexpr double atmosphere_temperature = 1e4;
constexpr double atmosphere_mu = 0.62;

void initialize_units(physcon::Unit& units, const config::SimulationConfig&)
{
    units.init_unitsystem_MLT(solar_gm / physcon::G_, solar_radius,
                             std::sqrt(std::pow(solar_radius, 3) / solar_gm), 1.0);
}

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (!cfg.init_from_mesa_profile || mesh.geometry != hd_core::Geometry::Spherical ||
        !cfg.spherical_gravity_enabled || !cfg.spherical_self_gravity_enabled ||
        cfg.spherical_central_mass <= 0 || cfg.constant_heating_enabled || cfg.constant_gravity_enabled)
        throw std::invalid_argument("MESA CCSN requires a loaded spherical profile and central/self gravity");
    const auto original_faces = mesh.face; // Assert the initializer leaves the imported grid intact.
    const auto& units = physcon::units;
    const double atmosphere_pressure_cgs = atmosphere_density_cgs * physcon::k_B_
        * atmosphere_temperature / (atmosphere_mu * physcon::m_u_);
    int native_cells = 0;
    bool exterior = false;
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        auto& cell = mesh.at_active({i});
        if (cell.U[hd_core::RHO] == 0.0)
        {
            exterior = true;
            hd_core::Primitives ambient{};
            ambient[hd_core::RHO] = atmosphere_density_cgs / units.unit_density;
            ambient[hd_core::PRE] = atmosphere_pressure_cgs / units.unit_pressure;
            cell = hd_core::Cell(ambient, cfg.gamma);
        }
        else
        {
            if (exterior)
                throw std::runtime_error("Uninitialized cells inside the MESA envelope");
            ++native_cells;
        }
    }
    if (native_cells == 0 || native_cells == mesh.nactive[0])
        throw std::runtime_error("This example expects a retained MESA envelope and exterior padding");

    std::filesystem::create_directories(cfg.output_directory);
    std::ofstream profile(std::filesystem::path(cfg.output_directory) / "initial_envelope.csv");
    profile << std::setprecision(17)
            << "r_left,r_right,volume,rho,vx,p,gamma,energy_density_before,is_atmosphere\n";
    long double envelope_mass = 0, atmosphere_mass = 0;
    double minimum_width = std::numeric_limits<double>::infinity();
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        const auto flat = mesh.flat_index({i}, true);
        const auto [a, b] = mesh.cell_face_pos(0, i + mesh.nghost[0]);
        const auto& cell = mesh.cells[flat];
        const auto w = cell.U.to_Primtivies(cell.gamma);
        const double volume = mesh.volume_cell[flat];
        (i < native_cells ? envelope_mass : atmosphere_mass) += w[hd_core::RHO] * volume;
        minimum_width = std::min(minimum_width, b - a);
        profile << a << ',' << b << ',' << volume << ',' << w[hd_core::RHO] << ','
                << w[hd_core::VEL(0)] << ',' << w[hd_core::PRE] << ',' << cell.gamma << ','
                << cell.U[hd_core::ENERGY] << ',' << (i >= native_cells) << '\n';
    }
    profile.close();
    if (!profile) throw std::runtime_error("Cannot write MESA initial envelope");

    // Same thermal bomb as CCSN_Simple: add 1.6e51 erg to the first active shell.
    // No mesh rebuild, no Lane-Emden model, and no change to imported density/momentum.
    const auto first = mesh.flat_index({0}, true);
    mesh.cells[first].U[hd_core::ENERGY] +=
        injection_energy_erg / units.unit_energy / mesh.volume_cell[first];
    mesh.validate_active();
    if (mesh.face != original_faces)
        throw std::logic_error("Thermal bomb changed the MESA grid");
    const double initial_dt = hydro::compute_cfl_dt(mesh, cfg.cfl);

    std::ofstream info(std::filesystem::path(cfg.output_directory) / "model.json");
    info << std::setprecision(17)
         << "{\n  \"mass_unit_g\": " << units.unit_mass
         << ",\n  \"length_unit_cm\": " << units.unit_length
         << ",\n  \"time_unit_s\": " << units.unit_time
         << ",\n  \"density_unit_g_cm3\": " << units.unit_density
         << ",\n  \"velocity_unit_cm_s\": " << units.unit_velocity
         << ",\n  \"pressure_unit_dyn_cm2\": " << units.unit_pressure
         << ",\n  \"energy_unit_erg\": " << units.unit_energy
         << ",\n  \"G_code\": " << physcon::G
         << ",\n  \"core_mass_code\": " << cfg.spherical_central_mass
         << ",\n  \"cut_radius_code\": " << cfg.lower[0]
         << ",\n  \"surface_radius_code\": " << mesh.face[0][mesh.nghost[0] + native_cells]
         << ",\n  \"outer_radius_code\": " << cfg.upper[0]
         << ",\n  \"native_cells\": " << native_cells
         << ",\n  \"atmosphere_cells\": " << mesh.nactive[0] - native_cells
         << ",\n  \"envelope_mass_code\": " << envelope_mass
         << ",\n  \"atmosphere_mass_code\": " << atmosphere_mass
         << ",\n  \"atmosphere_density_cgs\": " << atmosphere_density_cgs
         << ",\n  \"atmosphere_pressure_cgs\": " << atmosphere_pressure_cgs
         << ",\n  \"injection_energy_erg\": " << injection_energy_erg
         << ",\n  \"injection_cell_count\": 1"
         << ",\n  \"injection_mass_code\": " << mesh.cells[first].U[hd_core::RHO] * mesh.volume_cell[first]
         << ",\n  \"minimum_cell_width_cm\": " << minimum_width * units.unit_length
         << ",\n  \"initial_cfl_dt_s\": " << initial_dt * units.unit_time
         << ",\n  \"target_time_s\": " << cfg.t_end * units.unit_time
         << ",\n  \"grid_preserved\": true\n}\n";
    info.close();
    if (!info) throw std::runtime_error("Cannot write MESA model metadata");
    std::cout << "MESA CCSN: retained " << native_cells << " native shells, "
              << mesh.nactive[0] - native_cells << " atmosphere shells; initial dt="
              << initial_dt * units.unit_time << " s\n";
}

void outer_boundary(hd_core::MeshBlock& mesh, double, int direction,
                    boundary::BoundarySide side)
{
    if (direction != 0 || side != boundary::BoundarySide::Upper)
        throw std::invalid_argument("Only the outer radial boundary uses userbc");
    const int end = mesh.nghost[0] + mesh.nactive[0];
    auto state = mesh.at({end - 1}).U.to_Primtivies(mesh.at({end - 1}).gamma);
    state[hd_core::VEL(0)] = std::max(0.0, state[hd_core::VEL(0)]);
    for (int i = end; i < mesh.ntotal(0); ++i)
        mesh.at({i}) = hd_core::Cell(state, mesh.at({end - 1}).gamma);
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize_units = initialize_units;
    hooks.initialize = initialize;
    hooks.userbc = outer_boundary;
    return hooks;
}
