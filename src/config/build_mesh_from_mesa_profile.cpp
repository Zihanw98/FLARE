/**
 * @file build_mesh_from_mesa_profile.cpp
 * @brief Construct spherical FLARE cell geometry and optional hydrodynamic states from a MESA profile.
 */
#include "read_mesa_profile.h"
#include "physcon/unit.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace config
{
namespace
{
    /**
     * @brief Initialize the active cells that have a MESA zone row.
     * @param[in,out] mesh Allocated radial mesh. Row zero maps to the first active storage cell.
     * @param[in] zones MESA table in mesh order, including any copied inner padding rows.
     * @param[in] gamma Ideal-gas adiabatic index passed to the FLARE cell constructor.
     * @pre The number of rows is at most the active cell count.
     * @throws std::runtime_error If a required column is missing.
     * @throws std::invalid_argument If a unit scale or converted gas state is invalid.
     */
    void fill_mesa_states(hd_core::MeshBlock& mesh, const Table& zones, double gamma)
    {
        // CGS units per code unit of density, pressure and velocity.
        const auto& units = physcon::units;
        if (!(std::isfinite(units.unit_density) && units.unit_density > 0.0 &&
              std::isfinite(units.unit_pressure) && units.unit_pressure > 0.0 &&
              std::isfinite(units.unit_velocity) && units.unit_velocity > 0.0))
            throw std::invalid_argument("MESA state units must be finite and positive");

        // logRho and logP are base-10 logarithms of CGS values; velocity is signed cm/s.
        const auto log_rho = zones.column("logRho");
        const auto log_p = zones.column("logP");
        const auto velocity = zones.column("velocity");
        const auto first = static_cast<std::size_t>(mesh.nghost[0]); // First active storage cell.
        for (std::size_t i = 0; i < zones.data.size(); ++i)
        {
            hd_core::Primitives state{}; // FLARE density, radial velocity and pressure for row i.
            state[hd_core::RHO] = std::pow(10.0, log_rho[i]) / units.unit_density;
            state[hd_core::PRE] = std::pow(10.0, log_p[i]) / units.unit_pressure;
            state[hd_core::VEL(0)] = velocity[i] / units.unit_velocity;
            if (!(state[hd_core::RHO] > 0.0))
                throw std::invalid_argument("MESA density is not positive in code units");
            mesh.cells[first + i] = hd_core::Cell(state, gamma);
        }
    }
}

hd_core::MeshBlock build_mesh_from_mesa_profile(const MesaProfile& profile,
                                                SimulationConfig& cfg,
                                                bool initialize_state)
{
    if (hd_core::NDIM != 1)
        throw std::invalid_argument("Building a MESA radial mesh requires HD_NDIM=1");
    const auto& radii = profile.radius_faces; // Active faces from the inner cut to the outer boundary.
    if (radii.size() < 2)
        throw std::invalid_argument("A MESA mesh requires at least two active faces");
    const auto active = radii.size() - 1; // Includes exterior padding cells without MESA rows.
    if (active > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::overflow_error("MESA active cell count exceeds the mesh index range");
    if (initialize_state && profile.zones.data.size() > active)
        throw std::invalid_argument("MESA profile has more state rows than active cells");

    hd_core::MeshBlock mesh; // Build locally and publish the derived config after successful construction.
    mesh.geometry = hd_core::Geometry::Spherical;
    mesh.nactive[0] = static_cast<int>(active);
    mesh.nghost = cfg.nghost;
    mesh.bc_lower[0] = boundary::parse_bc(cfg.boundary_lower[0]);
    mesh.bc_upper[0] = boundary::parse_bc(cfg.boundary_upper[0]);
    const int total = mesh.ntotal(0); // Includes both ghost regions; checks size and overflow.
    if (total == std::numeric_limits<int>::max())
        throw std::overflow_error("MESA face count exceeds the mesh index range");
    auto& faces = mesh.face[0]; // Full storage coordinates including ghost faces.
    faces.resize(static_cast<std::size_t>(total) + 1);
    const int first = mesh.nghost[0]; // Inner active face in storage coordinates.
    const int last = first + mesh.nactive[0]; // Outer active face.
    std::copy(radii.begin(), radii.end(), faces.begin() + first);

    // Continue each boundary shell's linear width through its own ghost region.
    const double inner_width = faces[first + 1] - faces[first];
    const double outer_width = faces[last] - faces[last - 1];
    for (int g = 1; g <= mesh.nghost[0]; ++g) // g counts ghost layers away from the active face.
    {
        faces[first - g] = std::fma(-static_cast<double>(g), inner_width, faces[first]);
        faces[last + g] = std::fma(static_cast<double>(g), outer_width, faces[last]);
    }
    mesh.allocate();
    mesh.build_geometry();
    if (initialize_state && !profile.zones.data.empty())
        fill_mesa_states(mesh, profile.zones, cfg.gamma);

    cfg.geometry = "spherical";
    cfg.nactive = mesh.nactive;
    cfg.lower[0] = radii.front();
    cfg.upper[0] = radii.back();
    return mesh;
}
}
