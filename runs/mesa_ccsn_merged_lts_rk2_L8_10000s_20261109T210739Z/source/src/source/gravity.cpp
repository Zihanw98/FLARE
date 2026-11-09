/**
 * @file gravity.cpp
 * @brief Physical source terms.
 */
#include "source.h"
#include "utils/performance.h"
#include "core/core.h"
#include "core/mesh.h"
#include "physcon/unit.h"
#include "physcon/constants.h"
#include <array>
namespace hd_source
{
    SourceTerm make_constant_gravity(ConstantGravityParams params)
    {
        return {"constant_gravity", SourceKind::Force,
            [params](const hd_core::MeshBlock& mesh, double time,
                     const config::SimulationConfig& cfg,
                     std::vector<hd_core::Conserved>& rhs, double& dt_limit) {
                constant_gravity(mesh, time, cfg, params, rhs, dt_limit);
            }};
    }

    void constant_gravity(
        const hd_core::MeshBlock& mesh,
        double /*time*/,
        const config::SimulationConfig& /*cfg*/,
        const ConstantGravityParams& params,
        std::vector<hd_core::Conserved>& rhs,
        double& /*dt_limit*/
    )
    {
        performance::Scope timer(performance::Section::Gravity);
        //Don't change dt_limit
        for(int d=0; d < hd_core::NDIM; d++)
        {
            for(size_t flat = 0; flat < mesh.cells.size(); flat++)
            {
                rhs[flat][hd_core::MOMDEN(d)] += mesh.cells[flat].U[hd_core::RHO] * params.acceleration[d];
                rhs[flat][hd_core::ENERGY] += mesh.cells[flat].U[hd_core::MOMDEN(d)] * params.acceleration[d];
            }
        }
    }
    SourceTerm make_spherical_self_gravity(SphericalGravityParams params)
    {
        return {"spherical_gravity", SourceKind::Force,
            [params](const hd_core::MeshBlock& mesh, double time,
                     const config::SimulationConfig& cfg,
                     std::vector<hd_core::Conserved>& rhs, double& dt_limit) {
                spherical_self_gravity(mesh, time, cfg, params, rhs, dt_limit);
        }};
    }
    void spherical_self_gravity(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& /*cfg*/,
        const SphericalGravityParams& params,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    )
    {
        performance::Scope timer(performance::Section::Gravity);
        //Central Gravity + Self Gravity

        const std::size_t r_total = mesh.nactive[0];
        std::vector<double> shell_mass(r_total, 0.0);
        std::vector<double> shell_g(r_total, 0.0);
        for(size_t flat=0; flat<mesh.cells.size(); flat++)
        {
            const auto& Index = mesh.coordinates(flat);
            if(!mesh.is_active(flat))
            {
                continue;
            }
            if(params.self_gravity_enabled)
            {
                auto ir = mesh.coordinates_active(flat)[0];
                auto volume = mesh.volume_cell[flat];
                auto rho = mesh.cells[flat].U[hd_core::RHO];
                shell_mass[ir] += volume * rho;
            }
        }
        double enclosed = params.central_mass;
        for(int ir = 0;ir<r_total;ir++)
        {
            const int nr = ir + mesh.nghost[0];
            const auto [rL, rR] = mesh.cell_face_pos(0, nr);

            const double D = rL*rL + rL*rR + rR*rR;
            const double alpha = (rR*rR + 2.0*rL*rR + 3.0*rL*rL) / (4.0*D);

            shell_g[ir] = -3.0 * physcon::G * (enclosed + alpha * shell_mass[ir])/D;

            enclosed += shell_mass[ir];
        }
        for(size_t flat=0; flat<mesh.cells.size(); flat++)
        {
            if (!mesh.is_active(flat))
                continue;
            const int ir = mesh.coordinates_active(flat)[0];
            const double g = shell_g[ir];
            const auto& U = mesh.cells[flat].U;

            rhs[flat][hd_core::MOMDEN(0)] += U[hd_core::RHO] * g;
            rhs[flat][hd_core::ENERGY]    += U[hd_core::MOMDEN(0)] * g;
        }
    }
}
