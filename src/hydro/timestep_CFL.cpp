#include "timestep_CFL.h"
#include "core.h"
#include "mesh.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
namespace hydro
{
    double compute_cfl_dt(
        const hd_core::MeshBlock& mesh,
        double cfl
    )
    {
        if (!std::isfinite(cfl) || cfl <= 0.0)
            throw std::invalid_argument("CFL must be finite and positive");
        double max_inverse_dt = 0.0;
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {

            const auto& index = mesh.coordinates_active(flat);
            
            bool is_active = true;

            for (int d = 0; d < hd_core::NDIM; ++d)
            {
                if (index[d] < 0 || index[d] >= mesh.nactive[d])
                {
                    is_active = false;
                    break;
                }
            }

            if (!is_active)
                continue;
            
            double inverse_dt = 0.0;
            const auto gamma = mesh.cells[flat].gamma;
            const auto& W = mesh.cells[flat].U.to_Primtivies(gamma);
            if(W[hd_core::RHO] == 0.0)
            {
                continue;
            }
            const double c_s = std::sqrt(gamma *W[hd_core::PRE]/W[hd_core::RHO]);
            for(int d = 0; d < hd_core::NDIM; d++)
            {
                inverse_dt += (std::abs(W[hd_core::VEL(d)]) + c_s)/mesh.cell_width(d,index[d]+mesh.nghost[d]);
            }  
            if (!std::isfinite(inverse_dt))
                throw std::runtime_error("Non-finite CFL rate");
            max_inverse_dt = std::max(max_inverse_dt,inverse_dt);
        }
        if (max_inverse_dt <= 0.0)
            throw std::runtime_error("No positive CFL rate found");
        return cfl / max_inverse_dt;
    }
}