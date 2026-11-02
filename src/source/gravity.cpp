#pragma once
#include "source.h"
#include "core/core.h"
#include "core/mesh.h"
namespace hd_source
{
    void constant_gravity(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    )
    {
        //Don't change dt_limit
        for(int d=0; d < hd_core::NDIM; d++)
        {
            for(size_t flat = 0; flat < mesh.cells.size(); flat++)
            {
                rhs[flat][hd_core::MOMDEN(d)] += mesh.cells[flat].U[hd_core::RHO] * cfg.gravity[d];
                rhs[flat][hd_core::ENERGY] += mesh.cells[flat].U[hd_core::MOMDEN(d)] * cfg.gravity[d];
            }
        }
    }
}