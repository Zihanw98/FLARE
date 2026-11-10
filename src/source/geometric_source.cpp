#include "source.h"
#include "utils/performance.h"
namespace hd_source
{
    SourceTerm make_geometric_source()
    {
        return {"geometric", SourceKind::General, geometric_source, true};
    }

    void geometric_source(
        const hd_core::MeshBlock& mesh,
        double /*time*/,
        const config::SimulationConfig& /*cfg*/,
        std::vector<hd_core::Conserved>& rhs,
        double& /*dt_limit*/
    )
    {
        performance::Scope timer(performance::Section::Geometry);
        if(mesh.geometry==hd_core::Geometry::Spherical)
        {
            
            if(hd_core::NDIM==1)
            {
                for(size_t flat = 0; flat < mesh.cells.size(); flat++)
                {
                    bool is_active = true;
                    const auto index = mesh.coordinates(flat);
                    auto V = mesh.volume_cell[flat];
                    for (int d = 0; d < hd_core::NDIM; ++d)
                    {
                        if (index[d] < mesh.nghost[d] || index[d] >= (mesh.nghost[d]+mesh.nactive[d]))
                        {
                            is_active = false;
                            break;
                        }
                    }

                    if (!is_active)
                        continue;  
                    const auto& cell = mesh.cells[flat];
                    const auto W = cell.U.to_Primtivies(cell.gamma);
                    const double p = W[hd_core::PRE];
                    
                    
                    auto index_right = index;
                    index_right[0]++;
                    const double A_R = mesh.face_area[0].at(index_right);
                    const double A_L = mesh.face_area[0].at(index);
                    rhs[flat][hd_core::MOMDEN(0)] += p * (A_R - A_L) / V;
                }

            }
        }
    }
}
