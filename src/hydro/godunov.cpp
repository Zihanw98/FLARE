#include "flux.h"
#include "godunov.h"
#include "riemann.h"
#include "facefield.h"
#include "core.h"
#include "mesh.h"
#include <stdexcept>
namespace hydro
{
    void compute_godunov_rhs(
        const hd_core::MeshBlock& mesh,
        const RiemannSolver& solver,
        std::vector<hd_core::Conserved>& rhs
    )
    {
        auto nactive = mesh.nactive;
        auto nghost = mesh.nghost;
        auto ntotal = mesh.ntotal();
        std::array<FaceField, hd_core::NDIM> flux;
        for(int n=0;n<hd_core::NDIM;n++)
        {
            auto shape = nactive;
            shape[n]++;
            flux[n].allocate(shape);
        }
        for (int n = 0; n < hd_core::NDIM; ++n)
        {
            if (nghost[n] < 1)
                throw std::invalid_argument(
                    "Godunov requires at least one ghost cell per direction"
                );

            auto& field = flux[n];

            for (std::size_t i = 0; i < field.faces.size(); ++i)
            {
                auto right_index = field.coordinates(i);

                for (int d = 0; d < hd_core::NDIM; ++d)
                    right_index[d] += nghost[d];

                auto left_index = right_index;
                left_index[n]--;

                const auto& cellLeft  = mesh.at(left_index);
                const auto& cellRight = mesh.at(right_index);

                const auto gamma_L = cellLeft.gamma;
                const auto gamma_R = cellRight.gamma;

                const auto WL = cellLeft.U.to_Primtivies(gamma_L);
                const auto WR = cellRight.U.to_Primtivies(gamma_R);

                field.faces[i] =
                    solver.flux(WL, WR, n, gamma_L, gamma_R);
            }
        }
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {

            const auto index = mesh.coordinates_active(flat);
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

            const auto& cell = mesh.cells[flat];

            rhs[flat] = {{}};
            auto face_left_index = index;
            for(int n=0;n<hd_core::NDIM;n++)
            {
                auto ds_n = mesh.cell_width(n,index[n]+nghost[n]);
                auto face_right_index = index;
                face_right_index[n]++;
                const auto& F_nL = flux[n].at(face_left_index);
                const auto& F_nR = flux[n].at(face_right_index);
                rhs[flat] += -(F_nR-F_nL) / ds_n;
            }
        }

    }
}