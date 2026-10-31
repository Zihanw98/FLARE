#include "godunov.h"
#include "wave_timestep.h"
#include "riemann/riemann.h"
#include "facefield.h"
#include "core.h"
#include "core/mesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace hydro
{
    void compute_godunov_rhs(
        const hd_core::MeshBlock& mesh,
        const RiemannSolver& solver,
        std::vector<hd_core::Conserved>& rhs
    )
    {
        // 每次 RHS 都重新统计，避免沿用上一个时间步的限制。
        wave_timestep = {};
        if (mesh.cells.empty() || rhs.size() != mesh.cells.size())
            throw std::invalid_argument("Godunov requires allocated mesh and matching RHS storage");
        mesh.validate_geometry();
        auto nactive = mesh.nactive;
        auto nghost = mesh.nghost;
        std::array<FaceField, hd_core::NDIM> flux;
        // 与面通量使用相同的 flat 下标，避免重新求解 Riemann 问题。
        std::array<std::vector<double>, hd_core::NDIM> face_speed;
        for(int n=0;n<hd_core::NDIM;n++)
        {
            auto shape = nactive;
            shape[n]++;
            flux[n].allocate(shape);
            face_speed[n].resize(flux[n].faces.size());
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

                // 哨兵值：求解器若忘记写入当前面波速，就报错而不沿用旧值。
                wave_timestep.face_speed = -1.0;
                field.faces[i] = solver.flux(WL, WR, n, gamma_L, gamma_R);
                if (!std::isfinite(wave_timestep.face_speed) || wave_timestep.face_speed < 0.0)
                    throw std::runtime_error("Invalid Riemann wave speed");
                face_speed[n][i] = wave_timestep.face_speed;
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

            rhs[flat] = {{}};
            double inverse_dt = 0.0;
            auto face_left_index = index;
            for(int n=0;n<hd_core::NDIM;n++)
            {
                auto ds_n = mesh.cell_width(n,index[n]+nghost[n]);
                auto face_right_index = index;
                face_right_index[n]++;
                const auto& F_nL = flux[n].at(face_left_index);
                const auto& F_nR = flux[n].at(face_right_index);
                rhs[flat] += -(F_nR-F_nL) / ds_n;
                const auto left_flat = flux[n].flat_index(face_left_index);
                const auto right_flat = flux[n].flat_index(face_right_index);
                inverse_dt += std::max(face_speed[n][left_flat],
                                       face_speed[n][right_flat]) / ds_n;
            }
            if (!std::isfinite(inverse_dt))
                throw std::runtime_error("Non-finite Godunov CFL rate");
            wave_timestep.max_inverse_dt = std::max(wave_timestep.max_inverse_dt, inverse_dt);
        }
    }
}
