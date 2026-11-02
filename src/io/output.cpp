#include "output.h"
#include "config/config.h"
#include "core/mesh.h"
#include "user/user_api.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace io
{
    void Output(
        hd_core::MeshBlock& mesh,
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks,
        std::size_t step,
        double time,
        double dt
    )
    {
        if (cfg.step_every == 0)
            throw std::invalid_argument("step_every must be positive");
        if (!std::isfinite(time) || !std::isfinite(dt) || time < 0.0 || dt < 0.0)
            throw std::invalid_argument("Output requires finite nonnegative time and dt");
        if (mesh.cells.empty())
            throw std::invalid_argument("Cannot output an unallocated mesh");
        mesh.validate_geometry();

        const bool write_snapshot = step <= 1 || step % cfg.step_every == 0
                                 || time >= cfg.t_end || step >= cfg.maxsteps;
        std::filesystem::path filename;
        std::ofstream csv;
        if (write_snapshot)
        {
            if (cfg.output_directory.empty())
                throw std::invalid_argument("Output directory must not be empty");
            std::filesystem::create_directories(cfg.output_directory);
            std::ostringstream name;
            name << "mesh_" << cfg.output_suffix << '_'
                 << std::setfill('0') << std::setw(8) << step << ".csv";
            filename = std::filesystem::path(cfg.output_directory) / name.str();
            csv.open(filename);
            if (!csv)
                throw std::runtime_error("Cannot open output file: " + filename.string());
            csv.imbue(std::locale::classic());
            csv << std::setprecision(std::numeric_limits<double>::max_digits10);
            // Snapshot metadata belongs to the whole file. Python readers can
            // skip these three comment lines, e.g. pandas.read_csv(comment="#").
            csv << "# step=" << step << '\n'
                << "# time=" << time << '\n'
                << "# dt=" << dt << '\n';
            csv << "x,y,z,gamma,rho,vx,vy,vz,p";
            for (int k = 0; k < hd_core::NCOMP; ++k)
                csv << ",comp_" << k;
            for (int k = 0; k < hd_core::NUSER; ++k)
                csv << ",user_" << k;
            csv << ",is_ghost\n";
        }

        double rho_min = std::numeric_limits<double>::infinity();
        double rho_max = 0.0;
        double p_min = std::numeric_limits<double>::infinity();
        double p_max = 0.0;
        long double mass = 0.0L;
        long double energy = 0.0L;
        std::size_t active_cells = 0;
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            const auto index = mesh.coordinates(flat);
            bool is_ghost = false;
            for (int d = 0; d < hd_core::NDIM; ++d)
                if (index[d] < mesh.nghost[d]
                    || index[d] >= mesh.nghost[d] + mesh.nactive[d])
                    is_ghost = true;

            if (is_ghost && !write_snapshot)
                continue;
            const auto& cell = mesh.cells[flat];
            cell.U.validate(cell.gamma);
            const auto W = cell.U.to_Primtivies(cell.gamma);
            if (!is_ghost)
            {
                ++active_cells;
                rho_min = std::min(rho_min, W[hd_core::RHO]);
                rho_max = std::max(rho_max, W[hd_core::RHO]);
                p_min = std::min(p_min, W[hd_core::PRE]);
                p_max = std::max(p_max, W[hd_core::PRE]);
                long double volume = 1.0L;
                for (int d = 0; d < hd_core::NDIM; ++d)
                    volume *= mesh.cell_width(d, index);
                mass += cell.U[hd_core::RHO] * volume;
                energy += cell.U[hd_core::ENERGY] * volume;
            }
            if (write_snapshot)
            {
                std::array<double, 3> position{};
                std::array<double, 3> velocity{};
                for (int d = 0; d < hd_core::NDIM; ++d)
                {
                    position[d] = mesh.cell_center(d, index);
                    velocity[d] = W[hd_core::VEL(d)];
                }
                csv << position[0] << ',' << position[1] << ',' << position[2]
                    << ',' << cell.gamma << ',' << W[hd_core::RHO]
                    << ',' << velocity[0] << ',' << velocity[1] << ',' << velocity[2]
                    << ',' << W[hd_core::PRE];
                for (int k = 0; k < hd_core::NCOMP; ++k)
                    csv << ',' << W[hd_core::COMP(k)];
                for (int k = 0; k < hd_core::NUSER; ++k)
                    csv << ',' << W[hd_core::CONUSR(k)];
                csv << ',' << (is_ghost ? 1 : 0) << '\n';
            }
        }
        if (write_snapshot)
        {
            csv.close();
            if (!csv)
                throw std::runtime_error("Failed to write output file: " + filename.string());
        }

        // 使用独立的流设置精度，避免修改调用者的 cout 格式。
        std::ostringstream message;
        message.imbue(std::locale::classic());
        message << std::setprecision(10)
                << "step = " << step << ", time = " << time << ", dt = " << dt
                << ", active_cells = " << active_cells
                << ", rho = [" << rho_min << ", " << rho_max << ']'
                << ", p = [" << p_min << ", " << p_max << ']'
                << ", mass = " << mass << ", energy = " << energy;
        if (write_snapshot)
            message << ", csv = " << filename.string();
        std::cout << message.str() << '\n';

        if (write_snapshot && hooks.useroutput)
            hooks.useroutput(mesh, cfg);
    }
}
