/**
 * @file constant_heating.cpp
 * @brief Constant specific heating with an internal-energy timescale limit.
 */
#include "source.h"
#include "utils/performance.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hd_source
{
    // 也可以不用工厂函数，直接在注册位置创建 lambda（以下示例位于 hd_source 命名空间内）：
    //
    // double heating_rate = 0.5;  // 在这里提供参数
    // SourceFunction source =
    //     [heating_rate](const hd_core::MeshBlock& mesh,
    //                    double time,
    //                    const config::SimulationConfig& cfg,
    //                    std::vector<hd_core::Conserved>& rhs,
    //                    double& dt_limit)
    //     {
    //         // 可以直接用 heating_rate 写计算，也可以传给已有的计算函数。
    //         constant_heating(mesh, time, cfg, heating_rate, rhs, dt_limit);
    //     };
    // sources.push_back({"constant_heating", SourceKind::Thermal, source});
    //
    // [heating_rate] 在创建 lambda 时保存值的副本；之后修改外面的变量不会改变它。
    // 此写法与调用 make_constant_heating(0.5) 二选一，避免重复注册同一份加热。
    SourceTerm make_constant_heating(double heating_rate)
    {
        if (!std::isfinite(heating_rate) || heating_rate < 0.0)
            throw std::invalid_argument("heating_rate must be finite and nonnegative");
        return {"constant_heating", SourceKind::Thermal,
            [heating_rate](const hd_core::MeshBlock& mesh, double time,
                           const config::SimulationConfig& cfg,
                           std::vector<hd_core::Conserved>& rhs, double& dt_limit) {
                constant_heating(mesh, time, cfg, heating_rate, rhs, dt_limit);
            }};
    }

    void constant_heating(
        const hd_core::MeshBlock& mesh,
        double /*time*/,
        const config::SimulationConfig& /*cfg*/,
        double heating_rate,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    )
    {
        performance::Scope timer(performance::Section::Heating);
        if (!std::isfinite(heating_rate) || heating_rate < 0.0)
            throw std::invalid_argument("heating_rate must be finite and nonnegative");
        if (heating_rate == 0.0) return;

        // Fixed safety factor: the only source-specific parameter is heating_rate.
        constexpr double timescale_fraction = 0.1;
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            const auto index = mesh.coordinates_active(flat);
            bool active = true;
            for (int d = 0; d < hd_core::NDIM; ++d)
                active = active && index[d] >= 0 && index[d] < mesh.nactive[d];
            if (!active) continue;

            const auto& cell = mesh.cells[flat];
            const auto primitive = cell.U.to_Primtivies(cell.gamma);
            const double rho = cell.U[hd_core::RHO];
            if (rho == 0.0) continue;
            const double internal_energy = primitive[hd_core::PRE] / (cell.gamma - 1.0);
            const double energy_rate = rho * heating_rate;
            rhs[flat][hd_core::ENERGY] += energy_rate;
            dt_limit = std::min(dt_limit, timescale_fraction * internal_energy / energy_rate);
        }
    }
}
