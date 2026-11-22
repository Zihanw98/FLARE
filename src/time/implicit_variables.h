#pragma once
#include "core/mesh.h"
#include "implicit_linear.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace hd_implicit::detail
{
// Local Newton coordinates: d(log rho), dv / velocity_scale, d(log(E-K)).
// E-K is the non-kinetic energy density, including ionization/radiation when used.
// Cold cells use a linear non-kinetic-energy increment in the ENERGY slot:
// log(E-K) cannot represent E-K < 0, which dual energy deliberately permits.
// The independent INTERNAL_ENERGY and composition coordinates remain linear;
// EoS validation of trial states enforces a positive thermal-energy budget.
// The equations remain conservative. Row scaling is held fixed throughout
// each Jacobian construction, linear solve and its line search.
class LogVariables
{
    struct CellScale
    {
        hd_core::Primitives primitive;
        hd_core::Conserved residual;
        double velocity;
    };
    std::vector<CellScale> cells_;

public:
    // Both states must have thermodynamics refreshed for their current active U.
    LogVariables(const hd_core::MeshBlock& old, const hd_core::MeshBlock& trial,
                 const std::vector<std::size_t>& active)
    {
        cells_.reserve(active.size());
        for (auto flat : active)
        {
            const auto& cell = trial.cells[flat];
            const auto w = cell.Get_Primitives_From_Updated_State();
            const auto previous = old.cells[flat].Get_Primitives_From_Updated_State();
            if (!(w[hd_core::RHO] > 0 && w[hd_core::PRE] > 0
                  && previous[hd_core::RHO] > 0 && previous[hd_core::PRE] > 0))
                throw std::invalid_argument("Log implicit variables require positive density and pressure, not vacuum");
            CellScale scale{};
            scale.primitive = w;
            scale.velocity = std::max(cell.get_sound_speed(),
                old.cells[flat].get_sound_speed());
            for (int d = 0; d < hd_core::NDIM; ++d)
                scale.velocity = std::max({scale.velocity, std::abs(w[hd_core::VEL(d)]),
                                          std::abs(previous[hd_core::VEL(d)])});
            const double rho = std::max(w[hd_core::RHO], previous[hd_core::RHO]);
            scale.residual[hd_core::RHO] = rho;
            for (int d = 0; d < hd_core::NDIM; ++d)
                scale.residual[hd_core::MOMDEN(d)] = rho*scale.velocity;
            scale.residual[hd_core::ENERGY] = std::max(cell.U[hd_core::ENERGY], old.cells[flat].U[hd_core::ENERGY]);
            scale.residual[hd_core::INTERNAL_ENERGY] = std::max(
                cell.U.internal_energy(), old.cells[flat].U.internal_energy());
            // Additional conserved fields stay linear, with a density-based
            // scale; never drop them through the hydro-only primitive conversion.
            for (int k = hd_core::NHYDRO; k < hd_core::NVAR; ++k)
                scale.residual[k] = std::max({rho, std::abs(cell.U[k]), std::abs(old.cells[flat].U[k])});
            for (double value : scale.residual.data)
                if (!std::isfinite(value) || value <= 0)
                    throw std::invalid_argument("Invalid log implicit variable scale");
            cells_.push_back(scale);
        }
    }

    double row_scale(std::size_t row) const
    {
        return cells_[row/hd_core::NVAR].residual[row%hd_core::NVAR];
    }

    double norm(const Vector& residual) const
    {
        double worst = 0;
        for (std::size_t row = 0; row < residual.size(); ++row)
        {
            const double value = std::abs(residual[row]/row_scale(row));
            if (!std::isfinite(value)) return std::numeric_limits<double>::infinity();
            worst = std::max(worst, value);
        }
        return worst;
    }

    hd_core::Conserved update(std::size_t cell, const hd_core::Cell& base,
                              const double* dq, double alpha) const
    {
        const auto& scale = cells_[cell];
        const auto& w = scale.primitive;
        const double drho = w[hd_core::RHO]*std::expm1(alpha*dq[hd_core::RHO]);
        const double rho = w[hd_core::RHO]*std::exp(alpha*dq[hd_core::RHO]);
        const bool hot = base.U.total_energy_reliable();
        const double internal = base.U.internal_energy();
        const double updated_internal = hot ? internal*std::exp(alpha*dq[hd_core::ENERGY]) : internal;
        if (!std::isfinite(rho) || !std::isfinite(updated_internal) || rho <= 0 || updated_internal <= 0)
            throw std::invalid_argument("Log implicit update overflow/underflow");
        auto result = base.U;
        result[hd_core::RHO] = rho;
        // Increment form preserves the exact base state for dq=0, avoiding
        // primitive/conserved roundtrip noise in small finite differences.
        double energy_change = hot
            ? internal*std::expm1(alpha*dq[hd_core::ENERGY])
            : alpha*scale.residual[hd_core::ENERGY]*dq[hd_core::ENERGY];
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            const double velocity = w[hd_core::VEL(d)];
            const double dv = alpha*scale.velocity*dq[hd_core::VEL(d)];
            result[hd_core::MOMDEN(d)] += drho*velocity + rho*dv;
            energy_change += 0.5*drho*velocity*velocity + rho*velocity*dv + 0.5*rho*dv*dv;
        }
        result[hd_core::ENERGY] += energy_change;
        result[hd_core::INTERNAL_ENERGY] += alpha*scale.residual[hd_core::INTERNAL_ENERGY]*dq[hd_core::INTERNAL_ENERGY];
        for (int k = hd_core::NHYDRO; k < hd_core::NVAR; ++k)
            result[k] += alpha*scale.residual[k]*dq[k];
        result.validate();
        return result;
    }

    // Analytic dU/dq action for the nonlocal source Jacobian (including gravity).
    hd_core::Conserved direction(std::size_t cell, const hd_core::Cell& base,
                                 const double* dq) const
    {
        const auto& scale = cells_[cell];
        const auto& w = scale.primitive;
        hd_core::Conserved du{};
        du[hd_core::RHO] = w[hd_core::RHO]*dq[hd_core::RHO];
        du[hd_core::INTERNAL_ENERGY] = scale.residual[hd_core::INTERNAL_ENERGY]*dq[hd_core::INTERNAL_ENERGY];
        du[hd_core::ENERGY] = base.U.total_energy_reliable()
            ? base.U.internal_energy()*dq[hd_core::ENERGY]
            : scale.residual[hd_core::ENERGY]*dq[hd_core::ENERGY];
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            const double velocity = w[hd_core::VEL(d)];
            const double dv = scale.velocity*dq[hd_core::VEL(d)];
            du[hd_core::MOMDEN(d)] = velocity*du[hd_core::RHO] + w[hd_core::RHO]*dv;
            du[hd_core::ENERGY] += 0.5*velocity*velocity*du[hd_core::RHO] + w[hd_core::RHO]*velocity*dv;
        }
        for (int k = hd_core::NHYDRO; k < hd_core::NVAR; ++k)
            du[k] = scale.residual[k]*dq[k];
        return du;
    }
};
}
