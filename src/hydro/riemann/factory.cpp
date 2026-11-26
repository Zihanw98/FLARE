/**
 * @file factory.cpp
 * @brief Hydrodynamic spatial operators.
 */
#include "factory.h"
#include "exact.h"
#include "hll.h"
#include <memory>
#include <string>
namespace hydro
{
    std::unique_ptr<hydro::RiemannSolver> make_riemann_solver(
        const std::string& name,
        double tolerance,
        int max_iterations
    )
    {
        if (name == "exact")
        {
            return std::make_unique<hydro::ExactSolver>(
                tolerance, max_iterations
            );
        }

        if (name == "exact_bisection")
            return std::make_unique<hydro::ExactSolver>(
                tolerance, max_iterations, hydro::ExactPressureMethod::Bisection);

        if (name == "hll")
            return std::make_unique<hydro::HLLSolver>();

        throw std::invalid_argument(
            "Unknown or unavailable Riemann solver: " + name
        );
    }
    }
