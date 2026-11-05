/**
 * @file factory.h
 * @brief Hydrodynamic spatial operators.
 */
#pragma once
#include "riemann.h"
#include <memory>
#include <string>
namespace hydro
{
    /**
     * @brief Create a Riemann solver by its configured name.
     * @param[in] name Case-sensitive solver name; currently only exact is supported.
     * @param[in] tolerance Finite positive relative pressure tolerance.
     * @param[in] max_iterations Positive maximum number of pressure iterations.
     * @return Exclusive ownership of the newly created solver.
     * @throws std::invalid_argument If the name or solver settings are invalid.
     */
    std::unique_ptr<hydro::RiemannSolver> make_riemann_solver(    const std::string& name,double tolerance,int max_iterations);
}