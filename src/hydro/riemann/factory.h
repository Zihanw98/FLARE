#pragma once
#include "riemann.h"
#include <memory>
#include <string>
std::unique_ptr<hydro::RiemannSolver> make_riemann_solver(    const std::string& name,double tolerance,int max_iterations);