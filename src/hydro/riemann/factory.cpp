#include "factory.h"
#include "exact.h"
#include <memory>
#include <string>
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

    // 将来实现 HLLSolver、HLLCSolver 后，在这里增加分支。

    throw std::invalid_argument(
        "Unknown or unavailable Riemann solver: " + name
    );
}