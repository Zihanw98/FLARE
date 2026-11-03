#pragma once
#include <string>

namespace hd_core { struct MeshBlock; }
namespace hd_user { struct UserHooks; }
namespace boundary
{
    enum class BoundarySide
    {
        Upper,
        Lower
    };
    enum class BoundaryCondition
    {
        Periodic,
        Outflow,
        Reflective,
        Cont,
        User
    };
    BoundaryCondition parse_bc(const std::string& bc);
    /**
     * @brief Fill the boundary in specified side of specified dimension with given boundary condition (which are set in config)
     * @param[in,out] mesh Copy the entire mesh as input, and modify it.
     * @param[in] direction Dimension of the boundary, should be in [0,NDIM)
     * @param[in] side Side of the boundary, lower or upper
     * @param[in] time Physical time. Some boundary conditions may require it.
     * @param[in] hooks User hooks, in this case it's for user boundary
     * @note Currently the built-in BC are not completed. Only Periodic and Reflective are finished.
     * @throws std::invalid_argument direction, side or time invalid
     * @throws std::runtime_error User hook not registered, or BC that are not implemented.
     */
    void ApplyBoundaryCondition(
        hd_core::MeshBlock& mesh,
        int direction,
        BoundarySide side,
        double time,
        const hd_user::UserHooks& hooks
    );
    void FillBoundary(
        hd_core::MeshBlock& mesh,
        double time,
        const hd_user::UserHooks& hooks
    );
}
