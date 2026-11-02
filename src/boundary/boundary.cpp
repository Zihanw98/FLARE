#include "boundary.h"
#include "core/mesh.h"
#include "user/user_api.h"
#include "utils.h"
#include <cmath>
#include <stdexcept>      // std::runtime_error
namespace boundary
{
    BoundaryCondition parse_bc(const std::string& bc)
    {
        if(utils::iequals(bc,"periodic"))
        {
            return BoundaryCondition::Periodic;
        }
        if(utils::iequals(bc,"reflective"))
        {
            return BoundaryCondition::Reflective;
        }
        if(utils::iequals(bc,"cont"))
        {
            return BoundaryCondition::Cont;
        }
        if(utils::iequals(bc,"outflow"))
        {
            return BoundaryCondition::Outflow;
        }
        if(utils::iequals(bc,"user"))
        {
            return BoundaryCondition::User;
        }
        throw std::runtime_error("Not a valid Boundary Condition, currently support BC: periodic, symmetric, antisymmetric, outflow, user");
    }
    void ApplyBoundaryCondition(
        hd_core::MeshBlock& mesh,
        int direction,
        BoundarySide side,
        double time,
        const hd_user::UserHooks& hooks
    )
    {
        if (direction < 0 || direction >= hd_core::NDIM)
            throw std::invalid_argument("Invalid boundary direction");
        if (side != BoundarySide::Lower && side != BoundarySide::Upper)
            throw std::invalid_argument("Invalid boundary side");
        if (!std::isfinite(time))
            throw std::invalid_argument("Boundary time must be finite");
        if (mesh.cells.empty())
            throw std::invalid_argument("Boundary conditions require an allocated mesh");

        if (hooks.userbc)
        {
            hooks.userbc(mesh, time, direction, side);
            return;
        }

        const auto condition = side == BoundarySide::Lower
            ? mesh.bc_lower[direction] : mesh.bc_upper[direction];
        if (condition == BoundaryCondition::User)
            throw std::runtime_error("User boundary selected, but userbc is not registered");
        if (condition != BoundaryCondition::Periodic)
            throw std::runtime_error(
                "Only periodic built-in boundaries are implemented"
            );
        if (mesh.bc_lower[direction] != BoundaryCondition::Periodic
            || mesh.bc_upper[direction] != BoundaryCondition::Periodic)
            throw std::invalid_argument("Periodic boundaries must be paired on the same axis");

        mesh.ntotal(direction); // 检查该方向的 nactive/nghost。
        const int count = mesh.nactive[direction];
        const int first = mesh.nghost[direction];
        const int end = first + count;
        for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
        {
            auto index = mesh.coordinates(flat);
            const int coordinate = index[direction];
            const bool on_side = side == BoundarySide::Lower
                ? coordinate < first : coordinate >= end;
            if (!on_side)
                continue;

            // 包裹到同方向的活跃区，允许 ghost 层数超过活跃单元数。
            int wrapped = (coordinate - first) % count;
            if (wrapped < 0)
                wrapped += count;
            index[direction] = first + wrapped;
            mesh.cells[flat] = mesh.at(index);
        }
    }
}
