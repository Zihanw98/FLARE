#include "boundary.h"
#include "utils.h"
#include <stdexcept>      // std::runtime_error
namespace boundary
{
    BoundaryCondition parse_bc(const std::string& bc)
    {
        if(utils::iequals(bc,"periodic"))
        {
            return BoundaryCondition::Periodic;
        }
        if(utils::iequals(bc,"symmetric"))
        {
            return BoundaryCondition::Symmetric;
        }
        if(utils::iequals(bc,"antisymmetric"))
        {
            return BoundaryCondition::Antisymmetric;
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
    void ApplyBoundaryCondition(const hd_core::MeshBlock& mesh)
    {

    }
}