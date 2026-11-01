#pragma once
#include <string>
#include "mesh.h"
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
        Symmetric,
        Antisymmetric,
        User
    };
    BoundaryCondition parse_bc(const std::string& bc);
    void ApplyBoundaryCondition(const hd_core::MeshBlock);
}