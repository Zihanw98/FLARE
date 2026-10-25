#pragma once
#include "core.h"
#include "cell.h"
#include <array>
#include <vector>
namespace hd_core
{
    struct MeshBlock
    {
        std::array<int,NDIM> nactive{};
        std::array<int,NDIM> nghost{};
        
        int total(int d) const;
        
        std::array<double,NDIM> dx{};
        std::array<double,NDIM> length{};

        double spacing(int d) const;
        
        void allocate();

        std::vector<Cell> cells;

    };
}
