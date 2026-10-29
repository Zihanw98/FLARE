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
        
        int ntotal(int d) const;
        
        // Tensor-product coordinates, including ghost faces.
        // face[d].size() == total(d) + 1.
        // Active faces run from nghost[d] through nghost[d] + nactive[d].
        std::array<std::vector<double>,NDIM> face;

        // Set nactive/nghost first. Bounds describe the ACTIVE region;
        // ghost faces are extended with the same uniform spacing.
        void set_uniform_faces(int d, double lower, double upper);

        // i is the storage index along direction d, including ghost cells.
        double cell_width(int d, int i) const;
        double cell_center(int d, int i) const; // Coordinate midpoint, not volume centroid.
        double active_length(int d) const;

        void validate_geometry() const;
        // Validates geometry, then allocates cells. Does not set physical states.
        // After allocation, do not change sizes/coordinates without reinitializing
        // or conservatively remapping the physical data.
        void allocate();

        std::vector<Cell> cells;
 
    };
}
