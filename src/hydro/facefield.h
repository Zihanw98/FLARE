#pragma once
#include "core.h"
#include <cstddef>
#include <initializer_list>
#include <vector>

namespace hydro
{
    struct FaceField
    {

        hd_core::Index shape{};

        std::vector<hd_core::Conserved> faces;

        void allocate(
            const hd_core::Index& dimensions
        );

        // Face-field local coordinates; axis 0 varies fastest.
        hd_core::Index coordinates(std::size_t flat) const;

        hd_core::Conserved& at(
            std::initializer_list<int> coordinates
        );

        const hd_core::Conserved& at(
            std::initializer_list<int> coordinates
        ) const;

        hd_core::Conserved& at(const hd_core::Index& coordinates);
        const hd_core::Conserved& at(const hd_core::Index& coordinates) const;

        // Convert face-field local coordinates to a faces storage index.
        std::size_t flat_index(
            std::initializer_list<int> coordinates
        ) const;

        std::size_t flat_index(const hd_core::Index& coordinates) const;
    };
}
