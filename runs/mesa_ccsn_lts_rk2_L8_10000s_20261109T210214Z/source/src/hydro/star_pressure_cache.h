/**
 * @file star_pressure_cache.h
 * @brief Caller-owned star pressures used as initial guesses for face solves.
 */
#pragma once
#include "core.h"
#include <array>
#include <cstddef>
#include <vector>

namespace hydro
{
    /** One value per directional face, with the same flat indexing as FaceFlux.
     * Zero means no usable guess. Copies are independent trial snapshots;
     * integrators decide when to copy, restore and commit them with swap().
     * This storage does not track mesh geometry or schedule face updates.
     */
    class StarPressureCache
    {
    public:
        std::vector<double>& operator[](std::size_t direction) { return faces_[direction]; }
        const std::vector<double>& operator[](std::size_t direction) const { return faces_[direction]; }

        /** Preserve existing guesses only when the directional face count matches. */
        void ensure_size(std::size_t direction, std::size_t face_count)
        {
            if (faces_[direction].size() != face_count)
                reset(direction, face_count);
        }

        /** Invalidate a direction, including geometry changes with the same face count. */
        void reset(std::size_t direction, std::size_t face_count)
        {
            faces_[direction].assign(face_count, 0.0);
        }

        /** Discard every direction before reusing the cache for a new mesh or run. */
        void clear() noexcept
        {
            for (auto& faces : faces_) faces.clear();
        }

        void swap(StarPressureCache& other) noexcept { faces_.swap(other.faces_); }

        bool operator==(const StarPressureCache& other) const { return faces_ == other.faces_; }

    private:
        std::array<std::vector<double>, hd_core::NDIM> faces_;
    };
}
