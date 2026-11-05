/**
 * @file boundary.h
 * @brief Boundary-condition handling.
 */
#pragma once
#include <string>

namespace hd_core { struct MeshBlock; }
namespace hd_user { struct UserHooks; }
namespace boundary
{
    /** @brief Lower or upper side of a coordinate axis. */
    enum class BoundarySide
    {
        Upper, ///< Upper side of an axis.
        Lower ///< Lower side of an axis.
    };
    /** @brief Boundary identifiers; not every parsed identifier is implemented. */
    enum class BoundaryCondition
    {
        Periodic, ///< Wrap ghost states to the opposite active boundary.
        Outflow, ///< Reserved outflow boundary; application is not implemented.
        Reflective, ///< Mirror states and reverse the normal momentum.
        Cont, ///< Reserved continuous boundary; application is not implemented.
        User ///< Delegate boundary filling to the registered user callback.
    };
    /**
     * @brief Parse a case-insensitive boundary-condition name.
     * @param[in] bc Name: periodic, reflective, cont, outflow, or user.
     * @return The corresponding boundary-condition enum.
     * @throws std::runtime_error If the name is unknown.
     * @note Parsing cont or outflow succeeds, but applying either is not implemented.
     */
    BoundaryCondition parse_bc(const std::string& bc);
    /**
     * @brief Fill ghost cells on one side of a mesh axis.
     * @param[in,out] mesh Allocated mesh whose selected ghost cells are updated in place.
     * @param[in] direction Axis index in [0, hd_core::NDIM).
     * @param[in] side Lower or upper side of the axis.
     * @param[in] time Finite physical time passed to a user boundary callback.
     * @param[in] hooks Callbacks; userbc must be registered for a user boundary.
     * @pre Active states and geometry are initialized; periodic boundaries are paired on the same axis.
     * @throws std::invalid_argument If direction, side, time, allocation, or boundary constraints are invalid.
     * @throws std::runtime_error If the callback is missing or the boundary type is not implemented.
     * @note Built-in periodic and reflective boundaries copy the full Cell; reflection reverses normal momentum. Reflective boundaries require nghost <= nactive.
     */
    void ApplyBoundaryCondition(
        hd_core::MeshBlock& mesh,
        int direction,
        BoundarySide side,
        double time,
        const hd_user::UserHooks& hooks
    );
    /**
     * @brief Fill both boundaries of every mesh axis in axis order.
     * @param[in,out] mesh Allocated mesh whose ghost states are filled.
     * @param[in] time Finite physical time of the boundary update.
     * @param[in] hooks User boundary callbacks.
     * @pre Active cell states are initialized and boundary conditions are configured.
     * @note Delegates to ApplyBoundaryCondition(); its errors and callback exceptions propagate.
     */
    void FillBoundary(
        hd_core::MeshBlock& mesh,
        double time,
        const hd_user::UserHooks& hooks
    );
}
