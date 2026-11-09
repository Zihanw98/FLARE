/**
 * @file facefield.h
 * @brief Conserved-flux specialization of the generic face field.
 */
#pragma once
#include "../core/facefield.h"

namespace hydro
{
    // Flux values use Conserved; geometry areas use FaceField<double> directly.
    // Godunov allocates nactive with one additional face in the normal direction.
    using FaceFlux = hd_core::FaceField<hd_core::Conserved>;
}
