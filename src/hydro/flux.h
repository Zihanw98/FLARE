#pragma once
#include "core.h"

namespace hydro
{
    hd_core::Conserved physical_flux(
        const hd_core::Primitives& W,
        int direction,
        double gamma
    );
}