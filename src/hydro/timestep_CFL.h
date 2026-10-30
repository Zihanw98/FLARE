#pragma once

namespace hd_core
{
    struct MeshBlock;
}

namespace hydro
{
    double compute_cfl_dt(
        const hd_core::MeshBlock& mesh,
        double cfl = 0.4
    );
}