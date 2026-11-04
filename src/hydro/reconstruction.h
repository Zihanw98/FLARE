#pragma once
#include "core.h"
#include "mesh.h"
#include "limiter.h"

namespace hydro
{
    enum struct Reconstructor
    {
        PCM,
        PLM,
    };
    Reconstructor parse_reconstructor(const std::string& name);
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct(const Reconstructor& reconstruct, const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter,double limiter_beta=1.5);
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PCM(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter,double limiter_beta=1.5);
    std::pair<hd_core::Primitives,hd_core::Primitives> Reconstruct_PLM(const hd_core::MeshBlock& mesh, int direction, const hd_core::Index& right_index, const Limiter& limiter,double limiter_beta=1.5);
}