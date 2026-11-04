#pragma once
#include <string>
#include "mesh.h"
#include "config.h"
namespace hd_source
{
    // 将未乘 dt 的源项加到 rhs；dt_limit 只能收紧，不能增大。
    // RK2 以 +infinity 初始化每个源项的独立限制；rhs 不应依赖 dt_limit。
    using SourceFunction = void (*)(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    );
    enum class SourceKind
    {
        Thermal,
        Force,
        General
    };

    struct SourceTerm
    {
        std::string name;
        SourceKind kind;
        SourceFunction func;
    };

    using SourceList = std::vector<SourceTerm>;
    SourceList make_sources(
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks
    );

    //Built-in Source Function
    void constant_gravity(
        const hd_core::MeshBlock& mesh,
        double time,
        const config::SimulationConfig& cfg,
        std::vector<hd_core::Conserved>& rhs,
        double& dt_limit
    );
}
