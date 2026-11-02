#include "source.h"

namespace hd_source
{
    SourceList make_sources(
        const config::SimulationConfig& cfg,
        const hd_user::UserHooks& hooks
    )
    {
        SourceList sources;
        if(cfg.constant_gravity_enabled)
        {
            sources.push_back({"constant_gravity",SourceKind::Force,constant_gravity});
        }
        return sources;
    }
}