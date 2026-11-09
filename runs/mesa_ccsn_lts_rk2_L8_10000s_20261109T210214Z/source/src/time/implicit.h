#pragma once
#include "time/integrator.h"
#include <memory>

// Backward Euler workspace used by the common time-integrator dispatcher.
namespace hd_implicit
{
class Integrator
{
public:
    Integrator(const hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg);
    ~Integrator();
    Integrator(const Integrator&) = delete;
    Integrator& operator=(const Integrator&) = delete;
    double suggested_dt() const;
    // Commits only a converged, boundary-filled state. dt <= remaining.
    double advance(hd_core::MeshBlock& mesh, double time, double remaining,
                   hydro::SpatialScheme scheme, const hydro::RiemannSolver& solver,
                   const hd_source::SourceList& sources, const config::SimulationConfig& cfg,
                   const hd_user::UserHooks& hooks);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
}
