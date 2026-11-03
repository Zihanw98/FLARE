#include "user_config.h"
#include "user/user_api.h"
#include "../hydro_case_helpers.h"
static_assert(hd_core::NDIM==2, "2D Richtmyer-Meshkov");
namespace {
void initialize(hd_core::MeshBlock& mesh,const config::SimulationConfig& cfg)
{
    test_case::no_gravity(cfg);
    constexpr double pi=3.14159265358979323846;
    const double lx=mesh.active_length(0),ly=mesh.active_length(1);
    for(int j=0;j<mesh.nactive[1];++j) for(int i=0;i<mesh.nactive[0];++i)
    {
        const hd_core::Index index{i+mesh.nghost[0],j+mesh.nghost[1]};
        const double x=mesh.cell_center(0,index)-cfg.lower[0];
        const double y=mesh.cell_center(1,index)-cfg.lower[1];
        const double interface=.35*lx+.02*ly*std::cos(2*pi*y/ly);
        hd_core::Primitives w{};
        w[hd_core::RHO]=2+std::tanh((x-interface)/(.015*ly));
        w[hd_core::PRE]=1;
        if(x<.175*lx) w=test_case::shock_downstream(2);
        mesh.at(index)=hd_core::Cell(w,cfg.gamma);
    }
}
void fill_boundary(hd_core::MeshBlock& mesh,double,int direction,boundary::BoundarySide side)
{ test_case::shock_boundary(mesh,direction,side,2); }
}
hd_user::UserHooks hd_user::make_user_hooks() { UserHooks h;h.initialize=initialize;h.userbc=fill_boundary;return h; }
