#include "user_config.h"
#include "user/user_api.h"
#include "../hydro_case_helpers.h"
static_assert(hd_core::NDIM==3, "3D shock-cloud interaction");
namespace {
void initialize(hd_core::MeshBlock& mesh,const config::SimulationConfig& cfg)
{
    test_case::no_gravity(cfg);
    const double lx=mesh.active_length(0),ly=mesh.active_length(1),lz=mesh.active_length(2);
    const double center[3]={cfg.lower[0]+.4*lx,cfg.lower[1]+.5*ly,cfg.lower[2]+.5*lz};
    const double radius=.18*std::min(ly,lz);
    for(int k=0;k<mesh.nactive[2];++k) for(int j=0;j<mesh.nactive[1];++j) for(int i=0;i<mesh.nactive[0];++i)
    {
        const hd_core::Index index{i+mesh.nghost[0],j+mesh.nghost[1],k+mesh.nghost[2]};
        double r2=0;
        for(int d=0;d<3;++d) { const double v=mesh.cell_center(d,index)-center[d];r2+=v*v; }
        hd_core::Primitives w{};
        w[hd_core::RHO]=r2<radius*radius ? 5:1;
        w[hd_core::PRE]=1;
        if(mesh.cell_center(0,index)<cfg.lower[0]+.075*lx) w=test_case::shock_downstream(2);
        mesh.at(index)=hd_core::Cell(w,cfg.gamma);
    }
}
void fill_boundary(hd_core::MeshBlock& mesh,double,int direction,boundary::BoundarySide side)
{ test_case::shock_boundary(mesh,direction,side,2); }
}
hd_user::UserHooks hd_user::make_user_hooks() { UserHooks h;h.initialize=initialize;h.userbc=fill_boundary;return h; }
