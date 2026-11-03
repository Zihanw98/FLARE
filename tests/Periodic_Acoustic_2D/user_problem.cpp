#include "user_config.h"
#include "user/user_api.h"
#include "../hydro_case_helpers.h"
static_assert(hd_core::NDIM==2, "2D acoustic wave");
namespace {
void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    test_case::no_gravity(cfg);
    constexpr double pi=3.14159265358979323846, amplitude=1e-4;
    const double kx=2*pi/mesh.active_length(0), ky=2*pi/mesh.active_length(1);
    const double kmag=std::hypot(kx,ky);
    for (int j=0;j<mesh.nactive[1];++j) for(int i=0;i<mesh.nactive[0];++i)
    {
        const hd_core::Index index{i+mesh.nghost[0],j+mesh.nghost[1]};
        const double phase=kx*(mesh.cell_center(0,index)-cfg.lower[0])+ky*(mesh.cell_center(1,index)-cfg.lower[1]);
        const double ax=.5*kx*mesh.cell_width(0,index), ay=.5*ky*mesh.cell_width(1,index);
        const double delta=amplitude*std::sin(phase)*std::sin(ax)/ax*std::sin(ay)/ay;
        hd_core::Primitives w{};
        w[hd_core::RHO]=1+delta;
        w[hd_core::PRE]=1/cfg.gamma+delta;
        w[hd_core::VEL(0)]=delta*kx/kmag;
        w[hd_core::VEL(1)]=delta*ky/kmag;
        mesh.at(index)=hd_core::Cell(w,cfg.gamma);
    }
}
}
hd_user::UserHooks hd_user::make_user_hooks() { UserHooks h; h.initialize=initialize; return h; }
