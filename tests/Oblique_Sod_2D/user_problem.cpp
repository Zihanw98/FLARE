#include "user_config.h"
#include "user/user_api.h"
#include "../hydro_case_helpers.h"
static_assert(hd_core::NDIM==2, "2D oblique Sod");
namespace {
double cx=.5,cy=.5;
hd_core::Cell sample(double x,double y,double time)
{
    const double nx=2/std::sqrt(5.),ny=1/std::sqrt(5.);
    const auto s=test_case::sod(nx*(x-cx)+ny*(y-cy),time);
    hd_core::Primitives w{};
    w[hd_core::RHO]=s.rho;w[hd_core::PRE]=s.p;
    w[hd_core::VEL(0)]=nx*s.u;w[hd_core::VEL(1)]=ny*s.u;
    return hd_core::Cell(w,1.4);
}
hd_core::Cell average(hd_core::MeshBlock& mesh,const hd_core::Index& index,double time)
{
    // Subcell quadrature averages CONSERVED variables across the oblique front.
    hd_core::Conserved sum{};
    constexpr int count=8;
    for(int b=0;b<count;++b) for(int a=0;a<count;++a)
    {
        const double x=mesh.cell_center(0,index)+((a+.5)/count-.5)*mesh.cell_width(0,index);
        const double y=mesh.cell_center(1,index)+((b+.5)/count-.5)*mesh.cell_width(1,index);
        sum+=sample(x,y,time).U;
    }
    return hd_core::Cell(sum/(count*count),1.4);
}
void initialize(hd_core::MeshBlock& mesh,const config::SimulationConfig& cfg)
{
    test_case::no_gravity(cfg);
    cx=.5*(cfg.lower[0]+cfg.upper[0]);cy=.5*(cfg.lower[1]+cfg.upper[1]);
    for(int j=0;j<mesh.nactive[1];++j) for(int i=0;i<mesh.nactive[0];++i)
    {
        const hd_core::Index index{i+mesh.nghost[0],j+mesh.nghost[1]};
        mesh.at(index)=average(mesh,index,0);
    }
}
void fill_boundary(hd_core::MeshBlock& mesh,double time,int direction,boundary::BoundarySide side)
{
    const int first=mesh.nghost[direction],end=first+mesh.nactive[direction];
    for(std::size_t flat=0;flat<mesh.cells.size();++flat)
    {
        const auto index=mesh.coordinates(flat);
        if(side==boundary::BoundarySide::Lower ? index[direction]<first : index[direction]>=end)
            mesh.cells[flat]=average(mesh,index,time);
    }
}
}
hd_user::UserHooks hd_user::make_user_hooks() { UserHooks h;h.initialize=initialize;h.userbc=fill_boundary;return h; }
