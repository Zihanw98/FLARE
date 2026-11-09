#include "config/config.h"
#include "config/read_mesa_profile.h"
#include "core/mesh.h"
#include "hydro/riemann/factory.h"
#include "time/integrator.h"
#include "user/user_api.h"
#include "boundary/boundary.h"
#include "physcon/constants.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

struct Snapshot { double time=0,dt=0; int step=0; };
Snapshot read(const std::filesystem::path& path, hd_core::MeshBlock& mesh) {
    std::ifstream input(path); std::string line; Snapshot meta; std::size_t flat=0;
    if(!input) throw std::runtime_error("Missing snapshot");
    while(std::getline(input,line)) {
        if(line.rfind("# time=",0)==0) meta.time=std::stod(line.substr(7));
        else if(line.rfind("# dt=",0)==0) meta.dt=std::stod(line.substr(5));
        else if(line.rfind("# step=",0)==0) meta.step=std::stoi(line.substr(7));
        else if(line.empty() || line[0]=='#' || line[0]=='x') continue;
        else {
            std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);
            double x,y,z,g,r,vx,vy,vz,p;int ghost;
            if(!(row>>x>>y>>z>>g>>r>>vx>>vy>>vz>>p>>ghost) || flat>=mesh.cells.size())
                throw std::runtime_error("Invalid snapshot row");
            hd_core::Primitives w{};w[hd_core::RHO]=r;w[hd_core::PRE]=p;
            const double velocity[]={vx,vy,vz};
            for(int d=0;d<hd_core::NDIM;++d) w[hd_core::VEL(d)]=velocity[d];
            mesh.cells[flat++]=hd_core::Cell(w,g);
        }
    }
    if(flat!=mesh.cells.size()) throw std::runtime_error("Wrong snapshot cell count");
    mesh.validate_active();return meta;
}
int main(int argc,char** argv) {
    try {
        if(argc!=2 && !(argc==3 && std::string(argv[2])=="--allow-partial")) throw std::runtime_error("Expected parameter path");
        auto cfg=config::ReadParameters(argv[1]);auto hooks=hd_user::make_user_hooks();
        if(hooks.initialize_units) hooks.initialize_units(physcon::units,cfg);
        auto profile = config::read_mesa_profile(cfg.init_mesa_profile_name,cfg.mesa_radius_inner,cfg.mesa_radius_outer,"linear",100,100);
        auto mesh = config::build_mesh_from_mesa_profile(profile,cfg,true);
        auto sources=hd_source::make_sources(cfg,mesh,hooks);
        auto solver=hydro::make_riemann_solver(cfg.riemann_solver,cfg.riemann_tolerance,cfg.riemann_max_iterations);
        std::vector<std::filesystem::path> paths;
        for(const auto& p:std::filesystem::directory_iterator(cfg.output_directory))
            if(p.path().extension()==".csv" && p.path().filename().string().rfind(cfg.output_suffix+"_",0)==0) paths.push_back(p.path());
        std::sort(paths.begin(),paths.end());if(paths.size()<2)throw std::runtime_error("Need consecutive snapshots");
        auto previous=read(paths.front(),mesh);double worst=0;int checked=0;
        for(std::size_t i=1;i<paths.size();++i) {
            auto old=mesh;auto current=read(paths[i],mesh);
            if(current.step!=previous.step+1 || std::abs(current.time-previous.time-current.dt)>1e-13
               )
                throw std::runtime_error("Nonconsecutive snapshots");
            boundary::FillBoundary(mesh,current.time,hooks);double hydro_limit,source_limit;
            const auto result=hd_time::evaluate_RHS(mesh,current.time,hydro::SpatialScheme::Godunov,sources,*solver,cfg,cfg.cfl,hydro_limit,source_limit);
            for(std::size_t f=0;f<mesh.cells.size();++f) if(mesh.is_active(f)) for(int k=0;k<hd_core::NVAR;++k) {
                const auto w=mesh.cells[f].U.to_Primtivies(mesh.cells[f].gamma);
                const auto wold=old.cells[f].U.to_Primtivies(old.cells[f].gamma);
                const double rho=std::max(w[hd_core::RHO],wold[hd_core::RHO]);
                const double speed=std::max({std::sqrt(mesh.cells[f].gamma*w[hd_core::PRE]/w[hd_core::RHO]),
                    std::sqrt(old.cells[f].gamma*wold[hd_core::PRE]/wold[hd_core::RHO]),
                    std::abs(w[hd_core::VEL(0)]),std::abs(wold[hd_core::VEL(0)])});
                const double residual=mesh.cells[f].U[k]-old.cells[f].U[k]-current.dt*result.second[f][k];
                const double scale=1e-6*(k==hd_core::RHO ? rho : k==hd_core::ENERGY
                    ? std::max(mesh.cells[f].U[k],old.cells[f].U[k]) : rho*speed);
                worst=std::max(worst,std::abs(residual)/scale);
            }
            previous=current;++checked;
        }
        if((argc!=3 && std::abs(previous.time-cfg.t_end)>1e-13) || worst>1.01)
            throw std::runtime_error("End time or conservative residual check failed");
        std::cout<<std::setprecision(17)<<"{\"steps_checked\":"<<checked<<",\"max_scaled_residual\":"<<worst
                 <<",\"final_time\":"<<previous.time<<",\"reached_target\":"<<(std::abs(previous.time-cfg.t_end)<=1e-13 ? "true" : "false")<<",\"saved_steps_passed\":true}\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
