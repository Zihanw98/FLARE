#pragma once
#include "core.h"
#include "mesh.h"
#include "boundary.h"
#include <array>
#include <cstddef>
#include <string>
#include <unordered_set>
namespace config
{
    struct SimulationConfig
    {
        hd_core::Index nactive{};
        hd_core::Index nghost{};

        //Currenly, we use an uniform grid
        //Next step, we will try to implement different grid, including user-defined

        std::array<double,hd_core::NDIM> lower{};
        std::array<double,hd_core::NDIM> upper{};
        double t_end = 1.0;
        double cfl = 0.4;
        std::size_t maxsteps = 100000;

        //Numerical method settings
        std::string riemann_solver = "exact";
        std::string spatial_scheme = "godunov";
        std::string time_integrator = "euler";

        double riemann_tolerance = 1e-8;
        int riemann_max_iterations = 100;

        //Default gas adiabatic index
        double gamma=5.0/3.0;
        
        //Boundary conditions, by default it's periodic
        std::array<std::string, hd_core::NDIM> boundary_lower{};
        std::array<std::string, hd_core::NDIM> boundary_upper{};

        //Output
        std::string output_directory = "output";
        std::string output_suffix = "test";
        // 每隔多少个完成的时间步输出一次；初始、第一步和末步额外输出。
        std::size_t step_every = 10;

    };
    SimulationConfig ReadParameters(const std::string& filename);
    std::string trim(const std::string& text);
    template<typename T>T parse_value(const std::string& text);
    template<> std::string parse_value<std::string>(const std::string& text);
    void validate_cfg(const SimulationConfig& cfg,const std::unordered_set<std::string>& seen);
    hd_core::MeshBlock build_mesh(const SimulationConfig& cfg);
}
