#include "config/config.h"
#include "hydro/riemann/factory.h"
#include "hydro/spatial.h"
#include "time/integrator.h"
#include "user/user_api.h"
#include "boundary/boundary.h"
#include "source/source.h"
#include "io/output.h"
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
/**
 * @brief The main entry of this code
 * @details What we do in this main code? This gives an outline of how the code works. 
 * First of all, we read the parameter file by the -par argument in argv
 * If no argv then 
 */
int main(int argc, char* argv[])
{
    try
    {
        std::string parameter_file;

        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];

            if (arg == "-par")
            {
                if (i + 1 >= argc)
                    throw std::runtime_error("-par requires a file path");

                if (!parameter_file.empty())
                    throw std::runtime_error("Repeated -par option");

                parameter_file = argv[++i];
            }
            else
            {
                throw std::runtime_error("Unknown argument: " + arg);
            }
        }

        if (parameter_file.empty())
        {
            std::cout<<"Parameter File -par <parameterfile> not given, using default parameters.txt instead\n";
            parameter_file = "parameters.txt";
        }
        const auto cfg = config::ReadParameters(parameter_file);
        std::cout << "Loaded parameters from: "
                  << parameter_file << '\n';
        auto mesh = config::build_mesh(cfg);
        
        auto solver = hydro::make_riemann_solver(cfg.riemann_solver,cfg.riemann_tolerance,cfg.riemann_max_iterations);

        //In the future, if the simulation doesn't require a riemann solver, use a if/switch branch to avoid creating a default solver

        auto spatial = hydro::parse_spatial_scheme("Godunov");

        auto integrator = hd_time::parse_time_integrator(cfg.time_integrator);

        auto hooks = hd_user::make_user_hooks();
        if(!hooks.initialize)
        {
            throw std::runtime_error("User initializer is not registered");
        }
        hooks.initialize(mesh,cfg);
        
        double t_current = 0;
        std::size_t step_current = 0;

        boundary::FillBoundary(mesh,t_current,hooks);
        io::Output(mesh, cfg, hooks, step_current, t_current, 0.0);
        auto source_list = hd_source::make_sources(cfg,hooks);
        double dt_next = std::numeric_limits<double>::infinity();
        while(t_current < cfg.t_end)
        {
            if (step_current >= cfg.maxsteps)
            throw std::runtime_error(
                "Maximum number of steps reached before t_end"
            );
            // 初始状态和上一步末尾已经填好当前时刻的边界。
            const double remaining = cfg.t_end - t_current;
            //Avoid too small remaining time
            const double dt = hd_time::advance_one_step(
                mesh,
                t_current,
                remaining,
                integrator,
                spatial,
                *solver,
                source_list,
                cfg,
                hooks,
                cfg.cfl,
                dt_next
            );
            if (!std::isfinite(dt) || dt <= 0.0 ||
                dt > remaining || t_current + dt <= t_current)
            {
                throw std::runtime_error("Invalid time step");
            }
            t_current = (dt == remaining)
                ? cfg.t_end
                : t_current + dt;

            step_current++;

            io::Output(mesh, cfg, hooks, step_current, t_current, dt);
        }

    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
