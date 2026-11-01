#include "config/config.h"
#include "hydro/riemann/factory.h"
#include "hydro/spatial.h"
#include "time/integrator.h"
#include "user/user_api.h"
#include "boundary/boundary.h"
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

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

        auto integrator = hd_time::parse_time_integrator("euler");

        auto hooks = hd_user::make_user_hooks();
        if(!hooks.initialize)
        {
            throw std::runtime_error("User initializer is not registered");
        }
        hooks.initialize(mesh,cfg);
        
        double t_current = 0;
        int step_current = 0;
        while(t_current < cfg.t_end)
        {
            if (step_current >= cfg.maxsteps)
            throw std::runtime_error(
                "Maximum number of steps reached before t_end"
            );
            boundary::ApplyBoundaryCondition(mesh);
            //Apply BC before every step
            const double remaining = cfg.t_end - t_current;
            //Avoid too small remaining time
            const double dt = hd_time::advance_one_step(
                mesh,
                t_current,
                remaining,
                integrator,
                spatial,
                *solver,
                cfg.cfl
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

            std::cout << "step = " << step_current<< ", time = " << t_current<< ", dt = " << dt << '\n';
        }

    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}