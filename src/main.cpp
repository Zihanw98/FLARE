#include "config/config.h"
#include "hydro/riemann/factory.h"
#include "hydro/spatial.h"
#include "time/integrator.h"
#include "user/user_api.h"
#include "boundary/boundary.h"
#include "source/source.h"
#include "io/output.h"
#include "config/read_mesa_profile.h"
#include "utils/performance.h"
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
/**
 * @brief The main entry of this code
 * @details What we do in this main code? This gives an outline of how the code works. 
 * 1.First of all, we read the parameter file by the -par argument in argv  
 * If no argv then the code will read a default parameters.txt
 * 2.Build the mesh we use in the simulation, and the riemann solver(if applicable---Right now we use exact riemann solver and godunov scheme, so it's required.)  
 * 
 */
int main(int argc, char* argv[])
{
    int exit_code = 0;
    std::string performance_file;
    try
    {
        std::string parameter_file;
        auto performance_mode = performance::Mode::Off;

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
            else if (arg == "--performance")
            {
                if (performance_mode == performance::Mode::Off)
                    performance_mode = performance::Mode::Summary;
            }
            else if (arg == "--performance-detail")
            {
                performance_mode = performance::Mode::Detail;
            }
            else if (arg == "--performance-file")
            {
                if (i + 1 >= argc)
                    throw std::runtime_error("--performance-file requires a JSON file path");
                if (!performance_file.empty())
                    throw std::runtime_error("Repeated --performance-file option");
                performance_file = argv[++i];
                if (performance_file.empty())
                    throw std::runtime_error("--performance-file requires a nonempty path");
                if (performance_mode == performance::Mode::Off)
                    performance_mode = performance::Mode::Summary;
            }
            else if (arg == "--help" || arg == "-h")
            {
                std::cout << "Usage: flare [-par parameters.txt] [--performance] "
                             "[--performance-detail] [--performance-file report.json]\n"
                             "  --performance         Report serial wall-clock timings to stderr.\n"
                             "  --performance-detail  Also time each face reconstruction/Riemann call (higher overhead).\n"
                             "  --performance-file    Also write JSON; implies --performance.\n";
                return 0;
            }
            else
            {
                throw std::runtime_error("Unknown argument: " + arg);
            }
        }

        performance::configure(performance_mode);
        performance::Scope run_timer(performance::Section::Run);
        performance::Scope initialization_timer(performance::Section::Initialization);
        if (parameter_file.empty())
        {
            std::cout<<"Parameter File -par <parameterfile> not given, using default parameters.txt instead\n";
            parameter_file = "parameters.txt";
        }
        auto cfg_temp = config::ReadParameters(parameter_file);
        std::cout << "Loaded parameters from: "
                  << parameter_file << '\n';
        auto hooks = hd_user::make_user_hooks();
        if (hooks.initialize_units)
            hooks.initialize_units(physcon::units, cfg_temp);
        hd_core::MeshBlock mesh;
        if(cfg_temp.init_from_mesa_profile)
        {
            config::MesaProfile mesaprofile = config::read_mesa_profile(cfg_temp.init_mesa_profile_name,cfg_temp.mesa_radius_inner,cfg_temp.mesa_radius_outer,"linear",100,100);
            mesh = config::build_mesh_from_mesa_profile(mesaprofile,cfg_temp,true);
        }
        else
        {
            mesh.build_mesh(cfg_temp);
        }

        const auto cfg = cfg_temp;
        auto solver = hydro::make_riemann_solver(cfg.riemann_solver,cfg.riemann_tolerance,cfg.riemann_max_iterations);

        //In the future, if the simulation doesn't require a riemann solver, use a if/switch branch to avoid creating a default solver

        auto spatial = hydro::parse_spatial_scheme("Godunov");

        auto integrator = hd_time::parse_time_integrator(cfg.time_integrator);

        if(!hooks.initialize)
        {
            throw std::runtime_error("User initializer is not registered");
        }
        hooks.initialize(mesh,cfg);
        
        double t_current = 0;
        std::size_t step_current = 0;

        boundary::FillBoundary(mesh,t_current,hooks);
        io::Output(mesh, cfg, hooks, step_current, t_current, 0.0);
        auto source_list = hd_source::make_sources(cfg,mesh,hooks);
        double dt_next = std::numeric_limits<double>::infinity();
        std::uint64_t active_cells = 1;
        for (int d = 0; d < hd_core::NDIM; ++d)
            active_cells *= mesh.nactive[d];
        initialization_timer.stop();
        performance::Scope evolution_timer(performance::Section::Evolution);
        while(t_current < cfg.t_end)
        {
            if (step_current >= cfg.maxsteps)
            throw std::runtime_error(
                "Maximum number of steps reached before t_end"
            );

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
            performance::count(performance::Counter::Steps);
            performance::count(performance::Counter::CellUpdates, active_cells);

            io::Output(mesh, cfg, hooks, step_current, t_current, dt);
        }

        evolution_timer.stop();
        if (performance::enabled()) std::cout.flush();
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        exit_code = 1;
    }
    try
    {
        performance::report(std::cerr, exit_code == 0);
        if (performance::enabled() && !performance_file.empty())
        {
            std::ofstream report_file(performance_file);
            if (!report_file)
                throw std::runtime_error("Cannot open performance report: " + performance_file);
            performance::write_json(report_file, exit_code == 0);
            report_file.close();
            if (!report_file)
                throw std::runtime_error("Cannot write performance report: " + performance_file);
        }
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        exit_code = 1;
    }
    return exit_code;
}
