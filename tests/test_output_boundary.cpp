#include <unordered_set>
#include "config/config.h"
#include "core/mesh.h"
#include "boundary/boundary.h"
#include "io/output.h"
#include "user/user_api.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{
    int user_outputs = 0;

    void user_output(hd_core::MeshBlock&, const config::SimulationConfig&)
    {
        ++user_outputs;
    }

    template<class Function>
    void expect_failure(Function function)
    {
        bool failed = false;
        try { function(); }
        catch (const std::exception&) { failed = true; }
        assert(failed);
    }

    bool near(double a, double b)
    {
        return std::abs(a - b) <= 1e-11 * std::max({1.0, std::abs(a), std::abs(b)});
    }

    bool active(const hd_core::MeshBlock& mesh, const hd_core::Index& index)
    {
        for (int d = 0; d < hd_core::NDIM; ++d)
            if (index[d] < mesh.nghost[d]
                || index[d] >= mesh.nghost[d] + mesh.nactive[d])
                return false;
        return true;
    }

    hd_core::MeshBlock make_mesh()
    {
        hd_core::MeshBlock mesh{};
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            mesh.nactive[d] = 2 + d;
            mesh.nghost[d] = 3 + d; // More ghost layers than active cells.
            mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
            mesh.face[d].push_back(-1.0);
            for (int i = 0; i < mesh.ntotal(d); ++i)
                mesh.face[d].push_back(mesh.face[d].back() + 0.1 * (d + 1) * (1 + 0.05 * i));
        }
        mesh.allocate();
        for (std::size_t f = 0; f < mesh.cells.size(); ++f)
        {
            if (!active(mesh, mesh.coordinates(f))) continue;
            hd_core::Primitives W{};
            W[hd_core::RHO] = 1.0 + 0.01 * f;
            W[hd_core::PRE] = 2.0 + 0.02 * f;
            for (int d = 0; d < hd_core::NDIM; ++d)
                W[hd_core::VEL(d)] = -0.1 * (d + 1);
            mesh.cells[f] = hd_core::Cell(W, 1.4 + 0.0001 * f);
        }
        return mesh;
    }

    void periodic_tests(hd_core::MeshBlock& mesh)
    {
        const auto before = mesh.cells;
        boundary::ApplyBoundaryCondition(mesh, 0, boundary::BoundarySide::Lower, 0.0, {});
        for (std::size_t f = 0; f < mesh.cells.size(); ++f)
            if (mesh.coordinates(f)[0] >= mesh.nghost[0])
            {
                assert(mesh.cells[f].U.data == before[f].U.data);
                assert(mesh.cells[f].gamma == before[f].gamma);
            }
        for (int d = 0; d < hd_core::NDIM; ++d)
            for (const auto side : {boundary::BoundarySide::Lower, boundary::BoundarySide::Upper})
                boundary::ApplyBoundaryCondition(mesh, d, side, 0.1, {});

        for (std::size_t f = 0; f < mesh.cells.size(); ++f)
        {
            auto source = mesh.coordinates(f);
            for (int d = 0; d < hd_core::NDIM; ++d)
            {
                while (source[d] < mesh.nghost[d]) source[d] += mesh.nactive[d];
                while (source[d] >= mesh.nghost[d] + mesh.nactive[d]) source[d] -= mesh.nactive[d];
            }
            const auto& expected = before[mesh.flat_index(source)];
            assert(mesh.cells[f].U.data == expected.U.data);
            assert(mesh.cells[f].gamma == expected.gamma);
        }
        expect_failure([&] {
            boundary::ApplyBoundaryCondition(mesh, hd_core::NDIM, boundary::BoundarySide::Lower, 0.0, {});
        });
        mesh.bc_upper[0] = boundary::BoundaryCondition::Outflow;
        expect_failure([&] {
            boundary::ApplyBoundaryCondition(mesh, 0, boundary::BoundarySide::Lower, 0.0, {});
        });
        expect_failure([&] {
            boundary::ApplyBoundaryCondition(mesh, 0, boundary::BoundarySide::Upper, 0.0, {});
        });
        mesh.bc_upper[0] = boundary::BoundaryCondition::Periodic;
    }

    struct BoundaryCall
    {
        hd_core::MeshBlock* mesh;
        double time;
        int direction;
        boundary::BoundarySide side;
    };
    std::vector<BoundaryCall> boundary_calls;

    void user_boundary(hd_core::MeshBlock& mesh, double time, int direction,
                       boundary::BoundarySide side)
    {
        boundary_calls.push_back({&mesh, time, direction, side});
        auto index = mesh.nghost;
        index[direction] = side == boundary::BoundarySide::Lower
            ? mesh.nghost[direction] - 1 : mesh.nghost[direction] + mesh.nactive[direction];
        hd_core::Primitives W{};
        W[hd_core::RHO] = 90.0 + direction;
        W[hd_core::PRE] = 1.0;
        mesh.at(index) = hd_core::Cell(W, 1.7);
    }

    void hook_tests(hd_core::MeshBlock mesh)
    {
        hd_user::UserHooks hooks;
        hooks.userbc = user_boundary;
        boundary_calls.clear();
        // Registering a hook must not override explicitly selected periodic boundaries.
        auto periodic_reference = mesh;
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::Periodic;
            periodic_reference.bc_lower[d] = periodic_reference.bc_upper[d]
                = boundary::BoundaryCondition::Periodic;
            for (const auto side : {boundary::BoundarySide::Lower, boundary::BoundarySide::Upper})
            {
                boundary::ApplyBoundaryCondition(mesh, d, side, 0.25, hooks);
                boundary::ApplyBoundaryCondition(periodic_reference, d, side, 0.25, {});
                assert(boundary_calls.empty());
                for (std::size_t f = 0; f < mesh.cells.size(); ++f)
                {
                    assert(mesh.cells[f].U.data == periodic_reference.cells[f].U.data);
                    assert(mesh.cells[f].gamma == periodic_reference.cells[f].gamma);
                }
            }
        }
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            mesh.bc_lower[d] = mesh.bc_upper[d] = boundary::BoundaryCondition::User;
            for (const auto side : {boundary::BoundarySide::Lower, boundary::BoundarySide::Upper})
            {
                const auto count = boundary_calls.size();
                boundary::ApplyBoundaryCondition(mesh, d, side, 0.25, hooks);
                assert(boundary_calls.size() == count + 1);
                const auto& call = boundary_calls.back();
                assert(call.mesh == &mesh && call.time == 0.25);
                assert(call.direction == d && call.side == side);
                auto index = mesh.nghost;
                index[d] = side == boundary::BoundarySide::Lower
                    ? mesh.nghost[d] - 1 : mesh.nghost[d] + mesh.nactive[d];
                // Explicit user boundaries must preserve the hook's result.
                assert(mesh.at(index).U[hd_core::RHO] == 90.0 + d);
                assert(mesh.at(index).gamma == 1.7);
            }
        }
        mesh.bc_lower[0] = boundary::BoundaryCondition::User;
        mesh.bc_upper[0] = boundary::BoundaryCondition::User;
        expect_failure([&] {
            boundary::ApplyBoundaryCondition(mesh, 0, boundary::BoundarySide::Lower, 0.0, {});
        });
        boundary::ApplyBoundaryCondition(mesh, 0, boundary::BoundarySide::Lower, 0.5, hooks);
        assert(boundary_calls.back().time == 0.5);
        // Unsupported built-in types must still fail even when a hook is registered.
        mesh.bc_upper[0] = boundary::BoundaryCondition::Outflow;
        const auto count = boundary_calls.size();
        expect_failure([&] {
            boundary::ApplyBoundaryCondition(mesh, 0, boundary::BoundarySide::Upper, 0.75, hooks);
        });
        assert(boundary_calls.size() == count);
    }

    std::filesystem::path snapshot(const config::SimulationConfig& cfg, std::size_t step)
    {
        std::ostringstream name;
        name << "mesh_" << cfg.output_suffix << '_' << std::setw(8) << std::setfill('0') << step << ".csv";
        return std::filesystem::path(cfg.output_directory) / name.str();
    }

    void check_csv(const std::filesystem::path& path, const hd_core::MeshBlock& mesh,
                   std::size_t step, double time, double dt)
    {
        std::ifstream file(path);
        std::string line;
        assert(std::getline(file, line));
        assert(line == "# step=" + std::to_string(step));
        assert(std::getline(file, line));
        assert(line.rfind("# time=", 0) == 0);
        assert(std::stod(line.substr(7)) == time);
        assert(std::getline(file, line));
        assert(line.rfind("# dt=", 0) == 0);
        assert(std::stod(line.substr(5)) == dt);
        assert(std::getline(file, line));
        assert(line == "x,y,z,gamma,rho,vx,vy,vz,p,is_ghost");
        std::size_t f = 0;
        while (std::getline(file, line))
        {
            assert(f < mesh.cells.size());
            std::istringstream row(line);
            std::vector<double> fields;
            std::string field;
            while (std::getline(row, field, ',')) fields.push_back(std::stod(field));
            assert(fields.size() == 10);
            const auto index = mesh.coordinates(f);
            const auto& cell = mesh.cells[f];
            const auto W = cell.U.to_Primtivies(cell.gamma);
            for (int d = 0; d < 3; ++d)
            {
                assert(near(fields[d], d < hd_core::NDIM ? mesh.cell_center(d, index) : 0.0));
                assert(near(fields[5 + d], d < hd_core::NDIM ? W[hd_core::VEL(d)] : 0.0));
            }
            assert(near(fields[3], cell.gamma));
            assert(near(fields[4], W[hd_core::RHO]));
            assert(near(fields[8], W[hd_core::PRE]));
            assert(fields[9] == (active(mesh, index) ? 0.0 : 1.0));
            ++f;
        }
        assert(f == mesh.cells.size());
    }

    void output_tests(hd_core::MeshBlock& mesh, const std::filesystem::path& directory)
    {
        config::SimulationConfig cfg;
        cfg.output_directory = (directory / "snapshots").string();
        cfg.output_suffix = "check";
        cfg.output_step_every = 3;
        cfg.t_end = 0.8;
        hd_user::UserHooks hooks;
        hooks.useroutput = user_output;
        const auto before = mesh.cells;
        std::ostringstream log;
        auto* previous = std::cout.rdbuf(log.rdbuf());
        for (std::size_t step = 0; step <= 4; ++step)
            io::Output(mesh, cfg, hooks, step, step == 4 ? 0.8 : step * 0.1, step == 0 ? 0.0 : 0.1);
        std::cout.rdbuf(previous);
        assert(user_outputs == 4);
        const auto text = log.str();
        assert(std::count(text.begin(), text.end(), '\n') == 5);
        assert(text.find("mass = ") != std::string::npos);
        assert(text.find("energy = ") != std::string::npos);
        double expected_mass = 0.0, expected_energy = 0.0;
        for (std::size_t f = 0; f < mesh.cells.size(); ++f)
        {
            const auto index = mesh.coordinates(f);
            if (!active(mesh, index)) continue;
            double volume = 1.0;
            for (int d = 0; d < hd_core::NDIM; ++d)
                volume *= mesh.face[d][index[d] + 1] - mesh.face[d][index[d]];
            expected_mass += volume * mesh.cells[f].U[hd_core::RHO];
            expected_energy += volume * mesh.cells[f].U[hd_core::ENERGY];
        }
        std::istringstream messages(text);
        std::string message;
        while (std::getline(messages, message))
        {
            const double mass = std::stod(message.substr(message.find("mass = ") + 7));
            const double energy = std::stod(message.substr(message.find("energy = ") + 9));
            assert(std::abs(mass - expected_mass) < 1e-8 * (1 + expected_mass));
            assert(std::abs(energy - expected_energy) < 1e-8 * (1 + expected_energy));
        }
        for (const auto step : {0, 1, 3, 4})
            check_csv(snapshot(cfg, step), mesh, step,
                      step == 4 ? 0.8 : step * 0.1, step == 0 ? 0.0 : 0.1);
        assert(!std::filesystem::exists(snapshot(cfg, 2)));
        for (std::size_t f = 0; f < mesh.cells.size(); ++f)
        {
            assert(mesh.cells[f].U.data == before[f].U.data);
            assert(mesh.cells[f].gamma == before[f].gamma);
        }
        cfg.t_end = 10.0;
        cfg.maxsteps = 5;
        io::Output(mesh, cfg, {}, 5, 0.5, 0.1);
        check_csv(snapshot(cfg, 5), mesh, 5, 0.5, 0.1);
        const auto large_step = std::numeric_limits<std::size_t>::max();
        const double precise_time = 1.2345678901234567;
        const double small_dt = 1.2345678901234567e-12;
        io::Output(mesh, cfg, {}, large_step, precise_time, small_dt);
        check_csv(snapshot(cfg, large_step), mesh, large_step, precise_time, small_dt);
        cfg.output_step_every = 0;
        expect_failure([&] { io::Output(mesh, cfg, hooks, 6, 0.6, 0.1); });
        cfg.output_step_every = 3;
        const auto blocked = directory / "regular_file";
        std::ofstream(blocked) << "not a directory";
        cfg.output_directory = blocked.string();
        expect_failure([&] { io::Output(mesh, cfg, hooks, 6, 0.6, 0.1); });
        assert(user_outputs == 4);
    }

    void parameter_tests(const std::filesystem::path& directory)
    {
        const auto path = directory / "parameters.txt";
        const auto write = [&](const std::string& extra) {
            std::ofstream file(path);
            const char* axes[] = {"x", "y", "z"};
            for (int d = 0; d < hd_core::NDIM; ++d)
                file << 'n' << axes[d] << "a=2\nn" << axes[d] << "g=3\n"
                     << axes[d] << "_lower=0\n" << axes[d] << "_upper=1\n";
            file << extra;
        };
        write("output_step_every=7\n");
        assert(config::ReadParameters(path.string()).output_step_every == 7);
        write("step_every=4\n"); // Supported compatibility alias for output_step_every.
        assert(config::ReadParameters(path.string()).output_step_every == 4);
        for (const auto& value : {"0", "-1", "2.5"})
        {
            write(std::string("output_step_every=") + value + '\n');
            expect_failure([&] { config::ReadParameters(path.string()); });
        }
        write("riemann_tolerance=1e-12\nriemann_max_iterations=240\nmaxsteps=500000\n");
        const auto cfg = config::ReadParameters(path.string());
        assert(cfg.riemann_tolerance == 1e-12);
        assert(cfg.riemann_max_iterations == 240);
        assert(cfg.maxsteps == 500000);
        const auto largest = std::numeric_limits<std::size_t>::max();
        write("maxsteps=" + std::to_string(largest) + '\n');
        assert(config::ReadParameters(path.string()).maxsteps == largest);
        for (const auto& value : {"-1", "-0", "0", "2.5"})
        {
            write(std::string("maxsteps=") + value + '\n');
            expect_failure([&] { config::ReadParameters(path.string()); });
        }
        write("maxsteps=" + std::to_string(largest) + "0\n");
        expect_failure([&] { config::ReadParameters(path.string()); });
        for (const auto& value : {"0", "-1", "nan", "inf", "1e999"})
        {
            write(std::string("riemann_tolerance=") + value + '\n');
            expect_failure([&] { config::ReadParameters(path.string()); });
        }
        for (const auto& value : {"0", "-1", "2.5", "999999999999999999999"})
        {
            write(std::string("riemann_max_iterations=") + value + '\n');
            expect_failure([&] { config::ReadParameters(path.string()); });
        }
    }
}

int main(int argc, char** argv)
{
    assert(argc == 2);
    const std::filesystem::path directory = argv[1];
    std::filesystem::create_directories(directory);
    auto mesh = make_mesh();
    periodic_tests(mesh);
    for (std::size_t f = 0; f < mesh.cells.size(); ++f)
        if (active(mesh, mesh.coordinates(f))) mesh.cells[f].U *= 1.1;
    periodic_tests(mesh); // Refresh stale ghosts after the active state changes.
    hook_tests(mesh);
    output_tests(mesh, directory);
    parameter_tests(directory);
    std::cout << "NDIM=" << hd_core::NDIM << ": config, boundary hooks, periodic boundaries and CSV passed\n";
}
