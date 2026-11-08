#include "test_support.h"
#include "config/config.h"
#include "config/read_profile.h"
#include "config/read_mesa_profile.h"
#include "physcon/unit.h"
#include "utils/hdmath.h"
#include "user/user_api.h"
#include "hydro/reconstruction.h"
#include "hydro/riemann/factory.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <unordered_set>

namespace
{
config::SimulationConfig configuration()
{
    config::SimulationConfig cfg;
    cfg.nactive.fill(2);
    cfg.nghost.fill(1);
    cfg.lower.fill(0);
    cfg.upper.fill(1);
    cfg.boundary_lower.fill("periodic");
    cfg.boundary_upper.fill("periodic");
    return cfg;
}

void file_profile(const std::filesystem::path& directory)
{
    auto cfg = configuration();
    cfg.init_profile_from_file = true;
    cfg.init_profile_file_name = (directory / "profile.dat").string();
    const int count = 1 << hd_core::NDIM;
    {
        std::ofstream out(cfg.init_profile_file_name);
        out << "pre vz rho gamma vy vx\n";
        for (int i = 0; i < count; ++i)
            out << "2D0 0.3 " << i + 1 << " 1.4 0.2 0.1\n";
    }
    hd_core::MeshBlock mesh;
    // Profile-only input must work on an automatically generated uniform mesh.
    // The const overload is also part of the pre-update source API.
    const auto old_style_cfg = cfg;
    mesh.build_mesh(old_style_cfg);
    int row = 0;
    for (std::size_t flat = 0; flat < mesh.cells.size(); ++flat)
    {
        if (!mesh.is_active(flat)) continue;
        const auto& cell = mesh.cells[flat];
        const auto w = cell.U.to_Primtivies(cell.gamma);
        CHECK(near(w[hd_core::RHO], ++row));
        CHECK(near(w[hd_core::PRE], 2));
        for (int d = 0; d < hd_core::NDIM; ++d)
            CHECK(near(w[hd_core::VEL(d)], 0.1 * (d + 1)));
    }
    CHECK(row == count);
    cfg.init_mesh_from_file = true;
    cfg.init_mesh_file_name = (directory / "mesh.dat").string();
    {
        std::ofstream out(cfg.init_mesh_file_name);
        out << "x y z\n0 0 0\n0.2 0.3 0.4\n1 2 3\n";
    }
    mesh.build_mesh(cfg);
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        CHECK(cfg.lower[d] == 0);
        CHECK(cfg.upper[d] == d + 1);
        CHECK(near(mesh.face[d][0], -(0.2 + 0.1 * d)));
        CHECK(near(mesh.face[d].back(), 2 * (d + 1) - (0.2 + 0.1 * d)));
    }
    auto reflected = cfg;
    reflected.nghost.fill(2);
    reflected.boundary_lower.fill("reflective");
    reflected.boundary_upper.fill("reflective");
    hd_core::MeshBlock mirrored;
    mirrored.build_mesh(reflected);
    for (int d = 0; d < hd_core::NDIM; ++d)
    {
        CHECK(near(mirrored.face[d][0], -(d + 1)));
        CHECK(near(mirrored.face[d].back(), 2 * (d + 1)));
    }
    // A later invalid row must not partially overwrite valid loaded states.
    const auto before = mesh.cells;
    {
        std::ofstream out(cfg.init_profile_file_name);
        out << "rho vx vy vz pre gamma\n";
        for (int i = 0; i < count; ++i)
            out << "99 0 0 0 " << (i == count - 1 ? -1 : 1) << " 1.4\n";
    }
    expect_throw<std::runtime_error>([&] { config::read_profile(cfg.init_profile_file_name, mesh); });
    for (std::size_t i = 0; i < before.size(); ++i)
        for (int k = 0; k < hd_core::NVAR; ++k)
            CHECK(mesh.cells[i].U[k] == before[i].U[k]);

    cfg.init_from_mesa_profile = true;
    expect_throw<std::runtime_error>([&] { config::validate_cfg(cfg, {}); });
    // A rectangular table has the same row count in every axis column. Reject
    // incompatible per-axis cell counts instead of indexing beyond the data.
    if constexpr (hd_core::NDIM >= 2)
    {
        cfg.init_from_mesa_profile = false;
        cfg.nactive[1] = 3;
        expect_throw<std::runtime_error>([&] { mesh.build_mesh(cfg); });
    }
}

void mesa(const std::filesystem::path& directory)
{
    const auto path = (directory / "mesa.data").string();
    {
        std::ofstream out(path);
        out << "ignored\nignored\nignored\n\nignored\nlogR dr logRho logP velocity\n";
        out << std::setprecision(17);
        for (int i = 3; i >= 1; --i)
            out << std::log10(std::pow(2.0, i)) << " 1 "
                << std::log10(std::pow(2.0, 5 - i)) << ' '
                << std::log10(std::pow(2.0, 7 - i)) << ' ' << 2 * i << '\n';
    }
    physcon::units.init_unitsystem_MLT(8, 2, 4, 1);
    // In these units the native faces are 0.5, 1, 2, 4.
    const auto cut = config::read_mesa_profile(path, 1.5, 5, "linear", 2, 2, 1);
    CHECK(cut.radius_faces.size() == 5);
    CHECK(near(cut.radius_faces[0], 1.5));
    CHECK(near(cut.radius_faces[1], 2));
    CHECK(near(cut.radius_faces[2], 4));
    CHECK(near(cut.radius_faces[3], 4.5));
    CHECK(near(cut.radius_faces[4], 5));
    CHECK(cut.zones.data.size() == 2);
    auto cfg = configuration();
    cfg.gamma = 1.4;
    if constexpr (hd_core::NDIM == 1)
    {
        auto mesh = config::build_mesh_from_mesa_profile(cut, cfg, true);
        CHECK(cfg.geometry == "spherical");
        CHECK(cfg.nactive[0] == 4);
        CHECK(cfg.lower[0] == 1.5 && cfg.upper[0] == 5);
        auto w = mesh.at_active({0}).U.to_Primtivies(cfg.gamma);
        CHECK(near(w[hd_core::RHO], 8));
        CHECK(near(w[hd_core::PRE], 128));
        CHECK(near(w[hd_core::VEL(0)], 8));
        CHECK(mesh.at_active({2}).U[hd_core::RHO] == 0);

        // Regression: a reflecting wall on an unequal-width grid must have
        // exactly mirrored reconstructed velocities and zero mass/energy flux.
        cfg.nghost[0] = 2;
        cfg.boundary_lower[0] = "reflective";
        cfg.boundary_upper[0] = "reflective";
        mesh = config::build_mesh_from_mesa_profile(cut, cfg, false);
        for (int i = 0; i < mesh.nactive[0]; ++i)
        {
            hd_core::Primitives w{};
            w[hd_core::RHO] = 1;
            w[hd_core::PRE] = 1;
            w[hd_core::VEL(0)] = 1 + i;
            mesh.at_active({i}) = hd_core::Cell(w, cfg.gamma);
        }
        boundary::FillBoundary(mesh, 0, {});
        hd_core::Index wall{};
        wall[0] = mesh.nghost[0];
        auto states = hydro::Reconstruct_PLM(mesh, 0, wall, hydro::Limiter::MC);
        CHECK(near(states.first[hd_core::VEL(0)], -states.second[hd_core::VEL(0)], 1e-14));
        auto solver = hydro::make_riemann_solver("exact", 1e-10, 200);
        auto flux = solver->flux(states.first, states.second, 0, cfg.gamma, cfg.gamma);
        CHECK(std::abs(flux[hd_core::RHO]) < 1e-13);
        CHECK(std::abs(flux[hd_core::ENERGY]) < 1e-13);
    }
    else
        expect_throw<std::invalid_argument>([&] { config::build_mesh_from_mesa_profile(cut, cfg, true); });
    const auto inner = config::read_mesa_profile(path, 0, 3, "linear", 2, 2, 1);
    CHECK(inner.zones.data.size() == 5);
    CHECK(inner.radius_faces.size() == 6);
    CHECK(inner.radius_faces.front() == 0 && inner.radius_faces.back() == 3);
    CHECK(inner.zones.data[0] == inner.zones.data[2]);
    const auto log = config::read_mesa_profile(path, 0.25, 8, "log", 2, 2, 1);
    CHECK(near(log.radius_faces[1], std::sqrt(0.25 * 0.5)));
    CHECK(near(log.radius_faces[log.radius_faces.size() - 2], std::sqrt(4.0 * 8)));
    expect_throw<std::invalid_argument>([&] { config::read_mesa_profile(path, 5, 6, "linear", 2, 2, 1); });
    expect_throw<std::invalid_argument>([&] { config::read_mesa_profile(path, 0, 8, "log", 2, 2, 1); });
    physcon::units.init_unitsystem_MLT(1, 1, 1, 1);
}
}

int main(int argc, char** argv)
{
    return run_test("file and MESA initialization", [&] {
        CHECK(argc == 2);
        const std::filesystem::path directory(argv[1]);
        std::filesystem::create_directories(directory);
        file_profile(directory);
        mesa(directory);
        CHECK(near(hdmath::interp_linear({0, 1, 2, 3}, {0, 1, 0, 1}, {0.25, 0.75}, 2), 1.75));
    });
}
