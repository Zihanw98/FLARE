#include "test_support.h"
#include "config/profile_composition.h"
#include "config/read_profile.h"
#include "config/read_mesa_profile.h"
#include "microphysics/composition/species.h"
#include "microphysics/eos/eos.h"
#include "physcon/unit.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>

namespace
{
config::SimulationConfig configuration(const std::filesystem::path& directory)
{
    CHECK(hd_core::NCOMP == 10);
    config::SimulationConfig cfg;
    cfg.nactive.fill(2);
    cfg.nghost.fill(1);
    cfg.lower.fill(0);
    cfg.upper.fill(1);
    cfg.boundary_lower.fill("outflow");
    cfg.boundary_upper.fill("outflow");
    cfg.equation_of_state = "IdealGas_Mixed";
    hd_microphysics::eos_type = hd_microphysics::Equation_of_State::IdealGas_Mixed;
    cfg.composition_activate_file = (directory / "active.txt").string();
    // Deliberately reverse the isotope order, including two isotopes at the same Z.
    std::ofstream(cfg.composition_activate_file) << "c12 0 0\nhe4 0 2\nh1 -1 1\nhe3 0 2\n";
    hd_microphysics::initialize_species_table(cfg);
    return cfg;
}

void check_neutral(const hd_core::Cell& cell, double hydrogen, double helium3,
                   double helium4, double carbon)
{
    const auto& table = hd_microphysics::species_table;
    double sum = 0;
    for (int k = 0; k < hd_core::NCOMP; ++k)
    {
        const auto& sp = table.species_active[k];
        const double expected = sp.charge != 0 ? 0.0
            : sp.Z == 1 ? hydrogen : sp.Z == 6 ? carbon : sp.A == 3 ? helium3 : helium4;
        CHECK(near(cell.U[hd_core::COMP(k)], cell.U[hd_core::RHO] * expected));
        sum += cell.U[hd_core::COMP(k)];
    }
    CHECK(near(sum, cell.U[hd_core::RHO]));
}

void profiles(const std::filesystem::path& directory, const config::SimulationConfig& cfg)
{
    hd_core::MeshBlock mesh;
    mesh.build_mesh(cfg);
    const int count = 1 << hd_core::NDIM;
    const auto path = (directory / "profile.dat").string();
    auto write = [&](const std::string& names, const std::string& values,
                     const std::string& last_values = "") {
        std::ofstream out(path);
        out << "rho vx vy vz pre gamma " << names << '\n';
        for (int i = 0; i < count; ++i)
            out << i + 2 << " 0.1 0.2 0.3 2 1.4 "
                << (i == count - 1 && !last_values.empty() ? last_values : values) << '\n';
    };

    write("He_I C12 HI He3_I", "0.2 0.3 0.4 0.1");
    config::read_profile(path, mesh);
    for (std::size_t f = 0; f < mesh.cells.size(); ++f)
    {
        if (mesh.is_active(f)) check_neutral(mesh.cells[f], .4, .1, .2, .3);
        else CHECK(mesh.cells[f].U[hd_core::RHO] == 0.0);
    }
    // Isotope-name aliases can be used directly, independent of column order/case.
    write("he3 H1 c12 HE4", "0.1 0.4 0.3 0.2");
    config::read_profile(path, mesh);
    for (std::size_t f = 0; f < mesh.cells.size(); ++f)
        if (mesh.is_active(f)) check_neutral(mesh.cells[f], .4, .1, .2, .3);

    // Existing explicit charge populations are preserved, including H-.
    write("HI He3_I He_I C12 H- HII He3_II He3_III He_II He_III",
          "0.15 0.1 0.2 0.3 0.05 0.2 0 0 0 0");
    config::read_profile(path, mesh);
    const auto& table = hd_microphysics::species_table;
    const int neutral = table.find_species_from_available_table_index(1, 1, 0);
    const int ionized = table.find_species_from_available_table_index(1, 1, 1);
    const int negative = table.find_species_from_available_table_index(1, 1, -1);
    for (std::size_t f = 0; f < mesh.cells.size(); ++f)
    {
        if (!mesh.is_active(f)) continue;
        const auto& u = mesh.cells[f].U;
        CHECK(near(u[hd_core::COMP(neutral)] / u[hd_core::RHO], .15));
        CHECK(near(u[hd_core::COMP(ionized)] / u[hd_core::RHO], .2));
        CHECK(near(u[hd_core::COMP(negative)] / u[hd_core::RHO], .05));
    }
    write("HI He3_I He_I C12", "0.4 0.1 0.2 0.3");
    config::read_profile(path, mesh);
    for (std::size_t f = 0; f < mesh.cells.size(); ++f)
        if (mesh.is_active(f)) check_neutral(mesh.cells[f], .4, .1, .2, .3);
    const auto before = mesh.cells;
    auto rejected = [&](const std::string& names, const std::string& values,
                        const std::string& last_values = "") {
        write(names, values, last_values);
        expect_throw<std::runtime_error>([&] {
            config::read_profile(path, mesh);
        });
        for (std::size_t f = 0; f < mesh.cells.size(); ++f)
        {
            CHECK(mesh.cells[f].gamma == before[f].gamma);
            for (int k = 0; k < hd_core::NVAR; ++k)
                CHECK(mesh.cells[f].U[k] == before[f].U[k]);
        }
    };
    rejected("HI He_I C12", "0.4 0.3 0.3"); // He-3 must not be mistaken for He-4.
    rejected("HI HI He3_I He_I C12", "0.2 0.2 0.1 0.2 0.3");
    rejected("HI h1 He3_I He_I C12", "0.2 0.2 0.1 0.2 0.3");
    rejected("HI He3_I He_I C12 HII HII", "0.4 0.1 0.2 0.3 0 0");
    rejected("HI He3_I He_I C12", "0.4 0.1 0.2 0.3", "0.3 0.1 0.2 0.3");
    rejected("HI He3_I He_I C12", "0.4 0.1 0.2 0.3", "-0.1 0.1 0.2 0.8");
    rejected("HI He3_I He_I C12", "0.4 0.1 0.2 0.3", "1.1 0 0 0");

    const config::Table input{{"h1", "he3", "he4", "c12"}, {}};
    const auto columns = config::profile_composition_columns(
        input, config::CompositionColumnFormat::Species, "test");
    hd_core::Primitives state{};
    for (double invalid : {std::numeric_limits<double>::quiet_NaN(),
                           std::numeric_limits<double>::infinity()})
        expect_throw<std::runtime_error>([&] {
            config::set_profile_composition({invalid, .1, .2, .3}, columns, state, "test row");
        });
    expect_throw<std::runtime_error>([&] {
        config::set_profile_composition({.4, .1}, columns, state, "test row");
    });
}

void profiles_require_composition(const std::filesystem::path& directory,
                                  config::SimulationConfig cfg)
{
    cfg.init_profile_from_file = true;
    cfg.init_profile_file_name = (directory / "hydro_profile.dat").string();
    {
        std::ofstream out(cfg.init_profile_file_name);
        out << "rho vx vy vz pre gamma\n";
        for (int i = 0; i < (1 << hd_core::NDIM); ++i)
            out << i + 2 << " 0.1 0.2 0.3 2 1.4\n";
    }
    // A composition build always requires isotope fractions in state input.
    hd_core::MeshBlock mesh;
    expect_throw<std::runtime_error>([&] { mesh.build_mesh(cfg); });
}

void mesa(const std::filesystem::path& directory, config::SimulationConfig cfg)
{
    if constexpr (hd_core::NDIM != 1) return;
    const auto path = (directory / "mesa.data").string();
    {
        std::ofstream out(path);
        out << "ignored\nignored\nignored\n\nignored\n"
               "logR dr logRho logP velocity he4 c12 prot he3 h1\n" << std::setprecision(17);
        for (int i = 3; i >= 1; --i)
            out << std::log10(std::pow(2.0, i)) << " 1 "
                << std::log10(std::pow(2.0, 5 - i)) << ' '
                << std::log10(std::pow(2.0, 7 - i)) << ' ' << 2 * i
                << " 0.2 " << .3 - .01 * i << " 0.05 0.1 " << .35 + .01 * i << '\n';
    }
    physcon::units.init_unitsystem_MLT(8, 2, 4, 1);
    const auto profile = config::read_mesa_profile(path, .25, 5, "linear", 2, 2, 1);
    const auto mesh = config::build_mesh_from_mesa_profile(profile, cfg, true);
    CHECK(profile.zones.data.size() == 5); // Includes copied inner padding rows.
    CHECK(mesh.nactive[0] == 7); // Two exterior cells have no source rows.
    for (std::size_t i = 0; i < profile.zones.data.size(); ++i)
    {
        const int source = std::max(1, static_cast<int>(i) - 1);
        check_neutral(mesh.at_active({static_cast<int>(i)}), .4 + .01 * source, .1, .2, .3 - .01 * source);
    }
    CHECK(near(mesh.at_active({0}).U[hd_core::RHO], 16));
    for (int i = 5; i < mesh.nactive[0]; ++i)
        CHECK(mesh.at_active({i}).U[hd_core::RHO] == 0);

    const auto original_cfg = cfg;
    auto rejected = [&](const config::MesaProfile& bad) {
        expect_throw<std::runtime_error>([&] { config::build_mesh_from_mesa_profile(bad, cfg, true); });
        CHECK(cfg.nactive == original_cfg.nactive && cfg.lower == original_cfg.lower &&
              cfg.upper == original_cfg.upper && cfg.geometry == original_cfg.geometry);
    };
    auto bad = profile;
    bad.zones.headers.back() = "missing_h1";
    rejected(bad); // prot alone cannot account for the missing h1 abundance.
    bad = profile;
    bad.zones.headers[8] = "other_he3";
    rejected(bad);
    bad = profile;
    bad.zones.data.back()[9] = -.1;
    rejected(bad);
    bad = profile;
    bad.zones.data.back()[9] = .25;
    rejected(bad);
    bad = profile;
    bad.zones.headers.push_back("HE4");
    for (auto& row : bad.zones.data) row.push_back(.2);
    rejected(bad);

    // A network that supplies prot without h1 still maps its hydrogen to HI.
    auto protons = profile;
    protons.zones.headers.pop_back();
    for (auto& row : protons.zones.data) { row[7] += row.back(); row.pop_back(); }
    const auto proton_mesh = config::build_mesh_from_mesa_profile(protons, cfg, true);
    check_neutral(proton_mesh.at_active({0}), .41, .1, .2, .29);

    // Geometry-only loading needs no isotopes; state initialization does.
    auto hydro = profile;
    hydro.zones.headers.resize(5);
    for (auto& row : hydro.zones.data) row.resize(5);
    config::build_mesh_from_mesa_profile(hydro, cfg, false);
    expect_throw<std::runtime_error>([&] {
        config::build_mesh_from_mesa_profile(hydro, cfg, true);
    });
    physcon::units.init_unitsystem_MLT(1, 1, 1, 1);
}
}

int main(int argc, char** argv)
{
    return run_test("neutral profile and MESA isotope initialization", [&] {
        CHECK(argc == 2);
        const auto directory = std::filesystem::absolute(argv[1]);
        std::filesystem::create_directories(directory);
        const auto cfg = configuration(directory);
        profiles(directory, cfg);
        mesa(directory, cfg);
        profiles_require_composition(directory, cfg);
    });
}
