#include "test_support.h"
#include "config/config.h"
#include "microphysics/composition/species.h"
#include <algorithm>
#include <filesystem>
#include <fstream>

namespace
{
void test_activation(const std::filesystem::path& directory,
                     const std::vector<hd_microphysics::Species>& catalog)
{
    using hd_microphysics::read_active_species;
    using hd_microphysics::load_run_species;
    std::filesystem::create_directories(directory);
    const auto file = directory / "active.txt";
    std::ofstream(file) << "# base  min  max\n\nhe4 0 +1\nh1 -1 +1 # include H-\nhe3 0 0\n";
    const auto active = read_active_species(file.string(), catalog);
    CHECK(active.size() == 6);
    CHECK(active.find_species(2, 2, 1) == 1);
    CHECK(active.find_species(1, 0, -1) == 2);
    CHECK(active.find_species(2, 1, 0) == 5);
    CHECK(!active.find_species(2, 1, 1) && !active.find_species(2, 2, 2));
    CHECK(active[3].name == "HI" && active[3].ionization_lower == 2 && active[3].ionization_upper == 4);
    CHECK(active[1].ionization_upper == hd_microphysics::NO_IONIZATION_NEIGHBOR);

    const hd_microphysics::SpeciesList full(catalog);
    for (std::size_t k = 0; k < catalog.size(); ++k)
        CHECK(full.find_species(catalog[k].Z, catalog[k].A - catalog[k].Z, catalog[k].charge) == static_cast<int>(k));

    const auto parameters = directory / "parameters.txt";
    {
        std::ofstream out(parameters);
        for (const char axis : std::string("xyz"))
            out << 'n' << axis << "a = 2\nn" << axis << "g = 1\n" << axis << "_lower = 0\n" << axis << "_upper = 1\n";
        out << "enable_composition = true\ncomposition_activate_file = active.txt\n";
    }
    const auto cfg = config::ReadParameters(parameters.string());
    CHECK(cfg.enable_composition && cfg.composition_activate_file == "active.txt");
    if constexpr (hd_core::NCOMP == 6)
    {
        const auto selected = load_run_species(cfg, parameters.string());
        CHECK(selected.size() == 6);
        CHECK(selected.find_species(2, 1, 0) == 5);
    }
    else
        expect_throw<std::runtime_error>([&] { load_run_species(cfg, parameters.string()); });

    for (const char* invalid : {"he4 1 2", "h1 -1 -1", "he4 0 3", "h1 -2 0", "o16 0 1",
                                "he9 0 0", "he3 0 1\nhe3 0 0", "h1 0+1", "h1 0 1.5", "h1 -1 1 extra"})
    {
        std::ofstream(file) << invalid;
        expect_throw<std::runtime_error>([&] { read_active_species(file.string(), catalog); });
    }
    std::filesystem::remove(file);
    CHECK(config::ReadParameters(parameters.string()).enable_composition);
    expect_throw<std::runtime_error>([&] { load_run_species(cfg, parameters.string()); });
    std::ofstream(parameters, std::ios::app) << "enable_composition = false\n";
    const auto disabled = config::ReadParameters(parameters.string());
    const auto empty = load_run_species(disabled, parameters.string());
    CHECK(empty.size() == 0);
}
}

int main(int argc, char** argv)
{
    using hd_microphysics::NO_IONIZATION_NEIGHBOR;
    using hd_microphysics::read_species_table;
    using hd_microphysics::select_species;
    return run_test("species catalog, activation and indexed lookup", [&] {
        CHECK(argc == 2);
        const auto directory = std::filesystem::absolute(argv[1]);
        const auto root = std::filesystem::absolute(__FILE__).parent_path().parent_path();
        std::filesystem::current_path(root / "tests" / "sod_tube");
        const auto all = read_species_table();
        test_activation(directory, all);
        const std::string path = "tables/species.json";
        for (const auto& sp : all)
        {
            if (sp.charge == 0)
                CHECK(sp.energy_eV == 0.0);
            if (sp.ionization_upper != NO_IONIZATION_NEIGHBOR)
            {
                const auto& upper = all.at(sp.ionization_upper);
                CHECK(upper.energy_eV > sp.energy_eV);
                CHECK(all.at(upper.ionization_lower).name == sp.name);
            }
        }

        const auto h = select_species(all, {"HII", "H-", "HI"});
        CHECK(h.size() == 3 && h[2].name == "HI" && h[2].name_base == "h1");
        CHECK(h[0].charge == +1 && h[1].charge == -1 && h[2].charge == 0);
        CHECK(h[2].ionization_lower == 1 && h[2].ionization_upper == 0);
        CHECK(h[0].ionization_lower == 2 && h[1].ionization_upper == 2);
        CHECK(near(h[0].energy_eV, 13.598435, 1e-7) && near(h[1].energy_eV, -0.7542, 3e-5));

        const auto he = select_species(all, {"He_III", "He3_I", "He_II", "He3_III", "O16"});
        CHECK(he[0].A == 4 && he[0].Z == 2 && he[0].charge == 2);
        CHECK(he[3].A == 3 && he[3].Z == 2 && he[3].charge == 2);
        CHECK(he[0].ionization_lower == 2 && he[2].ionization_upper == 0);
        CHECK(he[2].ionization_lower == NO_IONIZATION_NEIGHBOR);
        CHECK(he[1].ionization_upper == NO_IONIZATION_NEIGHBOR);
        CHECK(he[3].ionization_lower == NO_IONIZATION_NEIGHBOR);
        CHECK(near(he[0].energy_eV - he[2].energy_eV, 54.417766, 1e-7));
        CHECK(near(he[3].energy_eV, 79.001418, 1e-7));
        CHECK(he[4].A == 16 && he[4].Z == 8 && he[4].charge == 0);

        CHECK(select_species(all, {}).empty());
        const auto neutral_h = read_species_table(path, {"HI"});
        CHECK(neutral_h.size() == 1);
        CHECK(neutral_h[0].ionization_upper == NO_IONIZATION_NEIGHBOR);
        CHECK(neutral_h[0].ionization_lower == NO_IONIZATION_NEIGHBOR);
        const auto hi = std::find_if(all.begin(), all.end(), [](const auto& sp) { return sp.name == "HI"; });
        CHECK(hi != all.end());
        CHECK(all.at(hi->ionization_upper).name == "HII" && all.at(hi->ionization_lower).name == "H-");
        CHECK(read_species_table("../../tables/species.json").size() == all.size());
        expect_throw<std::runtime_error>([] { read_species_table("../../tables/species.dat"); });
        expect_throw<std::runtime_error>([&] { select_species(all, {"missing"}); });
        expect_throw<std::runtime_error>([&] { select_species(all, {"HI", "HI"}); });
    });
}
