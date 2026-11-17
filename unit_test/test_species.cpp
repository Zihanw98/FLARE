#include "test_support.h"
#include "microphysics/composition/species.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>

int main(int argc, char** argv)
{
    using hd_microphysics::read_species_table;
    return run_test("species table and ionization name maps", [&] {
        CHECK(argc == 2);
        const auto directory = std::filesystem::absolute(argv[1]);
        std::filesystem::create_directories(directory);
        const auto root = std::filesystem::absolute(__FILE__).parent_path().parent_path();
        std::filesystem::current_path(root / "tests" / "sod_tube");
        const auto table = read_species_table();
        const auto& all = table.species_full;
        std::ifstream source(root / "tables/species.json");
        const auto document = nlohmann::json::parse(source);
        const auto& rows = document.at("species");
        CHECK(all.size() == rows.size());
        CHECK(table.ionization_upper_names.size() == all.size());
        CHECK(table.ionization_lower_names.size() == all.size());
        for (std::size_t k = 0; k < all.size(); ++k)
        {
            const auto& sp = all[k];
            const auto& row = rows[k];
            CHECK(sp.name == row.at("name").get<std::string>());
            CHECK(sp.A == row.at("A").get<int>() && sp.Z == row.at("Z").get<int>());
            CHECK(sp.charge == row.at("charge").get<int>());
            CHECK(!sp.ionization_upper.has_value() && !sp.ionization_lower.has_value());
            CHECK(table.ionization_upper_names.at(k) ==
                  (row.at("ionization_upper").is_null() ? "" : row.at("ionization_upper").get<std::string>()));
            CHECK(table.ionization_lower_names.at(k) ==
                  (row.at("ionization_lower").is_null() ? "" : row.at("ionization_lower").get<std::string>()));
            CHECK(std::isfinite(sp.ground_weight) && sp.ground_weight > 0.0);
            if (sp.charge == 0)
                CHECK(sp.energy_eV == 0.0);
        }
        CHECK(all[0].name == "H-" && all[1].name == "HI" && all[2].name == "HII");
        CHECK(near(all[0].energy_eV, -0.7542, 3e-5) && near(all[2].energy_eV, 13.598435, 1e-7));
        CHECK(all[0].ground_weight == 1.0 && all[1].ground_weight == 2.0 && all[2].ground_weight == 1.0);
        CHECK(all[3].ground_weight == 1.0 && all[4].ground_weight == 2.0 && all[5].ground_weight == 1.0);
        CHECK(all[6].ground_weight == 1.0 && all[7].ground_weight == 2.0 && all[8].ground_weight == 1.0);
        CHECK(near(all[8].energy_eV - all[7].energy_eV, 54.417766, 1e-7));
        CHECK(near(all[5].energy_eV, 79.001418, 1e-7));

        hd_microphysics::SpeciesTable active;
        expect_throw<std::out_of_range>([&] { active.find_species_from_available_table_index("H-"); });
        expect_throw<std::out_of_range>([&] { active.find_species_from_available_table_index(1, 1, -1); });
        // Sorted by Z, with He-4 before He-3 and only some charge states active.
        active.species_active = {all[0], all[7], all[4], all[3]};
        for (std::size_t k = 0; k < active.species_active.size(); ++k)
        {
            const auto& sp = active.species_active[k];
            CHECK(active.find_species_from_available_table_index(sp.name) == static_cast<int>(k));
            CHECK(active.find_species_from_available_table_index(sp.A, sp.Z, sp.charge) == static_cast<int>(k));
        }
        expect_throw<std::out_of_range>([&] { active.find_species_from_available_table_index("HI"); });
        expect_throw<std::out_of_range>([&] { active.find_species_from_available_table_index(4, 2, 0); });
        expect_throw<std::out_of_range>([&] { active.find_species_from_available_table_index(12, 6, 0); });
        CHECK(!active.find_species_from_available_table(4, 2, 0));

        CHECK(active.find_ionization_state_indices(2, 3) == std::vector<int>({3, 2}));
        CHECK(active.find_ionization_state_indices(2, 4) == std::vector<int>({1}));
        CHECK(active.find_ionization_state_indices(6, 12).empty());
        CHECK(hd_microphysics::SpeciesTable{}.find_ionization_state_indices(1, 1).empty());
        // Charge ordering must not depend on table ordering; include the negative ion.
        active.species_active = {all[2], all[0], all[1], all[7], all[3], all[6], all[4], all[8], all[5]};
        CHECK(active.find_ionization_state_indices(1, 1) == std::vector<int>({1, 2, 0}));
        CHECK(active.find_ionization_state_indices(2, 3) == std::vector<int>({4, 6, 8}));
        CHECK(active.find_ionization_state_indices(2, 4) == std::vector<int>({5, 3, 7}));

        // 换一个文件顺序：map 的 key 跟随全表下标，邻居仍保存名称，整数索引仍为空。
        const auto custom_file = directory / "custom.json";
        auto reordered = document;
        reordered["species"] = nlohmann::json::array({rows[2], rows[0], rows[1]});
        std::ofstream(custom_file) << reordered;
        const auto other = read_species_table(custom_file.string());
        CHECK(other.species_full.size() == 3);
        CHECK(other.ionization_upper_names.size() == 3 && other.ionization_lower_names.size() == 3);
        CHECK(other.species_full[0].name == "HII" && other.species_full[1].name == "H-" && other.species_full[2].name == "HI");
        CHECK(other.ionization_upper_names.at(2) == "HII" && other.ionization_lower_names.at(2) == "H-");
        CHECK(other.ionization_upper_names.at(0).empty() && other.ionization_lower_names.at(1).empty());
        for (const auto& sp : other.species_full)
            CHECK(!sp.ionization_upper && !sp.ionization_lower);
        CHECK(table.ionization_upper_names.at(0) == "HI" && table.ionization_lower_names.at(0).empty());

        // 保留基态统计权重的读取、默认值和非法输入检查。
        auto custom = document;
        custom["species"][0]["ground_weight"] = 2.5;
        std::ofstream(custom_file) << custom;
        CHECK(read_species_table(custom_file.string()).species_full[0].ground_weight == 2.5);
        for (const auto& invalid : nlohmann::json::array({0, -1, nullptr, "2", true}))
        {
            custom["species"][0]["ground_weight"] = invalid;
            std::ofstream(custom_file) << custom;
            expect_throw<std::runtime_error>([&] { read_species_table(custom_file.string()); });
        }
        custom["species"][0].erase("ground_weight");
        std::ofstream(custom_file) << custom;
        CHECK(read_species_table(custom_file.string()).species_full[0].ground_weight == 1.0);

        // 仅保留名称也必须拒绝未知邻居、错误电离关系和非字符串名称。
        for (const auto& invalid : nlohmann::json::array({"missing", "HII", 42}))
        {
            custom = document;
            custom["species"][0]["ionization_upper"] = invalid;
            std::ofstream(custom_file) << custom;
            expect_throw<std::runtime_error>([&] { read_species_table(custom_file.string()); });
        }
        custom = document;
        custom["species"][1]["name"] = custom["species"][0]["name"];
        std::ofstream(custom_file) << custom;
        expect_throw<std::runtime_error>([&] { read_species_table(custom_file.string()); });
        CHECK(read_species_table("../../tables/species.json").species_full.size() == all.size());
        expect_throw<std::runtime_error>([] { read_species_table("../../tables/species.dat"); });
    });
}
