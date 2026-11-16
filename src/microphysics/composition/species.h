#pragma once
#include "config.h"
#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace hd_microphysics
{
    inline constexpr char DEFAULT_SPECIES_TABLE[] = "tables/species.json";

    struct Species
    {
        std::string name;
        std::string name_base;
        int A;
        int Z;
        int charge = 0;
        // 读取全表时只记录邻居名称，整数索引尚未建立，统一为 nullopt。
        std::optional<int> ionization_upper = std::nullopt;
        std::optional<int> ionization_lower = std::nullopt;
        double gamma = 5.0/3.0;
        double energy_eV = 0.0;  // Per particle, relative to the neutral isotope's ground state.
        // Electronic ground-state weight g = 2J+1 for Saha ionization; nuclear spin omitted.
        // Untabulated species retain the placeholder value 1.
        double ground_weight = 1.0;
    };

    // 全表读取结果，仅存放数据，不包含筛选、排序或查询方法。
    struct SpeciesTable
    {
        std::vector<Species> species_full;  // 保持 JSON 文件中的顺序。
        // key 是自身在 species 中的下标，value 是相邻电离态的物种名称。
        // 每个物种在两个 map 中都有一项；JSON null 对应空字符串。
        std::map<std::size_t, std::string> ionization_upper_names;
        std::map<std::size_t, std::string> ionization_lower_names;
    };

    // 读取全部物种信息及上下电离态名称，校验数据，但不建立上下电离态整数索引。
    // 默认表路径也支持从算例目录启动时定位。
    SpeciesTable read_species_table(const std::string& filename = DEFAULT_SPECIES_TABLE);
}
