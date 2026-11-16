#pragma once
#include <cstddef>
#include <string>
#include <vector>
namespace config
{
    struct Table
    {
        std::vector<std::string> headers;
        std::vector<std::vector<double>> data;

        std::vector<double> column(const std::string& name) const;
    };
    struct CompositionRow
    {
        std::string name_base;  // 核素名称，例如 h1、he4。
        int q_min = 0;          // 最小电荷。
        int q_max = 0;          // 最大电荷。
    };

    struct CompositionTable
    {
        std::vector<CompositionRow> data;  // 按文件顺序保存三列数据。
    };

    // 读取无表头的组分启用表，每行必须是 name_base q_min q_max（string、int、int）。
    // 忽略空行、# 注释行及行内第一个 # 后的内容；列数错误或非整数时报出文件和行号。
    CompositionTable read_composition_table(const std::string& filename);

    // 数值表：支持跳过物理行、可选表头及 D/d 指数。
    Table read_table(const std::string& filename,
                    std::size_t skip_rows = 0,
                    bool has_header = true);
}
