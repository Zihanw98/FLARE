#pragma once
#include "core.h"
#include "config.h"
#include <string>
namespace config
{
    struct Table
    {
        std::vector<std::string> headers;
        std::vector<std::vector<double>> data;

        std::vector<double> column(const std::string& name) const;
    };
    Table read_table(const std::string& filename,
                    std::size_t skip_rows = 0,
                    bool has_header = true);
}