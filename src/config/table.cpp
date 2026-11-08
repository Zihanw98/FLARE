#include "table.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

namespace config
{
    Table read_table(const std::string& filename,
                    std::size_t skip_rows,
                    bool has_header)
    {
        std::ifstream file(filename);
        if(!file.is_open())
        {
            throw std::runtime_error("Cannot open table:" + filename);
        }
        std::string line;
        for(size_t i = 0; i < skip_rows; i++)
        {
            if (!std::getline(file, line)) {
                throw std::runtime_error("Not enough lines to skip");
            }
        }
        Table table;
        if(has_header)
        {
            if(!std::getline(file,line))
            {
                throw std::runtime_error("Cannot read header");
            }

            std::istringstream iss(line);
            std::string name;
            
            while(iss >> name)
            {
                table.headers.push_back(name);
            }
        }
        while (std::getline(file, line)) 
        {
            std::istringstream fields(line);
            std::vector<double> row;
            std::string token;
            while (fields >> token)
            {
                const auto pos = token.find_first_of("Dd");
                if (pos != std::string::npos)
                    token[pos] = 'E';

                std::istringstream number(token);
                double value;
                if (!(number >> value) || !number.eof())
                    throw std::runtime_error("Invalid number '" + token + "' in " + filename);

                row.push_back(value);
            }

            if (row.empty()) {
                continue;
            }

            if (has_header && row.size() != table.headers.size()) {
                throw std::runtime_error("Column count mismatch: " + filename);
            }

            table.data.push_back(std::move(row));
        }
        return table;
    }
    std::vector<double> Table::column(const std::string& name) const {
        auto it = std::find(headers.begin(), headers.end(), name);
        if (it == headers.end()) {
            throw std::runtime_error("Column not found: " + name);
        }

        const auto index =
            static_cast<std::size_t>(std::distance(headers.begin(), it));

        std::vector<double> values;
        values.reserve(data.size());
        for (const auto& row : data) {
            values.push_back(row.at(index));
        }
        return values;
    }
}
