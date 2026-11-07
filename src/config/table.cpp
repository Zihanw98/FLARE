#include "table.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>

namespace config
{
    Table read_table(const std::string& filename,
                    std::size_t skip_rows = 0,
                    bool has_header = true)
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
            double value;

            while (fields >> value) {
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
}