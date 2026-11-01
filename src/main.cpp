#include "config.h"
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[])
{
    try
    {
        std::string parameter_file;

        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];

            if (arg == "-par")
            {
                if (i + 1 >= argc)
                    throw std::runtime_error("-par requires a file path");

                if (!parameter_file.empty())
                    throw std::runtime_error("Repeated -par option");

                parameter_file = argv[++i];
            }
            else
            {
                throw std::runtime_error("Unknown argument: " + arg);
            }
        }

        if (parameter_file.empty())
        {
            std::cout<<"Parameter File -par <parameterfile> not given, using default parameters.txt instead";
            parameter_file = "parameters.txt";
        }
        const auto cfg = config::ReadParameters(parameter_file);

        // 后续用 cfg 创建网格、初始化模拟……
        std::cout << "Loaded parameters from: "
                  << parameter_file << '\n';

        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}