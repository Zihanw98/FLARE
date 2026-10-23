#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include "initialize.h"
#include "simulation.h"
int main(int argc, char* argv[])
{
    if (argc > 2)
    {
        std::cerr << "Usage: " << argv[0] << " [parameter_file]\n";
        return 1;
    }
    try
    {
        TempInit(argc == 2 ? argv[1] : "");
        simulation_main();
    }
    catch (const std::exception& error)
    {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
