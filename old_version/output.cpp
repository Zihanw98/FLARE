#include "output.h"
#include "globalvar.h"
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

void Output(double time, int timestep, bool is_last_step)
{
    if (output_every <= 0)
        throw std::invalid_argument("output_every must be positive");
    if (timestep < 0)
        throw std::invalid_argument("timestep must be nonnegative");

    if (timestep > 1 && !is_last_step && timestep % output_every != 0)
        return;

    if (nx <= 0 || Domain.size() != static_cast<std::size_t>(nx))
        throw std::runtime_error("Output: inconsistent grid size");

    const std::string filename = output_name + "_" + std::to_string(timestep) + ".txt";
    std::ofstream file(filename);
    if (!file)
        throw std::runtime_error("Cannot open output file: " + filename);

    file << std::setprecision(std::numeric_limits<double>::max_digits10);
    file << "# time=" << time << " timestep=" << timestep
         << " nx=" << nx << " nxg=" << nxg << " nxa=" << nxa << '\n';

    // All nx cells, including both ghost regions; columns: rho vx p.
    for (int i = 0; i < nx; ++i)
    {
        const State prim = Domain[i].get_Prim();
        file << prim[0] << ' ' << prim[1] << ' ' << prim[2] << '\n';
    }
    file.close();
    if (!file)
        throw std::runtime_error("Failed to write output file: " + filename);
}
