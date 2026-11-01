#include "initialize.h"
#include "globalvar.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace
{
std::string trim(const std::string& text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return "";
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

template<typename T>
T parse_value(const std::string& text)
{
    std::istringstream input(text);
    T value{};
    if (!(input >> value))
        throw std::runtime_error("Cannot parse value: " + text);
    input >> std::ws;
    if (!input.eof())
        throw std::runtime_error("Unexpected trailing text: " + text);
    return value;
}

template<>
std::string parse_value<std::string>(const std::string& text)
{
    return text; // Keep internal spaces; no quoting/escaping syntax.
}
} // namespace

void ReadParameters(const std::string& filename)
{
    using Setter = std::function<void(const std::string&)>;
    // Only these keys are accepted. The key determines the value's type.
    const std::unordered_map<std::string, Setter> setters = {
        {"nxa", [](const std::string& v) { nxa = parse_value<int>(v); }},
        {"nxg", [](const std::string& v) { nxg = parse_value<int>(v); }},
        {"Xa", [](const std::string& v) { Xa = parse_value<double>(v); }},
        {"tmax", [](const std::string& v) { tmax = parse_value<double>(v); }},
        {"output_every", [](const std::string& v) { output_every = parse_value<int>(v); }},
        {"tout", [](const std::string& v) { tout = parse_value<int>(v); }},
        {"boundary_conditions", [](const std::string& v) {
            boundary_conditions = parse_value<std::string>(v);
        }},
        {"output_name", [](const std::string& v) {
            output_name = parse_value<std::string>(v);
        }}
    };

    std::ifstream file(filename);
    if (!file)
        throw std::runtime_error("Cannot open parameter file: " + filename);

    std::string line;
    int line_number = 0;
    while (std::getline(file, line))
    {
        ++line_number;
        const auto comment = line.find('#');
        if (comment != std::string::npos)
            line.erase(comment);
        line = trim(line);
        if (line.empty())
            continue;

        const std::string location = filename + ":" + std::to_string(line_number) + ": ";
        const auto equal = line.find('=');
        if (equal == std::string::npos)
            throw std::runtime_error(location + "Missing '='");

        const std::pair<std::string, std::string> kv = {
            trim(line.substr(0, equal)), trim(line.substr(equal + 1))
        };
        const auto& [key, value] = kv;
        if (key.empty())
            throw std::runtime_error(location + "Empty parameter name");
        const auto it = setters.find(key);
        if (it == setters.end())
            continue; // Ignore unknown keys without parsing their values.
        try
        {
            it->second(value);
        }
        catch (const std::exception& error)
        {
            throw std::runtime_error(location + "Parameter '" + key + "': " + error.what());
        }
    }
    if (file.bad())
        throw std::runtime_error("Failed to read parameter file: " + filename);
}

void TempInit(const std::string& parameter_file)
{
    // Temporary 1-D Sod shock tube on the active interval [0, Xa].
    Xa = 100.0;
    nxa = 100;
    nxg = 1;                   // Ghost-cell count on EACH side.
    ny = nya = 1;
    nyg = 0;

    tmax = 20.0;               // Sod time 0.2 scaled to a length-100 tube.
    dt = 0.0;                 // Computed from the CFL condition each step.
    tout = 10;
    output_every = 1;
    boundary_conditions = "zero_gradient";
    output_name = "test";

    // Defaults -> file overrides -> validation -> derived sizes -> initial state.
    if (!parameter_file.empty())
        ReadParameters(parameter_file);
    if (nxa <= 0 || nxg < 1 || nxg > (std::numeric_limits<int>::max() - nxa) / 2)
        throw std::invalid_argument("Require nxa > 0, nxg >= 1 and nx within int range");
    if (!std::isfinite(Xa) || Xa <= 0.0 || !std::isfinite(tmax) || tmax < 0.0)
        throw std::invalid_argument("Require finite Xa > 0 and finite tmax >= 0");
    if (output_every <= 0)
        throw std::invalid_argument("output_every must be positive");
    nx = nxa + 2 * nxg;
    dxa = Xa / nxa;
    if (!std::isfinite(dxa) || dxa <= 0.0)
        throw std::invalid_argument("Invalid derived cell width dxa");

    const double gamma = 1.4;
    Domain.resize(nx);
    const Grid left(1.0, 0.0, 1.0, gamma);
    const Grid mid(0.5,0.0,0.2);
    const Grid right(0.125, 0.0, 0.1, gamma);
    for (int i = nxg; i < nxg + nxa; ++i)
    {
        const double x = (i - nxg + 0.5) * dxa;
        Domain[i] = (x<0.5*Xa) ? left : right;
        //Domain[i] = (x < 0.33 * Xa) ? left : (x<0.67*Xa? mid : right);
    }
    // const Grid gmoving(1.0,1.0,1.0,gamma);
    // const Grid gstatic(1.0,0.0,1.0,gamma);
    // for(int i=nxg; i< nxg+nxa; ++i)
    // {
    //     const double x = (i - nxg ) * dxa;
    //     Domain[i] = (x > 0.4 * Xa && x<0.6 *Xa) ? gmoving : gstatic;
    // }
    ApplyBoundaryConditions();
}

void ApplyBoundaryConditions()
{
    // Only list boundary conditions that have an implementation below.
    static const std::array<std::string, 6> supported_conditions = {
        "zero_gradient","periodic","symmetric","antisymmetric","inflow","outflow"
    };
    const auto found = std::find(supported_conditions.begin(),
                                 supported_conditions.end(), boundary_conditions);
    if (found == supported_conditions.end())
    {
        std::string allowed;
        for (const auto& name : supported_conditions)
        {
            if (!allowed.empty())
                allowed += ", ";
            allowed += name;
        }
        throw std::invalid_argument("Unsupported boundary condition '" +
                                    boundary_conditions + "'. Supported: " + allowed);
    }
    if (nxa <= 0 || nxg < 1 || nx != nxa + 2 * nxg ||
        Domain.size() != static_cast<std::size_t>(nx))
        throw std::runtime_error("Invalid 1-D domain layout");

    if (boundary_conditions == "zero_gradient")
    {
        // Copy the whole cell, including gamma, from the nearest active cell.
        for (int j = 0; j < nxg; ++j)
        {
            Domain[j] = Domain[nxg];
            Domain[nxg + nxa + j] = Domain[nxg + nxa - 1];
        }
        return;
    }
    if(boundary_conditions == "periodic")
    {
        for(int j=0; j<nxg; j++)
        {
            Domain[j] = Domain[nxa+j];
        }
        for(int j=nxa+nxg; j<nx; j++)
        {
            Domain[j]= Domain[j-nxa];
        }
        return;
    }
    throw std::logic_error("Listed boundary condition has no implementation");
}
