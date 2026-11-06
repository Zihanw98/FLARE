#include "user_config.h"
#include "user/user_api.h"
#include "source/source.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

static_assert(hd_core::NDIM == 1, "User heating example requires 1D");

namespace
{
struct HeatingParams
{
    double specific_heating_rate = 0.5;
    double max_fractional_energy_change = 0.1;
};

// This parser belongs to the problem, not SimulationConfig. Parse once when
// registering the source, then capture the validated parameters by value.
HeatingParams read_heating_params()
{
    const std::string filename = "heating.txt";
    std::ifstream input(filename);
    if (!input)
        throw std::runtime_error("Cannot open user source parameters: " + filename);
    HeatingParams params;
    std::unordered_set<std::string> seen;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line))
    {
        ++line_number;
        line = config::trim(line.substr(0, line.find('#')));
        if (line.empty()) continue;
        const std::string location = filename + ":" + std::to_string(line_number) + ": ";
        const auto equal = line.find('=');
        if (equal == std::string::npos)
            throw std::runtime_error(location + "Expected key=value");
        const auto key = config::trim(line.substr(0, equal));
        if (!seen.insert(key).second)
            throw std::runtime_error(location + "Duplicate parameter: " + key);
        std::istringstream value(line.substr(equal + 1));
        double number;
        if (!(value >> number) || !std::isfinite(number))
            throw std::runtime_error(location + "Expected a finite number for " + key);
        value >> std::ws;
        if (!value.eof())
            throw std::runtime_error(location + "Unexpected trailing text for " + key);
        if (key == "specific_heating_rate")
            params.specific_heating_rate = number;
        else if (key == "max_fractional_energy_change")
            params.max_fractional_energy_change = number;
        else
            throw std::runtime_error(location + "Unknown parameter: " + key);
    }
    if (input.bad())
        throw std::runtime_error("Failed to read " + filename);
    if (params.specific_heating_rate < 0.0)
        throw std::runtime_error("specific_heating_rate must be nonnegative");
    if (params.max_fractional_energy_change <= 0.0
        || params.max_fractional_energy_change > 1.0)
        throw std::runtime_error("max_fractional_energy_change must be in (0, 1]");
    return params;
}

void initialize(hd_core::MeshBlock& mesh, const config::SimulationConfig& cfg)
{
    if (mesh.geometry != hd_core::Geometry::Cartesian || cfg.constant_gravity_enabled)
        throw std::invalid_argument("Uniform heating requires Cartesian geometry without gravity");
    if (mesh.bc_lower[0] != boundary::BoundaryCondition::Periodic
        || mesh.bc_upper[0] != boundary::BoundaryCondition::Periodic)
        throw std::invalid_argument("Uniform heating requires periodic boundaries");
    hd_core::Primitives state{};
    state[hd_core::RHO] = 1.0;
    state[hd_core::PRE] = 1.0;
    for (int i = 0; i < mesh.nactive[0]; ++i)
        mesh.at_active({i}) = hd_core::Cell(state, cfg.gamma);
}

void heating(const hd_core::MeshBlock& mesh, const HeatingParams& params,
             std::vector<hd_core::Conserved>& rhs, double& dt_limit)
{
    if (params.specific_heating_rate == 0.0) return;
    for (int i = 0; i < mesh.nactive[0]; ++i)
    {
        const auto flat = mesh.flat_index({i}, true);
        const auto& cell = mesh.cells[flat];
        const auto primitive = cell.U.to_Primtivies(cell.gamma);
        const double internal_energy = primitive[hd_core::PRE] / (cell.gamma - 1.0);
        const double energy_rate = cell.U[hd_core::RHO] * params.specific_heating_rate;
        rhs[flat][hd_core::ENERGY] += energy_rate; // A rate: do not multiply by dt.
        dt_limit = std::min(dt_limit,
            params.max_fractional_energy_change * internal_energy / energy_rate);
    }
}

hd_source::SourceTerm make_heating(HeatingParams params)
{
    return {"user_heating", hd_source::SourceKind::Thermal,
        [params](const hd_core::MeshBlock& mesh, double /*time*/,
                 const config::SimulationConfig& /*cfg*/,
                 std::vector<hd_core::Conserved>& rhs, double& dt_limit) {
            heating(mesh, params, rhs, dt_limit);
        }};
}

void register_sources(const config::SimulationConfig&, const hd_core::MeshBlock&,
                      hd_source::SourceList& sources)
{
    sources.push_back(make_heating(read_heating_params()));
}
}

hd_user::UserHooks hd_user::make_user_hooks()
{
    UserHooks hooks;
    hooks.initialize = initialize;
    hooks.register_sources = register_sources;
    return hooks;
}
