#include "config.h"
#include "utils.h"
#include <array>
#include <string>         // std::string、std::getline、std::to_string
#include <sstream>        // std::istringstream
#include <istream>        // std::ws
#include <fstream>        // std::ifstream
#include <functional>     // std::function
#include <unordered_map>  // std::unordered_map
#include <utility>        // std::pair
#include <stdexcept>      // std::runtime_error
#include <exception>      // std::exception
#include <unordered_set> // std::unordered_set

namespace config
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
    template<>
    bool parse_value<bool>(const std::string& text)
    {
        const auto value = trim(text);

        if (utils::iequals(value, "true") || value == "1")
            return true;

        if (utils::iequals(value, "false") || value == "0")
            return false;

        throw std::runtime_error(
            "Expected true, false, 1 or 0; got: " + text
        );
    }
    SimulationConfig ReadParameters(const std::string& filename)
    {
        using Setter = std::function<void(const std::string&)>;
        // Only these keys are accepted. The key determines the value's type.
        SimulationConfig cfg{};

        using Setter = std::function<void(const std::string&)>;
        for(int d=0;d<hd_core::NDIM;d++)
        {
            cfg.boundary_lower[d] = "periodic";
            cfg.boundary_upper[d] = "periodic";
        }
        std::unordered_set<std::string> seen;
        const Setter set_output_step_every = [&cfg](const std::string& v) {
            const int interval = parse_value<int>(v);
            if (interval <= 0)
                throw std::runtime_error("output_step_every must be a positive integer");
            cfg.output_step_every = static_cast<std::size_t>(interval);
        };
        std::unordered_map<std::string, Setter> setters = {
            {"t_end", [&cfg](const std::string& v) {
                cfg.t_end = parse_value<double>(v);
            }},
            {"cfl", [&cfg](const std::string& v) {
                cfg.cfl = parse_value<double>(v);
            }},
            {"maxsteps", [&cfg](const std::string& v) {
                // 无符号流提取也可能接受 -1，因此必须在转换前拒绝负号。
                if (v.empty() || v.front() == '-')
                    throw std::runtime_error("maxsteps must be a positive integer");
                cfg.maxsteps = parse_value<std::size_t>(v);
            }},
            {"gamma", [&cfg](const std::string& v) {
                cfg.gamma = parse_value<double>(v);
            }},
            {"riemann_solver", [&cfg](const std::string& v) {
                cfg.riemann_solver = v;
            }},
            {"riemann_tolerance", [&cfg](const std::string& v) {
                cfg.riemann_tolerance = parse_value<double>(v);
            }},
            {"riemann_max_iterations", [&cfg](const std::string& v) {
                cfg.riemann_max_iterations = parse_value<int>(v);
            }},
            {"spatial_scheme", [&cfg](const std::string& v) {
                cfg.spatial_scheme = v;
            }},
            {"time_integrator", [&cfg](const std::string& v) {
                cfg.time_integrator = v;
            }},
            {"output_directory", [&cfg](const std::string& v) {
                cfg.output_directory = v;
            }},
            {"output_suffix", [&cfg](const std::string& v) {
                cfg.output_suffix = v;
            }},
            {"step_every", set_output_step_every},
            // 兼容原先配置字段的名字。
            {"output_step_every", set_output_step_every},
            {"constant_gravity_enabled",[&cfg](const std::string& v){
                cfg.constant_gravity_enabled = parse_value<bool>(v);
            }},
        };

        constexpr std::array<const char*, 3> axes{"x", "y", "z"};
        for (int d = 0; d < hd_core::NDIM; ++d)
        {
            const std::string axis = axes[d];

            setters.emplace("n" + axis + "a",
                [&cfg, d](const std::string& v) {
                    cfg.nactive[d] = parse_value<int>(v);
                });

            setters.emplace("n" + axis + "g",
                [&cfg, d](const std::string& v) {
                    cfg.nghost[d] = parse_value<int>(v);
                });

            setters.emplace(axis + "_lower",
                [&cfg, d](const std::string& v) {
                    cfg.lower[d] = parse_value<double>(v);
                });

            setters.emplace(axis + "_upper",
                [&cfg, d](const std::string& v) {
                    cfg.upper[d] = parse_value<double>(v);
                });

            setters.emplace("boundary_" + axis + "_lower",
                [&cfg, d](const std::string& v) {
                    cfg.boundary_lower[d] = v;
                });

            setters.emplace("boundary_" + axis + "_upper",
                [&cfg, d](const std::string& v) {
                    cfg.boundary_upper[d] = v;
                });
            setters.emplace("gravity_" + axis,
                [&cfg,d](const std::string& v)
                {
                    cfg.gravity[d] = parse_value<double>(v);
                }
            );
        }
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
                seen.insert(key);
            }
            catch (const std::exception& error)
            {
                throw std::runtime_error(location + "Parameter '" + key + "': " + error.what());
            }
        }
        if (file.bad())
            throw std::runtime_error("Failed to read parameter file: " + filename);
        validate_cfg(cfg,seen);
        return cfg;
    }
    void validate_cfg(const SimulationConfig& cfg, const std::unordered_set<std::string>& seen)
    {
        const auto require = [&seen](const std::string& key)
        {
            if (seen.count(key) == 0)
                throw std::runtime_error(
                    "Missing required parameter: " + key
                );
        };
        constexpr std::array<const char*, 3> axes{"x", "y", "z"};
        for(int d=0; d< hd_core::NDIM; d++)
        {
            const std::string axis = axes[d];
            require("n" + axis + "a");       //active zone required
            require("n" + axis + "g");
            require(axis + "_upper" );
            require(axis + "_lower");
            if (cfg.nactive[d] <= 0)
                throw std::runtime_error(
                    "n" + axis + "a must be > 0"
                );

            if (cfg.nghost[d] < 1)
                throw std::runtime_error(
                    "n" + axis + "g must be >= 1"
                );

            if (!std::isfinite(cfg.lower[d]) ||
                !std::isfinite(cfg.upper[d]) ||
                cfg.upper[d] <= cfg.lower[d])
            {
                throw std::runtime_error(
                    axis + " bounds must be finite and upper > lower"
                );
            }
            const bool& bc_upper_periodic = utils::iequals(cfg.boundary_upper[d],"periodic");
            const bool& bc_lower_periodic = utils::iequals(cfg.boundary_lower[d],"periodic");
            if(bc_upper_periodic!=bc_lower_periodic)
            {
                throw std::runtime_error("Periodic boundary condition should be double-sided in axis" + d);
            }
            if(!std::isfinite(cfg.gravity[d]))
                throw std::runtime_error("gravity not finite");
        }
        if (!std::isfinite(cfg.t_end) || cfg.t_end <= 0.0)
            throw std::runtime_error("t_end must be finite and > 0");

        if(cfg.time_integrator == "euler")
        {
            if (!std::isfinite(cfg.cfl) ||
                cfg.cfl <= 0.0 || cfg.cfl > 0.5)
            {
                throw std::runtime_error("cfl must satisfy 0 < cfl <= 0.5");
            }
        }

        if (!std::isfinite(cfg.gamma) || cfg.gamma <= 1.0)
            throw std::runtime_error("gamma must be finite and > 1");

        if (!std::isfinite(cfg.riemann_tolerance) || cfg.riemann_tolerance <= 0.0)
            throw std::runtime_error("riemann_tolerance must be finite and > 0");
        if (cfg.riemann_max_iterations <= 0)
            throw std::runtime_error("riemann_max_iterations must be > 0");

        if (cfg.maxsteps <= 0)
            throw std::runtime_error("maxsteps must be > 0");
    }
    hd_core::MeshBlock build_mesh(const SimulationConfig& cfg)
    {
        hd_core::MeshBlock mesh{};
        mesh.nactive = cfg.nactive;
        mesh.nghost = cfg.nghost;
        for(int d=0;d<hd_core::NDIM;d++)
        {
            mesh.set_uniform_faces(d,cfg.lower[d],cfg.upper[d]);
            mesh.bc_upper[d] = boundary::parse_bc(cfg.boundary_upper[d]);
            mesh.bc_lower[d] = boundary::parse_bc(cfg.boundary_lower[d]);
        }
        mesh.allocate();
        return mesh;
    }
}
