/**
 * @file config.h
 * @brief Simulation configuration.
 */
#pragma once
#include "core.h"
#include "mesh.h"
#include "boundary.h"
#include <array>
#include <cstddef>
#include <string>
#include <unordered_set>
namespace config
{
    /** @brief Runtime mesh, numerical-method, source, and output settings. */
    struct SimulationConfig
    {
        hd_core::Index nactive{};
        hd_core::Index nghost{};

        //Currenly, we use an uniform grid
        //Next step, we will try to implement different grid, including user-defined

        std::array<double,hd_core::NDIM> lower{};
        std::array<double,hd_core::NDIM> upper{};
        double t_end = 1.0;
        // 保留原 cfl 成员作为目标值；参数文件也支持 cfl_desired 这一别名。
        double cfl = 0.4;
        double cfl_max = 0.5;
        std::size_t maxsteps = 100000;

        //Numerical method settings
        std::string riemann_solver = "exact";
        std::string spatial_scheme = "godunov";
        std::string time_integrator = "euler";
        std::string reconstructor = "PCM";
        // Slope limiter for PLM; PCM does not use a limiter.
        std::string limiter = "minmod";
        double limiter_beta = 1.5; // Used by mcbeta only; valid range [1, 2].

        double riemann_tolerance = 1e-8;
        int riemann_max_iterations = 100;

        //Default gas adiabatic index
        double gamma=5.0/3.0;
        
        //Boundary conditions, by default it's periodic
        std::array<std::string, hd_core::NDIM> boundary_lower{};
        std::array<std::string, hd_core::NDIM> boundary_upper{};

        //Output
        std::string output_directory = "output";
        std::string output_suffix = "test";
        // 每隔多少个完成的时间步输出一次；初始、第一步和末步额外输出。
        std::size_t output_step_every = 10;

        bool constant_gravity_enabled = false;

        std::array<double,hd_core::NDIM> gravity{};

        bool constant_heating_enabled = false;
        double heating_rate = 0.0; ///< Energy added per unit mass per unit time.

        int integrator_retry_max = 20;

        std::string geometry = "Cartesian";
        bool spherical_gravity_enabled = false;
        bool spherical_self_gravity_enabled = false;
        double spherical_central_mass = 0.0;

    };
    /**
     * @brief Read and validate a key=value simulation parameter file.
     * @param[in] filename Path to the parameter file.
     * @return Configuration populated from the file and built-in defaults.
     * @throws std::runtime_error If opening, reading, parsing, or configuration validation fails.
     * @note Strips # comments and ignores unknown keys. Repeated recognized keys overwrite earlier values; cfl and cfl_desired cannot both be specified.
     */
    SimulationConfig ReadParameters(const std::string& filename);
    /**
     * @brief Remove leading and trailing spaces, tabs, and line endings.
     * @param[in] text Input text; internal whitespace is preserved.
     * @return The trimmed string, or an empty string for whitespace-only input.
     */
    std::string trim(const std::string& text);
    /**
     * @brief Parse one stream-extractable value and reject trailing non-whitespace text.
     * @param[in] text Text containing exactly one value of type T.
     * @return The parsed value.
     * @tparam T Type supporting stream extraction and default construction.
     * @throws std::runtime_error If extraction fails or trailing text remains.
     * @note The generic definition lives in config.cpp; arbitrary external instantiations are not exposed by this header.
     */
    template<typename T>T parse_value(const std::string& text);
    /**
     * @brief Return a string parameter without interpreting quotes or escapes.
     * @param[in] text String value; any internal or surrounding spaces are preserved here.
     * @return A copy of text.
     */
    template<> std::string parse_value<std::string>(const std::string& text);
    /**
     * @brief Parse a boolean parameter after trimming whitespace.
     * @param[in] text Case-insensitive true/false or numeric 1/0.
     * @return The corresponding boolean value.
     * @throws std::runtime_error If the value is not recognized.
     */
    template<> bool parse_value<bool>(const std::string& text);
    /**
     * @brief Check required mesh parameters and supported numerical value ranges.
     * @param[in] cfg Configuration to inspect without modifying it.
     * @param[in] seen Recognized parameter names present in the input file.
     * @throws std::runtime_error If required keys or checked configuration constraints are violated.
     * @throws std::invalid_argument If the limiter name is unknown.
     * @note Requires 0 < cfl <= cfl_max <= 0.5. This is not a complete check of every method name; some names are validated later by their factories.
     */
    void validate_cfg(const SimulationConfig& cfg,const std::unordered_set<std::string>& seen);
}
