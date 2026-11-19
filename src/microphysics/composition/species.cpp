#include "species.h"
#include "config.h"
#include "table.h"
#include "core.h"
#include "cell.h"
#include "utils.h"
#include "physcon/constants.h"
#include <algorithm>
#include <optional>
#include <nlohmann/json.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace hd_microphysics
{
    namespace
    {
        SpeciesTable simulation_species;

        std::vector<Species>::const_iterator find_active_species(
            const std::vector<Species>& active, int A, int Z, int charge)
        {
            auto it = std::lower_bound(active.begin(), active.end(), Z,
                [](const Species& species, int target_Z) { return species.Z < target_Z; });
            for (; it != active.end() && it->Z == Z; ++it)
            {
                if (it->A == A && it->charge == charge)
                    return it;
            }
            return active.end();
        }

        std::filesystem::path resolve_species_table(const std::filesystem::path& requested)
        {
            if (requested != DEFAULT_SPECIES_TABLE || std::filesystem::exists(requested))
                return requested;

            // Both FLARE makefiles compile framework sources by absolute path.
            // For a manual build with relative source paths, start at the case instead.
            const std::filesystem::path source(__FILE__);
            auto dir = source.is_absolute() ? source.parent_path() : std::filesystem::current_path();
            for (;;)
            {
                const auto candidate = dir / requested;
                if (std::filesystem::is_regular_file(candidate))
                    return candidate;
                if (dir == dir.root_path())
                    return requested;
                dir = dir.parent_path();
            }
        }

        int read_integer(const nlohmann::json& item, const char* key, const std::string& location)
        {
            const auto& value = item.at(key);
            if (!value.is_number_integer() ||
                (value.is_number_unsigned() && value.get<nlohmann::json::number_unsigned_t>() >
                                            static_cast<unsigned>(std::numeric_limits<int>::max())))
                throw std::runtime_error(location + ": " + key + " must be an integer in range");
            const auto n = value.get<nlohmann::json::number_integer_t>();
            if (n < std::numeric_limits<int>::min() || n > std::numeric_limits<int>::max())
                throw std::runtime_error(location + ": " + key + " is outside integer range");
            return static_cast<int>(n);
        }
    }

    /**
     * @brief 共享物种表的只读引用。Read-only reference to the shared species table.
     *
     * @par 中文
     * 这行定义在初始化引用时将 `species_table` 绑定到 `simulation_species`。
     * C++ 引用一旦绑定就不能改绑，`const` 表示不能通过该引用修改表的内容。
     * 模拟启动时，initialize_species_table() 直接给可写对象 `simulation_species` 赋值。
     * 此后通过已经绑定的 `species_table` 就能读到填好的内容。
     *
     * @par English
     * This definition initializes `species_table` by binding it to `simulation_species`.
     * A C++ reference cannot be rebound, and `const` prevents writes through this reference.
     * At startup, initialize_species_table() assigns the table contents to the mutable
     * `simulation_species` object. Reads through the existing reference then see those contents.
     */
    const SpeciesTable& species_table = simulation_species;

    void initialize_species_table(const config::SimulationConfig& cfg)
    {
        SpeciesTable table;
        if (cfg.enable_composition)
        {
            table = read_species_table();
            table.build_active_species_table(cfg);
        }
        simulation_species = std::move(table);
    }

    std::optional<Species> SpeciesTable::find_species_from_full_table(const std::string& name, int charge) const
    {
        for (const auto& species : species_full)
        {
            if (utils::iequals(species.name_base,name) && species.charge == charge)
                return species;
        }
        return std::nullopt;
    }

    std::optional<Species> SpeciesTable::find_species_from_available_table(int A, int Z, int charge) const
    {
        const auto it = find_active_species(species_active, A, Z, charge);
        if (it == species_active.end())
            return std::nullopt;
        return *it;
    }
    std::optional<Species> SpeciesTable::find_species_from_available_table(const std::string& name, int charge) const
    {
        for (const auto& species : species_active)
        {
            if (utils::iequals(species.name_base,name) && species.charge == charge)
                return species;
        }
        return std::nullopt;
    }

    int SpeciesTable::find_species_from_available_table_index(const std::string& name) const
    {
        for (std::size_t i = 0; i < species_active.size(); ++i)
        {
            if (species_active[i].name == name)
                return static_cast<int>(i);
        }
        throw std::out_of_range("Species not found in active table: " + name);
    }

    int SpeciesTable::find_species_from_available_table_index(int A, int Z, int charge) const
    {
        const auto it = find_active_species(species_active, A, Z, charge);
        if (it == species_active.end())
            throw std::out_of_range("Species not found in active table: A=" + std::to_string(A)
                                    + ", Z=" + std::to_string(Z) + ", charge=" + std::to_string(charge));
        return static_cast<int>(it - species_active.begin());
    }

    std::vector<int> SpeciesTable::find_ionization_state_indices(int Z, int A) const
    {
        std::vector<int> indices;
        auto it = std::lower_bound(species_active.begin(), species_active.end(), Z,
            [](const Species& species, int target_Z) { return species.Z < target_Z; });
        for (; it != species_active.end() && it->Z == Z; ++it)
        {
            if (it->A == A)
                indices.push_back(static_cast<int>(it - species_active.begin()));
        }
        std::stable_sort(indices.begin(), indices.end(), [&](int left, int right) {
            return species_active[left].charge < species_active[right].charge;
        });
        return indices;
    }

    SpeciesTable read_species_table(const std::string& filename)
    {
        const auto path = resolve_species_table(filename);
        std::ifstream file(path);
        if (!file)
            throw std::runtime_error("Cannot open species table: " + path.string());

        std::string location = path.string();
        try
        {
            const auto document = nlohmann::json::parse(file);
            const auto& rows = document.at("species");
            if (!rows.is_array() || rows.empty())
                throw std::runtime_error(location + ": species must be a nonempty array");
            if (rows.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw std::length_error(location + ": too many species for integer IDs");

            SpeciesTable table;
            table.species_full.resize(rows.size());
            std::unordered_map<std::string, int> by_name;
            for (std::size_t k = 0; k < rows.size(); ++k)
            {
                location = path.string() + ": species[" + std::to_string(k) + "]";
                const auto& row = rows[k];
                auto& sp = table.species_full[k];
                row.at("name").get_to(sp.name);
                row.at("name_base").get_to(sp.name_base);
                sp.A = read_integer(row, "A", location);
                sp.Z = read_integer(row, "Z", location);
                sp.charge = read_integer(row, "charge", location);
                row.at("energy_eV").get_to(sp.energy_eV);
                sp.ground_weight = row.value("ground_weight", sp.ground_weight);
                if (sp.name.empty() || sp.name == "-" || sp.name_base.empty())
                    throw std::runtime_error(location + ": invalid species name or base");
                if (!by_name.emplace(sp.name, static_cast<int>(k)).second)
                    throw std::runtime_error(location + ": duplicate species " + sp.name);
                if (sp.A < 1 || sp.Z < 0 || sp.Z > sp.A || sp.charge > sp.Z || !std::isfinite(sp.energy_eV))
                    throw std::runtime_error(location + ": invalid A, Z, charge or energy_eV for " + sp.name);
                if (!std::isfinite(sp.ground_weight) || sp.ground_weight <= 0.0)
                    throw std::runtime_error(location + ": ground_weight must be finite and positive for " + sp.name);
            }

            // 所有物种读完后校验邻居名称；两个 map 保留名称，物种的整数索引保持 nullopt。
            for (std::size_t k = 0; k < rows.size(); ++k)
            {
                auto& sp = table.species_full[k];
                location = path.string() + ": " + sp.name;
                auto neighbor_name = [&](const char* key, const char* reverse, int step) -> std::string {
                    const auto& link = rows[k].at(key);
                    if (link.is_null())
                        return {};
                    const auto& name = link.get_ref<const std::string&>();
                    const auto it = by_name.find(name);
                    if (it == by_name.end())
                        throw std::runtime_error(location + ": unknown neighbor " + name);
                    const int id = it->second;
                    const auto& other = table.species_full[id];
                    if (other.A != sp.A || other.Z != sp.Z || other.name_base != sp.name_base ||
                        other.charge != static_cast<long long>(sp.charge) + step || rows[id].at(reverse) != sp.name)
                        throw std::runtime_error(location + ": invalid ionization link to " + name);
                    return name;
                };
                table.ionization_lower_names.emplace(k, neighbor_name("ionization_lower", "ionization_upper", -1));
                table.ionization_upper_names.emplace(k, neighbor_name("ionization_upper", "ionization_lower", +1));
            }
            return table;
        }
        catch (const nlohmann::json::exception& error)
        {
            throw std::runtime_error(location + ": " + error.what());
        }
    }
    void SpeciesTable::build_active_species_table(const config::SimulationConfig& cfg)
    {
        const auto composition_file = cfg.composition_activate_file;
        auto composition_active_list = config::read_composition_table(composition_file);

        for (const auto& comp_row : composition_active_list.data)
        {
            auto name = comp_row.name_base;
            auto lower = comp_row.q_min;
            auto upper = comp_row.q_max;
            if(lower>upper)
            {
                throw std::runtime_error(name + "name: Lower ionization stage is higher than Upper ionization stage");
            }
            if(lower>0 || upper <0)
            {
                throw std::runtime_error(name + ": Ground state must be contained");
            }
            for(auto q = lower; q<=upper; q++)
            {
                auto result = find_species_from_full_table(name,q);
                if(!result.has_value())
                {
                    throw std::runtime_error("No such species of " + name + "with" +  std::to_string(q) + "ionization stage");
                }
                species_active.push_back({*result});
            }
        }
        auto& species = species_active;
        std::stable_sort(species.begin(), species.end(),
            [](const Species& a, const Species& b) {
                return a.Z < b.Z;
            });
        if(HD_NCOMP != species_active.size())
        {
            throw std::runtime_error("HD_NCOMP parameter doesn't equal to the number of species!");
        }

    }
    void validate_cell_composition(const hd_core::Cell& cell)
    {
        if constexpr (hd_core::NCOMP == 0) return;

        const double rho = cell.U[hd_core::RHO];
        if (!std::isfinite(rho) || rho < 0.0)
            throw std::invalid_argument(
                "validate_cell_composition: rho must be finite and nonnegative");

        constexpr double fraction_tolerance = 1e-12;
        constexpr double sum_tolerance = 1e-6;
        double sum = 0.0;
        for (int k = 0; k < hd_core::NCOMP; ++k)
        {
            const double partial_density = cell.U[hd_core::COMP(k)];
            if (!std::isfinite(partial_density))
                throw std::invalid_argument(
                    "validate_cell_composition: nonfinite partial density at component "
                    + std::to_string(k));
            if (rho == 0.0)
            {
                if (partial_density != 0.0)
                    throw std::invalid_argument(
                        "validate_cell_composition: vacuum has nonzero partial density at component "
                        + std::to_string(k));
                continue;
            }

            const double fraction = partial_density / rho;
            if (!std::isfinite(fraction) || fraction < -fraction_tolerance ||
                fraction > 1.0 + fraction_tolerance)
                throw std::invalid_argument(
                    "validate_cell_composition: mass fraction outside [0, 1] at component "
                    + std::to_string(k));
            sum += fraction;
        }
        if (rho > 0.0 && (!std::isfinite(sum) || std::abs(sum - 1.0) >= sum_tolerance))
            throw std::invalid_argument(
                "validate_cell_composition: mass fractions sum to "
                + std::to_string(sum) + "; expected one");
    }

    void evaluate_cell_density(hd_core::Cell& cell)
    {
        const double rho = cell.U[hd_core::RHO];

        if (!std::isfinite(rho) || rho <= 0.0)
        {
            throw std::invalid_argument(
                "evaluate_cell_density: rho must be finite and positive"
            );
        }

        const auto& active = species_table.species_active;

        if (active.empty() ||
            active.size() != static_cast<std::size_t>(hd_core::NCOMP))
        {
            throw std::logic_error(
                "evaluate_cell_density: active species must be nonempty "
                "and match NCOMP"
            );
        }

        double mu_inv=0.0;
        //Calculate mu first for every three EoS left
        double Y_e_free = 0.0;
        double Y_Ion = 0.0;
        double n_baryon = rho / physcon::m_u;
        cell.n_baryon = n_baryon;
        for(int k=0; k<hd_core::NCOMP; k++)
        {
            const auto A_k = hd_microphysics::species_table.species_active[k].A;
            const auto Z_k = hd_microphysics::species_table.species_active[k].charge;
            const auto X_k = cell.U[hd_core::COMP(k)]/rho;
            Y_e_free += Z_k*X_k/A_k;
            mu_inv += (1+Z_k)*X_k/A_k;
            //1/mu = \Sigma (1+Z_k)*X_k/A_k
            Y_Ion += X_k/A_k;
        }
        cell.mean_molecular_weight = 1.0/mu_inv;
        cell.Y_e_free = Y_e_free;
        cell.Y_Ion = Y_Ion;
        cell.n_e_free = n_baryon * Y_e_free;
        cell.n_Ion = n_baryon * Y_Ion;
        return;
    }
}
