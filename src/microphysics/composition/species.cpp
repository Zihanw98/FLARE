#include "species.h"
#include "config/config.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace hd_microphysics
{
namespace
{
    std::array<int, 3> species_key(const Species& sp)
    {
        return {sp.Z, sp.A - sp.Z, sp.charge};
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

std::vector<Species> read_species_table(const std::string& filename)
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

        std::vector<Species> catalog(rows.size());
        std::unordered_map<std::string, int> by_name;
        for (std::size_t k = 0; k < rows.size(); ++k)
        {
            location = path.string() + ": species[" + std::to_string(k) + "]";
            const auto& row = rows[k];
            auto& sp = catalog[k];
            row.at("name").get_to(sp.name);
            row.at("name_base").get_to(sp.name_base);
            sp.A = read_integer(row, "A", location);
            sp.Z = read_integer(row, "Z", location);
            sp.charge = read_integer(row, "charge", location);
            row.at("energy_eV").get_to(sp.energy_eV);
            if (sp.name.empty() || sp.name == "-" || sp.name_base.empty())
                throw std::runtime_error(location + ": invalid species name or base");
            if (!by_name.emplace(sp.name, static_cast<int>(k)).second)
                throw std::runtime_error(location + ": duplicate species " + sp.name);
            if (sp.A < 1 || sp.Z < 0 || sp.Z > sp.A || sp.charge > sp.Z || !std::isfinite(sp.energy_eV))
                throw std::runtime_error(location + ": invalid A, Z, charge or energy_eV for " + sp.name);
        }

        // Names can refer to later rows. Resolve links after all catalog IDs exist.
        for (std::size_t k = 0; k < rows.size(); ++k)
        {
            auto& sp = catalog[k];
            location = path.string() + ": " + sp.name;
            auto neighbor_id = [&](const char* key, const char* reverse, int step) {
                const auto& link = rows[k].at(key);
                if (link.is_null())
                    return NO_IONIZATION_NEIGHBOR;
                const auto& name = link.get_ref<const std::string&>();
                const auto it = by_name.find(name);
                if (it == by_name.end())
                    throw std::runtime_error(location + ": unknown neighbor " + name);
                const int id = it->second;
                const auto& other = catalog[id];
                if (other.A != sp.A || other.Z != sp.Z || other.name_base != sp.name_base ||
                    other.charge != static_cast<long long>(sp.charge) + step || rows[id].at(reverse) != sp.name)
                    throw std::runtime_error(location + ": invalid ionization link to " + name);
                return id;
            };
            sp.ionization_lower = neighbor_id("ionization_lower", "ionization_upper", -1);
            sp.ionization_upper = neighbor_id("ionization_upper", "ionization_lower", +1);
        }
        return catalog;
    }
    catch (const nlohmann::json::exception& error)
    {
        throw std::runtime_error(location + ": " + error.what());
    }
}

std::vector<Species> select_species(const std::vector<Species>& catalog,
                                    const std::vector<std::string>& active_names)
{
    std::unordered_map<std::string, std::size_t> by_name;
    for (std::size_t k = 0; k < catalog.size(); ++k)
        by_name.emplace(catalog[k].name, k);

    std::vector<int> active_ids(catalog.size(), NO_IONIZATION_NEIGHBOR);
    std::vector<Species> active;
    for (const auto& name : active_names)
    {
        const auto it = by_name.find(name);
        if (it == by_name.end())
            throw std::runtime_error("Unknown active species: " + name);
        auto& id = active_ids[it->second];
        if (id != NO_IONIZATION_NEIGHBOR)
            throw std::runtime_error("Duplicate active species: " + name);
        id = static_cast<int>(active.size());
        active.push_back(catalog[it->second]);
    }
    for (auto& sp : active)
    {
        if (sp.ionization_lower != NO_IONIZATION_NEIGHBOR)
            sp.ionization_lower = active_ids.at(sp.ionization_lower);
        if (sp.ionization_upper != NO_IONIZATION_NEIGHBOR)
            sp.ionization_upper = active_ids.at(sp.ionization_upper);
    }
    return active;
}

std::vector<Species> read_species_table(const std::string& filename,
                                        const std::vector<std::string>& active_names)
{
    return select_species(read_species_table(filename), active_names);
}

SpeciesList::SpeciesList(std::vector<Species> species)
    : species_(std::move(species)), sorted_ids_(species_.size())
{
    if (species_.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::length_error("Too many species for integer IDs");
    std::iota(sorted_ids_.begin(), sorted_ids_.end(), 0);
    std::sort(sorted_ids_.begin(), sorted_ids_.end(), [&](int a, int b) {
        return species_key(species_[a]) < species_key(species_[b]);
    });
    for (std::size_t i = 1; i < sorted_ids_.size(); ++i)
    {
        const auto& a = species_[sorted_ids_[i - 1]];
        const auto& b = species_[sorted_ids_[i]];
        if (species_key(a) == species_key(b))
            throw std::invalid_argument("Same proton, neutron and charge numbers for " + a.name + " and " + b.name);
    }
}

std::optional<int> SpeciesList::find_species(int protons, int neutrons, int charge) const
{
    const std::array<int, 3> key{protons, neutrons, charge};
    const auto it = std::lower_bound(sorted_ids_.begin(), sorted_ids_.end(), key,
        [&](int id, const auto& target) { return species_key(species_[id]) < target; });
    if (it == sorted_ids_.end() || species_key(species_[*it]) != key)
        return std::nullopt;
    return *it;
}

SpeciesList read_active_species(const std::string& filename, const std::vector<Species>& catalog)
{
    std::ifstream file(filename);
    if (!file)
        throw std::runtime_error("Cannot open composition activation file: " + filename);
    std::unordered_map<std::string, std::size_t> neutral;
    for (std::size_t k = 0; k < catalog.size(); ++k)
    {
        const auto& sp = catalog[k];
        if (sp.charge == 0 && !neutral.emplace(sp.name_base, k).second)
            throw std::runtime_error("Ambiguous neutral species for " + sp.name_base);
    }

    std::unordered_set<std::string> seen;
    std::vector<std::string> names;
    std::string line;
    for (std::size_t line_number = 1; std::getline(file, line); ++line_number)
    {
        std::istringstream fields(line.substr(0, line.find('#')));
        std::string base, lower, upper, extra;
        if (!(fields >> base))
            continue;
        const std::string location = filename + ":" + std::to_string(line_number) + ": ";
        int q_min, q_max;
        auto parse_charge = [](const std::string& token, int& q) {
            std::istringstream input(token);
            return (input >> q) && input.eof();
        };
        if (!(fields >> lower >> upper) || (fields >> extra) ||
            !parse_charge(lower, q_min) || !parse_charge(upper, q_max))
            throw std::runtime_error(location + "expected name_base q_min q_max with integer charges");
        if (q_min > 0 || q_max < 0)
            throw std::runtime_error(location + base + " requires q_min <= 0 <= q_max");
        const auto it = neutral.find(base);
        if (it == neutral.end())
            throw std::runtime_error(location + "no neutral species in catalog for " + base);
        if (!seen.insert(base).second)
            throw std::runtime_error(location + "duplicate nuclide " + base);

        // The catalog links were validated as reciprocal, adjacent charge states.
        const Species* sp = &catalog[it->second];
        while (sp->charge > q_min)
        {
            if (sp->ionization_lower == NO_IONIZATION_NEIGHBOR)
                throw std::runtime_error(location + base + " has no lower state below charge " + std::to_string(sp->charge));
            sp = &catalog.at(sp->ionization_lower);
        }
        for (;;)
        {
            names.push_back(sp->name);
            if (sp->charge == q_max)
                break;
            if (sp->ionization_upper == NO_IONIZATION_NEIGHBOR)
                throw std::runtime_error(location + base + " has no upper state above charge " + std::to_string(sp->charge));
            sp = &catalog.at(sp->ionization_upper);
        }
    }
    if (file.bad())
        throw std::runtime_error("Failed to read composition activation file: " + filename);
    return SpeciesList(select_species(catalog, names));
}

SpeciesList load_run_species(const config::SimulationConfig& cfg, const std::string& parameter_file)
{
    if (!cfg.enable_composition)
        return {};

    const auto selection = std::filesystem::path(parameter_file).parent_path() / cfg.composition_activate_file;
    const auto catalog = read_species_table();
    auto active = read_active_species(selection.string(), catalog);
    if (active.size() != static_cast<std::size_t>(hd_core::NCOMP))
        throw std::runtime_error(selection.string() + ": selected " + std::to_string(active.size())
            + " species, but HD_NCOMP is " + std::to_string(hd_core::NCOMP));
    return active;
}
}
