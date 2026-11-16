#pragma once
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace config
{
    struct SimulationConfig;
}

namespace hd_microphysics
{
    // No tracked species at the adjacent charge state.
    inline constexpr int NO_IONIZATION_NEIGHBOR = std::numeric_limits<int>::min();
    inline constexpr char DEFAULT_SPECIES_TABLE[] = "tables/species.json";

    struct Species
    {
        std::string name;
        std::string name_base;
        int A;
        int Z;
        int charge = 0;
        int ionization_upper = NO_IONIZATION_NEIGHBOR;
        int ionization_lower = NO_IONIZATION_NEIGHBOR;
        double gamma = 5.0/3.0;
        double energy_eV = 0.0;  // Per particle, relative to the neutral isotope's ground state.
        double ground_weight = 1.0;
    };

    // Read the whole JSON catalog in file order. Ionization neighbors index the
    // catalog. The default table can also be read from a problem directory.
    std::vector<Species> read_species_table(const std::string& filename = DEFAULT_SPECIES_TABLE);

    // Select from an already loaded catalog. active_names[k] supplies COMP(k).
    // Neighbors index the returned vector; omitted neighbors use NO_IONIZATION_NEIGHBOR.
    std::vector<Species> select_species(const std::vector<Species>& catalog,
                                        const std::vector<std::string>& active_names);

    // Convenience overload for loading a catalog and immediately selecting species.
    std::vector<Species> read_species_table(const std::string& filename,
                                            const std::vector<std::string>& active_names);

    // Keep the original IDs (e.g. COMP slots) and index them separately for lookup.
    class SpeciesList
    {
    public:
        SpeciesList() = default;
        explicit SpeciesList(std::vector<Species> species);
        std::size_t size() const { return species_.size(); }
        const Species& operator[](std::size_t id) const { return species_.at(id); }
        // O(log n); nullopt if the requested (Z, A-Z, charge) is absent.
        std::optional<int> find_species(int protons, int neutrons, int charge) const;

    private:
        std::vector<Species> species_;
        std::vector<int> sorted_ids_;
    };

    // Each non-comment line contains name_base, q_min and q_max (inclusive, containing zero).
    // IDs follow file order, with each nuclide's charges in ascending order.
    SpeciesList read_active_species(const std::string& filename, const std::vector<Species>& catalog);

    // Load the run's COMP(k) mapping using the activation file relative to the parameter file.
    // Disabled composition returns an empty list without reading either file.
    SpeciesList load_run_species(const config::SimulationConfig& cfg, const std::string& parameter_file);
}
