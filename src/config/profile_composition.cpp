#include "profile_composition.h"
#include "microphysics/composition/species.h"
#include "utils/utils.h"
#include <cmath>
#include <optional>
#include <stdexcept>

namespace config
{
CompositionColumns profile_composition_columns(const Table& table,
    CompositionColumnFormat format, const std::string& context)
{
    CompositionColumns columns;
    if constexpr (hd_core::NCOMP == 0) return columns;
    const auto& species = hd_microphysics::species_table.species_active;
    if (species.size() != static_cast<std::size_t>(hd_core::NCOMP))
        throw std::runtime_error(context + ": active species count must equal HD_NCOMP; initialize the species table first");

    for (std::size_t k = 0; k < species.size(); ++k)
    {
        const auto& sp = species[k];
        if (format == CompositionColumnFormat::MesaIsotopes && sp.charge != 0)
            continue;
        std::optional<std::size_t> column;
        std::optional<std::size_t> proton_column;
        for (std::size_t c = 0; c < table.headers.size(); ++c)
        {
            const auto& name = table.headers[c];
            const bool matches = format == CompositionColumnFormat::MesaIsotopes
                ? utils::iequals(name, sp.name_base)
                : name == sp.name || (sp.charge == 0 && utils::iequals(name, sp.name_base));
            if (matches)
            {
                if (column)
                    throw std::runtime_error(context + ": duplicate or ambiguous composition columns for " + sp.name);
                column = c;
            }
            // MESA nuclear networks may store h1 and prot separately. Both belong
            // to (Z,A)=(1,1); their total starts in HI until Saha is applied.
            if (format == CompositionColumnFormat::MesaIsotopes && sp.Z == 1 && sp.A == 1 &&
                utils::iequals(name, "prot"))
            {
                if (proton_column)
                    throw std::runtime_error(context + ": duplicate composition column prot");
                proton_column = c;
            }
        }
        if (!column && !proton_column && sp.charge == 0)
            throw std::runtime_error(context + ": missing neutral isotope column " + sp.name_base + " (" + sp.name + ")");
        if (column) columns.emplace_back(*column, static_cast<int>(k));
        if (proton_column) columns.emplace_back(*proton_column, static_cast<int>(k));
    }
    return columns;
}

void set_profile_composition(const std::vector<double>& row,
    const CompositionColumns& columns, hd_core::Primitives& state,
    const std::string& context)
{
    if constexpr (hd_core::NCOMP == 0) return;
    for (int k = 0; k < hd_core::NCOMP; ++k)
        state[hd_core::COMP(k)] = 0.0;
    double sum = 0.0;
    for (const auto& entry : columns)
    {
        if (entry.first >= row.size())
            throw std::runtime_error(context + ": missing composition value");
        const double fraction = row[entry.first];
        if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0)
            throw std::runtime_error(context + ": mass fraction must be finite and in [0, 1]");
        state[hd_core::COMP(entry.second)] += fraction;
        sum += fraction;
    }
    for (int k = 0; k < hd_core::NCOMP; ++k)
        if (state[hd_core::COMP(k)] > 1.0)
            throw std::runtime_error(context + ": combined isotope mass fraction exceeds one");
    if (!std::isfinite(sum) || std::abs(sum - 1.0) >= 1e-6)
        throw std::runtime_error(context + ": mass fractions sum to " + std::to_string(sum)
            + "; expected one across the active species (check missing isotopes)");
}
}
