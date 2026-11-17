#pragma once
#include "core/core.h"
#include "table.h"
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace config
{
    enum class CompositionColumnFormat { Species, MesaIsotopes };
    using CompositionColumns = std::vector<std::pair<std::size_t, int>>;

    // Resolve columns once against the initialized active species table.
    // Species format requires neutral names (or name_base aliases), with optional ions.
    // MESA format requires name_base columns and deposits isotope totals in charge 0.
    CompositionColumns profile_composition_columns(const Table& table,
        CompositionColumnFormat format, const std::string& context);

    // Set omitted states to zero; reject nonfinite/out-of-range values and |sum(X)-1| >= 1e-6.
    // Does not normalize or solve ionization equilibrium. Context should include the row.
    void set_profile_composition(const std::vector<double>& row,
        const CompositionColumns& columns, hd_core::Primitives& state,
        const std::string& context);
}
