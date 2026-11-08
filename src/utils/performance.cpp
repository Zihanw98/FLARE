#include "performance.h"
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <locale>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace performance
{
Mode mode = Mode::Off;
std::array<std::uint64_t, counter_count> counters{};
namespace
{
std::array<Measurement, section_count> totals{};
Scope* current = nullptr;
constexpr std::array<const char*, section_count> names{{
    "run", "initialization", "evolution", "advance", "rhs", "hydro",
    "face_fluxes", "reconstruction", "riemann", "flux_divergence",
    "sources", "gravity", "geometry_source", "heating", "boundary",
    "state_update", "mesh_copy", "validation", "output",
    "output_snapshot_and_diagnostics", "output_diagnostics_only"
}};
constexpr std::array<const char*, counter_count> counter_names{{
    "accepted_steps", "active_cell_updates", "rk2_attempts", "rk2_rejected_attempts",
    "snapshots"
}};
void require_closed()
{
    if (current) throw std::logic_error("Performance scopes must be closed before reset/report");
}
double wall_seconds() { return totals[static_cast<std::size_t>(Section::Run)].inclusive_seconds; }
}

void configure(Mode requested)
{
    require_closed();
    totals = {};
    counters = {};
    mode = requested;
}
const std::array<Measurement, section_count>& measurements() noexcept { return totals; }
const char* name(Section section) noexcept { return names[static_cast<std::size_t>(section)]; }

void Scope::start() noexcept
{
    parent_ = current;
    current = this;
    start_ = std::chrono::steady_clock::now();
}
void Scope::finish() noexcept
{
    const auto end = std::chrono::steady_clock::now();
    assert(current == this);
    const double elapsed = std::chrono::duration<double>(end - start_).count();
    auto& result = totals[static_cast<std::size_t>(section_)];
    ++result.calls;
    result.inclusive_seconds += elapsed;
    result.self_seconds += std::max(0.0, elapsed - children_seconds_);
    current = parent_;
    if (parent_) parent_->children_seconds_ += elapsed;
    active_ = false;
}

void report(std::ostream& out, bool success)
{
    if (!enabled()) return;
    require_closed();
    const double wall = wall_seconds();
    const double evolution = totals[static_cast<std::size_t>(Section::Evolution)].inclusive_seconds;
    const auto steps = counters[static_cast<std::size_t>(Counter::Steps)];
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << "\nFLARE performance (" << (detailed() ? "detail" : "summary") << ", "
         << (success ? "completed" : "failed; partial timings") << ")\n"
         << std::fixed << std::setprecision(6) << "Wall time: " << wall << " s; evolution: "
         << evolution << " s\n";
    for (std::size_t i = 0; i < counter_count; ++i)
        text << counter_names[i] << "=" << counters[i] << (i + 1 == counter_count ? '\n' : ' ');
    if (steps && evolution > 0)
        text << "Evolution throughput: " << std::setprecision(2) << steps / evolution
             << " steps/s; " << 1000 * evolution / steps << " ms/step (including output)\n";
    text << "Inclusive includes nested scopes; self excludes them. Self % uses total wall time.\n"
         << "Do not sum inclusive times. Snapshot timing includes its diagnostics and user output.\n";
    if (detailed())
        text << "Detail mode times EVERY face reconstruction/Riemann call; overhead can be significant.\n";
    else
        text << "Use --performance-detail to split face_fluxes into reconstruction and Riemann.\n";
    text << std::left << std::setw(34) << "Section" << std::right << std::setw(12) << "Calls"
         << std::setw(16) << "Inclusive(s)" << std::setw(16) << "Self(s)" << std::setw(12) << "Self %" << '\n';
    for (std::size_t i = 0; i < section_count; ++i)
    {
        const auto& item = totals[i];
        if (!item.calls) continue;
        text << std::left << std::setw(34) << names[i] << std::right << std::setw(12) << item.calls
             << std::fixed << std::setprecision(6) << std::setw(16) << item.inclusive_seconds
             << std::setw(16) << item.self_seconds << std::setprecision(2) << std::setw(12)
             << (wall > 0 ? 100 * item.self_seconds / wall : 0) << '\n';
    }
    out << text.str();
}

void write_json(std::ostream& out, bool success)
{
    require_closed();
    std::ostringstream text;
    text.imbue(std::locale::classic());
    text << std::setprecision(17)
         << "{\n  \"schema_version\": 1,\n  \"status\": \"" << (success ? "completed" : "failed")
         << "\",\n  \"mode\": \"" << (detailed() ? "detail" : (enabled() ? "summary" : "off"))
         << "\",\n  \"clock\": \"steady_clock_wall\",\n  \"wall_seconds\": " << wall_seconds()
         << ",\n  \"counters\": {";
    for (std::size_t i = 0; i < counter_count; ++i)
        text << (i ? "," : "") << "\n    \"" << counter_names[i] << "\": " << counters[i];
    text << "\n  },\n  \"sections\": [";
    bool first = true;
    for (std::size_t i = 0; i < section_count; ++i)
    {
        const auto& item = totals[i];
        if (!item.calls) continue;
        text << (first ? "" : ",") << "\n    {\"name\": \"" << names[i] << "\", \"calls\": " << item.calls
             << ", \"inclusive_seconds\": " << item.inclusive_seconds << ", \"self_seconds\": " << item.self_seconds << "}";
        first = false;
    }
    text << "\n  ]\n}\n";
    out << text.str();
}
}
