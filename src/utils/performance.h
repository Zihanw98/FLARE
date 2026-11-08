#pragma once

#include <array>
#include <chrono>
#include <cstdint>
#include <iosfwd>
#include <utility>

/** Serial wall-clock profiling. Disabled scopes never read the clock or allocate.
 * Timings include waiting/I/O and instrumentation overhead, not just CPU time.
 * Like the current hydro wave-speed record, this API is not thread-safe.
 */
namespace performance
{
enum class Mode { Off, Summary, Detail };
enum class Section
{
    Run, Initialization, Evolution, Advance, RHS, Hydro, Faces,
    Reconstruction, Riemann, Divergence, Sources, Gravity, Geometry, Heating,
    Boundary, StateUpdate, MeshCopy, Validation, Output, SnapshotOutput,
    DiagnosticOutput, Count
};
enum class Counter { Steps, CellUpdates, RK2Attempts, RK2Rejected, Snapshots, Count };

struct Measurement
{
    std::uint64_t calls = 0;
    double inclusive_seconds = 0;
    double self_seconds = 0;
};
inline constexpr auto section_count = static_cast<std::size_t>(Section::Count);
inline constexpr auto counter_count = static_cast<std::size_t>(Counter::Count);

// Changing the mode/resetting is only permitted when no enabled scope is active.
extern Mode mode;
extern std::array<std::uint64_t, counter_count> counters;
void configure(Mode requested);
const std::array<Measurement, section_count>& measurements() noexcept;
const char* name(Section section) noexcept;
inline bool enabled() noexcept { return mode != Mode::Off; }
inline bool detailed() noexcept { return mode == Mode::Detail; }
inline void count(Counter counter, std::uint64_t amount = 1) noexcept
{
    if (enabled()) counters[static_cast<std::size_t>(counter)] += amount;
}

class Scope
{
public:
    explicit Scope(Section section, bool detail_only = false) noexcept
        : active_(enabled() && (!detail_only || detailed())), section_(section)
    {
        if (active_) start();
    }
    ~Scope() noexcept { stop(); }
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;
    Scope(Scope&&) = delete;
    Scope& operator=(Scope&&) = delete;
    // Idempotent; enabled scopes must stop in reverse construction order.
    void stop() noexcept
    {
        if (active_) finish();
    }
private:
    void start() noexcept;
    void finish() noexcept;
    bool active_;
    Section section_;
    Scope* parent_ = nullptr;
    std::chrono::steady_clock::time_point start_;
    double children_seconds_ = 0;
};

// No callable wrapping/allocation or clock reads when profiling is disabled.
template<class Function>
decltype(auto) measure(Section section, Function&& function, bool detail_only = false)
{
    Scope scope(section, detail_only);
    return std::forward<Function>(function)();
}

// Reports require all scopes to be closed. Use independent streams internally,
// preserving the caller's formatting. JSON status is "completed" or "failed".
void report(std::ostream& out, bool success);
void write_json(std::ostream& out, bool success);
}
